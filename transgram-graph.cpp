//#define GFF_DEBUG 1 //debugging guides loading
#include <iostream>
#include <fstream>
#include "rlink.h"
#include "tmerge.h"
//#define NOTHREADS=0 //yuting New LR
#ifndef NOTHREADS
#include "GThreads.h"
#endif

//#define GMEMTRACE 1  //debugging mem allocation

#ifdef GMEMTRACE
#include "proc_mem.h"
#endif

#define VERSION "2.1.7"

//#define DEBUGPRINT 1

#ifdef DEBUGPRINT
#define DBGPRINT(x) GMessage(x)
#define DBGPRINT2(a,b) GMessage(a,b)
#define DBGPRINT3(a,b,c) GMessage(a,b,c)
#define DBGPRINT4(a,b,c,d) GMessage(a,b,c,d)
#define DBGPRINT5(a,b,c,d,e) GMessage(a,b,c,d,e)
#else
#define DBGPRINT(x)
#define DBGPRINT2(a,b)
#define DBGPRINT3(a,b,c)
#define DBGPRINT4(a,b,c,d)
#define DBGPRINT5(a,b,c,d,e)
#endif
#define USAGE "./NewAssembler file.bam -L -o output_dir"
//---- globals

std::ofstream out;
std::ofstream out_raw;
FILE* f_out=NULL;
FILE* c_out=NULL;
//#define B_DEBUG 1
#ifdef B_DEBUG
 FILE* dbg_out=NULL;
#endif
int rg_index=0;
GStr outfname;
GStr out_dir;
GStr tmp_path;
GStr tmpfname;
GStr genefname;
GStr traindir; // training directory for CDS option
bool guided=false;
bool trim=true;
bool viral=false;
bool eonly=false; // parameter -e ; for mergeMode includes estimated coverage sum in the merged transcripts
bool longreads=false;
bool rawreads=false;
bool nomulti=false;
bool enableNames=false;
bool includecov=false;
bool fr_strand=false;
bool rf_strand=false;
//bool complete=true; // set by parameter -i the reference annotation contains partial transcripts
bool retained_intron=false; // set by parameter -i for merge option
bool geneabundance=false;
//bool partialcov=false;
int num_cpus=1;
int mintranscriptlen=200; // minimum length for a transcript to be printed
//int sensitivitylevel=1;
uint junctionsupport=10; // anchor length for junction to be considered well supported <- consider shorter??
uint sserror=25; // window arround splice sites that we use to generate consensus in case of long read data
int junctionthr=1; // number of reads needed to support a particular junction
float readthr=1;     // read coverage per bundle bp to accept it; // paper uses 3
float singlethr=4.75;
uint bundledist=50;  // reads at what distance should be considered part of separate bundles
uint runoffdist=200;
float mcov=1; // fraction of bundle allowed to be covered by multi-hit reads paper uses 1
int allowed_nodes=1000;
//bool adaptive=true; // adaptive read coverage -> depends on the overall gene coverage
//GPVec<CDSparam> cds;

int no_xs=0; // number of records without the xs tag

float fpkm_thr=1;
float tpm_thr=1;

// different options of implementation reflected with the next three options
bool includesource=true;
//bool EM=false;
//bool weight=false;

float isofrac=0.01;
bool isunitig=true;
GStr label("TRGR");
GStr ballgown_dir;

GFastaDb* gfasta=NULL;

GStr guidegff; // -G option
GStr ptff; // --ptf option (point features)

bool debugMode=false;
bool verbose=false;
bool ballgown=false;

//int maxReadCov=1000000; //max local read coverage (changed with -s option)
//no more reads will be considered for a bundle if the local coverage exceeds this value
//(each exon is checked for this)

bool forceBAM = false; //useful for stdin (piping alignments into StringTie)

bool mergeMode = false; //--merge option
bool keepTempFiles = false; //--keeptmp

bool mixedMode = false; // both short and long read data alignments are provided

int GeneNo=0; //-- global "gene" counter
double Num_Fragments=0; //global fragment counter (aligned pairs)
double Frag_Len=0;
double Cov_Sum=0;
//bool firstPrint=true; //just for writing the GFF header before the first transcript is printed

GffNames* gseqNames=NULL; //used as a dictionary for reference sequence names and ids

int refseqCount=0; // number of reference sequences found in the guides file

#ifdef GMEMTRACE
 double maxMemRS=0;
 double maxMemVM=0;
 GStr maxMemBundle;
#endif


#ifndef NOTHREADS
//single producer, multiple consumers
//main thread/program is always loading the producer
GMutex dataMutex; //manage availability of data records ready to be loaded by main thread
GVec<int> dataClear; //indexes of data bundles cleared for loading by main thread (clear data pool)
GConditionVar haveBundles; //will notify a thread that a bundle was loaded in the ready queue
                           //(or that no more bundles are coming)
int bundleWork=1; // bit 0 set if bundles are still being prepared (BAM file not exhausted yet)
                  // bit 1 set if there are Bundles ready in the queue

//GFastMutex waitMutex;
GMutex waitMutex; // controls threadsWaiting (idle threads counter)

int threadsWaiting; // idle worker threads
GConditionVar haveThreads; //will notify the bundle loader when a thread
                          //is available to process the currently loaded bundle

GConditionVar haveClear; //will notify when bundle buf space available

GMutex queueMutex; //controls bundleQueue and bundles access

GFastMutex printMutex; //for writing the output to file

GFastMutex logMutex; //only when verbose - to avoid mangling the log output

GFastMutex bamReadingMutex;

GFastMutex countMutex;

#endif

GStrSet<> excludeGseqs; //hash of chromosomes/contigs to exclude (e.g. chrM)

bool NoMoreBundles=false;
bool moreBundles(); //thread-safe retrieves NoMoreBundles
void noMoreBundles(); //sets NoMoreBundles to true
//--
void processOptions(GArgs& args);

int loadPtFeatures(FILE* f, GArray<GRefPtData>& refpts);

char* sprintTime();

void processBundle(BundleData* bundle);
//void processBundle1stPass(BundleData* bundle); //two-pass testing

// yuting void writeUnbundledGuides(GVec<GRefData>& refdata, FILE* fout, FILE* gout=NULL);

#ifndef NOTHREADS

bool noThreadsWaiting();

void workerThread(GThreadData& td); // Thread function

//prepare the next free bundle for loading
int waitForData(BundleData* bundles);
#endif


TInputFiles bamreader;

int main(int argc, char* argv[]) {

 // == Process arguments.
 GArgs args(argc, argv,
   "debug;help;version;viral;conservative;mix;cds=;keeptmp;rseq=;ptf=;bam;fr;rf;merge;"
   "exclude=zihvteuLRx:n:j:s:D:G:C:S:l:m:o:a:j:c:f:p:g:P:M:Bb:A:E:F:T:");
 //args.printError(USAGE, true);

 processOptions(args);
 
 outfname=out_dir+"/transgram.reads.raw.align";
 f_out=fopen(outfname.chars(), "w");

 GStr graph_file = out_dir+"/transgram.graph";
 out.open(graph_file.chars());

 //GStr rawreads_file = out_dir+"/transgram.reads.raw.align";
 //out_raw.open(rawreads_file.chars());


 //return 0;
 GVec<GRefData> refguides; // plain vector with transcripts for each chromosome

 GArray<GRefPtData> refpts(true, true); // sorted,unique array of refseq point-features data

 //table indexes for Ballgown Raw Counts data (-B/-b option)
 GPVec<RC_TData> guides_RC_tdata(true); //raw count data or other info for all guide transcripts
 GPVec<RC_Feature> guides_RC_exons(true); //raw count data for all guide exons
 GPVec<RC_Feature> guides_RC_introns(true);//raw count data for all guide introns

 GVec<int> alncounts(30); //keep track of the number of read alignments per chromosome [gseq_id]

 int bamcount=bamreader.start(); //setup and open input files
#ifndef GFF_DEBUG
 if (bamcount<1) {
	 GError("%sError: no input files provided!\n",USAGE);
 }
#endif

#ifdef DEBUGPRINT
  verbose=true;
#endif

const char* ERR_BAM_SORT="\nError: the input alignment file is not sorted!\n";

 if(guided) { // read guiding transcripts from input gff file
	 if (verbose) {
		 printTime(stderr);
		 GMessage(" Loading reference annotation (guides)..\n");
	 }
   FILE* f=fopen(guidegff.chars(),"r");
   if (f==NULL) GError("Error: could not open reference annotation file (%s)!\n",
       guidegff.chars());
   //                transcripts_only    sort by location?
   GffReader gffr(f,       true,             true); //loading only recognizable transcript features
   gffr.setRefAlphaSorted(); //alphabetical sorting of refseq IDs
   gffr.showWarnings(verbose);
   //        keepAttrs    mergeCloseExons   noExonAttrs
   gffr.readAll(false,          true,        true);
   //the list of GffObj is in gffr.gflst, sorted by chromosome and start-end coordinates
   //collect them in other data structures, if it's kept for later call gffobj->isUsed(true)
   // (otherwise it'll be deallocated when gffr is destroyed due to going out of scope)
   refseqCount=gffr.gseqtable.Count();
   if (refseqCount==0 || gffr.gflst.Count()==0) {
	   GError("Error: could not any valid reference transcripts in %s (invalid GTF/GFF file?)\n",
			   guidegff.chars());
   }
   refguides.setCount(refseqCount); //maximum gseqid
   uint c_tid=0;
   uint c_exon_id=0;
   uint c_intron_id=0;
   GList<RC_Feature> uexons(true, false, true); //sorted, free items, unique
   GList<RC_Feature> uintrons(true, false, true);
   //assign unique transcript IDs based on the sorted order
   int last_refid=-1;
   bool skipGseq=false;
   for (int i=0;i<gffr.gflst.Count();i++) {
	   GffObj* m=gffr.gflst[i];
	   if (last_refid!=m->gseq_id) {
		   //chromosome switch
		   if (ballgown) { //prepare memory storage/tables for all guides on this chromosome
			   uexons.Clear();
			   uintrons.Clear();
		   }
		   last_refid=m->gseq_id;
		   skipGseq=excludeGseqs.hasKey(m->getGSeqName());
	   }
	   //sanity check: make sure there are no exonless "genes" or other
	   if (skipGseq) continue;
	   if (m->exons.Count()==0) {
		    if (verbose)
		    	GMessage("Warning: exonless GFF %s feature with ID %s found, added implicit exon %d-%d.\n",
		    			m->getFeatureName(), m->getID(), m->start, m->end);
		    m->addExon(m->start, m->end); //should never happen!
	   }
	   //DONE: always keep a RC_TData pointer around, with additional info about guides
	   RC_TData* tdata=new RC_TData(*m, ++c_tid);
	   m->uptr=tdata;
	   guides_RC_tdata.Add(tdata);
	   if (ballgown) { //already gather exon & intron info for all ref transcripts
		   tdata->rc_addFeatures(c_exon_id, uexons, guides_RC_exons,
		          c_intron_id, uintrons, guides_RC_introns);
	   }
	   GRefData& grefdata = refguides[m->gseq_id];
	   grefdata.add(&gffr, m); //transcripts already sorted by location
   }
	 if (verbose) {
		 printTime(stderr);
		 GMessage(" %d reference transcripts loaded.\n", gffr.gflst.Count());
	 }
 }

 gseqNames=GffObj::names; //might have been populated already by gff data
 gffnames_ref(gseqNames);  //initialize the names collection if not guided
 bool havePtFeatures=false;

 // -- loading point-feature data
 if (!ptff.is_empty()) {
   FILE* f=fopen(ptff.chars(),"r");
   if (f==NULL) GError("Error: could not open reference annotation file (%s)!\n",
       ptff.chars());
   //                transcripts_only    sort by location?
   int numptf=loadPtFeatures(f, refpts); //adds to gseqNames->gseqs accordingly, populates refpts
   havePtFeatures=(numptf>0);
   fclose(f);
 }


#ifdef GFF_DEBUG
  for (int r=0;r<refguides.Count();++r) {
	  GRefData& grefdata = refguides[r];
      for (int k=0;k<grefdata.rnas.Count();++k) {
    	  GMessage("#transcript #%d : %s (%d exons)\n", k, grefdata.rnas[k]->getID(), grefdata.rnas[k]->exons.Count());
    	  grefdata.rnas[k]->printGff(stderr);
      }
  }
  GMessage("GFF Debug mode, exiting...\n");
  exit(0);
#endif

 // --- input processing


 GHash<int> hashread;      //read_name:pos:hit_index => readlist index
 GList<GffObj>* guides=NULL; //list of transcripts on a specific reference
 GList<GPtFeature>* refptfs=NULL; //list of point-features on a specific reference
 int currentstart=0, currentend=0;
 int ng_start=0;
 int ng_end=-1;
 int ptf_idx=0; //point-feature current index in the current (*refptfs)[]
 int ng=0;
 GStr lastref;
 bool no_ref_used=true;
 int lastref_id=-1; //last seen gseq_id
 // int ncluster=0; used it for debug purposes only

 //Ballgown files
 FILE* f_tdata=NULL;
 FILE* f_edata=NULL;
 FILE* f_idata=NULL;
 FILE* f_e2t=NULL;
 FILE* f_i2t=NULL;
if (ballgown)
 Ballgown_setupFiles(f_tdata, f_edata, f_idata, f_e2t, f_i2t);
#ifndef NOTHREADS
//model: one producer, multiple consumers
#define DEF_TSTACK_SIZE 8388608
 size_t defStackSize=DEF_TSTACK_SIZE;
#ifdef _GTHREADS_POSIX_
 int tstackSize=GThread::defaultStackSize();
 if (tstackSize<DEF_TSTACK_SIZE) defStackSize=DEF_TSTACK_SIZE;
 if (verbose) {
   if (defStackSize>0){
    int ssize=defStackSize;
    GMessage("Default stack size for threads: %d (increased to %d)\n", tstackSize, ssize);
   }
   else GMessage("Default stack size for threads: %d\n", tstackSize);
 }
#endif
 GThread* threads=new GThread[num_cpus]; //bundle processing threads

 GPVec<BundleData> bundleQueue(false); //queue of loaded bundles
 //the consumers take (pop) bundles out of this queue for processing
 //the producer populates this queue with bundles built from reading the BAM input

 BundleData* bundles=new BundleData[num_cpus+1];
 //bundles[0..num_cpus-1] are processed by threads, loading bundles[num_cpus] first

 dataClear.setCapacity(num_cpus+1);
 for (int b=0;b<num_cpus;b++) {
	 threads[b].kickStart(workerThread, (void*) &bundleQueue, defStackSize);
	 bundles[b+1].idx=b+1;
	 dataClear.Push(b);
   }
 BundleData* bundle = &(bundles[num_cpus]);
#else
 BundleData bundles[1];
 BundleData* bundle = &(bundles[0]);
#endif
 GBamRecord* brec=NULL;
 bool more_alns=true;
 TAlnInfo* tinfo=NULL; // for --merge
 int prev_pos=0;
 bool skipGseq=false;
 while (more_alns) {
	 bool chr_changed=false;
	 int pos=0;
	 const char* refseqName=NULL;
	 char xstrand=0;
	 int nh=1;
	 int hi=0;
	 int gseq_id=lastref_id;  //current chr id
	 bool new_bundle=false;
	 //delete brec;
	 if ((brec=bamreader.next())!=NULL) {
		 if (brec->isUnmapped()) continue;
		 if (brec->start<1 || brec->mapped_len<10) {
			 if (verbose) GMessage("Warning: invalid mapping found for read %s (position=%d, mapped length=%d)\n",
					 brec->name(), brec->start, brec->mapped_len);
			 continue;
		 }
		 refseqName=brec->refName();
		 xstrand=brec->spliceStrand(); // tagged strand gets priority
		 if(xstrand=='.' && (fr_strand || rf_strand)) { // set strand if stranded library
			 if(brec->isPaired()) { // read is paired
				 if(brec->pairOrder()==1) { // first read in pair
					 if((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) xstrand='+';
					 else xstrand='-';
				 }
				 else {
					 if((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) xstrand='-';
					 else xstrand='+';
				 }
			 }
			 else {
				 if((rf_strand && brec->revStrand())||(fr_strand && !brec->revStrand())) xstrand='+';
				 else xstrand='-';
			 }
		 }

		 /*
		 if (xstrand=='.' && brec->exons.Count()>1) {
			 no_xs++;
			 continue; //skip spliced alignments lacking XS tag (e.g. HISAT alignments)
		 }
		 // I might still infer strand later */

		 if (refseqName==NULL) GError("Error: cannot retrieve target seq name from BAM record!\n");
		 pos=brec->start; //BAM is 0 based, but GBamRecord makes it 1-based
		 chr_changed=(lastref.is_empty() || lastref!=refseqName);
		 if (chr_changed) {
			 skipGseq=excludeGseqs.hasKey(refseqName);
			 gseq_id=gseqNames->gseqs.addName(refseqName);
			 if (guided) 
			 {
			      //
			 }

			 if (alncounts.Count()<=gseq_id) {
				 alncounts.Resize(gseq_id+1);
			 }
			 else if (alncounts[gseq_id]>0)
			           GError("%s\nAlignments (%d) already found for %s !\n",
			             ERR_BAM_SORT, alncounts[gseq_id], refseqName);
			 prev_pos=0;
		 }
		 if (pos<prev_pos) GError("%s\nread %s (start %d) found at position %d on %s when prev_pos=%d\n",
		       ERR_BAM_SORT, brec->name(), brec->start,  pos, refseqName, prev_pos);
		 prev_pos=pos;
		 if (skipGseq) continue;
		 alncounts[gseq_id]++;
		 nh=brec->tag_int("NH");
		 if (nh==0) nh=1;
		 hi=brec->tag_int("HI");
		 if (mergeMode) 
		 {
		    //
		 }

		 if (!chr_changed && currentend>0 && pos>currentend+(int)runoffdist) {
			 new_bundle=true;
		 }
	 }
	 else { //no more alignments
		 more_alns=false;
		 new_bundle=true; //fake a new start (end of last bundle)
	 }

	 if (new_bundle || chr_changed) 
	 {
		 //cerr<<endl<<".....New bundle......"<<endl;//yuting New LR
		 hashread.Clear();
		 if (bundle->readlist.Count()>0) { // process reads in previous bundle
			// (readthr, junctionthr, mintranscriptlen are globals)
			if (refptfs) { //point-features defined for this reference
				while (ptf_idx<refptfs->Count() && (int)(refptfs->Get(ptf_idx)->coord)<currentstart)
					ptf_idx++;
				//TODO: what if a PtFeature is nearby, just outside the bundle?
				while (ptf_idx<refptfs->Count() && (int)(refptfs->Get(ptf_idx)->coord)<=currentend) {
					bundle->ptfs.Add(refptfs->Get(ptf_idx)); //keep this PtFeature
					ptf_idx++;
				}
			}
			bundle->getReady(currentstart, currentend);
			if (gfasta!=NULL) { //genomic sequence data requested
				GFaSeqGet* faseq=gfasta->fetch(bundle->refseq.chars());
				if (faseq==NULL) {
					GError("Error: could not retrieve sequence data for %s!\n", bundle->refseq.chars());
				}
				bundle->gseq=faseq->copyRange(bundle->start, bundle->end, false, true);
			}
#ifndef NOTHREADS
			
			//push this in the bundle queue where it'll be picked up by the threads
			DBGPRINT2("##> Locking queueMutex to push loaded bundle into the queue (bundle.start=%d)\n", bundle->start);
			int qCount=0;
			queueMutex.lock();
			bundleQueue.Push(bundle);
			bundleWork |= 0x02; //set bit 1
			qCount=bundleQueue.Count();
			queueMutex.unlock();
			DBGPRINT2("##> bundleQueue.Count()=%d)\n", qCount);
			//wait for a thread to pop this bundle from the queue
			waitMutex.lock();
			DBGPRINT("##> waiting for a thread to become available..\n");
			while (threadsWaiting==0) {
				haveThreads.wait(waitMutex);
			}
			waitMutex.unlock();
			haveBundles.notify_one();
			DBGPRINT("##> waitMutex unlocked, haveBundles notified, current thread yielding\n");
			current_thread::yield();
			queueMutex.lock();
			DBGPRINT("##> queueMutex locked until bundleQueue.Count()==qCount\n");
			while (bundleQueue.Count()==qCount) {
				queueMutex.unlock();
				DBGPRINT2("##> queueMutex unlocked as bundleQueue.Count()==%d\n", qCount);
				haveBundles.notify_one();
				current_thread::yield();
				queueMutex.lock();
				DBGPRINT("##> queueMutex locked again within while loop\n");
			}
			queueMutex.unlock();
			

#else //no threads
			//Num_Fragments+=bundle->num_fragments;
			//Frag_Len+=bundle->frag_len;
			cerr<<"process bundle here..."<<endl;
			processBundle(bundle);
#endif
			// ncluster++; used it for debug purposes only
		 } //have alignments to process
		 else { //no read alignments in this bundle?
#ifndef NOTHREADS
			dataMutex.lock();
			DBGPRINT2("##> dataMutex locked for bundle #%d clearing..\n", bundle->idx);
#endif
			bundle->Clear();
#ifndef NOTHREADS
			dataClear.Push(bundle->idx);
			DBGPRINT2("##> dataMutex unlocking as dataClear got pushed idx #%d\n", bundle->idx);
			dataMutex.unlock();
#endif
		 } //nothing to do with this bundle

		 if (chr_changed) {
			 if (guided) {
				 ng=0;
				 guides=NULL;
				 ng_start=0;
				 ng_end=-1;
				 if (refguides.Count()>gseq_id && refguides[gseq_id].rnas.Count()>0) {
					 guides=&(refguides[gseq_id].rnas);
					 ng=guides->Count();
				 }
			 }
			 if (havePtFeatures) {
				 ptf_idx=-1;
				 //setup refptf
				 refptfs=NULL;
				 GRefPtData rd(gseq_id);
				 int ridx=refpts.IndexOf(rd);
				 if (ridx>=0) {
					refptfs=&(refpts[ridx].pfs);
					ptf_idx=0;
				 }
			 }
			 lastref=refseqName;
			 lastref_id=gseq_id;
			 currentend=0;
		 }

		 if (!more_alns) {
				if (verbose) {
#ifndef NOTHREADS
					GLockGuard<GFastMutex> lock(logMutex);
#endif
					if (Num_Fragments) {
					   printTime(stderr);
					   GMessage(" %g aligned fragments found.\n", Num_Fragments);
					}
					//GMessage(" Done reading alignments.\n");
				}
			 noMoreBundles();
			 break;
		 }
#ifndef NOTHREADS

		 int new_bidx=waitForData(bundles);
		 if (new_bidx<0) {
			 //should never happen!
			 GError("Error: waitForData() returned invalid bundle index(%d)!\n",new_bidx);
			 break;
		 }
		 bundle=&(bundles[new_bidx]);
#endif
		 currentstart=pos;
		 currentend=brec->end;
		 if (guides) 
		 { //guided and guides!=NULL
			//
		 } //guides present on the current chromosome
		bundle->refseq=lastref;
		bundle->start=currentstart;
		bundle->end=currentend;
		//cerr<<"-------Finish processing bunddle-----"<<endl;
	 } //<---- new bundle started

	 if (currentend<(int)brec->end) {
		 //current read extends the bundle
		 //this might not happen if a longer guide had already been added to the bundle
		 currentend=brec->end;
		 if (guides) { //add any newly overlapping guides to bundle
			 bool cend_changed;
			 do {
				 cend_changed=false;
				 while (ng_end+1<ng && (int)(*guides)[ng_end+1]->start<=currentend) {
					 ++ng_end;
					 //more transcripts overlapping this bundle?
					 if ((int)(*guides)[ng_end]->end>=currentstart) {
						 //it should really overlap the bundle
						 bundle->keepGuide((*guides)[ng_end],
								  &guides_RC_tdata, &guides_RC_exons, &guides_RC_introns);
						 if(currentend<(int)(*guides)[ng_end]->end) {
							 currentend=(*guides)[ng_end]->end;
							 cend_changed=true;
						 }
				 	 }
				 }
			 } while (cend_changed);
		 }
	 } //adjusted currentend and checked for overlapping reference transcripts
	 GReadAlnData alndata(brec, 0, nh, hi, tinfo);
     bool ovlpguide=bundle->evalReadAln(alndata, xstrand);
     if(!eonly || ovlpguide) { // in eonly case consider read only if it overlaps guide
    	 //check for overlaps with ref transcripts which may set xstrand
    	 if (xstrand=='+') alndata.strand=1;
    	 else if (xstrand=='-') alndata.strand=-1;
    	 //GMessage("%s\t%c\t%d\thi=%d\n",brec->name(), xstrand, alndata.strand,hi);
    	 //countFragment(*bundle, *brec, hi,nh); // we count this in build_graphs to only include mapped fragments that we consider correctly mapped
    	 //fprintf(stderr,"fragno=%d fraglen=%lu\n",bundle->num_fragments,bundle->frag_len);if(bundle->num_fragments==100) exit(0);
    	   processRead(currentstart, currentend, *bundle, hashread, alndata);
     }
 } //for each read alignment

 //cleaning up
 delete brec;

 bamreader.stop(); //close all BAM files

 if (guided && no_ref_used) {
	    GMessage("WARNING: no reference transcripts were found for the genomic sequences where reads were mapped!\n"
	    		"Please make sure the -G annotation file uses the same naming convention for the genome sequences.\n");
 }

 delete gfasta;

#ifndef NOTHREADS
 for (int t=0;t<num_cpus;t++)
	 threads[t].join();
 if (verbose) {
   printTime(stderr);
   GMessage(" All threads finished.\n");
 }
 delete[] threads;
 delete[] bundles;
#else
 if (verbose) {
    printTime(stderr);
    GMessage(" Done.\n");
 }
#endif

#ifdef B_DEBUG
 fclose(dbg_out);
#endif
// if (mergeMode && guided )
//	 writeUnbundledGuides(refguides, f_out);


 // clear refpts data, if loaded
  if (refpts.Count()>0)
	  for (int i=0;i<refpts.Count();i++) {
		  refpts[i].pfs.setFreeItem(true);
	  }

 fclose(f_out);
 if (c_out && c_out!=stdout) fclose(c_out);

 if(verbose && no_xs>0)
	 GMessage("Number spliced alignments missing the XS tag (skipped): %d\n",no_xs);

if(0 && !mergeMode) { //yuting New LR
	//yuting final output
	if(verbose) {
		GMessage("Total count of aligned fragments: %g\n", Num_Fragments);
		if (Num_Fragments)
		  GMessage("Fragment coverage length: %g\n", Frag_Len/Num_Fragments);
	}

	f_out=stdout;
	if(outfname!="stdout") { //outfname final gtf name //yuting
		f_out=fopen(outfname.chars(), "w");
		if (f_out==NULL) GError("Error creating output file %s\n", outfname.chars());
	}

	fprintf(f_out,"# ");
	args.printCmdLine(f_out);
	fprintf(f_out,"# StringTie version %s\n",VERSION);


	FILE *g_out=NULL;
	if(geneabundance) {
		g_out=fopen(genefname.chars(),"w");
		if (g_out==NULL)
			GError("Error creating gene abundance output file %s\n", genefname.chars());
		fprintf(g_out,"Gene ID\tGene Name\tReference\tStrand\tStart\tEnd\tCoverage\tFPKM\tTPM\n");
	}

	FILE* ftmp_in=fopen(tmpfname.chars(),"rt");
	if (ftmp_in!=NULL) {
		char* linebuf=NULL;
		int linebuflen=5000;
		GMALLOC(linebuf, linebuflen);
		int nl;
		int istr;
		int tlen;
		float tcov; //do we need to increase precision here ? (double)
		float calc_fpkm;
		float calc_tpm;
		int t_id;
		while(fgetline(linebuf,linebuflen,ftmp_in)) {
			sscanf(linebuf,"%d %d %d %d %g", &istr, &nl, &tlen, &t_id, &tcov);
			if (tcov<0) tcov=0;
			if (Frag_Len>0.001) calc_fpkm=tcov*1000000000/Frag_Len;
				else calc_fpkm=0.0;
			if (Cov_Sum>0.00001) calc_tpm=tcov*1000000/Cov_Sum;
				else calc_tpm=0.0;
			if(istr) { // this is a transcript
				if (ballgown && t_id>0) {
					guides_RC_tdata[t_id-1]->fpkm=calc_fpkm;
					guides_RC_tdata[t_id-1]->cov=tcov;
				}
				for(int i=0;i<nl;i++) {
					fgetline(linebuf,linebuflen,ftmp_in);
					if(!i) {
						//linebuf[strlen(line)-1]='\0';
						fprintf(f_out,"%s",linebuf);
						fprintf(f_out," FPKM \"%.6f\";",calc_fpkm);
						fprintf(f_out," TPM \"%.6f\";",calc_tpm);
						fprintf(f_out,"\n");
					}
					else fprintf(f_out,"%s\n",linebuf);
				}
			}
			else { // this is a gene -> different file pointer
				fgetline(linebuf, linebuflen, ftmp_in);
				fprintf(g_out, "%s\t%.6f\t%.6f\n", linebuf, calc_fpkm, calc_tpm);
			}
		}

		//if (guided) {
		//	writeUnbundledGuides(refguides, f_out, g_out);
		//}
		fclose(f_out);
		fclose(ftmp_in);
		if(geneabundance) fclose(g_out);
		GFREE(linebuf);
		if (!keepTempFiles) {
			remove(tmpfname.chars());
		}
	}
	else {
		fclose(f_out);
		GError("No temporary file %s present!\n",tmpfname.chars());
	}

	//lastly, for ballgown, rewrite the tdata file with updated cov and fpkm
	if (ballgown) {
		rc_writeRC(guides_RC_tdata, guides_RC_exons, guides_RC_introns,
				f_tdata, f_edata, f_idata, f_e2t, f_i2t);
	}
}

 if (!keepTempFiles) {
   tmp_path.chomp('/');
   remove(tmp_path);
 }


 gffnames_unref(gseqNames); //deallocate names collection


#ifdef GMEMTRACE
 if(verbose) GMessage(" Max bundle memory: %6.1fMB for bundle %s\n", maxMemRS/1024, maxMemBundle.chars());
#endif
} // -- END main

//----------------------------------------
char* sprintTime() {
	static char sbuf[32];
	time_t ltime; /* calendar time */
	ltime=time(NULL);
	struct tm *t=localtime(&ltime);
	sprintf(sbuf, "%02d_%02d_%02d:%02d:%02d",t->tm_mon+1, t->tm_mday,
			t->tm_hour, t->tm_min, t->tm_sec);
	return(sbuf);
}

void processOptions(GArgs& args) {


	if (args.getOpt('h') || args.getOpt("help")) {
		fprintf(stdout,"%s",USAGE);
	    exit(0);
	}
	if (args.getOpt("version")) {
	   fprintf(stdout,"%s\n",VERSION);
	   exit(0);
	}


	if (args.getOpt("viral")) {
		viral=true;
	}

	 longreads=(args.getOpt('L')!=NULL);
	 if(longreads) {
		 bundledist=0;
		 singlethr=1.5;
	 }
	 mixedMode=(args.getOpt("mix")!=NULL);
	 if(mixedMode) {
		 bundledist=0;
		 //isofrac=0.02; // allow mixedMode to be more conservative
	 }

	if (args.getOpt("conservative")) {
	  isofrac=0.05;
	  singlethr=4.75;
	  readthr=1.5;
	  trim=false;
	}

	if (args.getOpt('t')) {
		trim=false;
	}

	if (args.getOpt("fr")) {
		fr_strand=true;
	}
	if (args.getOpt("rf")) {
		rf_strand=true;
		if(fr_strand) GError("Error: --fr and --rf options are incompatible.\n");
	}

	 debugMode=(args.getOpt("debug")!=NULL || args.getOpt('D')!=NULL);
	 forceBAM=(args.getOpt("bam")!=NULL); //assume the stdin stream is BAM instead of text SAM
	 mergeMode=(args.getOpt("merge")!=NULL);
	 if(mergeMode) {
		 longreads=false; // these are not longreads
	 }
	 keepTempFiles=(args.getOpt("keeptmp")!=NULL);
	 //adaptive=!(args.getOpt('d')!=NULL);
	 verbose=(args.getOpt('v')!=NULL);
	 if (verbose) {
	     fprintf(stderr, "Running StringTie " VERSION ". Command line:\n");
	     args.printCmdLine(stderr);
	 }
	 //complete=!(args.getOpt('i')!=NULL);
	 // trim=!(args.getOpt('t')!=NULL);
	 includesource=!(args.getOpt('z')!=NULL);
	 //EM=(args.getOpt('y')!=NULL);
	 //weight=(args.getOpt('w')!=NULL);

	 GStr s=args.getOpt('m');
	 if (!s.is_empty()) {
	   mintranscriptlen=s.asInt();
	   if (!mergeMode) {
		   if (mintranscriptlen<30)
			   GError("Error: invalid -m value, must be >=30)\n");
	   }
	   else if (mintranscriptlen<0) GError("Error: invalid -m value, must be >=0)\n");
	 }
	 else if(mergeMode) mintranscriptlen=50;

	 /*
	 if (args.getOpt('S')) {
		 // sensitivitylevel=2; no longer supported from version 1.0.3
		 sensitivitylevel=1;
	 }
	*/

	 s=args.getOpt("rseq");
	 if (s.is_empty())
		 s=args.getOpt('S');
	 if (!s.is_empty()) {
		 gfasta=new GFastaDb(s.chars());
	 }

	 /*traindir=args.getOpt("cds");
	 if(!traindir.is_empty()) {
		 if(gfasta==NULL) GError("Genomic sequence file is required for --cds option.\n");
		 load_cds_param(traindir,cds);
	 }*/

     s=args.getOpt('x');
     if (!s.is_empty()) {
    	 //split by comma and populate excludeGSeqs
    	 s.startTokenize(" ,\t");
    	 GStr chrname;
    	 while (s.nextToken(chrname)) {
    		 excludeGseqs.Add(chrname.chars());
    	 }
     }

     /*
	 s=args.getOpt('n');
	 if (!s.is_empty()) {
		 sensitivitylevel=s.asInt();
		 if(sensitivitylevel<0) {
			 sensitivitylevel=0;
			 GMessage("sensitivity level out of range: setting sensitivity level at 0\n");
		 }
		 if(sensitivitylevel>3) {
			 sensitivitylevel=3;
			 GMessage("sensitivity level out of range: setting sensitivity level at 2\n");
		 }
	 }
	*/


	 s=args.getOpt('p');
	 if (!s.is_empty()) {
		   num_cpus=s.asInt();
		   if (num_cpus<=0) num_cpus=1;
	 }

	 s=args.getOpt('a');
	 if (!s.is_empty()) {
		 junctionsupport=(uint)s.asInt();
	 }

	 s=args.getOpt('j');
	 if (!s.is_empty()) junctionthr=s.asInt();

	 s=args.getOpt('E');
	 if (!s.is_empty()) sserror=s.asInt();

	 rawreads=(args.getOpt('R')!=NULL);
	 if(rawreads) {
		 if(mixedMode) {
			 GError("Mixed mode and rawreads options are incompatible!\n");
		 }

		 if(!longreads) {
			 if(verbose) GMessage("Enable longreads processing\n");
			 longreads=true;
			 bundledist=0;
		 }
		 readthr=0;

	 }

	 s=args.getOpt('c');
	 if (!s.is_empty()) {
		 readthr=(float)s.asDouble();
		 if (readthr<0.001 && !mergeMode) {
			 GError("Error: invalid -c value, must be >=0.001)\n");
		 }
	 }
	 else if(mergeMode) readthr=0;


	 s=args.getOpt('g');
	 if (!s.is_empty()) {
		 bundledist=s.asInt();
		 if(bundledist>runoffdist) runoffdist=bundledist;
	 }
	 else if(mergeMode) bundledist=250; // should figure out here a reasonable parameter for merge

	 s=args.getOpt('F');
	 if (!s.is_empty()) {
		 fpkm_thr=(float)s.asDouble();
	 }
	 //else if(mergeMode) fpkm_thr=0;

	 s=args.getOpt('T');
	 if (!s.is_empty()) {
		 tpm_thr=(float)s.asDouble();
	 }
	 //else if(mergeMode) tpm_thr=0;

	 s=args.getOpt('l');
	 if (!s.is_empty()) label=s;
	 else if(mergeMode) label="MSTRG";

	 s=args.getOpt('f');
	 if (!s.is_empty()) {
		 isofrac=(float)s.asDouble();
		 if(isofrac>=1) GError("Miminum isoform fraction (-f coefficient: %f) needs to be less than 1\n",isofrac);
	 }
	 else if(mergeMode) isofrac=0.01;
	 s=args.getOpt('M');
	 if (!s.is_empty()) {
		 mcov=(float)s.asDouble();
	 }

	 genefname=args.getOpt('A');
	 if(!genefname.is_empty()) {
		 geneabundance=true;
	 }

	 tmpfname=args.getOpt('o');

	 // coverage saturation no longer used after version 1.0.4; left here for compatibility with previous versions
	 s=args.getOpt('s');
	 if (!s.is_empty()) {
		 singlethr=(float)s.asDouble();
		 if (readthr<0.001 && !mergeMode) {
			 GError("Error: invalid -s value, must be >=0.001)\n");
		 }
	 }

	 if (args.getOpt('G')) {
	   guidegff=args.getOpt('G');
	   if (fileExists(guidegff.chars())>1)
	        guided=true;
	   else GError("Error: reference annotation file (%s) not found.\n",
	             guidegff.chars());
	 }
	 s=args.getOpt("ptf");
	 if (!s.is_empty()) {
	   ptff=s;
	   if (fileExists(ptff.chars())<=1)
		   GError("Error: point features data file (%s) not found.\n",
	             ptff.chars());
	 }

	 //enableNames=(args.getOpt('E')!=NULL);

	 retained_intron=(args.getOpt('i')!=NULL);

	 nomulti=(args.getOpt('u')!=NULL);

	 //isunitig=(args.getOpt('U')!=NULL);

	 eonly=(args.getOpt('e')!=NULL);
	 if(eonly && rawreads) {
		 if(verbose) GMessage("Error: can not use -e and -R at the same time; parameter -e will be ignored\n");
	 }
	 else if(eonly && mergeMode) {
		 eonly=false;
		 includecov=true;
	 }
	 else if(eonly && !guided)
		 GError("Error: invalid -e usage, GFF reference not given (-G option required).\n");


	 ballgown_dir=args.getOpt('b');
	 ballgown=(args.getOpt('B')!=NULL);
	 if (ballgown && !ballgown_dir.is_empty()) {
		 GError("Error: please use either -B or -b <path> options, not both.");
	 }
	 if ((ballgown || !ballgown_dir.is_empty()) && !guided)
		 GError("Error: invalid -B/-b usage, GFF reference not given (-G option required).\n");

	 /* s=args->getOpt('P');
	 if (!s.is_empty()) {
		 if(!guided) GError("Error: option -G with reference annotation file has to be specified.\n");
		 c_out=fopen(s.chars(), "w");
		 if (c_out==NULL) GError("Error creating output file %s\n", s.chars());
		 partialcov=true;
	 }
	 else { */
		 s=args.getOpt('C');
		 if (!s.is_empty()) {
			 if(!guided) GError("Error: invalid -C usage, GFF reference not given (-G option required).\n");
			 c_out=fopen(s.chars(), "w");
			 if (c_out==NULL) GError("Error creating output file %s\n", s.chars());
		 }
	 //}
	int numbam=args.startNonOpt();
#ifndef GFF_DEBUG
	if (numbam < 1 ) {
	 	 GMessage("%s\nError: no input file provided!\n",USAGE);
	 	 exit(1);
	}
#endif
	const char* ifn=NULL;
	while ( (ifn=args.nextNonOpt())!=NULL) {
		//input alignment files
		bamreader.Add(ifn);
	}
	//deferred creation of output path
	outfname="stdout";
	out_dir="./";
	//cout<<outfname<<endl; //yuting
	 if (!tmpfname.is_empty() && tmpfname!="-") { //yuting tmpfname is from -o option!!!!!
		 if (tmpfname[0]=='.' && tmpfname[1]=='/')
			 tmpfname.cut(0,2);
		 outfname=tmpfname;
		 int pidx=outfname.rindex('/');
		 if (pidx>=0) {//path given
			 out_dir=outfname.substr(0,pidx+1);
			 tmpfname=outfname.substr(pidx+1);
		 }
	 }
	 else { // stdout
		tmpfname=outfname;
		char *stime=sprintTime();
		tmpfname.tr(":","-");
		tmpfname+='.';
		tmpfname+=stime;
	 }
	 out_dir = outfname;
	 if (out_dir!="./") {
		 if (fileExists(out_dir.chars())==0) {
			//directory does not exist, create it
			if (Gmkdir(out_dir.chars()) && !fileExists(out_dir.chars())) {
				GError("Error: cannot create directory %s!\n", out_dir.chars());
			}
		 }
	 }
	 //cout<<"after creation: "<<outfname<<endl;
	 //cout<<"out_dir: "<<out_dir<<endl; //yuting
	 /* yuting 
	 { //prepare temp path
		 GStr stempl(out_dir);
		 cout<<"a "<<stempl<<endl; // cout: ./
		 stempl.chomp('/');
		 cout<<"b "<<stempl<<endl; // cout: .
		 stempl+="/tmp_XXXXXX";
		 cout<<"c "<<stempl<<endl; // cout: ./tmp_XXXXXX

		 //char* ctempl=Gstrdup(stempl.chars());
		 const char* ctempl=outfname.text();
		 Gmkdir(ctempl);
	     //Gmktempdir(ctempl);
	     tmp_path=ctempl;
	     tmp_path+='/';
	     //GFREE(ctempl);
	 }
	 */
	 tmpfname=out_dir+"/"+tmpfname;
	 //cout<<tmpfname<<endl;

#ifdef B_DEBUG
	 GStr dbgfname(tmpfname);
	 dbgfname+=".dbg";
	 dbg_out=fopen(dbgfname.chars(), "w");
	 if (dbg_out==NULL) GError("Error creating debug output file %s\n", dbgfname.chars());
#endif

	 {
		 /*yuting
		 tmpfname+=".tmp";
		 f_out=fopen(tmpfname.chars(), "w");
		 if (f_out==NULL) GError("Error creating output file %s\n", tmpfname.chars());
		 */
	 }
}

//---------------
bool moreBundles() { //getter (interogation)
	bool v=true;
#ifndef NOTHREADS
  GLockGuard<GFastMutex> lock(bamReadingMutex);
#endif
  v = ! NoMoreBundles;
  return v;
}

void noMoreBundles() {
#ifndef NOTHREADS
		bamReadingMutex.lock();
		NoMoreBundles=true;
		bamReadingMutex.unlock();
		queueMutex.lock();
		bundleWork &= ~(int)0x01; //clear bit 0;
		queueMutex.unlock();
		bool areThreadsWaiting=true;
		do {
		  waitMutex.lock();
		   areThreadsWaiting=(threadsWaiting>0);
		  waitMutex.unlock();
		  if (areThreadsWaiting) {
		    DBGPRINT("##> NOTIFY ALL workers: no more data!\n");
		    haveBundles.notify_all();
		    current_thread::sleep_for(1);
		    waitMutex.lock();
		     areThreadsWaiting=(threadsWaiting>0);
		    waitMutex.unlock();
		    current_thread::sleep_for(1);
		  }
		} while (areThreadsWaiting); //paranoid check that all threads stopped waiting
#else
	  NoMoreBundles=true;
#endif
}

void processBundle(BundleData* bundle) {
	if (verbose) {
	#ifndef NOTHREADS
		GLockGuard<GFastMutex> lock(logMutex);
	#endif
		printTime(stderr);
		GMessage(">bundle %s:%d-%d [%lu alignments (%d distinct), %d junctions, %d guides] begins processing...\n",
				bundle->refseq.chars(), bundle->start, bundle->end, bundle->numreads, bundle->readlist.Count(), bundle->junction.Count(),
                bundle->keepguides.Count());
	#ifdef GMEMTRACE
			double vm,rsm;
			get_mem_usage(vm, rsm);
			GMessage("\t\tstart memory usage: %6.1fMB\n",rsm/1024);
			if (rsm>maxMemRS) {
				maxMemRS=rsm;
				maxMemVM=vm;
				maxMemBundle.format("%s:%d-%d(%d)", bundle->refseq.chars(), bundle->start, bundle->end, bundle->readlist.Count());
			}
	#endif
	}
#ifdef B_DEBUG
	for (int i=0;i<bundle->keepguides.Count();++i) {
		GffObj& t=*(bundle->keepguides[i]);
		RC_TData* tdata=(RC_TData*)(t.uptr);
		fprintf(dbg_out, ">%s (t_id=%d) %s%c %d %d\n", t.getID(), tdata->t_id, t.getGSeqName(), t.strand, t.start, t.end );
		for (int fe=0;fe < tdata->t_exons.Count(); ++fe) {
			RC_Feature& exoninfo = *(tdata->t_exons[fe]);
			fprintf(dbg_out, "%d\texon\t%d\t%d\t%c\t%d\t%d\n", exoninfo.id, exoninfo.l, exoninfo.r,
					    exoninfo.strand, exoninfo.rcount, exoninfo.ucount);
			if (! (exoninfo==*(bundle->rc_data->guides_RC_exons->Get(exoninfo.id-1))))
				 GError("exoninfo with id (%d) not matching!\n", exoninfo.id);
		}
		for (int fi=0;fi < tdata->t_introns.Count(); ++fi) {
			RC_Feature& introninfo = *(tdata->t_introns[fi]);
			fprintf(dbg_out, "%d\tintron\t%d\t%d\t%c\t%d\t%d\n", introninfo.id, introninfo.l, introninfo.r,
					introninfo.strand, introninfo.rcount, introninfo.ucount);
			if (! (introninfo==*(bundle->rc_data->guides_RC_introns->Get(introninfo.id-1))))
				 GError("introninfo with id (%d) not matching!\n", introninfo.id);
		}
		//check that IDs are properly assigned
		if (tdata->t_id!=(uint)t.udata) GError("tdata->t_id(%d) not matching t.udata(%d)!\n",tdata->t_id, t.udata);
		if (tdata->t_id!=bundle->rc_data->guides_RC_tdata->Get(tdata->t_id-1)->t_id)
			 GError("tdata->t_id(%d) not matching rc_data[t_id-1]->t_id (%d)\n", tdata->t_id, bundle->rc_data->g_tdata[tdata->t_id-1]->t_id);

	}
#endif

	infer_transcripts(bundle); //yuting New LR

	if (ballgown && bundle->rc_data) {
		rc_update_exons(*(bundle->rc_data));
	}
	if (bundle->pred.Count()>0 || ((eonly || geneabundance) && bundle->keepguides.Count()>0)) {
#ifndef NOTHREADS
		GLockGuard<GFastMutex> lock(printMutex);
#endif
		//yuting New LR
		//if(mergeMode) GeneNo=printMergeResults(bundle, GeneNo,bundle->refseq);
		//else GeneNo=printResults(bundle, GeneNo, bundle->refseq);
	}

	if (bundle->num_fragments) {
		#ifndef NOTHREADS
				GLockGuard<GFastMutex> lock(countMutex);
		#endif
		Num_Fragments+=bundle->num_fragments;
		Frag_Len+=bundle->frag_len;
		Cov_Sum+=bundle->sum_cov;
	}

	if (verbose) {
		#ifndef NOTHREADS
				GLockGuard<GFastMutex> lock(logMutex);
		#endif
	  /*
	  SumReads+=bundle->sumreads;
	  SumFrag+=bundle->sumfrag;
	  NumCov+=bundle->num_cov;
	  NumReads+=bundle->num_reads;
	  NumFrag+=bundle->num_frag;
	  NumFrag3+=bundle->num_fragments3;
	  SumFrag3+=bundle->sum_fragments3;
	  fprintf(stderr,"Number of fragments in bundle: %g with length %g\n",bundle->num_fragments,bundle->frag_len);
	  fprintf(stderr,"Number of fragments in bundle: %g with sum %g\n",bundle->num_fragments,bundle->frag_len);
	  */
	  printTime(stderr);
	  GMessage("^bundle %s:%d-%d done (%d processed potential transcripts).\n",bundle->refseq.chars(),
	  		bundle->start, bundle->end, bundle->pred.Count());
	#ifdef GMEMTRACE
		    double vm,rsm;
		    get_mem_usage(vm, rsm);
		    GMessage("\t\tfinal memory usage: %6.1fMB\n",rsm/1024);
		    if (rsm>maxMemRS) {
			    maxMemRS=rsm;
			    maxMemVM=vm;
			    maxMemBundle.format("%s:%d-%d(%d)", bundle->refseq.chars(), bundle->start, bundle->end, bundle->readlist.Count());
		    }
	#endif
	    }
	bundle->Clear();
}

#ifndef NOTHREADS

bool noThreadsWaiting() {
	waitMutex.lock();
	int v=threadsWaiting;
	waitMutex.unlock();
	return (v<1);
}

void workerThread(GThreadData& td) {
	cerr<<"--workerThread--"<<endl;
	GPVec<BundleData>* bundleQueue = (GPVec<BundleData>*)td.udata;
	//wait for a ready bundle in the queue, until there is no hope for incoming bundles
	DBGPRINT2("---->> Thread%d starting..\n",td.thread->get_id());
	DBGPRINT2("---->> Thread%d locking queueMutex..\n",td.thread->get_id());
	queueMutex.lock(); //enter wait-for-notification loop
	while (bundleWork) {
		DBGPRINT3("---->> Thread%d: waiting.. (queue len=%d)\n",td.thread->get_id(), bundleQueue->Count());
		waitMutex.lock();
		 threadsWaiting++;
		queueMutex.unlock();
		waitMutex.unlock();
		haveThreads.notify_one(); //in case main thread is waiting
		current_thread::yield();
		queueMutex.lock();
		while (bundleWork && bundleQueue->Count()==0) {
		    haveBundles.wait(queueMutex);//unlocks queueMutex and wait until notified
		               //when notified, locks queueMutex and resume
		}
		waitMutex.lock();
		if (threadsWaiting>0) threadsWaiting--;
		waitMutex.unlock();
		DBGPRINT3("---->> Thread%d: awakened! (queue len=%d)\n",td.thread->get_id(),bundleQueue->Count());
		BundleData* readyBundle=NULL;
		if ((bundleWork & 0x02)!=0 && (readyBundle=bundleQueue->Pop())!=NULL) { //is bit 1 set?
				if (bundleQueue->Count()==0)
					 bundleWork &= ~(int)0x02; //clear bit 1 (queue is empty)
				//Num_Fragments+=readyBundle->num_fragments;
				//Frag_Len+=readyBundle->frag_len;
				queueMutex.unlock();
				processBundle(readyBundle);
				DBGPRINT3("---->> Thread%d processed bundle #%d, now locking back dataMutex and queueMutex\n",
						td.thread->get_id(), readyBundle->idx);
				dataMutex.lock();
				dataClear.Push(readyBundle->idx);
				DBGPRINT3("---->> Thread%d pushed bundle #%d into dataClear",
										td.thread->get_id(), readyBundle->idx);
				dataMutex.unlock();
				DBGPRINT2("---->> Thread%d informing main thread and yielding", td.thread->get_id());
				haveClear.notify_one(); //inform main thread
				current_thread::yield();
				DBGPRINT2("---->> Thread%d processed bundle, now locking back queueMutex\n", td.thread->get_id());
				queueMutex.lock();
				DBGPRINT2("---->> Thread%d locked back queueMutex\n", td.thread->get_id());

		}
		//haveThreads.notify_one();
	} //while there is reason to live
	queueMutex.unlock();
	DBGPRINT2("---->> Thread%d DONE.\n", td.thread->get_id());
}

//prepare the next available bundle slot for loading
int waitForData(BundleData* bundles) {
	int bidx=-1;
	dataMutex.lock();
	DBGPRINT("  #waitForData: locking dataMutex");
	while (dataClear.Count()==0) {
		DBGPRINT("  #waitForData: dataClear.Count is 0, waiting for dataMutex");
		haveClear.wait(dataMutex);
	}
	bidx=dataClear.Pop();
	if (bidx>=0) {
	  bundles[bidx].status=BUNDLE_STATUS_LOADING;
	}

	DBGPRINT("  #waitForData: unlocking dataMutex");
	dataMutex.unlock();
	return bidx;
}

#endif

void addPtFeature(const char* refname, GPtFeature* pf, GArray<GRefPtData>& refpts) {
  //expects gseqNames to be set to GffObj::names and initialized/populated already!
  //MUST be called AFTER the guides file has been loaded (if given)
  int gseq_id=gseqNames->gseqs.addName(refname);
  pf->ref_id=gseq_id;
  int ridx=-1;
  GRefPtData* rpd=NULL;
  GRefPtData rd(gseq_id);
  if (refpts.Count()==0) {
	  ridx=refpts.Add(rd);
  } else {
	  ridx=refpts.IndexOf(rd);
	  if (ridx<0) {
		  ridx=refpts.Add(rd);
	  }
  }
  if (ridx<0) GError("Error adding GRefPtData entry (bug!)\n");
  rpd = & refpts.Get(ridx);
  rpd->add(pf);
}

int loadPtFeatures(FILE* f, GArray<GRefPtData>& refpts) {
  //expected format:
  //<chromosome> <coordinate> <strand> <feature_type>
  int num=0;
  GLineReader lr(f);
  char* line=NULL;
  GDynArray<char*> tokens;
  while ((line=lr.nextLine())!=NULL) {
    strsplit(line, tokens);
    if (tokens.Count()<4)
    	GError("Error parsing point-feature line (not enough columns):\n%s\n",line);
    int start;
    if (!strToInt(tokens[1], start))
    	GError("Error parsing point-feature line (invalid coordinate):\n%s\n",line);
    int8_t strand=-2;
    if (strlen(tokens[2])==1) {
    	if (tokens[2][0]=='+')
    		strand=1;
    	else if (tokens[2][0]=='-')
    		strand=-1;
    	else if (tokens[2][0]=='.')
    		strand=0;
    }
    if (strand==-2)
    	GError("Error parsing point-feature line (invalid strand):\n%s\n",line);
    GPFType pftype=GPFT_NONE;
    if (strcmp(tokens[3], "TSS")==0)
			pftype=GPFT_TSS;
	else if (strcmp(tokens[3], "CPAS")==0)
			pftype=GPFT_CPAS;
    if (pftype==0)
    	GError("Error parsing point-feature line (unrecognized type):\n%s\n",line);
    GPtFeature* ptf=new GPtFeature(pftype, -1, strand, start);
    addPtFeature(tokens[0], ptf, refpts);
    num++;
  } //while line
  return num;
}





