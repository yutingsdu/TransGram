#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<algorithm>
#include<vector>
#include<map>
#include "simplify-graph.h"

using namespace std;
typedef vector<int> path_t;
string outdir;
ofstream out_info;
map<string, map<pair<int,int>,bool> > Chr_Junc_map;
map<string, vector<vector<int> > > Chr_Trans_map;
typedef map<string, map<pair<int,int>,bool> >::iterator Chr_Junc_map_iter;

ofstream out_gtf;
ofstream out_graph;
int pack_graph_num;
int unpack_graph_num;
bool PackingFlag;
double SampleSize;
bool unstranded;
bool SFlag;
string out_name;
int Path_length = 300;
double AVERAGE_REMOVE_RATE = 0.05;
double UNBALANCE_RATE = 0.03;
double SEED_Filter = 2;//1.01;
int rg_index = 0;
int trans_id;
bool MyFlag=false;
bool GFlag_=false;//true for beta type;
vector<int> edges_of_graph, edgesInLRP_of_graph;
vector< vector<int> > AllTrans;
vector< pair<string,string> > AllTransChr; 
vector< pair<int,int> > datatype_info;//edges; edges_NotInLRP;
/*
void get_graph_info(vector<Node> & node_set,string chr, string strand)
{
    
    out_info<<"Edges: "<<endl;
    int Graph_size = node_set.size();
    for(int i=0;i<Graph_size;i++)
    {
        Node n = node_set[i];
        for(size_t j=0;j<n.children.size();j++)
        {
            out_info<<i<<"->"<<n.children[j].first<<": "<<n.children[j].second<<endl;
        }
    }
    out_info<<"Nodes: "<<endl;
    for(int i=0;i<Graph_size;i++)
    {
        out_info<<i<<": "<<chr<<" "<<strand<<" "
                 <<node_set[i].str()<<": "<<node_set[i].coverage()
                 <<" ("<<node_set[i].JS_cov.first<<","<<node_set[i].JS_cov.second<<")"<<endl;
    }

    map<string, vector<vector<int> > >::iterator it = Chr_Trans_map.find(chr+strand);
    if(it == Chr_Trans_map.end()) return;
    vector< vector<int> > &Trans = it->second;
    int gl = node_set.front().sequence.front(), gr = node_set.back().sequence.back();
    //out_info<<"hh"<<Trans.size()<<endl;
    for(size_t i=0;i<Trans.size();i++)
    {
	//out_info<<"Trans: ";
	//for(size_t j=0;j<Trans[i].size();j++) out_info<<Trans[i][j]<<" ";
	//out_info<<endl;
	if( (Trans[i].front() >= gl && Trans[i].front()<=gr)
	    || (Trans[i].back() >= gl && Trans[i].back()<=gr)
	    || Trans[i].front() >= gl && Trans[i].back()<=gr
	    || Trans[i].front() <= gl && Trans[i].back() >= gr)
	{
	    out_info<<"Trans in graph: ";
	    out_info<<gl<<" "<<gr<<": ";
	    for(size_t j=0;j<Trans[i].size()-1;){ 
		    out_info<<Trans[i][j]<<"->"<<Trans[i][j+1]<<" ";
		    j += 2;
	    
	    }
	    out_info<<endl;
	}
	    
    }


}
*/
void load_graph_and_assemble(char* file,bool dt_data_type_flag) //dt_data_type_flag: true for determine data type
{
    if(!dt_data_type_flag)cerr<<"Loading graph "<<file<<endl;
    ifstream in(file);
    istringstream istr;
    string s,temp;
    bool edge_flag = false, node_flag = false, pair_flag = false;
    string strd = "", chr = "";
    vector<double> gene_vec;
    int start_pos=0;
    vector<int> exon_l,exon_r,junc_l,junc_r,edge_out,edge_in;
    vector<double> exon_cov, junc_cov, edge_weight;
    vector<pair<int,int> > junction, exon;
    int I = 0;
    vector<path_t> ppaths;
    vector<double> ppaths_cov;
    while(getline(in,s))
    {
	I++;
        if(I % 100000 == 0) cerr<<"Loading "<<I<<" lines"<<endl;
        if( s == "Edges"){ edge_flag = true; node_flag = false;pair_flag = false;continue;}
        if( s == "Nodes") { edge_flag = false; node_flag = true;pair_flag = false;continue;}
        if( s == "Pair") {edge_flag = false; node_flag = false;pair_flag = true;continue;}
        if( s[0] == '#') {
            istr.str(s);
	    bool invalid_flag = false;

	    //out_info<<"# Graph "<<chr<<" "<<strd<<" >>>"<<endl;
	    //out_info<<"Edges"<<endl;
	    string s_ = chr + strd;
	    Chr_Junc_map_iter it = Chr_Junc_map.find(s_);

	    if(it != Chr_Junc_map.end())
	    {
	    	map<pair<int,int>,bool> &Juncs = it->second;
	    	for(int i=0;i<edge_out.size();i++)
	    	{
		    //junction.push_back(make_pair(exon_r[edge_out[i]], exon_l[edge_in[i]]));
		    pair<int,int> junc = make_pair(exon_r[edge_out[i]], exon_l[edge_in[i]]);
		    //if(junc.second - junc.first == 1) continue;
		    //cout<<strd<<" "<<chr<<" "<<junc.first<<" "<<junc.second<<endl;
		    char Flag = 'y';
		    if(Juncs.find(junc) == Juncs.end()) Flag = 'x';
		    if(junc.second - junc.first == 1) Flag = 'p';
		    //out_info<<edge_out[i]<<" -> "<<edge_in[i]<<" : "<<edge_weight[i]<<" "<<Flag<<endl;
	    	}
	    }
	    else
	    {
	        for(int i=0;i<edge_out.size();i++)
		{
		    //out_info<<edge_out[i]<<" -> "<<edge_in[i]<<" : "<<edge_weight[i]<<" "<<"xx"<<endl;
		}
	    }
	    //out_info<<"Nodes"<<endl;
	    /*
	    for(size_t i=0;i<exon_l.size();i++)
	        out_info<<i<<" "<<exon_l[i]<<" "<<exon_r[i]<<": "<<exon_cov[i]<<endl;
	    */
	    //out_info<<"LRP"<<endl;
	    /*
	    for(size_t i=0;i<ppaths.size();i++)
	    {
	        path_t p = ppaths[i];
		for(size_t j=0;j<p.size();j++) out_info<<p[j]<<" ";
		out_info<<": "<<ppaths_cov[i]<<endl;
	    }
	    */
	    {
	    //build graph structure
	    /*
	      vector<Node> node_set;
	      for(int i=0;i < exon_l.size() ;i++)
	      {
	        Node node;
        	node.coverage_ = exon_cov[i];
        	node.length_ = exon_r[i] - exon_l[i] + 1;
		vector<int> vec;
		vec.push_back(exon_l[i]);
		vec.push_back(exon_r[i]);
		node.addseq(vec);
		node_set.push_back(node);
	      }
	      for(size_t i=0;i<edge_out.size();i++)
	      {
	          node_set[edge_out[i]].add_child(edge_in[i],edge_weight[i]);
		  node_set[edge_in[i]].add_parent(edge_out[i],edge_weight[i]);
	      }
	      get_graph_info(node_set,chr,strd);
	      */
	      SimplifyGraph simplify;
    	      simplify.build_node_set_and_simplify(exon_l,exon_r,exon_cov,edge_out,edge_in,edge_weight,ppaths,ppaths_cov,strd, chr,dt_data_type_flag);
	    }

	    istr.clear();
	    junc_l.clear();junc_r.clear();junc_cov.clear(); gene_vec.clear();
	    junction.clear();
	    exon_l.clear();exon_r.clear();exon_cov.clear();
	    exon.clear();
	    edge_out.clear();edge_in.clear();edge_weight.clear();
	    ppaths.clear(); ppaths_cov.clear();
            edge_flag = false; node_flag = false;pair_flag = false;continue;
        }

        if(edge_flag)
        {
            istr.str(s);
            int out,in;
            vector<int> v;
            double cov;
            istr>>out>>temp>>in>>temp>>cov;
            istr.clear();
            edge_out.push_back(out); edge_in.push_back(in);
            edge_weight.push_back(cov);
        }
        if(node_flag)
        {
 	    istr.str(s);
            int left, right;
            double cov;
            istr>>strd>>chr>>left>>right>>cov;
	    if( cov == 0) cov = 0.01;
            istr.clear();
            exon_l.push_back(left); exon_r.push_back(right);exon_cov.push_back(cov);
	    exon.push_back(make_pair(left,right));
	    if(gene_vec.empty()) start_pos = left;
	    int S = gene_vec.size();
	    for(int i=gene_vec.size(); i < left-start_pos;i++) gene_vec.push_back(0); //intron 
	    for(int i=gene_vec.size(); i <= right-start_pos;i++) gene_vec.push_back(cov); //exon 

        }  
	if(pair_flag)
        {
            istr.str(s);
            double cov;
            vector<int> v;
            while(istr>>temp)
            {
                if(temp == ":"){
                    istr>>cov;
                    break;
                }
                else
                    v.push_back(atoi(temp.c_str()) );
            }
            istr.clear();
	    ppaths.push_back(v);
            ppaths_cov.push_back(cov);
        }
    }
    return;
}
void load_raw_reads_align(char*file)
{
    cerr<<"Loading raw reads..."<<endl;
    ifstream in(file);
    istringstream istr;
    string s;

    //vector<int> vecExon;

    while(getline(in,s))
    {
        istr.str(s);
	string temp,chr,strand;
	double cov;
	int exon_bd;
	vector<int> vecExon;
	istr>>temp>>chr>>strand>>cov>>temp;
	while(istr>>exon_bd) vecExon.push_back(exon_bd);
	istr.clear();
	chr+=strand;
	if(1)
	{

	    sort(vecExon.begin(),vecExon.end());

	    if(Chr_Junc_map.find(chr) == Chr_Junc_map.end())
            {
                map<pair<int,int>,bool> m;
                for(size_t i=1;i<vecExon.size()-1;)
                {
                    pair<int,int> junc = make_pair(vecExon[i],vecExon[i+1]);
                    m[junc] = true;
                    i += 2;
                }
                Chr_Junc_map[chr] = m;
              } 
	    else 
	    {
                for(size_t i=1;i<vecExon.size()-1;)
                {
                    pair<int,int> junc = make_pair(vecExon[i],vecExon[i+1]);
                    Chr_Junc_map[chr][junc] = true;
                    i += 2;
                }
            }

	    if(Chr_Trans_map.find(chr) == Chr_Trans_map.end())
            {
                   vector<vector<int> > vec(1,vecExon);
                   Chr_Trans_map[chr] = vec;
            }
            else Chr_Trans_map[chr].push_back(vecExon);

	}

    }
    return;
}
void load_annotation(char* file)
{
    cerr<<"Loading annotation "<<file<<"..."<<endl;
    ifstream in(file);
    istringstream istr;
    string s;

    string chr, strand, lable;
    string tranid;
    string temp;
    int exon_l,exon_r;
    vector<int> vecExon;

    getline(in,s);
    istr.str(s);

    istr>>chr>>temp>>lable>>exon_l>>exon_r>>temp>>strand;
    while(istr>>temp) if( temp == "transcript_id") istr>>tranid;
    if(lable == "exon"){ vecExon.push_back(exon_l); vecExon.push_back(exon_r);}
    istr.clear();


    while(getline(in,s))
    {
	//cout<<s<<endl;
	istr.str(s);
 	string current_id;
	string curr_chr, curr_strand;
	int curr_el, curr_er;

	istr>>curr_chr>>temp>>lable>>curr_el>>curr_er>>temp>>curr_strand;
	if(lable != "exon") continue;

	while(istr>>temp) if( temp == "transcript_id") istr>>current_id;
	istr.clear();
	if(current_id == tranid)//commom transcript
	{
	    vecExon.push_back(curr_el); vecExon.push_back(curr_er);
	}
	else 
	{
	    vector<int> vecExon_= vecExon;
	    chr += strand;
	    if(vecExon.size() != 2)
	    {
	      sort(vecExon.begin(),vecExon.end());
	      vecExon.erase(vecExon.begin()); vecExon.pop_back();
	      //for(int i=0;i<vecExon.size();i++) cerr<<vecExon[i]<<" ";
	      //cerr<<endl;
	      for(size_t i=0;i<vecExon.size();)
	      {
		//cout<<strand<<" "<<chr<<" "<<vecExon[i]<<" "<<vecExon[i+1]<<endl;
	  	i += 2;
	      }
	      if(Chr_Junc_map.find(chr) == Chr_Junc_map.end())
	      {
		map<pair<int,int>,bool> m;
		for(size_t i=0;i<vecExon.size();)
	 	{
		    pair<int,int> junc = make_pair(vecExon[i],vecExon[i+1]);
		    m[junc] = true;
		    i += 2;
		}
		Chr_Junc_map[chr] = m;
	      } else {
		for(size_t i=0;i<vecExon.size();)
		{
		    pair<int,int> junc = make_pair(vecExon[i],vecExon[i+1]);
		    Chr_Junc_map[chr][junc] = true;
		    i += 2;
		}
	      }
	    }
	    vecExon = vecExon_;
	    //cout<<"Size: "<<vecExon.size()<<endl;
	    if(1)
	    {
	        sort(vecExon.begin(),vecExon.end());
		//cout<<vecExon.size()<<endl;
		for(size_t i=0;i<vecExon.size();)
		{
		    //cout<<strand<<" "<<chr<<" "<<vecExon[i]<<" "<<vecExon[i+1]<<endl;
		    i += 2;
		}
		if(Chr_Trans_map.find(chr) == Chr_Trans_map.end())
		{
		   //cout<<"HHH: "<<vecExon.size()<<endl;
		   vector<vector<int> > vec(1,vecExon);
		   Chr_Trans_map[chr] = vec;
		}
		else Chr_Trans_map[chr].push_back(vecExon);
	    }
	    vecExon.clear();
	    vecExon.push_back(curr_el); vecExon.push_back(curr_er);
	    chr = curr_chr; strand = curr_strand;
	    tranid = current_id;
	    
	}
    }
    chr += strand;
    if(vecExon.size() != 2) {
     vecExon.erase(vecExon.begin()); vecExon.pop_back();
     sort(vecExon.begin(),vecExon.end());
     //for(int i=0;i<vecExon.size();i++) cerr<<vecExon[i]<<" ";
     //cerr<<endl;
     for(size_t i=0;i<vecExon.size();)
     {
	break;
	cout<<strand<<" "<<chr<<" "<<vecExon[i]<<" "<<vecExon[i+1]<<endl;
	i += 2;	
     }
	      if(Chr_Junc_map.find(chr) == Chr_Junc_map.end())
              {
                map<pair<int,int>,bool> m;
                for(size_t i=0;i<vecExon.size();)
                {
                    pair<int,int> junc = make_pair(vecExon[i],vecExon[i+1]);
                    m[junc] = true;
                    i += 2;
                }
                Chr_Junc_map[chr] = m;
              } else {
                for(size_t i=0;i<vecExon.size();)
                {
                    pair<int,int> junc = make_pair(vecExon[i],vecExon[i+1]);
                    Chr_Junc_map[chr][junc] = true;
                    i += 2;
                }
              }
    }
    {
        sort(vecExon.begin(),vecExon.end());
        if(Chr_Trans_map.find(chr) == Chr_Trans_map.end())
        {
	    vector<vector<int> > vec(1,vecExon);
	    Chr_Trans_map[chr] = vec;
        }
        else Chr_Trans_map[chr].push_back(vecExon);
     }

    return;
}
 
void determine()
{
  
    int good=0, bad = 0;
    for(size_t i=0;i<datatype_info.size();i++)
    {
        double ratio = 1.0*datatype_info[i].second/(1.0*datatype_info[i].first);
	if(ratio >= 0.3) bad ++;
	else good ++;
    }
    out_info<<good<<" "<<bad<<" "<<good+bad<<endl;
    int sum = good+bad;
    if(1.0*good/(1.0*sum) > 0.9) {
	    out_info<<"beta"<<endl;
	    GFlag_ = true; //beta
    
    }
    else{
	    out_info<<"alpha"<<endl;
	    GFlag_ = false;//alpha
    }
}

int main(int argc,char* argv[])
{
    if(argc == 1)
    {
        cout<<"This is a program to searching paths in a splicing graphs with a GTF file as input(RAW reads for 3rd)."<<endl;
	cout<<"Run: ./exe raw-reads.gtf splicing.graph graph.info"<<endl;
	return 0;
    }
    outdir = argv[4];

    //out_graph.open(outdir+"/A-myGraph");
    out_gtf.open(outdir+"/transgram-temp.gtf");
    out_info.open(outdir+"/data.info");

    /*
    string datatype = argv[3];
    if(datatype == "alpha") GFlag_ = false;
    else if(datatype == "beta") GFlag_ = true;
    */

    //cerr<<"Data Type: "<<datatype<<endl;
    //out_info.open(argv[4]);
    //load_annotation(argv[1]);

    load_raw_reads_align(argv[1]);
    
    //dertemine data type
    out_info<<"Gene_id edgesNum edges_NotIn_LRP P-edgesNum ratio"<<endl;
    load_graph_and_assemble(argv[2],true);
    determine();

    //path search
    load_graph_and_assemble(argv[2],false);

    return 0;
}
