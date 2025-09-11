
#include <iostream> 
#include <map>
#include <string>
#include <vector>
#include <vector>
#include <algorithm>
#include <set>
#include <fstream>
#include "trans2path.h"

using namespace std;
extern map<string, map<pair<int,int>,bool> > Chr_Junc_map;
extern map<string, vector<vector<int> > > Chr_RawAlign_Trans_map;
extern map<string, vector<double > > Chr_RawAlign_TransCov_map;
extern map<string, vector<vector<int> > > Chr_Trans_map;
extern map<string, vector<string> > Chr_TransIDs_map;
typedef map<string, map<pair<int,int>,bool> >::iterator Chr_Junc_map_iter;
extern double SampleSize;
extern bool unstranded;
extern double SEED_Filter;
extern string mode;//R|I|U  reference/ipac/combine
extern int pack_graph_num;
extern int unpack_graph_num;
extern bool PackingFlag;
extern bool MyFlag;
extern double AVERAGE_REMOVE_RATE;
extern double UNBALANCE_RATE ;
extern bool SFlag;
extern int rg_index;
extern int trans_id;
extern string out_name;
extern ofstream out_gtf;
extern ofstream out_info;
extern ofstream out_graph;
extern string DT;

typedef vector<int> triple_t;

void SimplifyGraph::get_raw_align_info()
{
    /*
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
    */
    map<string, vector<vector<int> > >::iterator it = Chr_RawAlign_Trans_map.find(chr+strand);
    map<string, vector<double > >::iterator it2 = Chr_RawAlign_TransCov_map.find(chr+strand);
    if(it == Chr_RawAlign_Trans_map.end()) return;

    vector< vector<int> > &Trans = it->second;
    vector< double > &TransCov = it2->second;

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
	    
            //out_info<<"Trans in graph: ";
            //out_info<<gl<<" "<<gr<<": ";
	    /*
	    for(size_t j=0;j<Trans[i].size()-1;){
		    out_info<<Trans[i][j]<<"->"<<Trans[i][j+1]<<" ";
		    j += 2;
	    }
	    out_info<<endl;
	    */
	    //find path for the transcript in the graph
	    path_t p;
	    bool flag = false;
            for(size_t k=0;k<Trans[i].size()-1;)
	    { 
		
		    node_idx_t el = Trans[i][k], er = Trans[i][k+1];
		    //for each exon in the transcript, Do:
		    //out_info<<" exon: "<<el<<" - "<<er<<endl;
		    int J0 = 0;
		    if(k == 0) //first exon
		    {
			//out_info<<" * k=0: "<<k<<endl;
		        int J = -1;
			for(size_t j=0;j<nodes.size();j++)
			{
			    //out_info<<" a: "<<nodes[j].second<<endl;
			    if(nodes[j].second == er)
			    {
			        J = j;
				break;
			    }
			}
			//out_info<<" here:"<<J<<endl;
			if(J == -1) break; //not good trans
			bool partial_flag = true;//current support partial
			if(J == 0) {
			    p.push_back(J);
			    flag = true;//this exon is good;
			}
			else if(J > 0) {
			    for(size_t j=0;j<J;j++)
			    {
				//out_info<<" exon-graph: "<<nodes[j].first<<" - "<<nodes[j].second<<endl;
			        if(nodes[j].second < el) continue;
				if(nodes[j+1].first - nodes[j].second == 1){
				    p.push_back(j);
				}
				else{
				    p.clear();
				    /* new ren
				    partial_flag = false;
				    break;
				    */
				}
			    }
			    if(partial_flag) p.push_back(J);
			}
			J0=J;
			if(partial_flag) flag = true;//a good exon
			//out_info<<" 0## "<<p.size()<<endl;
		    }
		    else if(flag && k < Trans[i].size()-2)
		    {
			//out_info<<" * k<Trans[i].size()-2: "<<k<<endl;
		        flag = false;//I am not sure if this seg is good;
			int J1 = -1, J2 = -1;
			for(size_t j=J0;j<nodes.size();j++)
			{
			    if(nodes[j].first == el) J1 = j;
			    if(nodes[j].second == er){
			        J2 = j;
				break;
			    }
			}
			if(J1 == -1 || J2 == -1 || J1>J2) break;//bad exon
			bool partial_flag = true;//current support partial
			if( J1 == J2)
			{
			    p.push_back(J1);
			    flag = true;//this exon is good;
			}
			else
			{
			     for(size_t j=J1;j<J2;j++)
			     {
			         if(nodes[j+1].first - nodes[j].second == 1) p.push_back(j);
				 else{
				     partial_flag = false; // bad exon
				     break;
				 }
			     }
			     if(partial_flag)  p.push_back(J2);
			     
			}
			J0=J2;
			if(partial_flag) flag = true;//this exon is a good
			//out_info<<" 1## "<<p.size()<<endl;

		    } else if(flag && k == Trans[i].size()-2){
			//out_info<<" * k==last exon "<<k<<endl;
		        flag = false;
			int J = -1;
			for(size_t j=J0;j<nodes.size();j++)
			{
			    if(nodes[j].first == el) {
			        J = j;
				break;
			    }
			}
			if(J == -1) break;
			bool partial_flag = true;
			if(J == nodes.size() - 1){
			    p.push_back(J);
			    flag = true;
			}
			else
			{
			    int j=J;
			    for(;j<nodes.size()-1;j++)
			    {
			        if(nodes[j].second >= er) break;
				
				if(nodes[j+1].first - nodes[j].second == 1)
					 p.push_back(j);
				else{//not partial

				    if(nodes[j].second<er)//new ren
				    {
				        p.push_back(j);
					flag = true;
					partial_flag = false;
				    }

				    partial_flag = false; //bad exon
				    break;
				}
			    }
			    if(partial_flag) p.push_back(j);
			}
			if(partial_flag) flag = true;//this seg is a good seg
			 //out_info<<" n## "<<p.size()<<endl;
		    }

		    k += 2;
            
            }//for(size_t k=0;k<Trans[i].size()-1;)
	
	    if(!flag) p.clear();

	    //out_info<<"Path: ";
	    if(p.empty()){
		    //out_info<<"EMPTY"<<endl;
	    }
	    else{
		  for(size_t i=0;i<p.size();i++) {
			  //out_info<<p[i]<<"->";
			  //
		  }
		  raw_read_paths.push_back(p);
		  path_t pp = p;
		  if(pp.size() >= 2 )
		  {
            	    for(size_t j=0;j<pp.size()-1;)
              	    {
                	if(partial(pp[j],pp[j+1]))
                    	    pp.erase(pp.begin() + j);
                	else break;
            	    }
            	    if(pp.size() >= 2 ) 
		    {
            	      for(size_t j = pp.size()-1;j>0;j--)
            	      {
                	if(partial(pp[j-1],pp[j]))
                        	pp.erase(pp.begin() + j);
                	else break;
            	      }
		    }
		  }
		  if(!pp.empty());
		      raw_read_paths_remove_end_partial.push_back(pp);
		      raw_read_paths_cov.push_back(TransCov[i]);
		     
		  //final_paths.push_back(p);
	    }
	    //out_info<<endl;
 	    

        }//for trans in current graph
            
    }//for all Trans
}
void SimplifyGraph::get_graph_info()
{
    /*
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
    */
    map<string, vector<vector<int> > >::iterator it = Chr_Trans_map.find(chr+strand);
    map<string, vector<string> >::iterator it_ = Chr_TransIDs_map.find(chr+strand);
    if(it == Chr_Trans_map.end()) return;
    vector< vector<int> > &Trans = it->second;
    vector<string>& IDs = it_->second;

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
	    
            //out_info<<"Trans in graph: ";
            //out_info<<gl<<" "<<gr<<": ";
	    /*
	    for(size_t j=0;j<Trans[i].size()-1;){
		    out_info<<Trans[i][j]<<"->"<<Trans[i][j+1]<<" ";
		    j += 2;
	    }
	    out_info<<endl;
	    */
	    //find path for the transcript in the graph
	    path_t p;
	    bool flag = false;
            for(size_t k=0;k<Trans[i].size()-1;)
	    { 
		
		    node_idx_t el = Trans[i][k], er = Trans[i][k+1];
		    //for each exon in the transcript, Do:
		    //out_info<<" exon: "<<el<<" - "<<er<<endl;
		    int J0 = 0;
		    if(k == 0) //first exon
		    {
			//out_info<<" * k=0: "<<k<<endl;
		        int J = -1;
			for(size_t j=0;j<nodes.size();j++)
			{
			    //out_info<<" a: "<<nodes[j].second<<endl;
			    if(nodes[j].second == er)
			    {
			        J = j;
				break;
			    }
			}
			//out_info<<" here:"<<J<<endl;
			if(J == -1) break; //not good trans
			bool partial_flag = true;//current support partial
			if(J == 0) {
			    p.push_back(J);
			    flag = true;//this exon is good;
			}
			else if(J > 0) {
			    for(size_t j=0;j<J;j++)
			    {
				//out_info<<" exon-graph: "<<nodes[j].first<<" - "<<nodes[j].second<<endl;
			        if(nodes[j].second < el) continue;
				if(nodes[j+1].first - nodes[j].second == 1){
				    p.push_back(j);
				}
				else{
				    p.clear();
				    /* new ren
				    partial_flag = false;
				    break;
				    */
				}
			    }
			    if(partial_flag) p.push_back(J);
			}
			J0=J;
			if(partial_flag) flag = true;//a good exon
			//out_info<<" 0## "<<p.size()<<endl;
		    }
		    else if(flag && k < Trans[i].size()-2)
		    {
			//out_info<<" * k<Trans[i].size()-2: "<<k<<endl;
		        flag = false;//I am not sure if this seg is good;
			int J1 = -1, J2 = -1;
			for(size_t j=J0;j<nodes.size();j++)
			{
			    if(nodes[j].first == el) J1 = j;
			    if(nodes[j].second == er){
			        J2 = j;
				break;
			    }
			}
			if(J1 == -1 || J2 == -1 || J1>J2) break;//bad exon
			bool partial_flag = true;//current support partial
			if( J1 == J2)
			{
			    p.push_back(J1);
			    flag = true;//this exon is good;
			}
			else
			{
			     for(size_t j=J1;j<J2;j++)
			     {
			         if(nodes[j+1].first - nodes[j].second == 1) p.push_back(j);
				 else{
				     partial_flag = false; // bad exon
				     break;
				 }
			     }
			     if(partial_flag)  p.push_back(J2);
			     
			}
			J0=J2;
			if(partial_flag) flag = true;//this exon is a good
			//out_info<<" 1## "<<p.size()<<endl;

		    } else if(flag && k == Trans[i].size()-2){
			//out_info<<" * k==last exon "<<k<<endl;
		        flag = false;
			int J = -1;
			for(size_t j=J0;j<nodes.size();j++)
			{
			    if(nodes[j].first == el) {
			        J = j;
				break;
			    }
			}
			if(J == -1) break;
			bool partial_flag = true;
			if(J == nodes.size() - 1){
			    p.push_back(J);
			    flag = true;
			}
			else
			{
			    int j=J;
			    for(;j<nodes.size()-1;j++)
			    {
			        if(nodes[j].second >= er) break;
				
				if(nodes[j+1].first - nodes[j].second == 1)
					 p.push_back(j);
				else{//not partial

				    if(nodes[j].second<er)//new ren
				    {
				        p.push_back(j);
					flag = true;
					partial_flag = false;
				    }

				    partial_flag = false; //bad exon
				    break;
				}
			    }
			    if(partial_flag) p.push_back(j);
			}
			if(partial_flag) flag = true;//this seg is a good seg
			 //out_info<<" n## "<<p.size()<<endl;
		    }

		    k += 2;
            
            }//for(size_t k=0;k<Trans[i].size()-1;)
	
	    if(!flag) p.clear();

	    //out_info<<"Path: ";
	    if(p.empty()){
		    //out_info<<"EMPTY"<<endl;
	    }
	    else{
		  for(size_t i=0;i<p.size();i++) {
			  //out_info<<p[i]<<"->";
			  //
		  }
		  right_paths.push_back(p);
		  right_paths_ids.push_back(IDs[i]);
		  path_t pp = p;
		  if(pp.size() >= 2 )
		  {
            	    for(size_t j=0;j<pp.size()-1;)
              	    {
                	if(partial(pp[j],pp[j+1]))
                    	    pp.erase(pp.begin() + j);
                	else break;
            	    }
            	    if(pp.size() >= 2 ) 
		    {
            	      for(size_t j = pp.size()-1;j>0;j--)
            	      {
                	if(partial(pp[j-1],pp[j]))
                        	pp.erase(pp.begin() + j);
                	else break;
            	      }
		    }
		  }
		  if(!pp.empty());
		      right_paths_remove_end_partial.push_back(pp);
		  //final_paths.push_back(p);
	    }
	    //out_info<<endl;
 	    

        }//for trans in current graph
            
    }//for all Trans



}

bool sorter(pair<edge_t,double> p1, pair<edge_t,double>p2)
{
    return ((p1.second>p2.second) || (p1.second == p2.second && p1.first < p2.first));
}
  void SimplifyGraph::get_reserved_junc(vector<path_t> LongReadPath)
  {
	for(size_t i=0;i<LongReadPath.size();i++)
	{
	    path_t p = LongReadPath[i];

	    if(LongReadPath_cov[i] <= 2) continue;

	    for(size_t j=0;j<p.size()-1;j++)
	    {
		edge_t e = make_pair(p[j],p[j+1]);
		reserved_junc.push_back(e);
	    }
	}
	sort(reserved_junc.begin(),reserved_junc.end());
	reserved_junc.erase(unique(reserved_junc.begin(),reserved_junc.end()),reserved_junc.end());

	return;
  }
  bool SimplifyGraph::keep_edge(edge_t e)
  {
	if(find(reserved_junc.begin(),reserved_junc.end(),e) == reserved_junc.end()) 
		return false;

	return true;
  }
  void SimplifyGraph::show_block(block_info bi)
  {
      out_info<<"Block: "<<bi.start<<"--"<<bi.end<<endl;
      for(size_t i=0;i<bi.paths.size();i++){
          path_t p = bi.paths[i];
	  out_info<<" bpath: ";
	  for(size_t j=0;j<p.size();j++)
		  out_info<<p[j]<<" ";
	  out_info<<": "<<bi.paths_cov[i]<<endl;
      }
      return;

  }
  void SimplifyGraph::GetBlock()
  {
    //out_graph<<"Block info"<<endl;
    //out_graph<<"Graph_"<<rg_index<<" - "<<node_set.front().str()<<" - (GraphSize:"<<Graph_size<<") ";

    int bsum = 0;
    node_idx_t max_node = 0;
    vector<node_idx_t> blocks_intervals;

    //out_graph<<max_node<<"; ";
    blocks_intervals.push_back(max_node);

    for(int i=0;i<Graph_size;i++)
    {
	if(i<max_node) continue;;
	bool flag = true;//first node in a new block
	int last_node_of_last_block = max_node;
	for(int k=last_node_of_last_block; k<=max_node;k++)
	{
	        for(size_t j=0;j<node_set[k].children.size();j++)
		{
		    node_idx_t c = node_set[k].children[j].first;
		    if(c>max_node)
		    {
			    if(k<max_node || flag)
			    {
			        max_node = c;
				flag = false;
				
			    }
			    else{ //k == max_node, next block
				    break;
			    }
		    }
		}
	}

	if(last_node_of_last_block == max_node)//NO edge out from last_node_of_last_block(tail or one connected component?)
	{
	    if(max_node == Graph_size - 1) //tail of graph
	    {
		 break;
	    }
	    else //maybe next connected component e.g. graph:0->1, 2->3;no edge between 1 and 2
	    {
		//out_graph<<"(next cc) ";
	        max_node += 1;
	    } 
	}

	bsum++;
	//out_graph<<max_node<<"; ";
        blocks_intervals.push_back(max_node);

	//if(max_node == Graph_size - 1) break;
    }
    BlockNumber = bsum;
    //out_info<<"(BlockSize:"<<bsum<<")"<<endl;
    vector<path_t> FinalPaths;
    for(size_t i=0;i<blocks_intervals.size()-1;i++)
    {
	vector<path_t> Bps;
	vector<double> Bps_cov;

	//out_graph<<" Block:"<<blocks[i]<<" "<<blocks[i+1]<<endl;
	get_path_of_one_block(blocks_intervals[i], blocks_intervals[i+1],Bps,Bps_cov);
        block_info bi(blocks_intervals[i],blocks_intervals[i+1],Bps,Bps_cov);
	Blocks.push_back(bi);
	//out_info<<"Block: "<<i<<endl;
	//show_block(bi);
	/*
	if(bsum == 1 && !bi.paths.empty())
	{
	  final_paths.push_back(bi.paths.front());
	  //FinalPaths.push_back(bi.paths.front());
	}
	*/
    }
    trans_id = 1;
    //path_search(strand,chr);
    get_LRP();
    //get_full_length_path_from_blocks();
    path_search(strand,chr);
    //get_full_length_path_from_blocks();
    //get_LRP();
    output(strand,chr,final_paths);
return;
/*
    out_graph<<"HH"<<final_paths_temp.size()<<" "<<final_paths.size()<<endl;
    if(!final_paths_temp.empty() && !final_paths.empty())
    {
        if(final_paths_temp.front() != final_paths.front())
	{
	    out_graph<<"WARNING"<<endl;
	    path_t p1 = final_paths_temp.front(),p2 = final_paths.front();

	    out_graph<<"block: ";
	    for(size_t i=0;i<p1.size();i++) out_graph<<p1[i]<<" ";
            out_graph<<endl;

	    out_graph<<"edge: ";
	    for(size_t i=0;i<p2.size();i++) out_graph<<p2[i]<<" ";
            //out_graph<<endl;
	    output(strand,chr,final_paths);
	    output(strand,chr,final_paths_temp);
	}
	else output(strand,chr,final_paths);
    }
*/
/*
    if(bsum == 1 && !FinalPaths.empty() && !final_paths.empty() && FinalPaths.front() != final_paths.front())
    {
	path_t p1 = FinalPaths.front(),p2 = final_paths.front();
	for(size_t i=0;i<p1.size();i++) out_graph<<p1[i]<<" ";
	out_graph<<endl;

	for(size_t i=0;i<p2.size();i++) out_graph<<p2[i]<<" ";
	out_graph<<endl;

        out_graph<<"WARNING"<<endl;
    }
 */
    return;
  }
  void SimplifyGraph::get_full_length_path_from_blocks()
  {
    //size_t BI = -1;//for block
    //double cov = 0;

    int iteration = 1;
    while(1)
    {
	//out_graph<<"iteration: "<<iteration<<"..."<<endl;
	int BI = -1;//for block
	int BPI=-1; //for the path in a block
	double cov = 0;
        for(size_t i=0;i<Blocks.size();i++)
	{
	    if(Blocks[i].paths.empty()) continue;
	    int index = Blocks[i].get_index_of_maxcov_path(false);//used can't be used again
	    if(index == -1) continue;
	    if(Blocks[i].paths_cov[index] > cov)
	    {
	        cov = Blocks[i].paths_cov[index];
		BI = i;
		BPI = index;
	    }
	} 
	//out_graph<<"A:"<<BI<<endl;
	if(cov <4) break;
	if(BI == -1) break;
	Blocks[BI].paths_flag[BPI] = 1;
	path_t p = Blocks[BI].paths[BPI];
	for(int i=BI-1;i>=0;i--)
	{
	    if(Blocks[i].end == p.front() && !Blocks[i].paths.empty())//reverse extend
	    {
	        int index = Blocks[i].get_index_of_maxcov_path(true);
		if(index == -1) continue;

		path_t temp = p;
		p = Blocks[i].paths[index];
		for(size_t j=1;j<temp.size();j++) p.push_back(temp[j]);

		Blocks[i].paths_flag[index] = 1;//this path is used
	    }
	}
	for(int i=BI+1;i<Blocks.size();i++)
	{
	    if(Blocks[i].start == p.back() && !Blocks[i].paths.empty())//forward extend
	    {
	        int index = Blocks[i].get_index_of_maxcov_path(true);
		if(index == -1) continue;
		for(size_t j=1;j<Blocks[i].paths[index].size();j++)
			p.push_back(Blocks[i].paths[index][j]);

		Blocks[i].paths_flag[index] = 1;//this path is used
	    }
	}
	iteration++;
	final_paths.push_back(p);
	//break;

    }

    size_t I = -1;
    double cov = 0;
    while(0)
    {
      for(size_t i=0;i<Blocks.size();i++)
      {
	if(Blocks[i].paths.empty()) continue;
	if(Blocks[i].paths_cov.front() > cov)
	{
	    I = i;
	    cov = Blocks[i].paths_cov.front();
	}
    
      }
      if(I == -1) return;
      path_t p = Blocks[I].paths.front();
      for(int i=I-1;i>=0;i--)
      {
        //current_path.insert(current_path.begin(), e.first);
	if(Blocks[i].end == p.front() && !Blocks[i].paths.empty())
	{
	    path_t temp = p;
	    p = Blocks[i].paths.front();
	    for(size_t j=1;j<temp.size();j++) p.push_back(temp[j]);
	}
      }
      for(int i=I+1;i<Blocks.size();i++)
      {
         if(Blocks[i].start == p.back() && !Blocks[i].paths.empty())
	 {
	     for(size_t j=1;j<Blocks[i].paths.front().size();j++)
		     p.push_back(Blocks[i].paths.front()[j]);
	 }
      }

      final_paths.push_back(p);
      break;

    }//extend to get all paths
    /*
	while(1)
	{
	    for(int j=i+1;j<Blocks.size();j++)
	    {
	        if(Blocks[j].start = p.back() && !Blocks[j].paths.empty())
		{
		    for(size_t k=1;k<Blocks[j].paths.front().size();k++)
		    {
		        p.push_back(Blocks[j].paths.front()[k]);
		    }
		}
		else 
		{
		     break;
		}
	    }
	    break;
	}
    */
    //output(strand,chr,final_paths_temp);
  }
  void SimplifyGraph::get_path_of_one_block(node_idx_t start, node_idx_t end,vector<path_t>& Bpaths, vector<double>& Bpaths_cov)
  {

    vector<edge_t> seeds;
    vector<double> seeds_cov;
    //get seed from block B[start,end]
    for(int i = start;i<=end;i++)
    {
	//break;
        for(size_t j=0;j<node_set[i].children.size();j++)
	{
	    int c = node_set[i].children[j].first;
            double cov = node_set[i].children[j].second;
	    if(partial(i,c)) continue;
	    if(c > end) continue;
	    edge_t edge = make_pair(i,c);
	    seeds.push_back(edge);
	    seeds_cov.push_back(cov);
	}
    }
    if( seeds.empty())//only partial edge exist??
    {
        for(int i = start;i<=end;i++)
	{
	    for(size_t j=0;j<node_set[i].children.size();j++)
	    {
	        int c = node_set[i].children[j].first;
		double cov = node_set[i].children[j].second;
		if(c > end) continue;
		edge_t edge = make_pair(i,c);
		seeds.push_back(edge);
		seeds_cov.push_back(cov);
	    }
	}
    }
    if(seeds.empty()) return;
    double minSeed = 5.5; 
    map<edge_t, bool> used_edges;
    while(1)//edge as seed
    {
        vector<double>::iterator biggest = max_element(seeds_cov.begin(),seeds_cov.end());
	int k = biggest - seeds_cov.begin();
	double seed_cov = seeds_cov[k];
	edge_t seed = seeds[k];
	if(seed_cov < minSeed) break;
	else seeds_cov[k] = 0;

	if(used_edges.find(seed) != used_edges.end()) continue;

	vector<path_t> forward_paths, reverse_paths;
	forward_extend(seed,forward_paths,end);
	reverse_extend(seed,reverse_paths,start);
	path_t final_path = reverse_paths.front();

	for(size_t i=2;i<forward_paths.front().size();i++) final_path.push_back(forward_paths.front()[i]);
	
	double cov_ = 0;
	for(size_t i=0;i<final_path.size()-1;i++)
	{
	    edge_t edge = make_pair(final_path[i],final_path[i+1]);
	    used_edges[edge] = true;
	    cov_ += node_set[final_path[i]].get_child_coverage(final_path[i+1]);
	}
	cov_ = cov_/(1.0*final_path.size()-1);

	Bpaths.push_back(final_path);
	Bpaths_cov.push_back(seed_cov);
	//Bpaths_cov.push_back(cov_);
    }
    return;
  }
  void SimplifyGraph::forward_extend(edge_t e, vector<path_t>& paths, node_idx_t end)
  {
      //out_graph<<"forward extend..."<<endl;
      //out_graph<<e.first<<" "<<e.second<<endl;
      path_t path_temp;
      path_temp.push_back(e.first);
      path_temp.push_back(e.second);
      vector<path_t> paths_temp;
      paths_temp.push_back(path_temp);
      while(1)
      {
	  if(paths_temp.empty()) break;

	  path_t current_path = paths_temp.front();
	  paths_temp.erase(paths_temp.begin());

          edge_t current_edge = make_pair(current_path[current_path.size() - 2], current_path[current_path.size() - 1]);
	  //out_graph<<"* "<<current_edge.first<<" "<<current_edge.second<<endl;
	  //if(node_set[current_edge.second].children.size() == 0)
	  if(current_edge.second == end || node_set[current_edge.second].children.size() == 0) //for block
	  {
	      paths.push_back(current_path);    
	  }
	  else
	  {
	      vector<pair<edge_t,double> > edges_cov;
	      for(size_t i=0;i<node_set[current_edge.second].children.size() ;i++)
	      {
		  edge_t e = make_pair(current_edge.second, node_set[current_edge.second].children[i].first);
		  edges_cov.push_back(make_pair(e,node_set[current_edge.second].children[i].second));
	      }
	      sort(edges_cov.begin(),edges_cov.end(),sorter);
	      if(edges_cov.size() == 1)
	      {
	          edge_t e = edges_cov[0].first;
		  current_path.push_back(e.second);
		  paths_temp.push_back(current_path);
	      }
	      else {
	        for(size_t i=0;i<edges_cov.size();i++)
	        {
	          edge_t e = edges_cov[i].first;
		  //if(!partial(e.first,e.second)) continue;
	          pair<edge_t ,edge_t> ee = make_pair(current_edge, e);
		  if(1 || packing_map.find(ee) != packing_map.end())
		  {
		      current_path.push_back(e.second);
		      paths_temp.push_back(current_path);
		      break;
		  }
		}
	      }
	  }
      }
      return;

  }

  void SimplifyGraph::reverse_extend(edge_t e, vector<path_t>& paths,node_idx_t start)
  {
      path_t path_temp;
      path_temp.push_back(e.first);
      path_temp.push_back(e.second);
      vector<path_t> paths_temp;
      paths_temp.push_back(path_temp);
      while(1)
      {
	  if(paths_temp.empty()) break;

	  path_t current_path = paths_temp.front();
	  paths_temp.erase(paths_temp.begin());

          edge_t current_edge = make_pair(current_path[0], current_path[1]);
	  if(current_edge.first == start || node_set[current_edge.first].parents.size() == 0)//for block
	  {
	      paths.push_back(current_path);    
	  }
	  else
	  {
	      vector<pair<edge_t,double> > edges_cov;
	      for(size_t i=0;i<node_set[current_edge.first].parents.size() ;i++)
	      {
		  edge_t e = make_pair(node_set[current_edge.first].parents[i].first,current_edge.first);
		  edges_cov.push_back(make_pair(e,node_set[current_edge.first].parents[i].second));
	      }
	      sort(edges_cov.begin(),edges_cov.end(),sorter);
	      if(edges_cov.size() == 1)
	      {
	          edge_t e = edges_cov[0].first;
		  current_path.insert(current_path.begin(), e.first);
		  paths_temp.push_back(current_path);
	      }
	      else 
	      {
	        for(size_t i=0;i<edges_cov.size();i++)
	        {
	          edge_t e = edges_cov[i].first;
		  //if(!partial(e.first,e.second)) continue;
	          pair<edge_t ,edge_t> ee = make_pair(e,current_edge);
		  if(1|| packing_map.find(ee) != packing_map.end())
		  {
		      //current_path.push_back(e.second);
		      current_path.insert(current_path.begin(), e.first);
		      paths_temp.push_back(current_path);
		      break;
		  }
		}
	      }
	  }
      }
      return;
  }

  void SimplifyGraph::get_CovInfo()
  {
    //cout<<"Info: "<<Graph_size<<" ";
    int good=0,bad=0;
    for(int i=0;i<Graph_size;i++)
    {
	double Ccov=0,Pcov=0;
	if(node_set[i].children.size()>0 && node_set[i].parents.size()>0)
	{
	    if(node_set[i].children.size()>1 || node_set[i].parents.size()>1)
	    {
		for(size_t j=0;j<node_set[i].children.size();j++) Ccov += node_set[i].children[j].second;
		for(size_t j=0;j<node_set[i].parents.size();j++) Pcov += node_set[i].parents[j].second;
		double min = Ccov<Pcov?Ccov:Pcov;
		double max = Ccov<Pcov?Pcov:Ccov;
		double rate = min/max;
		//cout<<rate<<" ";
		if(rate>0.7) good++;
		else bad++;
	    }
	}
    }
    //cout<<"; "<<good<<"/"<<bad;
    //cout<<endl;
    if(good > 5*bad) PackingFlag = true;
    if(PackingFlag) pack_graph_num++;
    else unpack_graph_num++;
  }
  node_idx_t SimplifyGraph::find_partial(node_idx_t n, bool downstream)
  {
    if(n >= Graph_size) return false;
    if(downstream)
    {
    	for(size_t i=0;i< node_set[n].children.size();i++)
    	{
            if(partial(n,node_set[n].children[i].first))
	        return node_set[n].children[i].first;
    	}
    }
    else
    {
        for(size_t i=0;i< node_set[n].parents.size();i++)
	{
	    if(partial(node_set[n].parents[i].first,n))
		    return node_set[n].parents[i].first;
	}
    }

    return -1;

  }
  void SimplifyGraph::get_end_node(map<int,bool>& head, map<int,bool>& tail)
  {
    for(int i=0;i < Graph_size;i++)
    {
        if(node_set[i].parents.empty()) head[i] = true;
	if(node_set[i].children.empty()) tail[i] = true;
/*
	//new1
	node_idx_t p = find_partial(i,false);
	if(p != -1) head[i] = true;
	node_idx_t c = find_partial(i,true);
	if(c != -1) tail[i] = true;
*/	

	//new2 consider, 1-2-3-4...all partial
	node_idx_t current = i;
	while(1)
	{
	    node_idx_t p = find_partial(current,false);
	    if( p != -1)//partial
	    {
	        if( node_set[p].parents.empty() ) //p is head in graph, then i is head
		{
		    head[i] = true;
		    break;
		}
		else //p is not head
		    current = p;
	    }
	    else break;//i have no partial parent
	}
	while(1)
	{
	    node_idx_t c = find_partial(current,true);
	    if(c != -1) //partial
	    {
		if( node_set[c].children.empty() ) //c is tail in graph, then i is tail
		{
		    tail[i] = true;
		    break;
		}
		else //c is not tail
		    current = c;
	    }
	    else break; //i have no partial children
	}
	
	/*
	for(size_t j = 0;j<node_set[i].parents.size();j++)
	{
	    if(partial(node_set[i].parents[j].first,i))
	    {
	        head[i] = true;
		break;
	    }
	}
	node_idx_t current = i;
	for(size_t j = 0;j<node_set[current].children.size();j++)
	{
		node_idx_t c = node_set[i].children[j].first;
		if(partial(i,c)) 
		{
		    tail[i] = true;
		    break;
		}
	}
	*/
    }
    return;
  }
  void SimplifyGraph::contract_LRP()
  {
	//cout<<node_set.front().str()<<endl;
	for(size_t i=0;i<LongReadPath.size();i++)
	{
	    /*
	    cout<<"** "<<i<<endl;
	    for(size_t j=0;j<LongReadPath[i].size();j++) cout<<LongReadPath[i][j]<<" ";
	    cout<<endl;
	    */
	    /*
	    for(size_t j=0;j<LongReadPath[i].size();j++) 
		    LongReadPath[i][j] += 1;
	    */
	    if(LongReadPath[i].size() < 2 ) continue;
	    for(size_t j=0;j<LongReadPath[i].size()-1;)
	    {
	        if(partial(LongReadPath[i][j],LongReadPath[i][j+1]))
		    LongReadPath[i].erase(LongReadPath[i].begin() + j);
		else break;
	    }
	    if(LongReadPath[i].size() < 2 ) continue;
	    for(size_t j = LongReadPath[i].size()-1;j>0;j--)
	    {
	        if(partial(LongReadPath[i][j-1],LongReadPath[i][j]))
			LongReadPath[i].erase(LongReadPath[i].begin() + j);
		else break;
	    }
	}
  }
  void SimplifyGraph::add_source_and_sink()
  {
    for(int i = 1;i<=exon_number;i++)
    {
        if(node_set[i].parents.empty())
	{
	    double cov=0;
	    for(size_t j=0;j<node_set[i].children.size();j++) 
		    cov += node_set[i].children[j].second;
	    node_set.front().add_child(i,cov);
	    node_set[i].add_parent(0,cov);
	}
	if(node_set[i].children.empty())
	{
	    double cov=0;
	    for(size_t j=0;j<node_set[i].parents.size();j++)
		    cov += node_set[i].parents[j].second;

	    node_set.back().add_parent(i,cov);
	    node_set[i].add_child(Graph_size-1,cov);
	}
    }
    return;
  }
  void SimplifyGraph::build_node_set_and_simplify(vector<int>& exon_l, vector<int>& exon_r, vector<double>& exon_cov,
			       	vector<node_idx_t>&edge_out ,vector<node_idx_t>&edge_in, vector<double>& edge_weight,
				vector<path_t>& LRP, vector<double> LRP_cov,
				string strand_, string chr_)
  {
    SSFlag = false;//add source&sink or not

    exon_number = exon_l.size();

    if(SSFlag) Graph_size = exon_l.size()+2;// add source and sink
    else Graph_size = exon_l.size();

    size_ = Graph_size;
    RawSize = Graph_size;
    strand = strand_;
    chr = chr_;

    
    //source
    if(SSFlag)
    {
      Node source;
      source.coverage_ = 1;
      source.length_ = 1;
      vector<int> vec_;
      vec_.push_back(exon_l.front()-2);
      vec_.push_back(exon_l.front()-1);
      source.addseq(vec_);
      node_set.push_back(source);
    }
    
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
	nodes.push_back(make_pair(exon_l[i],exon_r[i]));
    }
    /*
    out_info<<"****"<<endl;
    for(size_t i=0;i<nodes.size();i++)
    {
        out_info<<"*"<<nodes[i].first<<" "<<nodes[i].second<<endl;
    }
    out_info<<"****"<<endl;
    */
    //sink
    if(SSFlag)
    {
      Node sink;
      sink.coverage_ = 1;
      sink.length_ = 1;
      vector<int> vec_;
      vec_.push_back(exon_r.back()+1);
      vec_.push_back(exon_r.back()+2);
      sink.addseq(vec_);
      node_set.push_back(sink);
    }
    
    for(size_t i=0;i<edge_out.size();i++)
    {
      if(edge_weight[i]>=1)
      {
	      
	//add_source_and_sink
	if(SSFlag)
	{
	  node_set[edge_out[i]+1].add_child(edge_in[i]+1,edge_weight[i]);
 	  node_set[edge_in[i]+1].add_parent(edge_out[i]+1,edge_weight[i]);
	}
	else 
	{
	  node_set[edge_out[i]].add_child(edge_in[i],edge_weight[i]);
	  node_set[edge_in[i]].add_parent(edge_out[i],edge_weight[i]);
	  edges.push_back(make_pair(edge_in[i],edge_out[i]));
	}
      }

    }

    if(SSFlag) add_source_and_sink();
    get_graph_info();
    get_raw_align_info();

    LongReadPath = LRP;
    LongReadPath_cov = LRP_cov;


    //process graph
    vector<edge_t> delete_edges_child,delete_edges_parent;

    remove_intron_contamination(delete_edges_child,delete_edges_parent);
    AVERAGE_REMOVE_RATE = 0.01;
    remove_edges_by_average_coverage(delete_edges_child,delete_edges_parent);
    

    delete_children_edges(delete_edges_child);
    delete_parents_edges(delete_edges_parent);

    process_LRP_after_simplify_graph();

    get_LRP();//get long read path to get data type

    vector<pair<int,int> > unused_junctions;
    int edge_number = 0,partial_edge_number = 0,edges_Notin_LRP=0;
    double max_cov = 0, ave_cov = 1; 
    for(int i = 0;i< Graph_size;i++) //add source and sink
    {   
        for(size_t j=0;j<node_set[i].children.size();j++)
        {   
            if(partial(i,node_set[i].children[j].first)){
                    partial_edge_number++;
                    continue;
            }
            edge_number++;
            int c = node_set[i].children[j].first;
            
            double cov = node_set[i].children[j].second;
            if(cov > max_cov) max_cov = cov;
            ave_cov += cov;
            pair<int,int> junc_ = make_pair(i,node_set[i].children[j].first);
            unused_junctions.push_back(junc_);
        }
    
    }
    if(edge_number >= 1) ave_cov = ave_cov/(1.0*edge_number);

    map<edge_t,bool> used_edges;
    for(size_t i=0;i<final_paths.size();i++)
    {
        path_t p = final_paths[i];
        for(size_t j=0;j<p.size()-1;j++)
        {
            edge_t edge = make_pair(p[j],p[j+1]);
            used_edges[edge] = true;
        }
    }

    for(size_t i=0;i<unused_junctions.size();i++)
    {
        pair<int,int> junc = unused_junctions[i];
        if(used_edges.find(junc) == used_edges.end())
        {
            edges_Notin_LRP++;
        }
    }
    if(edge_number >=2)
    {
	
        double flag_ratio = 1.0*edges_Notin_LRP/(1.0*edge_number);
	/*
        cout<<"Gene_"<<rg_index<<" "<<edge_number<<" "<<edges_Notin_LRP<<" "<<partial_edge_number<<" "
                <<1.0*edges_Notin_LRP/(1.0*edge_number)<<" "
                <<max_cov<<" "<<ave_cov
                <<endl;*/
		
    }
    double R=1.0;
    if(DT == "alpha") R=0.0;
    out_info<<"Assembled_Paths:"<<endl;
    for(size_t i=0;i<right_paths.size();i++)
    {
        path_t p = right_paths[i];
        remove_partial_end_for_a_path(p);//only consider intron chain
	vector<vector<int> >::iterator it = find(raw_read_paths_remove_end_partial.begin(),raw_read_paths_remove_end_partial.end(),p);
	double RawReadsCov = 0;
	if(it != raw_read_paths_remove_end_partial.end())
	{
	    int k = it - raw_read_paths_remove_end_partial.begin();
	    RawReadsCov = raw_read_paths_cov[k];
	}
        out_info<<right_paths_ids[i]<<" ";
        for(size_t j=0;j<p.size();j++) out_info<<p[j]<<"->";
        out_info<<" "<<R*RawReadsCov<<" "<<edges_Notin_LRP<<" "<<edge_number<<" "<<1.0*edges_Notin_LRP/(1.0*edge_number)<<endl;
    }
    return;


return;
    int partial_number = 0;
    for(int i=0;i<Graph_size-1;i++)
    {
	if(partial(i,i+1)) partial_number++;
    }
    if(partial_number < 10) get_reserved_junc(LRP);
    return;
  }
  void SimplifyGraph::remove_all_partial(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent)
  {
    for(int i=0;i<Graph_size;i++)
    {
	for(size_t j = 0;j<node_set[i].children.size();j++)
	{
	    node_idx_t c = node_set[i].children[j].first;
	    if(MultiJunction(i,c)) continue;
	    if(partial(i,c) && node_set[i].children[j].second < 0.05)
	    {
		delete_edges_child.push_back(make_pair(i,c));
	    }
	}
    }
    return;
  }
  double SimplifyGraph::compute_JScov_for_Node(node_idx_t n, bool donorFlag)//n,be donor or acceptor
  {
    double supCov = 0;
    if(donorFlag)
    {
        for(size_t i=0;i<node_set[n].children.size();i++)
	{
	    supCov += node_set[n].children[i].second;
	}
    }
    else 
    {
        for(size_t i=0;i<node_set[n].parents.size();i++)
	{
	    supCov += node_set[n].parents[i].second;
	}
    }
    return supCov;
  }
  void SimplifyGraph::compute_JScov_for_AllNodes()
  {

    for(size_t i = 0;i<Graph_size;i++)
    {
 	double c1=compute_JScov_for_Node(i,true);//donor
 	double c2=compute_JScov_for_Node(i,false);//acceptor
	node_set[i].JS_cov = make_pair(c2,c1);
 		
    }
  }
  void SimplifyGraph::remove_edges_basedon_JS_support()//suppose junction jl->jr, then junction site is jl or jr.
  {
    for(int i=0;i<Graph_size;i++)
    {
		
    }
  }
  void SimplifyGraph::get_LRP()
  {
    map<int,bool> head, tail;
    get_end_node(head, tail);

    if(1)
    {
	while(1 && !LongReadPath.empty())
        {
            vector<double>::iterator biggest = max_element(LongReadPath_cov.begin(), LongReadPath_cov.end());
            int k = biggest - LongReadPath_cov.begin();
            path_t pp = LongReadPath[k];

            final_paths.push_back(pp);

            if(LongReadPath_cov[k] < 2) break;
            else LongReadPath_cov[k] = 0;

        }
	return;
	//remove redundancy;
        for(size_t i=0;i<LongReadPath.size();i++)
	{
	    if(LongReadPath_cov[i] <= 1) continue;
	    for(size_t j=0;j<LongReadPath.size();j++)
	    {
	        if( i == j) continue;
		if(LongReadPath[i].size() > LongReadPath[j].size())
		{
		    if(check_overlap_of_2paths(LongReadPath[i],LongReadPath[j])){
			//if(LongReadPath_cov[i] > LongReadPath_cov[j])
			    LongReadPath_cov[j] = 0;
		    }
		}
	    }
	}
	//path end nodes == graph end nodes
	/*
	for(size_t i=0;i<LongReadPath.size();i++)
	{
	    break;
	    if(LongReadPath_cov[i] <= 1) continue;
	    path_t p = LongReadPath[i];
	    //if(!(node_set[p.front()].parents.empty() && node_set[p.back()].children.empty())) //good
	    //if((node_set[p.front()].parents.empty() && node_set[p.back()].children.empty()))
	    //if(node_set[p.front()].parents.empty())
	    //if(node_set[p.back()].children.empty())
	    if(node_set[p.front()].parents.empty() || node_set[p.back()].children.empty())
		    LongReadPath_cov[i] = 0;
	}
	*/
	
	for(size_t i=0;i<LongReadPath.size();i++)
	{
	    if(LongReadPath_cov[i] <= 1) continue;
	    path_t p = LongReadPath[i];

	    //LP.end-node == graph.end-node(head==head and tail==tail) //good
	    //if((node_set[p.front()].parents.empty() && node_set[p.back()].children.empty()))
	    if(LongReadPath_cov[i] >= 10)
	    {
	        final_paths.push_back(LongReadPath[i]);
		LongReadPath_cov[i] = 0;
	    }
	    /*
	    if(head.find(p.front()) != head.end() && tail.find(p.back()) != tail.end())
	    {    
		    final_paths.push_back(LongReadPath[i]);
		    LongReadPath_cov[i] = 0;
	    }
	    else if(LongReadPath_cov[i] >= 2)
	    {
	            final_paths.push_back(LongReadPath[i]);
		    LongReadPath_cov[i] = 0;
	    }
	    */

	    //!(head==head and tail==tail)
	    /*
	    if(!(node_set[p.front()].parents.empty() && node_set[p.back()].children.empty()))
		    final_paths.push_back(LongReadPath[i]);
		    */

	    //LP.tail==graph.tail
	    //if(node_set[p.front()].parents.empty()) final_paths.push_back(LongReadPath[i]);

	    //LP.head==graph.head
	    //if(node_set[p.back()].children.empty()) final_paths.push_back(LongReadPath[i]);

	    //head != head && tail != tail
	    /*
	    if(!node_set[p.front()].parents.empty() && !node_set[p.back()].children.empty()) 
		    final_paths.push_back(LongReadPath[i]);
	    */
	    //head==head or tail==tail
	    /*
	    if(node_set[p.front()].parents.empty() || node_set[p.back()].children.empty())
	        final_paths.push_back(LongReadPath[i]);
	    */

	}

	//sort(final_paths.begin(),final_paths.end());
	//final_paths.erase(unique(final_paths.begin(),final_paths.end()),final_paths.end());
    }
    return;

  }
  void SimplifyGraph::remove_edges_by_average_coverage(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent)
  {
	double all_coverage = 0;
	int edges_number = 0;
	map<node_idx_t, bool> used_nodes;
	vector<vector<int> > Comps;

	for(int i=0;i<Graph_size;i++)
	{
	    if(used_nodes[i]) continue;
	    used_nodes[i] = true;
	    vector<node_idx_t> comp_nodes;
	    comp_nodes.push_back(i);
	    for(size_t j=0;j<comp_nodes.size();j++)
	    {
		for(size_t k=0;k<node_set[comp_nodes[j]].children.size();k++)
		{
		    node_idx_t c = node_set[comp_nodes[j]].children[k].first;
		    if(find(comp_nodes.begin(),comp_nodes.end(), c) == comp_nodes.end())
		    {
			//cerr<<"  "<<c<<endl;
			comp_nodes.push_back(c);
			used_nodes[c] = true;
		    }
		}
		for(size_t k=0;k<node_set[comp_nodes[j]].parents.size();k++)
		{
		    node_idx_t p = node_set[comp_nodes[j]].parents[k].first;
		    if(find(comp_nodes.begin(),comp_nodes.end(), p) == comp_nodes.end())
		    {
		 	comp_nodes.push_back(p);
			used_nodes[p] = true;
		    }
		}
	    }
	    Comps.push_back(comp_nodes);
 	}
	map<node_idx_t,double> nodes_average_coverage_map;
	for(size_t i=0;i<Comps.size();i++)
	{
	    vector<int> comp = Comps[i];
	    double all_cov = 0;
	    double edges_number = 0;
	    for(size_t j=0;j<comp.size();j++)
	    {
		node_idx_t n = comp[j];
		for(size_t k=0;k<node_set[n].children.size();k++)
		{
		    if(partial(n,node_set[n].children[k].first)) continue;
		    all_cov += node_set[n].children[k].second;
		    edges_number++;
		}
	    }
	    double ave = 1;
	    if(edges_number != 0) ave = all_cov/edges_number;
	    if(ave < 2.0)//4.0/SampleSize) 
	    {
		for(size_t j=0;j<comp.size();j++) node_set[comp[j]].clear();
	    }
		
	    for(size_t j=0;j<comp.size();j++)
		nodes_average_coverage_map[comp[j]] = ave;
	}

	
	for(int i=0;i<Graph_size;i++)
	{
	    for(size_t j=0;j<node_set[i].children.size();j++)
	    {
		double rate = node_set[i].children[j].second / nodes_average_coverage_map[i];

		node_idx_t c = node_set[i].children[j].first;
		if(rate < AVERAGE_REMOVE_RATE)
	 	{
		    delete_edges_child.push_back(make_pair(i,node_set[i].children[j].first));
		}
	    }
 	}
	return;
  }
  bool SimplifyGraph::MultiJunction(node_idx_t i, node_idx_t c)
  {
	  return false;
      if(node_set[i].sequence.empty() || node_set[c].sequence.empty()) return false;
      pair<int,int> junc = make_pair(node_set[i].sequence.back(), node_set[c].sequence.front());
      vector<pair<int,int> >::iterator it = find(junctions.begin(),junctions.end(),junc);
      if(it == junctions.end()) return false;
      int k = it - junctions.begin();
      if(junctions_MappingInfo[k].size() <= 1) return false;
      if(junctions_MappingInfo[k].size() <= 0.05*SampleSize) return false;//1->0.05*SampleSize ipsc

      return true;
  }
  void SimplifyGraph::remove_partial_junction(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent)
  {
	double R=0.2; //raw 0.5
    	for(int i=1;i<Graph_size - 1;i++)
        {
	    if(partial(i,i+1))
	    {
	  	double partial_cov = node_set[i].get_child_coverage(i+1);
		if(partial_cov != -1)
		{
		    for(size_t j=0;j<node_set[i].children.size();j++)
	   	    {
			node_idx_t c = node_set[i].children[j].first;
			if(MultiJunction(i,c)) continue;
		 	if(node_set[i].children[j].second/partial_cov < R) 
			    delete_edges_child.push_back(make_pair(i,c));
		    }
		}
	    }
	    if(partial(i-1,i))
	    {
		double partial_cov = node_set[i-1].get_child_coverage(i);
		if(partial_cov != -1)
		{
		    for(size_t j=0;j<node_set[i].parents.size();j++)
		    {
			node_idx_t p = node_set[i].parents[j].first;
			if(MultiJunction(p,i)) continue;
			if(node_set[i].parents[j].second / partial_cov < R)
			    delete_edges_parent.push_back(make_pair(p,i));
		    }
		}
	    }
	}
	return;
  }
  void SimplifyGraph::remove_unbalance_edges( vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent)
  {
	for(int i=0;i<Graph_size;i++)
	{
	    if( !node_set[i].parents.empty() && node_set[i].children.size() == 1)
	    {
		double Incov = InCov(i);
		node_idx_t c = node_set[i].children.front().first;
		double rate = node_set[i].children.front().second/Incov;
		if(node_set[i].children.front().second > 5 && keep_edge(make_pair(i,c)) ) continue;
		if(rate<UNBALANCE_RATE) delete_edges_child.push_back(make_pair(i,c));
	    }
	    if(node_set[i].parents.size() == 1 && !node_set[i].children.empty())
	    {
		double outcov = OutCov(i);
		node_idx_t p = node_set[i].parents.front().first;
		double rate = node_set[i].parents.front().second/outcov;
		if(node_set[i].parents.front().second > 5 && keep_edge(make_pair(p,i)) ) continue;
		if(rate<UNBALANCE_RATE) delete_edges_parent.push_back(make_pair(p,i));
	    }
	}
	return;
  }
  void SimplifyGraph::remove_small_exons(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent)
  {
	for(int i=0;i<Graph_size;i++)
	{
	    //one out
	    if(node_set[i].parents.empty() && node_set[i].children.size() == 1)
	    {
		node_idx_t c =  node_set[i].children.front().first;
		double ecov = node_set[i].get_child_coverage(c);
		double r1 = double((1.0*ecov)/node_set[c].coverage());
		double r2 = double((1.0*node_set[i].coverage())/node_set[c].coverage());
		if(MultiJunction(i,c)) continue;
		//if(node_set[i].coverage() <= 2  && ecov < 3 && r1<0.2 &&  node_set[c].parents.size() > 1)
		if(node_set[i].coverage() < 5.0/SampleSize  && ecov < 5.0/SampleSize && r1<0.2 && r2<0.2 && node_set[c].parents.size() > 1 )
		//if(node_set[i].coverage() < 1 && ecov < 1 && r1<0.8 && r2<0.8 && node_set[c].parents.size() > 1 )
		{
		    delete_edges_child.push_back(make_pair(i,c));
		}
	    }
	    //one in
	    if(node_set[i].children.empty() && node_set[i].parents.size() == 1)
	    {
		node_idx_t p = node_set[i].parents.front().first;
		double ecov = node_set[p].get_child_coverage(i);

		double r1 = double(1.0*ecov/node_set[p].coverage());
		double r2 = double((1.0*node_set[i].coverage())/node_set[p].coverage());
		if(MultiJunction(p,i)) continue;
		//if(node_set[i].coverage() <=2 && ecov < 3 && r1<0.15 && node_set[p].children.size() > 1)
		if(node_set[i].coverage() < 5.0/SampleSize && ecov < 5.0/SampleSize && r1<0.2 && r2<0.2 && node_set[p].children.size() > 1)
		//if(node_set[i].coverage() < 1&& ecov < 1 && r1<0.8 && r2<0.8 && node_set[p].children.size() > 1)
		{
		    delete_edges_parent.push_back(make_pair(p,i));
		}
	    }
	}
        return;
  }
  bool SimplifyGraph::partial(node_idx_t n1, node_idx_t n2)
  {
	if(node_set[n1].sequence.empty() || node_set[n2].sequence.empty() ) return false;

	if(node_set[n1].sequence.back() + 1 == node_set[n2].sequence.front()) return true;
	return false;
  }
  void SimplifyGraph::remove_intron_contamination(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent)
  {
	for(int i=0;i<Graph_size;i++)
	{
	    if(node_set[i].children.size() == 1 && node_set[i].parents.size() == 1)
	    {
		bool partial_flag = false;
		node_idx_t p = node_set[i].parents.front().first;
		node_idx_t c = node_set[i].children.front().first;
		double pc_cov = node_set[p].get_child_coverage(c);
		if(partial(p,i) && partial(i,c) && pc_cov != -1)
		{
		    //if( (node_set[i].coverage() < 10.0/SampleSize && node_set[i].coverage()<pc_cov) )
		    if( node_set[i].coverage()<0.1*pc_cov) 		      
		    {
			delete_edges_child.push_back(make_pair(i,c));
			delete_edges_parent.push_back(make_pair(p,i));
		    } 
		}

	    }
	}
	return;
  }

  void SimplifyGraph::remove_partial_end_by_edge_coverage(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent)
  {
	double R=100;//1
	size_t L = 20000;
	for(int i=0;i<Graph_size;i++)
	{
	    if(node_set[i].children.size() > 1 || node_set[i].parents.size() > 1)
	    {
		bool partial_flag = false;
		node_idx_t partial_node;
		double max_cov = 0,partial_cov = 0;

		for(size_t j = 0;j<node_set[i].children.size();j++)
		{
		    if(partial(i,node_set[i].children[j].first)){
			partial_node = node_set[i].children[j].first;
			partial_cov = node_set[i].children[j].second;
			partial_flag = true;
		    }
		    if(node_set[i].children[j].second > max_cov) max_cov = node_set[i].children[j].second;
		}
	    	if(partial_flag)
		{
		    double r = double(1.0*partial_cov)/max_cov; //0.2
		    //if(node_set[partial_node].children.empty() && r <= 0.2 && node_set[partial_node].length() < 10000 && node_set[partial_node].parents.size() == 1)
		    //delete_edges_child.push_back(make_pair(i,partial_node));
		    size_t len = node_set[partial_node].length();
                    if( (node_set[partial_node].children.empty() && node_set[partial_node].parents.size() == 1)
                           && r<=R && len < L) 
                           //10.11-noref
                            delete_edges_child.push_back(make_pair(i,partial_node));
		}

		partial_flag = false; max_cov = 0;
		for(size_t j=0;j<node_set[i].parents.size();j++)
		{
		    if(partial(node_set[i].parents[j].first,i)){
			partial_node = node_set[i].parents[j].first;
			partial_cov = node_set[i].parents[j].second;
			partial_flag = true;
		    }
		    if(node_set[i].parents[j].second > max_cov) max_cov = node_set[i].parents[j].second;
		}
	  	if(partial_flag)
		{
		    double r = double(1.0*partial_cov)/max_cov; //0.2
		    //if(node_set[partial_node].parents.empty() && r <= 0.2&& node_set[partial_node].length() < 10000 && node_set[partial_node].children.size() == 1)
			//delete_edges_parent.push_back(make_pair(partial_node,i));

		    size_t len = node_set[partial_node].length();
                    if( (node_set[partial_node].parents.empty() && node_set[partial_node].children.size() == 1)
                           && r<=R && len < L) //1,20
                         //10.11-noref
                        delete_edges_parent.push_back(make_pair(partial_node,i));
		}
	    } 
	}
	return;
  }
  void SimplifyGraph::remove_edges_onemulti(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent)
  {
    for(int i=0;i<Graph_size;i++)
    {
	if(node_set[i].children.size() > 1 && node_set[i].parents.size() <=1)//one to multi
	{
	    double max_cov = 0;
	    for(size_t j = 0;j<node_set[i].children.size();j++)
		if(node_set[i].children[j].second > max_cov) max_cov = node_set[i].children[j].second;
	    for(size_t j = 0;j<node_set[i].children.size();j++)
	    {
		node_idx_t c = node_set[i].children[j].first;
		double cov = node_set[i].children[j].second;
		if(MultiJunction(i,c)) continue;
//		if(cov <= 1.0/SampleSize  ||(cov <= 10.0/SampleSize && double(cov/max_cov) <= 0.5))
		if(double(cov/max_cov) <= 0.01)
	 	{
		  if( node_set[c].children.empty())// || node_set[c].parents.size() > 1)
		  {
		    node_set[c].length_ = -1;
		    delete_edges_child.push_back(make_pair(i,c));
		  }
		}
	    }
	}
	if(node_set[i].children.size() <= 1 && node_set[i].parents.size() > 1)//multi to one
	{
	    double max_cov = 0;
	    for(size_t j=0;j<node_set[i].parents.size();j++)
	    	if(node_set[i].parents[j].second > max_cov) max_cov = node_set[i].parents[j].second;
	    for(size_t j=0;j<node_set[i].parents.size();j++)
	    {
	    	node_idx_t p = node_set[i].parents[j].first;
	    	double cov = node_set[i].parents[j].second;
		if(MultiJunction(p,i)) continue;
	    	//if(cov <= 1.0/SampleSize  ||(cov <= 10.0/SampleSize && double(cov/max_cov) <= 0.5))
		if(double(cov/max_cov) <= 0.01)
	    	{
		  if(node_set[p].parents.empty())// || node_set[p].children.size() > 1)
		  {
		    node_set[p].length_ = -1;
		    delete_edges_parent.push_back(make_pair(p,i));
		  }
	    	}
	    }
	}
    }
    return;
  }
  void SimplifyGraph::remove_lowcov_edges_of_bifurcation_nodes2(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent)
  {
    for(int i=0;i<Graph_size;i++)
    {
	if(node_set[i].children.size() > 1 )
	{
	    double max_cov = 0;
	    for(size_t j = 0;j<node_set[i].children.size();j++) 
	    {
		if( !partial(i,node_set[i].children[j].first)) continue;

		if(node_set[i].children[j].second > max_cov) max_cov = node_set[i].children[j].second;
	    }

	    for(size_t j = 0;j<node_set[i].children.size();j++)
	    {
		node_idx_t c = node_set[i].children[j].first;
	   	double cov = node_set[i].children[j].second;
		double r = double(cov/max_cov);
		if(MultiJunction(i,c)) continue;
		if( !partial(i,c) && cov < max_cov)
		{
			node_set[c].length_ = -1;
			delete_edges_child.push_back(make_pair(i,c));
		}
	    }
	    
	    max_cov = 0;
	    for(size_t j=0;j<node_set[i].parents.size();j++)
	    {
		if( !partial(node_set[i].parents[j].first,i)) continue;
		if(node_set[i].parents[j].second > max_cov) max_cov = node_set[i].parents[j].second;
	    }

	    for(size_t j=0;j<node_set[i].parents.size();j++)
	    {
		node_idx_t p = node_set[i].parents[j].first;
		double cov = node_set[i].parents[j].second;
		double r = double(cov/max_cov);
		if(MultiJunction(p,i)) continue;
		if(!partial(p,i) &&  cov < max_cov)//5 0.18
		{
			node_set[p].length_ = -1;
			delete_edges_parent.push_back(make_pair(p,i));
		}
	    }
	    //cerr<<endl;
	}	
    }
    return;
  }
  void SimplifyGraph::remove_lowcov_edges_of_bifurcation_nodes(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent)
  {
    double COV1 = 3.0;
    double R = 0.2;
    double COV2 = 5.0;
    //if(PackingFlag) R = 0.1;
    for(int i=0;i<Graph_size;i++)
    {
	//cerr<<"node: "<<i<<endl;
	if(node_set[i].children.size() > 1 && node_set[i].parents.size() > 1)
	{
	    //cerr<<"child"<<endl;
	    double max_cov = 0;
	    for(size_t j = 0;j<node_set[i].children.size();j++) 
	    {
		if(node_set[i].children[j].second > max_cov) max_cov = node_set[i].children[j].second;
	    }

	    for(size_t j = 0;j<node_set[i].children.size();j++)
	    {
		node_idx_t c = node_set[i].children[j].first;
	   	double cov = node_set[i].children[j].second;
		//if(MultiJunction(i,c)) continue;
		//if(keep_edge(make_pair(i,c))) continue;
		double r = double(cov/max_cov);
		if(cov < COV1 && (cov <= COV2 && r <= R)) //5,0.18
		{
		    if(1 || node_set[c].children.empty() || node_set[c].parents.size() > 1)
		    {
			node_set[c].length_ = -1;
			delete_edges_child.push_back(make_pair(i,c));
		    }
		}
	    }
	    
	    max_cov = 0;
	    for(size_t j=0;j<node_set[i].parents.size();j++)
	    {
		if(node_set[i].parents[j].second > max_cov) max_cov = node_set[i].parents[j].second;
	    }

	    for(size_t j=0;j<node_set[i].parents.size();j++)
	    {
		node_idx_t p = node_set[i].parents[j].first;
		double cov = node_set[i].parents[j].second;
		//if(MultiJunction(p,i)) continue;
		//if(keep_edge(make_pair(p,i))) continue;
		double r = double(cov/max_cov);
		if(cov < COV1 && (cov <= COV2 && r <= R))//5 0.18
		{
		    if(1 || node_set[p].parents.empty() || node_set[p].children.size() > 1)
		    {
			node_set[p].length_ = -1;
			delete_edges_parent.push_back(make_pair(p,i));
		    }
		}
	    }
	    //cerr<<endl;
	}	
    }
    return;
  }
  void SimplifyGraph::delete_children_edges( vector<edge_t> delete_edges_child)
  {
    sort(delete_edges_child.begin(),delete_edges_child.end());
    delete_edges_child.erase(unique(delete_edges_child.begin(),delete_edges_child.end()),delete_edges_child.end());
    for(size_t i = 0;i<delete_edges_child.size();i++)
    {
	node_idx_t n1 = delete_edges_child[i].first, n2 = delete_edges_child[i].second;
    	if( !node_set[n1].sequence.empty() && !node_set[n2].sequence.empty() )
        {
            double cov = node_set[n1].get_child_coverage(n2);

            pair<int,int> junction = make_pair(node_set[n1].sequence.back(),node_set[n2].sequence.front());
	}
	node_set[n1].delete_child(n2);
	node_set[n2].delete_parent(n1);
    }
    /*
    for(size_t i = 0;i<delete_edges_child.size();i++)
    {
	node_idx_t c = delete_edges_child[i].second;
	vector< pair<node_idx_t ,double> > c_children = node_set[c].children;
	while(1)
	{
	    if(c_children.empty() ) break;
	    if(c_children.size() >= 2) break;

	    node_idx_t c_child = c_children.front().first;
	    if(c_children.front().second > 1.5) break;
	    node_set[c].delete_child(c_child);
	    node_set[c_child].delete_parent(c);

	    if(!node_set[c_child].parents.empty()) break;//other in_edges

	    c_children = node_set[c_child].children;
	}
    }
    */
    
    return;
  }    
  void SimplifyGraph::delete_parents_edges( vector<edge_t> delete_edges_parent)
  {
    sort(delete_edges_parent.begin(),delete_edges_parent.end());
    delete_edges_parent.erase(unique(delete_edges_parent.begin(),delete_edges_parent.end()), delete_edges_parent.end());
    for(size_t i = 0;i<delete_edges_parent.size();i++)
    {
	node_idx_t n1 = delete_edges_parent[i].first, n2 = delete_edges_parent[i].second;
	if( !node_set[n1].sequence.empty() && !node_set[n2].sequence.empty() )
        {
             double cov = node_set[n1].get_child_coverage(n2);
 
             pair<int,int> junction = make_pair(node_set[n1].sequence.back(),node_set[n2].sequence.front());
        }
	node_set[n1].delete_child(n2);
	node_set[n2].delete_parent(n1);
    }
    /*
    for(size_t i = 0;i<delete_edges_parent.size();i++)
    {
	node_idx_t p = delete_edges_parent[i].first;
	vector< pair<node_idx_t ,double> > p_parents = node_set[p].parents;
	while(1)
	{
	    if(p_parents.empty()) break;
	    if(p_parents.size() >= 2) break;
	    node_idx_t p_parent = p_parents.front().first;
	    if(p_parents.front().second > 1.5) break;
	    node_set[p_parent].delete_child(p);
	    node_set[p].delete_parent(p_parent);

	    if(!node_set[p_parent].children.empty()) break; //other out_edges;
	    p_parents = node_set[p_parent].parents;
	}
    }
    */
    return;
  }
  void SimplifyGraph::show_graph()
  {
    out_graph<<endl<<"#Graph "<<rg_index<<" >>>"<<endl;
    out_graph<<"Edges: "<<endl;

    for(int i=0;i<Graph_size;i++)
    {
        Node n = node_set[i];
        for(size_t j=0;j<n.children.size();j++)
        {
            out_graph<<i<<"->"<<n.children[j].first<<": "<<n.children[j].second<<endl;
        }
    }
    out_graph<<"Nodes: "<<endl;
    for(int i=0;i<Graph_size;i++)
    {
	out_graph<<i<<": "<<chr<<" "<<strand<<" "
		 <<node_set[i].str()<<": "<<node_set[i].coverage()
		 <<" ("<<node_set[i].JS_cov.first<<","<<node_set[i].JS_cov.second<<")"<<endl;
    }
    out_graph<<"LRP: "<<endl;
    for(size_t i=0;i<LongReadPath.size();i++) {
         vector<int> pp = LongReadPath[i];
         for(size_t j=0;j<pp.size();j++) out_graph<<pp[j]<<" ";
         out_graph<<": "<<LongReadPath_cov[i]<<endl;
    }
    //out_graph<<"#Graph "<<rg_index<<endl;
    return;
  }
  void SimplifyGraph::show_junction()
  {
    for(size_t i=0;i<Graph_size;i++)
    {
        for(size_t j=0;j<node_set[i].children.size();j++)
	{
	    node_idx_t c = node_set[i].children[j].first;
	    if(partial(i,c)) continue;
	    cout<<strand<<" "<<chr<<" "<<node_set[i].sequence.back()<<" "<<node_set[c].sequence.front()<<endl;
	}
    }
  }
  void SimplifyGraph::contract_graph()
  {
    map<node_idx_t,bool> used_nodes;
    int S = Graph_size;
    for(int i=0;i<S;i++)
    {
	if(used_nodes[i]) continue;
	used_nodes[i] = true;
	path_t P(1,i);
	while(1)
	{
	    node_idx_t n = P.back();
	    bool eflag=false;
	    if(node_set[n].children.size() == 1) 
	    {
		node_idx_t c = node_set[n].children.front().first;
		if(node_set[c].parents.size() == 1){
		   P.push_back(c);
		   used_nodes[c] = true;
		   eflag = true;
		}
	    }
	    if(!eflag) break;
	}
	while(1)
	{
	    node_idx_t n = P.front();
	    bool eflag=false;
	    if(node_set[n].parents.size() == 1)
	    {
		node_idx_t p = node_set[n].parents.front().first;
		if(node_set[p].children.size() == 1){
		    P.insert(P.begin(),p);
		    used_nodes[p] = true;
		    eflag = true;
		}
	    }
	    if(!eflag) break;
	}
	if(P.size() == 1){
	    rawNode_eNode_map[i] = i;
	    continue;
	}
	Node N;
	N.coverage_ = 1;
	N.length_ = 1;
	N.path = P;
	node_set.push_back(N);
	Graph_size++;
	node_idx_t Nid = Graph_size - 1;
	node_set[Nid].children = node_set[P.back()].children;
	node_set[Nid].parents = node_set[P.front()].parents;
	for(size_t i=0;i<node_set[P.back()].children.size();i++){
	    node_idx_t c = node_set[P.back()].children[i].first;
	    double cov = node_set[P.back()].children[i].second;
	    node_set[ c ].delete_parent(P.back());
	    node_set[ c ].add_parent(Nid,cov);
	}
	for(size_t i=0;i<node_set[P.front()].parents.size();i++){
	    node_idx_t p = node_set[P.front()].parents[i].first;
	    double cov = node_set[P.front()].parents[i].second;
	    node_set[ p ].delete_child(P.front());
	    node_set[ p ].add_child(Nid,cov);
	}
	for(size_t i=0;i<P.size();i++){ 
	    rawNode_eNode_map[P[i]] = Nid;
	    node_set[P[i]].clear();
	}
    }
    return;
  }

  void SimplifyGraph::process_LRP_after_simplify_graph()
  {
	double R = 1;
	for(size_t i=0;i<LongReadPath.size();)
	{
	    path_t p = LongReadPath[i];
	    bool flag = false;
	    if(p.size() < 3) flag = true;//NEW
	    int normal_number = 0, partial_number = 0;
	    double cov = 10000000000;
	    for(size_t j=0;j<p.size() - 1;j++)
	    {
		double c_ = node_set[p[j]].get_child_coverage(p[j+1]);
	        if(c_ == -1)
		{
	 	    flag = true;
		    break;
		}
		else if(c_ < cov) cov = c_;
		
		if(partial(p[j],p[j+1])){
		    partial_number++;
		}
		else normal_number++;
		
	    }


	    if(normal_number == 0) flag = true;
	    if(LongReadPath_cov[i]<R) {//3.1 spk
		flag = true;
	    }
	    //flag = false;
	    if(flag)
	    {
	    	LongReadPath.erase(LongReadPath.begin() + i);
	    	LongReadPath_cov.erase(LongReadPath_cov.begin() + i);
	    }
	    else i++;
	}
	//out_graph<<"Here2"<<endl;
	return;
  }
  void SimplifyGraph::get_graph(vector<vector<int> >& vecEdges,vector<double>& vecCov)
  {
	//out_graph<<"Here:"<<endl;
	for(int i=0;i<Graph_size;i++)
	{
	    for(size_t j=0;j<node_set[i].children.size();j++)
	    {
		vector<int> edge;
		double cov;
		edge.push_back(i);edge.push_back(node_set[i].children[j].first);
	   	cov = node_set[i].children[j].second;

		vecEdges.push_back(edge); 
		if(cov == 0) cov = 0.1;
		vecCov.push_back(cov);
	    }
	}
	//out_graph<<"Here1"<<endl;
	double R = 1;//3.1
	for(size_t i=0;i<LongReadPath.size();)
	{
	    path_t p = LongReadPath[i];
	    bool flag = false;
	    if(p.size() < 3) flag = true;//NEW
	    int normal_number = 0, partial_number = 0;
	    double cov = 10000000000;
	    for(size_t j=0;j<p.size() - 1;j++)
	    {
		double c_ = node_set[p[j]].get_child_coverage(p[j+1]);
	        if(c_ == -1)
		{
	 	    flag = true;
		    break;
		}
		else if(c_ < cov) cov = c_;
		
		if(partial(p[j],p[j+1])){
		    partial_number++;
		}
		else normal_number++;
		
	    }


	    if(normal_number == 0) flag = true;
	    if(LongReadPath_cov[i]<R) {//3.1 spk
		flag = true;
	    }
	    //flag = false;
	    if(flag)
	    {
	    	LongReadPath.erase(LongReadPath.begin() + i);
	    	LongReadPath_cov.erase(LongReadPath_cov.begin() + i);
	    }
	    else i++;
	}
	//out_graph<<"Here2"<<endl;
	return;
  }
  double SimplifyGraph::InCov(node_idx_t node)
  {
    double cov = 0;
    for(size_t i=0;i<node_set[node].parents.size();i++)
    {
	cov += node_set[node].parents[i].second;
    }
    return cov;
  }

  double SimplifyGraph::OutCov(node_idx_t node)
  {
    double cov = 0;
    for(size_t i=0;i<node_set[node].children.size();i++)
    {
	cov += node_set[node].children[i].second;
    }
    return cov;
  }
  double SimplifyGraph::average_coverage()
  {
	double cov = 0;
	int number = 0;
	for(int i=0;i<Graph_size;i++)
	{
	    for(size_t j=0;j<node_set[i].children.size();j++)
	    {
		cov += node_set[i].children[j].second;
		number++;
	    }
	}
	if(number == 0) return 1;
	return (cov/number);
  }
  void SimplifyGraph::get_packing_result_new(vector<vector<int> >Vec_edges, vector<int> Edges_legt, vector<int>Edges_right, vector<double> Weights)
  {
//      packing_map
	packing_map.clear();      
	for(size_t i=0;i<Weights.size();i++)
	{
	    if(Weights[i] == 1) continue;
	    int i1 = Edges_legt[i], i2 = Edges_right[i];
	    vector<int> E1 = Vec_edges[i1], E2 = Vec_edges[i2];

	    pair<int,int> e1 = make_pair(E1[0],E1[1]), e2 = make_pair(E2[0],E2[1]);
	    packing_map[make_pair(e1,e2)] = 0;
	}
	return;
  }
  void SimplifyGraph::forward_extend(edge_t e, vector<path_t>& paths)
  {
      path_t path_temp;
      path_temp.push_back(e.first);
      path_temp.push_back(e.second);
      vector<path_t> paths_temp;
      paths_temp.push_back(path_temp);
      while(1)
      {
	  if(paths_temp.empty()) break;

	  path_t current_path = paths_temp.front();
	  paths_temp.erase(paths_temp.begin());

          edge_t current_edge = make_pair(current_path[current_path.size() - 2], current_path[current_path.size() - 1]);
	  if(node_set[current_edge.second].children.size() == 0)
	  {
	      paths.push_back(current_path);    
	  }
	  else
	  {
	      vector<pair<edge_t,double> > edges_cov;
	      for(size_t i=0;i<node_set[current_edge.second].children.size() ;i++)
	      {
		  edge_t e = make_pair(current_edge.second, node_set[current_edge.second].children[i].first);
		  edges_cov.push_back(make_pair(e,node_set[current_edge.second].children[i].second));
	      }
	      sort(edges_cov.begin(),edges_cov.end(),sorter);
	      if(edges_cov.size() == 1)
	      {
	          edge_t e = edges_cov[0].first;
		  current_path.push_back(e.second);
		  paths_temp.push_back(current_path);
	      }
	      else {
	        for(size_t i=0;i<edges_cov.size();i++)
	        {
	          edge_t e = edges_cov[i].first;
		  //if(!partial(e.first,e.second)) continue;
	          pair<edge_t ,edge_t> ee = make_pair(current_edge, e);
		  if(1 || packing_map.find(ee) != packing_map.end())
		  {
		      current_path.push_back(e.second);
		      paths_temp.push_back(current_path);
		      break;
		  }
		}
	      }
	  }
      }
      return;

  }
  void SimplifyGraph::reverse_extend(edge_t e, vector<path_t>& paths)
  {
      path_t path_temp;
      path_temp.push_back(e.first);
      path_temp.push_back(e.second);
      vector<path_t> paths_temp;
      paths_temp.push_back(path_temp);
      while(1)
      {
	  if(paths_temp.empty()) break;

	  path_t current_path = paths_temp.front();
	  paths_temp.erase(paths_temp.begin());

          edge_t current_edge = make_pair(current_path[0], current_path[1]);
	  if(node_set[current_edge.first].parents.size() == 0)
	  {
	      paths.push_back(current_path);    
	  }
	  else
	  {
	      vector<pair<edge_t,double> > edges_cov;
	      for(size_t i=0;i<node_set[current_edge.first].parents.size() ;i++)
	      {
		  edge_t e = make_pair(node_set[current_edge.first].parents[i].first,current_edge.first);
		  edges_cov.push_back(make_pair(e,node_set[current_edge.first].parents[i].second));
	      }
	      sort(edges_cov.begin(),edges_cov.end(),sorter);
	      if(edges_cov.size() == 1)
	      {
	          edge_t e = edges_cov[0].first;
		  current_path.insert(current_path.begin(), e.first);
		  paths_temp.push_back(current_path);
	      }
	      else 
	      {
	        for(size_t i=0;i<edges_cov.size();i++)
	        {
	          edge_t e = edges_cov[i].first;
		  //if(!partial(e.first,e.second)) continue;
	          pair<edge_t ,edge_t> ee = make_pair(e,current_edge);
		  if(1|| packing_map.find(ee) != packing_map.end())
		  {
		      //current_path.push_back(e.second);
		      current_path.insert(current_path.begin(), e.first);
		      paths_temp.push_back(current_path);
		      break;
		  }
		}
	      }
	  }
      }
      return;
  }

  bool SimplifyGraph::check_overlap_of_2paths(path_t p1_, path_t p2_)
  {
    if(p1_.empty() || p2_.empty()) return true;
    if(p1_ == p2_) return true;

    path_t p1,p2;
    if(p1_.size() <= p2_.size()) { p1 = p1_; p2= p2_;}
    else {p1 = p2_; p2 = p1_;}

    //p2.size() > p1.size()
    //out_graph<<" size: "<<p2.size()<<" "<<p1.size()<<endl;
    int I = -1;
    for(size_t i=0;i<p2.size();i++)
    {
        if(p2[i] == p1.front())
	{
	    I = i;
	    break;
	}
    }
    //out_graph<<" c: "<<I<<endl;
    if(I == -1) return false;

    for(size_t i=0;i<p1.size();i++)
    {
	//out_graph<<" "<<p1[i]<<" "<<p2[I+i]<<endl;
	if(I+i >= p2.size()) return false;
        if(p1[i] != p2[I+i]) return false;
    }
    return true;

    stringstream ss1,ss2;
    string s_p1, s_p2;
    for(size_t i=0;i<p1.size();i++) ss1<<p1[i]<<"_";
    s_p1 = ss1.str();
    for(size_t i=0;i<p2.size();i++) ss2<<p2[i]<<"_";
    s_p2 = ss2.str();

    if(s_p2.find(s_p1) != string::npos) return true;
    return false;   
  }
  void SimplifyGraph::get_components()
  {
        map<node_idx_t, bool> used_nodes;
	for(int i=0;i<Graph_size;i++)
	{
	    if(used_nodes[i]) continue;
	    used_nodes[i] = true;
	    vector<node_idx_t> comp_nodes;
	    comp_nodes.push_back(i);
	    for(size_t j=0;j<comp_nodes.size();j++)
	    {
		for(size_t k=0;k<node_set[comp_nodes[j]].children.size();k++)
		{
		    node_idx_t c = node_set[comp_nodes[j]].children[k].first;
		    if(find(comp_nodes.begin(),comp_nodes.end(), c) == comp_nodes.end())
		    {
			//cerr<<"  "<<c<<endl;
			comp_nodes.push_back(c);
			used_nodes[c] = true;
		    }
		}
		for(size_t k=0;k<node_set[comp_nodes[j]].parents.size();k++)
		{
		    node_idx_t p = node_set[comp_nodes[j]].parents[k].first;
		    if(find(comp_nodes.begin(),comp_nodes.end(), p) == comp_nodes.end())
		    {
		 	comp_nodes.push_back(p);
			used_nodes[p] = true;
		    }
		}
	    }
	    //Comps.push_back(comp_nodes);
	    nodes_of_components.push_back(comp_nodes);
 	}
  }
  void SimplifyGraph::update_graph(path_t path)
  {
      double mincov = 10000000;
      for(size_t i=0;i<path.size() -  1;i++)
      {
          int n1=path[i],n2=path[i+1];
	  double cov = node_set[n1].get_child_coverage(n2);
	  if(cov < mincov) mincov=cov;
      }
      for(size_t i=0;i<path.size() -  1;i++)
      {
          int n1=path[i],n2=path[i+1];
	  node_set[n1].reduced_child_coverage(n2,mincov/2);
	  node_set[n2].reduced_parent_coverage(n1,mincov/2);
      }
      return;
  }
  void SimplifyGraph::path_search(string strand, string chr)
  {
  //out_graph<<"raw pathsearch..."<<endl;
  //get_components(); //Liger
  //for(size_t I=0;I<nodes_of_components.size();I++)//Liger
  {
    //vector<int> nodes = nodes_of_components[I];//Liger
    //vector<int> nodes;
    //for(int i=0;i<Graph_size;i++) nodes.push_back(i);
    vector<pair<int,int> > unused_junctions;
    vector<double> unused_junctions_cov;
    vector<int> unused_junctions_sample_number;

    for(int i = 0;i< Graph_size;i++) //add source and sink
    {
	if(SSFlag)
	{
	    if(i == 0) continue;
	}
        for(size_t j=0;j<node_set[i].children.size();j++)
	{
	    if(partial(i,node_set[i].children[j].first)) continue; 
	    int c = node_set[i].children[j].first;

	    if(SSFlag){
	      if( c == Graph_size - 1) continue;//add source and sink
	    }

	    double cov = node_set[i].children[j].second;

	    pair<int,int> junc_ = make_pair(i,node_set[i].children[j].first);
	    unused_junctions.push_back(junc_);
	    unused_junctions_cov.push_back(node_set[i].children[j].second);
	}

    }
    /*
    if(unused_junctions.empty())
    {
        for(int k=0;k<nodes.size();k++)
	{
	    int i = nodes[k];
	    path_t p(1,i);
	    continue;
	    final_paths.push_back(p);
	}
    }
    */
    //get pair of current components
    vector<path_t> LongReadPath_current = LongReadPath;
    vector<double> LongReadPath_current_cov = LongReadPath_cov;

    //final_paths.clear();
    map<edge_t,bool> used_edges;
    for(size_t i=0;i<final_paths.size();i++)
    {
	path_t p = final_paths[i];
	for(size_t j=0;j<p.size()-1;j++)
	{
            edge_t edge = make_pair(p[j],p[j+1]);
            used_edges[edge] = true;
	}
    }


    while(1 && !LongReadPath_current.empty()) //pair as seed
    {
	vector<double>::iterator biggest = max_element(LongReadPath_current_cov.begin(), LongReadPath_current_cov.end());
	int k = biggest - LongReadPath_current_cov.begin();
	path_t pp = LongReadPath_current[k];

	if(LongReadPath_current_cov[k] < 2) break;
	else LongReadPath_current_cov[k] = 0;

        vector<path_t> forward_paths, reverse_paths;

	edge_t e = make_pair(pp[pp.size() - 2], pp[pp.size() - 1]);
	forward_extend(e,forward_paths);
	e = make_pair(pp[0],pp[1]);
	reverse_extend(e,reverse_paths);
	path_t final_path = reverse_paths.front();
	for(size_t i=2;i<pp.size();i++) final_path.push_back(pp[i]);
	for(size_t i=2;i<forward_paths.front().size();i++) final_path.push_back(forward_paths.front()[i]);
	for(size_t i=0;i<final_path.size()-1;i++){
	    edge_t edge = make_pair(final_path[i],final_path[i+1]);
	    used_edges[edge] = true;
	}

	final_paths.push_back(final_path);
	//final_paths.push_back(pp);
	for(size_t j=0;j<LongReadPath_current.size();){
	    if(check_overlap_of_2paths(LongReadPath_current[j],final_path))
	    {
	        LongReadPath_current.erase(LongReadPath_current.begin() + j);
		LongReadPath_current_cov.erase(LongReadPath_current_cov.begin() + j);
	    }
	    else j++;
	}
	//break;
	
    }
    double seed_filter = 3;
    int path_number = 0;
    while(1 && !unused_junctions.empty())//junction as seed
    {
        vector<double>::iterator biggest = max_element(unused_junctions_cov.begin(),unused_junctions_cov.end());
	int k = biggest - unused_junctions_cov.begin();
	double seed_cov = unused_junctions_cov[k];
	edge_t e = unused_junctions[k];
	//int sample_number = unused_junctions_sample_number[k];
	if(unused_junctions_cov[k] < seed_filter) break;
	else unused_junctions_cov[k] = 0;	
	
	if(used_edges.find(e) != used_edges.end()) continue;

	vector<path_t> forward_paths, reverse_paths;
	forward_extend(e,forward_paths);
	reverse_extend(e,reverse_paths);
	path_t final_path = reverse_paths.front();
	for(size_t i=2;i<forward_paths.front().size();i++) final_path.push_back(forward_paths.front()[i]);
	for(size_t i=0;i<final_path.size()-1;i++)
	{
	    edge_t edge = make_pair(final_path[i],final_path[i+1]);
	    used_edges[edge] = true;
	}
	final_paths.push_back(final_path);
	path_number++;
	//if(path_number>=6) break;
    }
    //output(strand,chr,final_paths);
    //get_LRP();
    //final_paths.clear();
    //final_paths_info.clear();
 }
    return;
}
  void SimplifyGraph::remove_partial_end_for_a_path(path_t& pp)
  {
        if(pp.size() >= 2 )
        {
                    for(size_t j=0;j<pp.size()-1;)
                    {
                        if(partial(pp[j],pp[j+1]))
                            pp.erase(pp.begin() + j);
                        else break;
                    }
                    if(pp.size() >= 2 )
                    {
                      for(size_t j = pp.size()-1;j>0;j--)
                      {
                        if(partial(pp[j-1],pp[j]))
                                pp.erase(pp.begin() + j);
                        else break;
                      }
                    }
          }
        return;
  }
  void SimplifyGraph::output(string strand, string chr, vector<path_t> paths)
  {
    //trans_id = 1;
    //if(paths.empty()) return;	
    /*
    if(!paths.empty())
    {
      vector<path_t> FP = paths;
      sort(paths.begin(),paths.end());
      paths.erase(unique(paths.begin(),paths.end()),paths.end());
    }
    */
    //out_graph<<"output..."<<endl;
    vector<int>paths_flag(paths.size(),1);
    //remove redundancey
    
    for(size_t i=0;i<paths.size();i++)
    {
            if(paths_flag[i] == 0) continue;
            for(size_t j=0;j<paths.size();j++)
            {
                if( i == j) continue;
                if(paths[i].size() >= paths[j].size())
                {
		    path_t p = paths[i];
		    //out_graph<<" here"<<endl;
		    //out_graph<<" ";
		    //for(size_t k=0;k<p.size();k++) out_graph<<p[k]<<"->";
		    //out_graph<<endl;

		    p = paths[j];
		    /*
		    out_graph<<" ";
		    for(size_t k=0;k<p.size();k++) out_graph<<p[k]<<"->";
		    out_graph<<endl;
		    */
                    if(check_overlap_of_2paths(paths[i],paths[j]))
		    {

			    /*
			    path_t p = paths[i];
			    out_graph<<"here"<<endl;
        		    for(size_t k=0;k<p.size();k++) out_graph<<p[k]<<"->";
        		    out_graph<<endl;

			    p = paths[j];
        		    for(size_t k=0;k<p.size();k++) out_graph<<p[k]<<"->";
        		    out_graph<<endl;
			    */
                            paths_flag[j] = 0;
			    //out_graph<<" YES"<<endl;
                    }
                }
            }
    }
    //out_graph<<"A"<<endl;
    out_info<<"Assembled_Paths:"<<endl;
    for(size_t i=0;i<right_paths.size();i++)
    {
	path_t p = right_paths[i];
	remove_partial_end_for_a_path(p);//only consider intron chain
	out_info<<right_paths_ids[i]<<" ";
        for(size_t j=0;j<p.size();j++) out_info<<p[j]<<"->";
	out_info<<endl;
    }
    return;
  }
