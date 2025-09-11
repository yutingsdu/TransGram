#ifndef SIMPLFYGRAPH_H
#define SIMPLFYGRAPH_H

#include <iostream> 
#include <map>
#include <string>
#include <vector>
#include <vector>
#include <algorithm>
#include <sstream>
using namespace std;
extern map<string, map<pair<int,int>,bool> > Chr_Junc_map;
extern map<string, vector<vector<int> > > Chr_Trans_map;
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
extern vector< vector<int> > AllTrans;
extern vector< pair<string,string> > AllTransChr;
extern bool GFlag_;

typedef int node_idx_t;
typedef vector<node_idx_t> path_t;
typedef pair<node_idx_t,node_idx_t> edge_t;
typedef pair<edge_t,double> edge_with_cov_t;


class SimplifyGraph{
public:
  
   class Node{
     public:
	double coverage_;
	int length_;
	vector<int> sequence;
	vector< pair<node_idx_t,double> > children;
	vector< pair<node_idx_t,double> > parents;
	pair<int,int> JS_cov; //acceptor and donor


	path_t path;//some node merge to one
     public:
  	void add_child(node_idx_t child,double weight){

	  pair<node_idx_t,double> p = make_pair(child,weight);
	  this->children.push_back(p);
	  return;
	}
	void add_parent(node_idx_t parent,double weight) 
	{
	  pair<node_idx_t,double> p = make_pair(parent,weight);
          this->parents.push_back(p);
          return;
	}
	void reduced_child_coverage(node_idx_t child,double weight)
	{
	    for(size_t i=0;i<children.size();i++)
	    {
	        if(child == children[i].first)
		{
			children[i].second -= weight;
			if(children[i].second < 0) children[i].second = 0;
			break;
		}
	    }
	    return;
	}
	void reduced_parent_coverage(node_idx_t parent,double weight)
	{
	    for(size_t i=0;i<parents.size();i++)
	    {
		if(parent == parents[i].first)
		{
		    parents[i].second -= weight;
		    if(parents[i].second < 0) parents[i].second = 0;
		    break;
		}	    
	    }
	}
	void delete_parent(node_idx_t n)
	{
	    for(size_t i=0;i<parents.size();)
            {
                if(n == parents[i].first) parents.erase(parents.begin() + i);
                else i++;
            }
            return;
	}
	void delete_child(node_idx_t n)
	{
	    for(size_t i=0;i<children.size();)
	    {
		if(n == children[i].first) children.erase(children.begin() + i);
		else i++;
	    }
	    return;
	}
	void add_path_node(node_idx_t node) {
	  path.push_back(node);
	  return;
	}
	double get_child_coverage(node_idx_t child)
     	{
	    for(size_t i=0;i<children.size();i++)
	    {
		if(child == children[i].first) return children[i].second;
	    }
	    return -1;
    	}
	string path_str()
	{
	    if(path.empty()) return "*; ";
	    stringstream ss;
	    for(size_t i=0;i<path.size();i++) ss<<path[i]<<"-";
	    return ss.str();
	}
	string str()
	{
	    if(sequence.empty()) return "*; ";
	    stringstream ss;
	  
	    for(size_t i=0;i<sequence.size()-1;){
		ss<<sequence[i]<<" "<<sequence[i+1]<<"; ";
		i = i+2;
	    }
	    return ss.str();
	}
	void addseq(vector<int> vec)
	{
	    for(size_t i=0;i<vec.size();i++) sequence.push_back(vec[i]);
	    return;
	}
	double coverage()
	{
	    return coverage_;
	}
	int length()
	{
	    return length_;
	}
	void clear()
	{
	    children.clear();
	    parents.clear();
	    path.clear();
	    sequence.clear();
	    coverage_ = 0;
	    length_ = -1;
	}
	
   };

   class block_info
   {
    public:
       node_idx_t start;
       node_idx_t end;
       vector<path_t> paths;
       vector<double> paths_cov;
       vector<int> paths_flag;//path used or not
       block_info(node_idx_t s, node_idx_t e, vector<path_t> paths_,vector<double> paths_cov_)
       {
           start = s;
	   end = e;
	   paths = paths_;
	   paths_cov = paths_cov_;
	   for(size_t i=0;i<paths.size();i++) paths_flag.push_back(0);
       }
       int get_index_of_maxcov_path(bool flag)//flag = true -> used paths can be used again, 
       {
	   int index= -1;
	   double cov = 0;
           for(size_t i=0;i<paths_flag.size();i++)
	   {
	       	if(!flag && paths_flag[i] == 1) continue;//has been used
		//vector<double>::iterator biggest = max_element(paths_cov.begin(),paths_cov.end());
		//index = biggest - paths_cov.begin();
		if(paths_cov[i] > cov)
		{
		    index = i;
		    cov = paths_cov[i];
		}
	   }

	   return index;
       }
   };

public:
   struct path_info
   {
       double seed_cov;
       int normal_edge_number;
       int partial_edge_number;
       double max_cov;
       double min_cov;
       double ave_cov;
       int mj_number;
       int Nmj_number;
       int seed_sample_number;
   };

public:
   vector<Node> node_set;
   bool SSFlag;//add source&sink or not;
   bool GFlag; //most edges covered by LRP or not -> true: yes 
   int Graph_size;
   int exon_number;
   int RawSize;
   int size_;
   string strand, chr;
   map<pair<int,int>,bool> guided_as_map;

   vector< pair<int,int> > nodes,edges;
   vector<path_t> right_paths;
   vector<path_t> right_paths_remove_end_partial; 

   vector<path_t> anno_paths;
   vector<path_t> anno_paths_remove_end_partial; 

   map<node_idx_t, node_idx_t> rawNode_eNode_map;

   map<edge_t,edge_t > pair_Left_to_Right_edge_map;

   map<edge_t,edge_t > pair_Right_to_Left_edge_map;

   map<edge_t,vector<edge_with_cov_t> > Left_to_MultiRight_edge_map;
   map<edge_t,vector<edge_with_cov_t> > Right_to_MultiLeft_edge_map; 
   typedef map<edge_t,edge_t >::iterator iter;
   vector< pair<edge_t,double> > edge_descendcov;

   map<edge_t,double> EdgeCov_map;
   vector< edge_t >  reserved_junc;
   vector<path_t> LongReadPath;
   vector<double> LongReadPath_cov;
   vector<path_t> final_paths; //LRP
   vector<double> final_paths_cov_onlyforLRP; //LRP
   vector<path_t> final_paths_extended;
   vector<path_t> final_paths_temp;
   vector<path_info> final_paths_info;
   vector<path_t> final_paths_single_exon;


   map<pair<edge_t,edge_t>,int> packing_map;
   vector<vector<int> > nodes_of_components;
   vector< vector< pair<int,double> > >junctions_MappingInfo;
   vector< pair<int,int> > junctions;

   vector<block_info> Blocks;
   int BlockNumber;

public:
  void get_reserved_junc(vector<path_t> PairPath);
  void get_CovInfo();
  void add_source_and_sink();
  void get_raw_align_info();
  void get_annotation_info();
  void build_node_set_and_simplify(vector<int>& exon_l, vector<int>& exon_r, vector<double>& exon_cov,
		      vector<node_idx_t>& ,vector<node_idx_t>&, vector<double>&, 
		      vector<path_t>& LRP,vector<double> LRP_cov,
		      string strand_, string chr_,
		      bool dt_data_type);
  void show_block(block_info bi);
  void GetBlock();
  void get_full_length_path_from_blocks();
  void get_path_of_one_block(node_idx_t start, node_idx_t end,vector<path_t>& Bpaths, vector<double>& Bpaths_cov);
  void forward_extend(edge_t e, vector<path_t>& paths, node_idx_t end);
  void reverse_extend(edge_t e, vector<path_t>& paths,node_idx_t start);
  void get_end_node(map<int,bool>& head, map<int,bool>& tail);
  node_idx_t find_partial(node_idx_t n, bool downstream);
  void contract_LRP();
  void get_LRP();
  void process_LRP_after_simplify_graph();
  void get_components();
  double compute_JScov_for_Node(node_idx_t n, bool donorFlag);
  void compute_JScov_for_AllNodes();
  void remove_edges_basedon_JS_support();
  void remove_all_partial(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent);
  bool keep_edge(edge_t e);
  bool keep_edge2(edge_t e);


  void remove_edges_by_average_coverage(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent);
  void remove_unbalance_edges( vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent);
  bool MultiJunction(node_idx_t i, node_idx_t c);
  void remove_partial_junction(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent);



  void remove_small_exons(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent);
  bool partial(node_idx_t n1, node_idx_t n2);
  void remove_intron_contamination(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent);
  void remove_partial_end_by_edge_coverage(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent);

  void remove_edges_onemulti(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent);
  void remove_lowcov_edges_of_bifurcation_nodes2(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent);
  void remove_lowcov_edges_of_bifurcation_nodes(vector<edge_t>& delete_edges_child,vector<edge_t>& delete_edges_parent);

  void delete_children_edges(vector<edge_t> );
  void delete_parents_edges(vector<edge_t> );
  
  void show_graph(); 
  void show_junction();
  void contract_graph();
  void get_graph(vector<vector<int> >& vecEdges,vector<double>& vecCov);

  double InCov(node_idx_t);
  double OutCov(node_idx_t);

  double average_coverage();
  void get_packing_result_new(vector<vector<int> >Vec_edges, vector<int> Edges_legt, vector<int>Edges_right, vector<double> Weights);
  void forward_extend(edge_t e, vector<path_t>& paths);
  void reverse_extend(edge_t e, vector<path_t>& paths);
  bool check_overlap_of_2paths(path_t p1_, path_t p2_);
  void update_graph(path_t path);
  void path_search(string strand, string chr,bool dt_data_type);
  void remove_partial_end_for_a_path(path_t& pp);
  void remove_redundancy(vector<path_t>& Paths_,vector<int>& Paths_flag_);
  void output(string strand, string chr,vector<path_t> paths);
};

#endif
