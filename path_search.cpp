//#include"expression_level.h"

//#include"get_junction_graph_new.h"
//#include"junction_paths_new.h"
//#include"recover_paths.h"
//#include"PairPacker.h"
#include "path_search.h"
using namespace std;

string out_name;
double SampleSize;
bool unstranded;
bool SFlag;
int Path_length = 300;
double AVERAGE_REMOVE_RATE = 0.05;
double UNBALANCE_RATE = 0.03;
float SEED_Filter = 2;//1.01;
int rg_index = 0;
int trans_id;
bool MyFlag=false;
int main (int argc, char* argv[])
{
    out_name = argv[2];
    load_graph(argv[1]);
    return 0;
}
