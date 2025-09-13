#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<map>
#include<vector>

using namespace std;
bool MyFlag = false;
string Flag = "isoQuant";
map<string,double> trans_count_map;
void load_Trans_count(char* file)
{
    ifstream in(file);
    string s;
    istringstream istr;

    if(Flag == "isoQuant")
    {
      getline(in,s);
      while(getline(in,s))
      {
        string id;
	double count;
	//cout<<s<<endl;
	istr.str(s);
	istr>>id>>count;
	istr.clear();
	trans_count_map[id] = count;
      }
    }
    else if(Flag == "Bambu")
    {
        getline(in,s);
	while(getline(in,s))
	{
	    string id,temp;
	    double count;
	    istr.str(s);
	    istr>>id>>temp>>count;
	    istr.clear();
	    trans_count_map[id] = count;
	}
    }
    else if(Flag == "transgram")
    {
        getline(in,s);
	while(getline(in,s))
	{
	    string id,count_;
	
	    size_t pos = s.find(":");
	    id = s.substr(0,pos);
	    count_ = s.substr(pos+1);
	    //cout<<id<<" "<<count_<<endl;
	    double count = atof(count_.c_str());
	    trans_count_map[id] = count;
	}
    }
    return;
}
int main(int argc,char* argv[])
{
    if(argc == 1)
    {
	cout<<""<<endl;
        cout<<"A simple program to get and filter transcripts based on the kallisto abundance!"<<endl;
	cout<<"./exe abundance.tsv input.gtf output_filtered.gtf filtercov(e.g.,1) flag(Bambu or isoQuant or transgram)"<<endl;
	cout<<endl;
	return 0;
    }

    Flag = argv[5];//Flag: isoquant or bambu or transgram
    load_Trans_count(argv[1]);//count.txt

    string s;
    ifstream in(argv[2]);//in.gtf

    ofstream out(argv[3]);//out.gtf

    string s_ = argv[4];//filtered cov
    double filter_cov = atof(s_.c_str());


    istringstream istr;
    while(getline(in,s))
    {
	if(s[0] == '#') continue;
        istr.str(s);
	string temp,flag, id;
	istr>>temp>>temp>>flag;
	if(flag == "gene") continue;
	while(istr>>temp)
	{
	    if(temp == "transcript_id") {
		    istr>>id;
	    	    break;
	    }
	}
	istr.clear();
	id = id.substr(1,id.length() - 3);//chr21_ML143377v1_fix.novel.38623.1 an example for novel id
	size_t first = id.find("novel");
	size_t last = id.rfind("novel");

	bool novel = false;
	if(first != std::string::npos && first == last)
	{
		novel = true;
	}
	if(novel)
	{
		if(trans_count_map.find(id) == trans_count_map.end()) {

			//cerr<<"WARRNING: "<<id<<endl;dd
		}
		else 
		{
	   		if(trans_count_map[id] > filter_cov)
		   		out<<s<<'\n';
		}
	}
	else{
	       	out<<s<<'\n';
	}


    }
    return 0;

}
