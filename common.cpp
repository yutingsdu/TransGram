/*
 *  common.h
 */
#include "common.h"
//#include <sys/types.h>
//#include <sys/stat.h>
//#include <bits/types.h>
//#include <bits/typesizes.h>
#include <cstring>
#include <cstdlib>
#include <errno.h>
#include <libgen.h>

// general
vector<double> g_plus_gene;
vector<double> g_minus_gene;
vector<double> g_zero_gene;
vector<pair<int,int> > g_junction;
vector<double> g_junction_coverage;
vector<double> g_junction_coverage_good;
vector<char> g_junction_strand;
bool g_help = false;

