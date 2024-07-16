#include <iostream>
#include "rlink.h"
#include "tmerge.h"
#ifndef COMMON_H
#define COMMON_H

/*
 *  common.h
 */

#include <vector>
#include <string>
#include <map>
#include <sys/stat.h>
#include <sys/types.h>
#include <bits/typesizes.h>
#include <bits/types.h>


// general options

extern vector<double> g_plus_gene;
extern vector<double> g_minus_gene;
extern vector<double> g_zero_gene;
extern vector<pair<int,int> > g_junction;
extern vector<double> g_junction_coverage;
extern vector<double> g_junction_coverage_good;
extern vector<char> g_junction_strand;
extern bool g_help;


#endif
