#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <iomanip>
#include <vector>
using namespace std;
void ReadMatrix(vector<vector<double> > &M, string filename);
void ReadVector(vector<double> &V, string filename);
void WriteMatrix(const vector< vector<double > > &M, string filename, bool flag_append=false);
void WriteVector(const vector<double> &V, string filename, bool flag_append=false);
ostream &operator<<(ostream &os, const vector< vector<double> > &M);
ostream &operator<<(ostream &os, const vector<double> &V);
#endif
