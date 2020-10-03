#include <bits/stdc++.h>
#include "Tools.h"
/*
   Realiza la lectura de una matriz, se hace paso por referencia.
   M: es la matriz que tiene los valores.
   datafile: es el fichero donde está almacenada la matriz
*/
void ReadMatrix(vector<vector<double> > &M, string filename)
{
    std::ifstream readf(filename);
    int n,m;
    readf >>n>>m; 
    M.assign(n, vector<double> (m, 0.0));
    for(int i = 0; i < n; i++)
    {
       for(int j = 0; j < m; j++)
       {
          readf >> M[i][j];
       }
    }
   readf.close();
}
/*
   Realiza la lectura del vector, se hace paso por referencia.
   V: es el vector que tiene los valores.
   datafile: es el fichero donde está almacenada el vector
*/
void ReadVector(vector<double> &V, string filename)
{
    std::ifstream readf(filename);
    int n, m;
    readf >>n>>m; 
    V.assign(n, 0.0);
    for(int i = 0; i < n; i++)
    {
          readf >> V[i];
    }
   readf.close();
}
void WriteMatrix(const vector< vector<double > > &M, string filename, bool flag_append)
{
    std::fstream fout;
    
    if(flag_append)
       fout.open(filename,fstream::app|fstream::out );
    else
    fout.open(filename,std::ios::out);
    fout << (int)M.size() << " " << (int)M[0].size() <<endl;
    for(int i = 0; i < M.size(); i++)
    {
       for(int j = 0; j < M[i].size(); j++)
	fout <<setw(12) << M[i][j]<<  setw(12);
	fout << endl;
    }
    fout.close();
}
void WriteVector(const vector<double> &V, string filename, bool flag_append)
{
    std::fstream fout;
    if(flag_append)
       fout.open(filename,fstream::app|fstream::out );
    else
    fout.open(filename,std::ios::out);
    fout << (int)V.size() << " 1 "<<endl;
    for(int i = 0; i < V.size(); i++)
	fout <<V[i]<< endl;
    fout.close();

}

ostream &operator<<(ostream &os, const vector< vector<double> > &M)
{
   for(int i = 0; i < M.size(); i++)
   {
      for(auto x:M[i])
      {
	 os << setw(12) << x << setw(12);
      }
      os <<endl;
   }
   return os;
}
ostream &operator<<(ostream &os, const vector<double> &V)
{
      for(auto x:V)
	 os << x << " ";
      os <<endl;
   return os;
}


