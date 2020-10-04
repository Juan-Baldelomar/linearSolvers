#ifndef SOLVERS_HPP
#define SOLVERS_HPP
#include <bits/stdc++.h>
#include <math.h>
using namespace std;
/*
   A: vector de vectores nxn
   b: es un vector de n
   x: es un vector de n
   Nota: 
     Todos los parámetros son paso por referencia, esto quiere decir que si se modifican dentro de la función también serán modificados fuera de la función, por lo tanto no es necesrio regresar ningun valor.
*/
void Diagonal(vector<vector<double> > &A, vector<double> &b, vector<double> &x);
void Triangular_Superior(vector<vector<double> > &U, vector<double> &b, vector<double> &x);
void Triangular_Inferior(vector<vector<double> > &L, vector<double> &b, vector<double> &x);
void Eliminacion_Gaussiana(vector<vector<double> > &A, vector<double> &b, vector<double> &x);
void Eliminacion_Gaussiana_Pivoteo(vector<vector<double> > &A, vector<double> &b, vector<double> &x);
void LU_Solve(vector<vector<double> > &A, vector<double>&b, vector<double>&x);
void Descomposicion_LU(vector<vector<double> > &A, vector<vector<double> > &L, vector<vector<double> > &U);
void Descomposicion_LU_Dolittle(vector<vector<double> > &A, vector<vector<double> > &L, vector<vector<double> > &U);
void LU_Solve(vector<vector<double> > &A, vector<double> &b, vector<double> &x);
void Matrix_Mult(vector<vector<double>>&A, vector<vector<double>>&B, vector<vector<double>>&C);
void VectorConstant_Mult(vector< double> &v, double c);
void MatrixLU_Mult(vector<vector<double>>&M, vector<vector<double>>&R);
void VectorConstant_Mult(vector< double> &v, double c);
void grand_schmidt(vector<vector<double>> &lista_vectores, vector<double> &v, int nValidVectors = 0, bool isNormalized = false);
void Try_Sol(vector< vector<double > > &A, vector<double>&b, vector<double>&x);
void Jacobi(vector<vector <double> > &A, vector < double> &b, vector< double> &x, double error, int max_it=1000);
void Gauss_Seidel(vector<vector <double> > &A, vector < double> &b, vector< double> &x, double error, int max_it =1000);
void metodoPotencia(vector < vector <double>> &A, vector<double> &v, double &lambda, double error, int max_it = 1000);
void metodoPotenciaDeflacion(vector < vector <double>> &A, vector<vector<double>> &eigenvectors, vector<double> &eigenvalues, double error, int range, int max_it=1000);
void PotenciaInversa(vector < vector < double >> &A, vector <double> &v, double &lambda, double tol, int max_it=1000) ;
void PotenciaInversaDeflacion(vector < vector <double>> &A, vector<vector<double>> &eigenvectors, vector<double> &eigenvalues, double tol, int range, int max_it);
void JacobiEigenValues(vector< vector < double>> &A, vector<double> &eigenvalues, double tol, int max_it=1000);
void cocienteRayleigh(vector<vector<double>> &A, vector<double> &v, double &lambda, double tol, int max_it=1000);
void metodoSubespacio(vector<vector< double>> &A, vector<vector<double>> &eigenvectors, vector<double> &eigenvalues, double tol, bool powerMethodType = true, int max_it=1000);
#endif
