/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: juan
 *
 * Created on October 1, 2020, 7:04 PM
 */

#include <cstdlib>
#include "Tools.h"
#include "Solvers.h"

using namespace std;

/*
 * 
 */

void ejercicio1(){
    vector<vector<double>> A;
    vector<double> v;
    
    ReadMatrix(A, "Insumo_Tarea07/Eigen_3x3.txt");
    v.assign(A.size(), 0.0);
    double lambda;

    cocienteRayleigh(A, v, lambda, 0.000001, 1000);

    WriteVector(v, "Out/Eigenvector_3x3.txt");
    ReadMatrix(A, "Insumo_Tarea07/Eigen_3x3.txt");
    cout << "Valor propio mas pequeño "<< lambda << endl;
    
    cout << "Vector Propio asociado " <<  v << endl;
    printf("%.10f", lambda);

    vector<double> b;
    b.assign(v.size(), 0);

    for (int i = 0; i < v.size(); i++)
        b[i] = lambda * v[i];
    
    Try_Sol(A, b, v);
}

void ejercicio2(){
    vector<vector<double>> A;
    vector<double> v;

    ReadMatrix(A, "Insumo_Tarea07/Eigen_50x50.txt");
    v.assign(A.size(), 0.0);
    double lambda;

    cocienteRayleigh(A, v, lambda, 0.000001, 1000);

    WriteVector(v, "Out/Eigenvector_50x50.txt");
    ReadMatrix(A, "Insumo_Tarea07/Eigen_50x50.txt");
    cout << "Valor propio mas pequeño "<< lambda << endl;

    cout << "Vector Propio asociado " <<  v << endl;
    printf("%.10f", lambda);

    vector<double> b;
    b.assign(v.size(), 0);

    for (int i = 0; i < v.size(); i++)
        b[i] = lambda * v[i];

    Try_Sol(A, b, v);
}

int main() {
    ejercicio2();
    return 0;
}

