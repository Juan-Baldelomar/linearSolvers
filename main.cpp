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

void ejercicio3(){
    vector<vector<double>> A;
    vector<double> v;
    vector<vector<double>> eigenvectors;
    vector<double> eigenvalues;



    ReadMatrix(A, "Insumo_Tarea07/Eigen_3x3.txt");

    eigenvalues.assign(2, 0.0);
    eigenvectors.assign(2, vector<double>(A.size(), 0.0));


    metodoSubespacio(A, eigenvectors, eigenvalues, 0.000001, false, 1000);

    WriteVector(v, "Out/Eigenvector_3x3.txt");
    ReadMatrix(A, "Insumo_Tarea07/Eigen_3x3.txt");
    cout << "Valor propio mas pequeño "<< eigenvalues << endl;

    cout << "Vector Propio asociado " <<  eigenvectors << endl;

    vector<double> b;
    b.assign(eigenvectors[0].size(), 0);

    for (int i = 0; i < eigenvectors[1].size(); i++)
        b[i] = eigenvalues[1] * eigenvectors[1][i];

    Try_Sol(A, b, eigenvectors[1]);
}

void ejercicio4(){
    vector<vector<double>> A;
    vector<double> v;
    vector<vector<double>> eigenvectors;
    vector<double> eigenvalues;



    ReadMatrix(A, "Insumo_Tarea07/Eigen_50x50.txt");

    eigenvalues.assign(10, 0.0);
    eigenvectors.assign(10, vector<double>(A.size(), 0.0));


    metodoSubespacio(A, eigenvectors, eigenvalues, 0.000001, false, 1000);

    WriteVector(v, "Out/Eigenvector_50x50.txt");
    ReadMatrix(A, "Insumo_Tarea07/Eigen_50x50.txt");
    cout << "Valor propio mas pequeño "<< eigenvalues << endl;

    cout << "Vector Propio asociado " <<  eigenvectors << endl;

    vector<double> b;
    b.assign(eigenvectors[8].size(), 0);

    for (int i = 0; i < eigenvectors[8].size(); i++)
        b[i] = eigenvalues[8] * eigenvectors[8][i];

    Try_Sol(A, b, eigenvectors[8]);
}



void ejercicioP(){

    //variables necesarias
    vector<vector<double>> A;
    vector<vector<double>> eigenvectors;
    vector<double> eigenvalues;

    //entrada
    ReadMatrix(A, "Insumo_Tarea07/Eigen_50x50.txt");
    int nEigen = 10;

    //procesamiento
    PotenciaInversaDeflacion(A, eigenvectors, eigenvalues, 0.000001, nEigen, 3000);

    //salida
    WriteMatrix(eigenvectors, "Out/N_eigenvectors_50x50.txt");
    WriteVector(eigenvalues, "Out/N_eigenvalues_50x50.txt");
    cout << "los eigenvalores son "<< eigenvalues << endl;
    cout << "los eigenvectores (almacenados cada uno en una fila) son "<< endl << eigenvectors << endl;

    for (int n = 0; n < nEigen; n++) {
        ReadMatrix(A, "Insumo_Tarea07/Eigen_50x50.txt");
        vector<double> b;
        b.assign(50, 0.0);

        for (int i = 0; i < 50 ; i++)
            b[i] = eigenvalues[n] * eigenvectors[n][i];


        Try_Sol(A, b, eigenvectors[n]);
    }
}


int main() {
    ejercicio3();
    return 0;
}

