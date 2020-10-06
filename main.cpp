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

void ejercicio1() {
    vector<vector<double>> A;
    vector<double> v;

    ReadMatrix(A, "Insumo_Tarea07/Eigen_3x3.txt");
    v.assign(A.size(), 0.0);
    double lambda;

    cocienteRayleigh(A, v, lambda, 0.000001, 1000);

    WriteVector(v, "Out/Eigenvector_3x3.txt");
    ReadMatrix(A, "Insumo_Tarea07/Eigen_3x3.txt");
    cout << "Valor propio  " << lambda << endl;

    cout << "Vector Propio asociado " << v << endl;
    printf("%.10f", lambda);

    vector<double> b;
    b.assign(v.size(), 0);

    for (int i = 0; i < v.size(); i++)
        b[i] = lambda * v[i];

    Try_Sol(A, b, v);
}

void ejercicio2() {
    vector<vector<double>> A;
    vector<double> v;

    ReadMatrix(A, "Insumo_Tarea07/Eigen_50x50.txt");
    v.assign(A.size(), 0.0);
    double lambda;

    cocienteRayleigh(A, v, lambda, 0.000001, 1000);

    WriteVector(v, "Out/Eigenvector_50x50.txt");
    ReadMatrix(A, "Insumo_Tarea07/Eigen_50x50.txt");
    cout << "Valor propio  " << lambda << endl;

    cout << "Vector Propio asociado " << v << endl;
    printf("%.10f", lambda);

    vector<double> b;
    b.assign(v.size(), 0);

    for (int i = 0; i < v.size(); i++)
        b[i] = lambda * v[i];

    Try_Sol(A, b, v);
}

void ejercicio3() {
    int nEigen = 2;
    vector<vector<double>> A;
    vector<vector<double>> eigenvectors;
    vector<double> eigenvalues;

    ReadMatrix(A, "Insumo_Tarea07/Eigen_3x3.txt");

    eigenvalues.assign(2, 0.0);
    eigenvectors.assign(2, vector<double>(A.size(), 0.0));

    metodoSubespacio(A, eigenvectors, eigenvalues, 0.000001, true, 1000);

    WriteVector(eigenvalues, "Out/Eigenvector_3x3_SUB.txt");
    ReadMatrix(A, "Insumo_Tarea07/Eigen_3x3.txt");
    cout << "Valor propio mas pequeño " << eigenvalues << endl;

    cout << "Vector Propio asociado " << eigenvectors << endl;

    vector<double> b;
    b.assign(eigenvectors[0].size(), 0);
    for (int n = 0; n < nEigen; n++) {
        for (int i = 0; i < eigenvectors[1].size(); i++)
            b[i] = eigenvalues[n] * eigenvectors[n][i];

        Try_Sol(A, b, eigenvectors[n]);
    }
}

void ejercicio4() {
    vector<vector<double>> A;
    vector<vector<double>> eigenvectors;
    vector<double> eigenvalues;


    ReadMatrix(A, "Insumo_Tarea07/Eigen_50x50.txt");

    eigenvalues.assign(10, 0.0);
    eigenvectors.assign(10, vector<double>(A.size(), 0.0));


    metodoSubespacio(A, eigenvectors, eigenvalues, 0.000001, true, 1000);

    WriteVector(eigenvalues, "Out/Eigenvector_50x50_SUB.txt");
    ReadMatrix(A, "Insumo_Tarea07/Eigen_50x50.txt");
    cout << "Valor propio " << eigenvalues << endl;

    cout << "Vector Propio asociado " << eigenvectors << endl;

    vector<double> b;
    b.assign(eigenvectors[8].size(), 0);

    for (int i = 0; i < eigenvectors[8].size(); i++)
        b[i] = eigenvalues[8] * eigenvectors[8][i];

    Try_Sol(A, b, eigenvectors[8]);
}


void ejercicioP() {

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
    cout << "los eigenvalores son " << eigenvalues << endl;
    cout << "los eigenvectores (almacenados cada uno en una fila) son " << endl << eigenvectors << endl;

    for (int n = 0; n < nEigen; n++) {
        ReadMatrix(A, "Insumo_Tarea07/Eigen_50x50.txt");
        vector<double> b;
        b.assign(50, 0.0);

        for (int i = 0; i < 50; i++)
            b[i] = eigenvalues[n] * eigenvectors[n][i];


        Try_Sol(A, b, eigenvectors[n]);
    }
}


void produceTriD_Sim_Entry(int n, vector<double> &D, vector<double> &L) {
    D.assign(n, 0.0);
    L.assign(n, 0.0);

    for (int i = 0; i < n; i++) {
        D[i] = 20;
        if (i>0)
            L[i] = -i -1;
    }
    cout << D<< endl;
    cout << L << endl;
}

void pruevaTRID(){
    vector<double> D, L, U;
    U.assign(30, 0.0);
    produceTriD_Sim_Entry(30, D, L);
    FactorizarCholeskyTriD(D, L, U);

    vector<vector<double>> MLD, MLU, A;
    MLD.assign(30, vector<double>(30, 0.0));
    MLU.assign(30, vector<double>(30, 0.0));
    A.assign(30, vector<double>(30, 0.0));

    for (int i = 0; i<30; i++){
        for(int j = 0; j<30; j++){
            if (i==j){
                MLD[i][j] = D[i];
                MLU[i][j] = D[i];
            }
            else if (j==i-1)
                MLD[i][j] = L[i];
            else if (j==i+1)
                MLU[i][j]=U[i];
        }
    }

    Matrix_Mult(MLD, MLU, A);
    cout << A << endl;

}

void pruevaTRIDLU(){
    vector<double> D, L;
    produceTriD_Sim_Entry(30, D, L);
    vector<vector<double>> MLD, MLU, A, AF, U;
    U.assign(30, vector<double>(2, 0.0));
    AF.assign(30, vector<double>(3, 0.0));
    copyVectorToCol(D, AF, 1);
    copyVectorToCol(L, AF, 0);
    copyVectorToCol(L, AF,2);
    Descomposicion_LU_Dolittle_TriD(AF, L, U);


    MLD.assign(30, vector<double>(30, 0.0));
    MLU.assign(30, vector<double>(30, 0.0));
    A.assign(30, vector<double>(30, 0.0));

    for (int i = 0; i<30; i++){
        for(int j = 0; j<30; j++){
            if (i==j){
                MLD[i][j] = 1;
                MLU[i][j] = U[i][0];
            }
            else if (j==i-1)
                MLD[i][j] = L[i];
            else if (j==i+1)
                MLU[i][j]=U[i][1];
        }
    }

    Matrix_Mult(MLD, MLU, A);
    cout << A << endl;

}

int main() {
    //ejercicio4();
    pruevaTRIDLU();
    return 0;
}

