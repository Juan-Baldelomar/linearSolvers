

#include "Solvers.h"
#include "Tools.h"

//sobrecarga metodos
void LU_pivot(vector<vector<double> > &A, int k);

void LU_pivot(vector<vector<double> > &A, int k, vector<double> &b, vector<int> &index);

void Partial_pivot(vector<vector<double> > &A, int k, vector<double> &b);

void pivot(vector<vector<double> > &A, int k, vector<double> &b, vector<int> &index);

void Diagonal(vector<vector<double> > &A, vector<double> &b, vector<double> &x) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        if (A[i][i] == 0)
            cout << "Diag WARNING: Elemento nulo en posicion [" << i << "][" << i << "]" << endl;

        x[i] = b[i] / A[i][i];
    }

}

// ******************************************************************************** LINEAR SOLVERS ********************************************************************************
// ******************************************************************************** LINEAR SOLVERS ********************************************************************************
// ******************************************************************************** LINEAR SOLVERS ********************************************************************************
// ******************************************************************************** LINEAR SOLVERS ********************************************************************************
// ******************************************************************************** LINEAR SOLVERS ********************************************************************************


void Triangular_Superior(vector<vector<double> > &U, vector<double> &b, vector<double> &x) {
    int n = U.size();
    double X[n];

    for (int i = n - 1; i >= 0; i--) {
        double acc = 0;
        for (int j = n - 1; j > i; j--) {
            acc += U[i][j] * x[j];
        }
        if (U[i][i] == 0)
            cout << "Triangular Sup  Msg WARNING: A_[" << i << ", " << i << "] es nulo, se dividira por 0" << endl;

        x[i] = (b[i] - acc) / U[i][i];
    }
}

void Triangular_Inferior(vector<vector<double> > &L, vector<double> &b, vector<double> &x) {
    int n = L.size();

    for (int i = 0; i < n; i++) {
        double acc = 0;
        for (int j = 0; j < i; j++) {
            acc += L[i][j] * x[j];
        }
        if (L[i][i] == 0)
            cout << "Triangular Inf  Msg WARNING: A_[" << i << ", " << i << "] es nulo, se dividira por 0" << endl;
        x[i] = (b[i] - acc) / L[i][i];
    }
}

/*
 * U[i][2] valores de diagonal superior
 * U[i][1] valores de diagonal principal
 * index d indicates col position of de diagonal
 * */
void Triangular_Superior_TriD(vector<vector<double>>&U, vector<double>&b, vector<double>&x, int d){
    int n = U.size();
    for (int i =n-1; i>=0; i--){
        double acc = 0;
        if (i<n-1){
            x[i] = (b[i]- x[i+1] * U[i][d+1]) / U[i][d];
        }else{
            x[i] = b[i] / U[i][d];
        }
    }
}


/*
 * L[i][0] valores de diagonal inferior
 * L[i][1] valores de diagonal principal
 * index d indicates col position of de diagonal
 * */
void Triangular_Inferior_TriD(vector<vector<double>>&L, vector<double>&b, vector<double>&x, int d){
    int n = L.size();
    for (int i =0; i<n; i++){
        double acc = 0;
        if (i>0){
            x[i] = (b[i]-x[i-1]*L[i][d-1])/L[i][d];
        }else{
            x[i] = b[i]/L[i][d];
        }
    }
}

void Eliminacion_Gaussiana(vector<vector<double> > &A, vector<double> &b, vector<double> &x) {

    int n = A.size();
    //indice k indica la fila con la cual estamos transformando a una Triangular Superior
    for (int k = 0; k < n - 1; k++) {

        if (A[k][k] == 0)
            Partial_pivot(A, k, b);

        if (A[k][k] == 0)
            cout << "GEM  Msg WARNING: A_[" << k << ", " << k << "] es nulo, se dividira por 0" << endl;

        for (int i = k + 1; i < n; i++) {
            double m_ik = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) {
                A[i][j] = A[i][j] - m_ik * A[k][j];
            }
            b[i] = b[i] - m_ik * b[k];
        }
    }
    Triangular_Superior(A, b, x);
}

void Eliminacion_Gaussiana_Pivoteo(vector<vector<double> > &A, vector<double> &b, vector<double> &x) {
    int n = A.size();
    vector<int> index;
    vector<double> orderedX;

    //inicializar
    orderedX.assign(n, 0);
    index.assign(n, 0);

    // llenar de las posiciones correspondientes de cada variable
    for (int i = 0; i < n; i++) {
        index[i] = i;
    }

    //indice k indica la fila con la cual estamos transformando a una Triangular Superior
    for (int k = 0; k < n - 1; k++) {
        pivot(A, k, b, index);
        if (A[k][k] == 0)
            cout << "GEM  Msg WARNING: A_[" << k << ", " << k << "] es nulo, se dividira por 0" << endl;

        for (int i = k + 1; i < n; i++) {
            double m_ik = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) {
                A[i][j] = A[i][j] - m_ik * A[k][j];
            }
            b[i] = b[i] - m_ik * b[k];
        }

    }
    Triangular_Superior(A, b, x);

    //ordenar arreglo x para dar la solucion
    for (int i = 0; i < n; i++) {
        orderedX[index[i]] = x[i];
    }

    //copiar data a arreglo X
    for (int i = 0; i < n; i++) {
        x[i] = orderedX[i];
    }

}

// LU Fact in same Matrix used to solve Equation System

void LU_Solve(vector<vector<double> > &A, vector<double> &b, vector<double> &x) {

    //variables necesarias
    int n = A.size();
    vector<int> index;
    vector<double> orderedX;

    //inicializar vectores
    orderedX.assign(n, 0);
    index.assign(n, 0);

    // llenar de las posiciones correspondientes de cada variable
    for (int i = 0; i < n; i++) {
        index[i] = i;
    }

    for (int k = 0; k < n; k++) {
        double acc; // variable que contiene la suma

        LU_pivot(A, k, b, index);

        //l_ik 
        for (int i = k; i < n; i++) {
            acc = 0;

            for (int r = 0; r <= k - 1; r++) {
                acc += A[i][r] * A[r][k];
            }
            A[i][k] = A[i][k] - acc;
        }

        // verificar si se dividira por 0 y lanzar un warning
        if (A[k][k] == 0)
            cout << "LU Solve Msg WARNING: A_[" << k << ", " << k << "] es nulo, se dividira por 0" << endl;

        //u_kj 
        for (int j = k + 1; j < n; j++) {
            acc = 0;
            for (int r = 0; r <= k - 1; r++) {
                acc += A[k][r] * A[r][j];
            }

            A[k][j] = (A[k][j] - acc) / A[k][k];
        }
    }

    //RESOLVER SISTEMA

    Triangular_Inferior(A, b, x);

    //Transformar para obtener U
    for (int i = 0; i < n; i++) {
        A[i][i] = 1;
    }

    Triangular_Superior(A, x, x);

    //ordenar arreglo x para dar la solucion
    for (int i = 0; i < n; i++) {
        orderedX[index[i]] = x[i];
    }

    //copiar data a arreglo X
    for (int i = 0; i < n; i++) {
        x[i] = orderedX[i];
    }
}

//M es una matriz de nx4
void TriDLU_Solve(vector<vector<double>> &M, vector<double>&b, vector<double>&x){
    Triangular_Inferior_TriD(M, b, x, 1);
    Triangular_Superior_TriD(M, x, x, 2);
}



// ******************************************************************************** DESCOMPONER MATRIZ ********************************************************************************
// ******************************************************************************** DESCOMPONER MATRIZ ********************************************************************************
// ******************************************************************************** DESCOMPONER MATRIZ ********************************************************************************
// ******************************************************************************** DESCOMPONER MATRIZ ********************************************************************************
// ******************************************************************************** DESCOMPONER MATRIZ ********************************************************************************
// ******************************************************************************** DESCOMPONER MATRIZ ********************************************************************************


// Metodo para descomponer una matriz unicamente

void Descomposicion_LU(vector<vector<double> > &A, vector<vector<double> > &L, vector<vector<double> > &U) {
    int n = A.size();

    // inicializar posiciones de U
    for (int i = 0; i < n; i++) {
        U[i][i] = 1;
    }

    for (int k = 0; k < n; k++) {
        double acc; // variable que contiene la suma

        //LU_pivot(A, k);

        //l_ik 
        for (int i = k; i < n; i++) {
            acc = 0;

            for (int r = 0; r <= k - 1; r++) {
                acc += L[i][r] * U[r][k];
            }
            L[i][k] = A[i][k] - acc;
        }

        // verificar division por 0
        if (L[k][k] == 0)
            cout << "LU Decomposition  Msg WARNING: L_[" << k << ", " << k << "] es nulo, se dividira por 0" << endl;

        //u_kj 
        for (int j = k + 1; j < n; j++) {
            acc = 0;
            for (int r = 0; r <= k - 1; r++) {
                acc += L[k][r] * U[r][j];
            }
            U[k][j] = (A[k][j] - acc) / L[k][k];
        }
    }

    //    //PRUEBA LU FUNCIONA
    //    vector<vector<double>> C;
    //    C.assign(n, vector<double>(n, 0.0));
    //    Matrix_Mult(L, U, C);
    //    cout << "MULT R" << endl;
    //    cout << C << endl;
}


// EXTRA

void Descomposicion_LU_Dolittle(vector<vector<double> > &A, vector<vector<double> > &L, vector<vector<double> > &U) {
    int n = A.size();
    // inicializar posiciones de U
    for (int i = 0; i < n; i++) {
        L[i][i] = 1;
    }

    for (int k = 0; k < n; k++) {
        double acc; // variable que contiene la suma
        //u_kj
        for (int j = k; j < n; j++) {
            acc = 0;

            for (int r = 0; r <= k - 1; r++) {
                acc += L[k][r] * U[r][j];
            }
            U[k][j] = A[k][j] - acc;
        }

        //l_ik
        for (int i = k + 1; i < n; i++) {
            acc = 0;
            for (int r = 0; r <= k - 1; r++) {
                acc += L[i][r] * U[r][k];
            }
            L[i][k] = (A[i][k] - acc) / U[k][k];
        }
    }

    //PRUEBA LU FUNCIONA
    vector<vector<double>> C;
    C.assign(n, vector<double>(n, 0.0));
    Matrix_Mult(L, U, C);
    cout << "MULT R" << endl;
    cout << C << endl;
}


/**
 * A[i][0] diagonal inferior
 * A[i][1] diagonal principal
 * A[i][2] diagonal superior
 *
 * LU[i][0] Diagonal Inferior L
 * LU[i][1] Diagonal Principal L
 * LU[i][2] diagonal principal U
 * LU[i][3] diagonal superior U
 */
void Descomposicion_LU_Dolittle_TriD(vector<vector<double> > &A, vector<vector<double>> &LU) {
    int n = A.size();
    // inicializar posiciones de U
    LU[0][2] = A[0][1]; //Diagonal Principal U
    LU[0][3] = A[0][2]; //Diagonal Superior U
    LU[0][1] = 1; //Diagonal Principal L

    for (int i = 1; i<n; i++){
        // Diagonal Inferior L
        LU[i][0] = A[i][0]/LU[i-1][2];
        // Diagonal Principal L
        LU[i][1] = 1;
        // Diagonal Superior U
        if (i<n-1)
            LU[i][3] = A[i][2];

        //Diagonal Principal U
        LU[i][2] = A[i][1]-LU[i][0]*LU[i-1][3];
    }
}

void FactorizarCholesky(vector<vector<double> > &A, vector<vector<double> > &H) {

    int n = A.size();

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            double acc = 0;

            if (i != j) {
                //h_ij
                for (int k = 0; k < j; k++) {
                    acc += H[i][k] * H[j][k];
                }
                H[i][j] = (A[i][j] - acc) / H[j][j];
            }
        }

        double acc = 0;
        // h_ii
        for (int k = 0; k < i; k++) {
            acc += H[i][k] * H[i][k];
        }

        H[i][i] = sqrt(A[i][i] - acc);

    }
    //calcular transpuesta
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            H[i][j] = H[j][i];
        }
    }
}

void factorizarLDLT(vector<vector<double> > &A, vector<vector<double> > &L, vector<vector<double> > &DLT) {
    int n = A.size();

    vector<double> D;
    D.assign(n, 0.0);

    for (int i = 0; i < n; i++) {


        L[i][i] = 1;
        //L_ij
        for (int j = 0; j <= i; j++) {

            double acc = 0;
            // D_jj
            for (int k = 0; k < j; k++) {
                acc += L[j][k] * L[j][k] * D[k];
            }
            D[j] = A[j][j] - acc;


            acc = 0;

            if (i > j) {
                //h_ij
                for (int k = 0; k < j; k++) {
                    acc += L[i][k] * L[j][k] * D[k];
                }
                L[i][j] = (A[i][j] - acc) / D[j];
            }
        }
    }
    //calcular DLT
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            DLT[i][j] = L[j][i];
            L[j][i] = L[j][i] * D[i];
        }
    }
}

/*void FactorizarCholeskyTriD(vector<vector<double> > &A, vector<vector<double> > &H) {
    int n = A.size();

    for (int i = 0; i < n; i++) {
        // no need for acc, hir*hjr = 0

        // h_ij
        if (i != 0)
            H[i][i - 1] = A[i][i - 1] / H[i - 1][i - 1];
        // h_ii
        H[i][i] = sqrt(A[i][i] - H[i][i - 1] * H[i][i - 1]);
    }

    //calcular transpuesta
    for (int i = 0; i < n - 1; i++) {
        H[i][i + 1] = H[i + 1][i];
    }
}*/


void FactorizarCholeskyTriD(vector<vector<double>> &D, vector<vector<double>>&LLT) {
    int n = D.size();
    for (int i = 0; i < n; i++) {

        // h_ij
        if (i != 0)
            LLT[i][0] = D[i][0]/ LLT[i - 1][1];
        // h_ii
        LLT[i][1] = sqrt(D[i][1] - LLT[i][0] * LLT[i][0]);
    }
    //calcular transpuesta
    for (int i = 0; i < n - 1; i++) {
        LLT[i][2] = LLT[i+1][0];
    }
}



// ******************************************************************************** ALGORITMOS PARA PIVOTEO ********************************************************************************
// ******************************************************************************** ALGORITMOS PARA PIVOTEO ********************************************************************************
// ******************************************************************************** ALGORITMOS PARA PIVOTEO ********************************************************************************
// ******************************************************************************** ALGORITMOS PARA PIVOTEO ********************************************************************************
// ******************************************************************************** ALGORITMOS PARA PIVOTEO ********************************************************************************
// ******************************************************************************** ALGORITMOS PARA PIVOTEO ********************************************************************************
// ******************************************************************************** ALGORITMOS PARA PIVOTEO ********************************************************************************



//GEM Pivot

void pivot(vector<vector<double> > &A, int k, vector<double> &b, vector<int> &index) {
    int n = A.size();
    int i_max = k, j_max = k;
    double max = fabs(A[k][k]);

    //find MAX
    for (int i = k; i < n; i++) {
        for (int j = k; j < n; j++) {
            if (fabs(A[i][j]) > max) {
                i_max = i;
                j_max = j;
                max = fabs(A[i][j]);
            }
        }
    }
    double temp;
    // swap file
    if (i_max != k) {
        for (int j = 0; j < n; j++) {
            temp = A[k][j];
            A[k][j] = A[i_max][j];
            A[i_max][j] = temp;
        }
        temp = b[k];
        b[k] = b[i_max];
        b[i_max] = temp;
    }
    // swap column
    if (j_max != k) {
        for (int i = 0; i < n; i++) {
            temp = A[i][k];
            A[i][k] = A[i][j_max];
            A[i][j_max] = temp;
        }
        //swap variable order
        int at_j = index[j_max];
        int at_k = index[k];
        index[j_max] = at_k;
        index[k] = at_j;
    }
}

//GEM Partial Pivot

void Partial_pivot(vector<vector<double> > &A, int k, vector<double> &b) {
    int n = A.size();
    int i_max = k;
    double max = fabs(A[k][k]);

    //find MAX
    for (int i = k + 1; i < n; i++) {
        if (fabs(A[i][k]) > max) {
            i_max = i;
            max = fabs(A[i][k]);
        }
    }

    double temp;
    // swap file
    if (i_max != k) {
        for (int j = 0; j < n; j++) {
            temp = A[k][j];
            A[k][j] = A[i_max][j];
            A[i_max][j] = temp;
        }
        temp = b[k];
        b[k] = b[i_max];
        b[i_max] = temp;
    }
}

// Metodo para pivotear cuando se descompone una matriz en LU 

void LU_pivot(vector<vector<double> > &A, int k) {
    int n = A.size();
    int i_max = k, j_max = k;
    double max = fabs(A[k][k]);

    //find MAX
    for (int i = k; i < n; i++) {
        for (int j = k; j < n; j++) {
            if (fabs(A[i][j]) > max) {
                i_max = i;
                j_max = j;
                max = fabs(A[i][j]);
            }
        }
    }
    double temp;
    // swap file
    if (i_max != k) {
        for (int j = 0; j < n; j++) {
            temp = A[k][j];
            A[k][j] = A[i_max][j];
            A[i_max][j] = temp;
        }
    }
    // swap column
    if (j_max != k) {
        for (int i = 0; i < n; i++) {
            temp = A[i][k];
            A[i][k] = A[i][j_max];
            A[i][j_max] = temp;
        }
    }
}

// Metodo para pivotear con LU_Solve. Es otro metodo aparte de LU_Decomposition_Pivot porque 

void LU_pivot(vector<vector<double> > &A, int k, vector<double> &b, vector<int> &index) {
    int n = A.size();
    int i_max = k, j_max = k;
    double max = fabs(A[k][k]);

    //find MAX
    for (int i = k; i < n; i++) {
        for (int j = k; j < n; j++) {
            if (fabs(A[i][j]) > max) {
                i_max = i;
                j_max = j;
                max = fabs(A[i][j]);
            }
        }
    }

    double temp;

    // swap file
    if (i_max != k) {
        for (int j = 0; j < n; j++) {
            temp = A[k][j];
            A[k][j] = A[i_max][j];
            A[i_max][j] = temp;
        }
        temp = b[k];
        b[k] = b[i_max];
        b[i_max] = temp;
    }

    // swap column
    if (j_max != k) {
        for (int i = 0; i < n; i++) {
            temp = A[i][k];
            A[i][k] = A[i][j_max];
            A[i][j_max] = temp;
        }
        //swap variable order
        int at_j = index[j_max];
        int at_k = index[k];
        index[j_max] = at_k;
        index[k] = at_j;
    }
}


// ******************************************************************************** MATRIX OPERATIONS ********************************************************************************
// ******************************************************************************** MATRIX OPERATIONS ********************************************************************************
// ******************************************************************************** MATRIX OPERATIONS ********************************************************************************
// ******************************************************************************** MATRIX OPERATIONS ********************************************************************************
// ******************************************************************************** MATRIX OPERATIONS ********************************************************************************
// ******************************************************************************** MATRIX OPERATIONS ********************************************************************************
// ******************************************************************************** MATRIX OPERATIONS ********************************************************************************
// ******************************************************************************** MATRIX OPERATIONS ********************************************************************************
// ******************************************************************************** MATRIX OPERATIONS ********************************************************************************




// multiplicacion de matrices

void Matrix_Mult(vector<vector<double> > &A, vector<vector<double> > &B, vector<vector<double>> &C) {
    int n = A.size();
    int m = B[0].size();
    int p = A[0].size();
    omp_set_num_threads(omp_get_num_procs());
    int i, j, k;
    #pragma omp parallel for private(i,j,k) shared(A,B,C,n, m, p) default(none)
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            //cout << i << " " << j << endl;
            double acc = 0;
            for (k = 0; k < p; k++) {
                acc += A[i][k] * B[k][j];
            }
            C[i][j] = acc;
        }
    }
}

//recibe dos matrices y realiza la multiplicacion transponiendo la primera
void Transpose_Matrix_Mult(vector<vector<double> > &AT, vector<vector<double> > &B, vector<vector<double>> &C) {
    int n = AT[0].size();
    int m = B[0].size();
    int p = AT.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double acc = 0;
            for (int k = 0; k < p; k++) {
                acc += AT[k][i] * B[k][j];
            }
            C[i][j] = acc;
        }
    }
}


// multiplicacion de matrices

void MatrixVector_Mult(vector<vector<double> > &A, vector<double> &v, vector<double> &v_r) {
    int n = A.size();
    int m = v.size();
    int p = A[0].size();
    if (p != m) {
        cout << "ERROR: size of matrix and vector dont match" << endl;
        return;
    }
    for (int i = 0; i < n; i++) {
        double acc = 0;
        for (int j = 0; j < m; j++)
            acc += A[i][j] * v[j];

        v_r[i] = acc;
    }
}

// stores the result in same vector, also paralelizable
void MatrixVector_Mult(vector<vector<double> > &A, vector<double> &v) {
    int n = A.size();
    int m = v.size();
    vector<double> tmp(n, 0.0);
    if (n != m) {
        cout << "ERROR: size of matrix and vector dont match" << endl;
        return;
    }
    for (int i = 0; i < n; i++) {
        double acc = 0;
        for (int j = 0; j < m; j++)
            acc += A[i][j] * v[j];

        tmp[i] = acc;
    }

    for (int i = 0; i < n; i++)
        v[i] = tmp[i];
}


double TransposeVectorMult(vector<double> &vT, vector<double> &w) {
    int n = vT.size();
    double acc = 0;

    for (int i = 0; i < n; i++)
        acc += vT[i] * w[i];
    return acc;
}

void VectorConstant_Mult(vector<double> &v, double c) {
    int n = v.size();
    for (int i = 0; i < n; i++)
        v[i] = v[i] * c;
}

// normalizar un vector
void normalize(vector<double> &v) {
    int n = v.size();
    double suma = 0;
    for (int i = 0; i < n; i++)
        suma += v[i] * v[i];

    suma = sqrt(suma);

    for (int i = 0; i < n; i++)
        v[i] = v[i] / suma;
}

// obtener la norma de un vector
double normaVector(vector<double> &v) {
    int n = v.size();
    double suma = 0;
    for (int i = 0; i < n; i++)
        suma += v[i] * v[i];

    return sqrt(suma);
}

double normaCuadradaVector(vector<double> &v) {
    int n = v.size();
    double suma = 0;
    for (int i = 0; i < n; i++)
        suma += v[i] * v[i];

    return suma;
}


void transposeMatrix(vector<vector<double>> &M, vector<vector<double>> &MT) {
    int n = M.size();
    int m = M[0].size();
    if (MT.empty())
        MT.assign(m, vector<double>(n, 0.0));

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            MT[i][j] = M[j][i];
        }
    }
}



// multiplicar matriz Triangular Superior e Inferior almacenadas en misma matriz

void MatrixLU_Mult(vector<vector<double> > &M, vector<vector<double>> &R) {
    int n = M.size();

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double acc = 0;
            for (int k = 0; k <= i && k <= j; k++) {
                acc += M[i][k] * M[k][j];
            }
            R[i][j] = acc;
        }
    }
}


void TriD_Matrix_Mult(vector<double>&D, vector<double>&L, vector<double>&U){
   /* int n = D.size();
    for (int i =0; i<n; i++){
        if (i == 0)
            cout << L[i]*L[i] << endl;
        else{
            double l = L[i]*D[i];
            double d = D[i]*D[i]+l*l;
            cout << l << "  " << d << " "
        }
    }*/
}


// ortogonalizar un vector respecto a una lista de  vectores
// nValidVectors = 0 por defecto. Este parametro indica si se debe tomar todos los vectores de la lista si tiene el valor de 0 o de lo contrario considera el valor enviado
// isNormalized = false por defecto
void grand_schmidt(vector<vector<double>> &lista_vectores, vector<double> &v, int nValidVectors, bool isNormalized) {
    int nVectors;
    int n = v.size();
    if (nValidVectors < 0)
        nVectors = lista_vectores.size();
    else
        nVectors = nValidVectors;

    for (int j = 0; j < nVectors; j++) {
        double a_k = TransposeVectorMult(lista_vectores[j], v);
        if (!isNormalized)
            a_k = a_k / normaCuadradaVector(lista_vectores[j]);

        for (int k = 0; k < n; k++) {
            v[k] = v[k] - a_k * lista_vectores[j][k];
        }
    }
}

void Identity(vector<vector<double>> &M) {
    int n = M.size();
    int m = M[0].size();

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j)
                M[i][i] = 1;
            else
                M[i][j] = 0;
        }
    }
}

void copyMatrix(vector<vector<double>> &Origin, vector<vector<double>> &Dest) {
    int n = Origin.size();
    int m = Origin[0].size();

    if (n != Dest.size() || m != Dest[0].size()) {
        cout << "ERR COPY MATRIX: SIZES DONT MATCH" << endl;
        return;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Dest[i][j] = Origin[i][j];
        }
    }
}

void copyVector(vector<double> &origin, vector<double> &dest){
    int n = origin.size();
    if (n!= dest.size()){
        cout << "ERR CPY VECTOR: Sizes dont match" << endl;
        return;
    }
    for (int i =0; i<n; i++)
        dest[i] = origin[i];
}

void copyVectorToCol(vector<double> &origin, vector<vector<double>> &dest, int col){
    int n = origin.size();
    if (n!= dest.size()){
        cout << "ERR CPY VECTOR TO COL: Sizes dont match" << endl;
        return;
    }
    for (int i =0; i<n; i++)
        dest[i][col] = origin[i];
}

void copyColToCol(vector<vector<double>> &origin, vector<vector<double>> &dest, int colOrigin, int colDest){
    int n = origin.size();
    if (n!= dest.size()){
        cout << "ERR CPY COL TO COL: Sizes dont match" << endl;
        return;
    }

    if (colOrigin>=origin[0].size() || colDest>= dest[0].size()){
        cout << "ERR: CPY COL TO COL: COLS DONT EXIST" <<endl;
        return;
    }

    for (int i =0; i<n; i++)
        dest[i][colDest] = origin[i][colOrigin];
}

double findMax(vector<vector<double>>&A, int &i_max, int &j_max){
    int n = A.size();
    int i, j;
    vector<double>max(n, 0.0);
    vector<int>max_j(n, 0);
    omp_set_num_threads(omp_get_num_procs());
    #pragma omp parallel for private(i, j) shared(A, max, max_j, n) default(none)
    for (i = 0; i<n; i++){
        for (j=0; j<n; j++){
            if (i == j)
                continue;

            if(max[i]<fabs(A[i][j])){
                max[i] = fabs(A[i][j]);
                max_j[i] = j;
            }
        }
    }

    double realMax = 0;
    for (int k=0; k<n;k++){
        if (max[k]>realMax){
            realMax = max[k];
            i_max = k;
            j_max = max_j[k];
        }
    }
    return realMax;
}


// Metodo para probar que el vector solucion de la igualdad esperada

void Try_Sol(vector<vector<double> > &A, vector<double> &b, vector<double> &x) {
    cout
            << "-------------------------------------------------- PRUEBA VECTOR SOLUCION -------------------------------------------------- "
            << endl;
    int n = A.size();

    for (int i = 0; i < n; i++) {
        double acc = 0;
        for (int j = 0; j < n; j++) {
            acc += A[i][j] * x[j];
        }
        cout << "valor b_" << i << ":" << setw(5) << b[i] << setw(15) << " Fila A_" << i << " * x:" << setw(5) << acc
             << endl;
    }
}



// ******************************************************************************** METODOS ITERATIVOS ********************************************************************************
// ******************************************************************************** METODOS ITERATIVOS ********************************************************************************
// ******************************************************************************** METODOS ITERATIVOS ********************************************************************************
// ******************************************************************************** METODOS ITERATIVOS ********************************************************************************
// ******************************************************************************** METODOS ITERATIVOS ********************************************************************************
// ******************************************************************************** METODOS ITERATIVOS ********************************************************************************
// ******************************************************************************** METODOS ITERATIVOS ********************************************************************************




void Jacobi(vector<vector<double> > &A, vector<double> &b, vector<double> &x, double error, int max_it) {

    int n = b.size();
    vector<double> xn;
    xn.assign(n, 0);

    //inicializacion
    for (int i = 0; i < n; i++) {
        x[i] = b[i] / A[i][i];
    }

    // calculo de nuevo vector 
    for (int it = 0; it < max_it; it++) {
        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < n; j++) {
                if (i == j)
                    continue;

                sum += A[i][j] * x[j];
            }
            xn[i] = (b[i] - sum) / A[i][i];
        }

        double num = 0;
        double den = 0;
        double fact;
        // actualizacion de vector y calculo de error
        for (int i = 0; i < n; i++) {
            double dif_v_num = (xn[i] - x[i]);
            double dif_v_den = xn[i] * xn[i];
            num += dif_v_num * dif_v_num;
            den += dif_v_den;
            x[i] = xn[i];
        }
        fact = sqrt(num);

        if (fact < error) {
            cout << "iteracion: " << it << endl;
            return;
        }
    }

    cout << "Sistema no pudo converger con " << max_it << " iteraciones " << endl;

}

void Gauss_Seidel(vector<vector<double> > &A, vector<double> &b, vector<double> &x, double error, int max_it) {

    int n = b.size();
    vector<double> xn;
    xn.assign(n, 0);

    // inicializacion
    for (int i = 0; i < n; i++) {
        x[i] = b[i] / A[i][i];
    }

    // calculo de nuevo vector
    for (int it = 0; it < max_it; it++) {
        for (int i = 0; i < n; i++) {
            double sum = 0;
            //valores nuevos
            for (int j = 0; j < i; j++) {
                sum += A[i][j] * xn[j];
            }
            // valores viejos
            for (int j = i + 1; j < n; j++) {
                sum += A[i][j] * x[j];
            }

            xn[i] = (b[i] - sum) / A[i][i];
        }


        double num = 0;
        double den = 0;
        double fact;

        //actualizacion del vectors y calculo de error
        for (int i = 0; i < n; i++) {
            double dif_v_num = (xn[i] - x[i]);
            double dif_v_den = xn[i] * xn[i];
            num += dif_v_num * dif_v_num;
            den += dif_v_den;
            x[i] = xn[i];
        }
        fact = sqrt(num);

        if (fact < error) {
            cout << "iteracion: " << it << endl;
            return;
        }

    }

    cout << "Sistema no pudo converger con " << max_it << " iteraciones " << endl;

}

/*
 * para el siguiente metodo la primera iteracion es realizada afuera del ciclo para poder optimizar las iteraciones siguientes
 * la idea consiste en que la multiplicacion A*v es utilizada para el calculo de lambda, pero tambien en la siguiente iteracion
 * entonces para no repetir el calculo, en el ciclo posterior unicamente actualizamos el valor del vector nuevo previo a entrar a la iteracion
 */
void metodoPotencia(vector<vector<double>> &A, vector<double> &v, double &lambda, double error, int max_it) {
    int n = v.size();
    vector<double> v0;
    lambda = 0;
    v0.assign(n, 1);

    // primera iteracion realizada afuera del ciclo
    MatrixVector_Mult(A, v0, v);
    normalize(v);

    for (int it = 0; it < max_it; it++) {
        // calculo A*v utilizado posteriormente
        MatrixVector_Mult(A, v, v0);
        lambda = TransposeVectorMult(v, v0);

        double num = 0;

        // actualizacion del nuevo vector
        for (int i = 0; i < n; i++) {
            double dif = v0[i] - lambda * v[i];
            num += dif * dif;
            v[i] = v0[i];
        }

        // la norma es calculada aca para que al momento de salir el eigenvector este normalizado
        normalize(v);
        if (sqrt(num) < error)
            return;
    }
    cout << "WARNING: Metodo de potencia no pudo converger con " << max_it << " iteaciones" << endl;


}

void
metodoPotenciaDeflacion(vector<vector<double>> &A, vector<vector<double>> &eigenvectors, vector<double> &eigenvalues,
                        double error, int range, int max_it) {
    int n = A.size();
    double a_k; // componente del vector propio a eliminar
    double lambda = 0;
    vector<double> v0;
    eigenvalues.assign(range, 0.0);
    eigenvectors.assign(range, vector<double>(n, 0.0));


    for (int e_i = 0; e_i < range; e_i++) {
        v0.assign(n, 1);

        // correcion del vector para eliminar componentes de vectores propios ya encontrados
        for (int k = 0; k < e_i; k++) {
            a_k = TransposeVectorMult(eigenvectors[k], v0);
            for (int j = 0; j < n; j++) {
                v0[j] = v0[j] - eigenvectors[k][j] * a_k;
            }
        }

        // al igual que antes, primera iteracion afuera del ciclo
        MatrixVector_Mult(A, v0, eigenvectors[e_i]);
        normalize(eigenvectors[e_i]);

        for (int it = 0; it < max_it; it++) {

            MatrixVector_Mult(A, eigenvectors[e_i], v0);
            lambda = TransposeVectorMult(eigenvectors[e_i], v0);

            // correcion del vector para eliminar componentes de vectores propios ya encontrados
            for (int k = 0; k < e_i; k++) {
                a_k = TransposeVectorMult(eigenvectors[k], v0);
                for (int j = 0; j < n; j++) {
                    v0[j] = v0[j] - eigenvectors[k][j] * a_k;
                }
            }

            // actualizacion del vector para la nueva iteracion
            double num = 0;
            for (int i = 0; i < n; i++) {
                double dif = v0[i] - lambda * eigenvectors[e_i][i];
                num += dif * dif;
                eigenvectors[e_i][i] = v0[i];
            }

            normalize(eigenvectors[e_i]);

            if (sqrt(num) < error)
                break;
        }
        eigenvalues[e_i] = lambda;
    }
}

void PotenciaInversa(vector<vector<double >> &A, vector<double> &v, double &lambda, double tol, int max_it) {
    int n = A.size();

    //declaracion e inicializacion de variables necesarias
    vector<vector<double>> L, U;
    vector<double> v0, Av;
    L.assign(n, vector<double>(n, 0.0));
    U.assign(n, vector<double>(n, 0.0));
    v0.assign(n, 1.0);
    Av.assign(n, 0.0);

    // descomposicion de matriz
    Descomposicion_LU(A, L, U);

    // norma vector inicial
    normalize(v0);

    for (int it = 0; it < max_it; it++) {

        //RESOLVER SISTEMA
        Triangular_Inferior(L, v0, v);
        Triangular_Superior(U, v, v);

        // calculo de eigenvalor
        normalize(v);
        MatrixVector_Mult(A, v, Av);
        lambda = TransposeVectorMult(v, Av);

        //actualizar iteracion
        double error = 0;
        for (int i = 0; i < n; i++) {
            double dif = Av[i] - lambda * v[i];
            error += dif * dif;
            v0[i] = v[i];
        }// i

        if (sqrt(error) < tol)
            return;
    }// it

    cout << "Metodo no pudo converger " << endl;

}

void
PotenciaInversaDeflacion(vector<vector<double>> &A, vector<vector<double>> &eigenvectors, vector<double> &eigenvalues,
                         double tol, int range, int max_it) {
    int n = A.size();

    //declaracion e inicializacion de variables necesarias
    vector<vector<double>> L, U;
    vector<double> v0, Av;
    L.assign(n, vector<double>(n, 0.0));
    U.assign(n, vector<double>(n, 0.0));
    v0.assign(n, 1.0);
    Av.assign(n, 0.0);

    double a_k; //componente del vector propio que se desea eliminar
    double lambda = 0;
    eigenvalues.assign(range, 0.0);
    eigenvectors.assign(range, vector<double>(n, 0.0));

    // Descomposicion Matriz
    Descomposicion_LU(A, L, U);

    int total_iterations = 0;

    for (int e_i = 0; e_i < range; e_i++) {
        v0.assign(n, 1);

        // correcion del vector para eliminar componentes de vectores propios ya encontrados
        for (int k = 0; k < e_i; k++) {
            a_k = TransposeVectorMult(eigenvectors[k], v0);
            for (int j = 0; j < n; j++) {
                v0[j] = v0[j] - eigenvectors[k][j] * a_k;
            }
        }

        // norma vector inicial
        normalize(v0);

        for (int it = 0; it < max_it; it++) {

            //RESOLVER SISTEMA
            Triangular_Inferior(L, v0, eigenvectors[e_i]);
            Triangular_Superior(U, eigenvectors[e_i], eigenvectors[e_i]);

            // calculo de eigenvalor
            normalize(eigenvectors[e_i]);
            MatrixVector_Mult(A, eigenvectors[e_i], Av);
            lambda = TransposeVectorMult(eigenvectors[e_i], Av);

            //actualizar iteracion
            double error = 0;
            for (int i = 0; i < n; i++) {
                double dif = Av[i] - lambda * eigenvectors[e_i][i];
                error += dif * dif;
                v0[i] = eigenvectors[e_i][i];
            }// i

            if (sqrt(error) < tol) {
                total_iterations += it + 1;
                break;
            }


            // correcion del vector para eliminar componentes de vectores propios ya encontrados
            // la correccion se realiza despues de verificar el error en caso de no ser necesario calculara
            for (int k = 0; k < e_i; k++) {
                a_k = TransposeVectorMult(eigenvectors[k], v0);
                for (int j = 0; j < n; j++) {
                    v0[j] = v0[j] - eigenvectors[k][j] * a_k;
                }
            }
        }
        eigenvalues[e_i] = lambda;
        cout << e_i << ":" << lambda << endl;
    }

    cout << "convergencia en iteracion: " << total_iterations << endl;
}

void JacobiEigenValues(vector<vector<double>> &A, vector<double> &eigenvalues, double tol, int max_it) {
    int n = A.size();
    int i_max = 0, j_max = 0;

    vector<vector<double>> G, A_tmp;

    //inicializar G
    G.assign(n, vector<double>(n, 0.0));
    A_tmp.assign(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; i++)
        G[i][i] = 1;


    for (int it = 0; it < max_it; it++) {
        //find MAX 
      //  double max = 0;
        double max = findMax(A, i_max, j_max);

       /* for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j)
                    continue;
                if (fabs(A[i][j]) > max) {
                    i_max = i;
                    j_max = j;
                    max = fabs(A[i][j]);
                }
            }
        }*/
        cout << max << " " << i_max << " " << j_max << endl;

        // criterio de paro
        if (max < tol) {
            for (int i = 0; i < n; i++)
                eigenvalues[i] = A[i][i];
            return;
        }

        // calculo de valores para matriz G
        double delta = (A[j_max][j_max] - A[i_max][i_max]) / (2 * A[i_max][j_max]);
        double sign = delta < 0 ? -1 : delta > 0 ? 1 : 0;
        double t = sign / (fabs(delta) + sqrt(1 + delta * delta));
        double c = 1 / (sqrt(1 + t * t));
        double s = c * t;

        // calculo de matriz G
        G[i_max][i_max] = G[j_max][j_max] = c;
        G[i_max][j_max] = s;
        G[j_max][i_max] = -s;

        // A*G
        Matrix_Mult(A, G, A_tmp);

        // G^T
        double tmp = G[i_max][j_max];
        G[i_max][j_max] = G[j_max][i_max];
        G[j_max][i_max] = tmp;

        // G^T(AG)
        Matrix_Mult(G, A_tmp, A);

        // Restaurar G
        G[i_max][i_max] = G[j_max][j_max] = 1;
        G[i_max][j_max] = G[j_max][i_max] = 0;
        //cout << it << endl;
    } //it

    cout << "Metodo no pudo converger " << endl;
}

int JacobiKrylov(vector<vector<double>> &A, vector<double> &eigenvalues, vector<vector<double>> &G, double tol,
                 int max_it) {
    int n = A.size();
    int i_max = 0, j_max = 0;
    int i_max_prev = 1, j_max_prev = 0;

    vector<vector<double>> A_tmp;

    //inicializar G
    G.assign(n, vector<double>(n, 0.0));
    A_tmp.assign(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; i++)
        G[i][i] = 1;


    for (int it = 0; it < max_it; it++) {

        //find MAX
        double max = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j)
                    continue;
                if (fabs(A[i][j]) > max) {
                    i_max = i;
                    j_max = j;
                    max = fabs(A[i][j]);
                }
            }
        }

        // criterio de paro
        if (max < tol) {
            // G[i_max_prev][j_max_prev] = -G[i_max_prev][j_max_prev];
            // G[j_max_prev][i_max_prev] = -G[j_max_prev][i_max_prev];

            for (int i = 0; i < n; i++)
                eigenvalues[i] = A[i][i];
            return it;
        }

        // Restaurar G
        G[i_max_prev][i_max_prev] = G[j_max_prev][j_max_prev] = 1;
        G[i_max_prev][j_max_prev] = G[j_max_prev][i_max_prev] = 0;

        // actualizar valores previos
        i_max_prev = i_max;
        j_max_prev = j_max;

        // calculo de valores para matriz G
        double delta = (A[j_max][j_max] - A[i_max][i_max]) / (2 * A[i_max][j_max]);
        double sign = delta < 0 ? -1 : delta > 0 ? 1 : 0;
        double t = sign / (fabs(delta) + sqrt(1 + delta * delta));
        double c = 1 / (sqrt(1 + t * t));
        double s = c * t;


        // calculo de matriz G
        G[i_max][i_max] = G[j_max][j_max] = c;
        G[i_max][j_max] = s;
        G[j_max][i_max] = -s;

        // A*G
        Matrix_Mult(A, G, A_tmp);

        // G^T
        double tmp = G[i_max][j_max];
        G[i_max][j_max] = G[j_max][i_max];
        G[j_max][i_max] = tmp;

        // G^T(AG)
        Matrix_Mult(G, A_tmp, A);


    } //it

    cout << "WARNING JACOBI :Metodo no pudo converger " << endl;
    return max_it;
}

void cocienteRayleigh(vector<vector<double>> &A, vector<double> &v, double &lambda, double tol, int max_it) {
    int n = A.size();

    //inicializar
    double rho;
    vector<vector<double>> L, DLT, A_rho;
    vector<double> Av;

    L.assign(n, vector<double>(n, 0.0));
    DLT.assign(n, vector<double>(n, 0.0));
    A_rho.assign(n, vector<double>(n, 0.0));
    v.assign(n, 1.0);
    Av.assign(n, 0.0);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            A_rho[i][j] = A[i][j];

    normalize(v);

    // calcular rho_0
    MatrixVector_Mult(A, v, Av);
    rho = TransposeVectorMult(v, Av);

    // iterar
    for (int it = 0; it < max_it; it++) {

        //calcular A_rho
        for (int i = 0; i < n; i++)
            A_rho[i][i] = A[i][i] - rho;

        // factorizar matriz
        factorizarLDLT(A_rho, L, DLT);


        //resolver sistema
        Triangular_Inferior(L, v, v);
        Triangular_Superior(DLT, v, v);

        //calcular nuevo vector de iteracion normalizado
        normalize(v);

        // calcular rho_k
        MatrixVector_Mult(A, v, Av);
        rho = TransposeVectorMult(v, Av);

        lambda = rho;

        // verificar tolerancia
        double err = 0;
        for (int i = 0; i < n; i++) {
            double dif = Av[i] - lambda * v[i];
            err += dif * dif;
        }
        if (sqrt(err) < tol){
            cout << "Convergencia Iteracion: " << it + 1 << endl;
            return;
        }
    }
    cout << "WARNING RAYLEIGH: Metodo no pudo converger " << endl;

}

/*
 * eigenvectors en realidad sera la matriz de vectores propios pero transpuesta para mayor facilidad de poder accesar a ellos y operar con ellos
 * powerMethodType = 1 inverse power methdo
 * powerMethodType = 0 power method
 * */
void metodoSubespacio(vector<vector<double>> &A, vector<vector<double>> &eigenvectors, vector<double> &eigenvalues,
                      double tol, bool powerMethodType, int max_it) {
    int n = A.size();
    int s = eigenvectors.size();

    //inicializacion
    vector<double> Av;

    // matrices para factorizar, phi^T*A, matriz de vectores propios transpuesta, rotacion de Givens, Matriz de metodo de Jacobi
    vector<vector<double>> L, U, phiT_A, phi_A, eigenvectors_T, G, J;

    L.assign(n, vector<double>(n, 0.0));
    U.assign(n, vector<double>(n, 0.0));

    //matriz de eigenvectores y eigenvectores transpuestos
    phiT_A.assign(s, vector<double>(n, 0.0));
    phi_A.assign(n, vector<double>(s, 0.0));

    // Matriz de Givens y de iteracion de Jacobi
    G.assign(s, vector<double>(s, 0.0));
    J.assign(s, vector<double>(s, 0.0));

    // inicializar identidad
    Identity(eigenvectors);

    //descomposicion
    if (powerMethodType)
        factorizarLDLT(A, L, U);


    for (int it = 0; it < max_it; it++) {

        //iteration of inverse power or power method with deflation
        for (int ev_pos = 0; ev_pos < s; ev_pos++) {
            grand_schmidt(eigenvectors, eigenvectors[ev_pos], ev_pos, true);
            if (powerMethodType) {
                Triangular_Inferior(L, eigenvectors[ev_pos], eigenvectors[ev_pos]);
                Triangular_Superior(U, eigenvectors[ev_pos], eigenvectors[ev_pos]);
            } else {
                MatrixVector_Mult(A, eigenvectors[ev_pos]);
            }
            normalize(eigenvectors[ev_pos]);
        }


        /* calcular transpuesta de matriz phi (eigenvectores)
         * recordemos que almacenamos los eigenvectores por fila y no por columna para facilidad de operarlos y accesarlos como eigenvector[pos]
         * entonces al querer multiplicar la matriz de eigenvectores por columnas como se hace en la teoria, debemos primero transponer la matriz de eigenvectores
         * */
        transposeMatrix(eigenvectors, eigenvectors_T);

        /* A_phi = A*phi (como los eigenvectores los accesamos por fila, eigenvectors_T = phi)
         * J =  phi_T* A_phi
         * */
        Matrix_Mult(A, eigenvectors_T, phi_A);
        Matrix_Mult(eigenvectors, phi_A, J);

        /*
         * otra forma de realizar la operacion superior, pero sin la variabla phi_A declarada. Optimizamos espacio
         * utilizar esta forma en caso de querer calcular el error viendo de la matriz J es diagonal o G es la identidad
         * de manera que el error en los eigenvectores sera mayor pero la aproximacion de los eigenvalores sera optima y el metodo va a converge mas rapido
         * */
        //Matrix_Mult(eigenvectors, A, phiT_A);
        //Matrix_Mult(phiT_A, eigenvectors_T, J);


        /* aplicar Jacobi
         * Obtener matriz de Givens transpuesta (G^T)
         * */
        it += JacobiKrylov(J, eigenvalues, G, tol, max_it) + 1;

        /* aplicar rotacion de givens a matriz phiT_A
         *  phi^(k) = phi^(k-1) * G
         *  phi^(k) = G^T * phi^{(k-1), T}   [recordando (AB)^T = B^T*A^T]
         * */
        Matrix_Mult(G, eigenvectors, phiT_A);

        // actualizar para nueva iteracion
        double err = 0;
        for (int i = 0; i < s; i++) {
            double v_err = 0;
            //v_err = (G[i][i] - 1)*(G[i][i]-1);
            for (int j = 0; j < n; j++) {
                eigenvectors[i][j] = phiT_A[i][j];
                double dif = eigenvalues[i] * eigenvectors[i][j] - phi_A[j][i];
                v_err += dif * dif;
            }

            if (sqrt(v_err) > err)
                err = sqrt(v_err);
        }

        if (err < tol) {
            cout << " Convergencia en iteracion: " << it << endl;
            return;
        }
    }

    cout << "WARNING METODO SUBESPACIO: metodo no pudo converger" << endl;

}

