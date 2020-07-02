
#include "cblas.h"
#include "lapacke.h"

#include "string.h"

typedef double *Matrix;
typedef double *Vector;

// M rows, N columns
Matrix zero_matrix(const int m, const int n) {
    double *result;

    result = (double *)calloc(n * m, sizeof(double));

    return result;
}

inline void set_element(Matrix a, const int m, const int n, const int i,
                        const int j, const double val) {
    a[i * m + j] = val;
}

// x should be of size n
inline void matmul(const Matrix a, const int m, const int n, Vector *x) {
    double *y = calloc(m, sizeof(double));
    cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, a, n, *x, 1, 0.0, y, 1);
    *x = y;
}

// Ku = v, solve for u, K is tridiagonal
int solve_tridiag_system(Vector KL, Vector K, Vector KU, const int n,
                         const Vector v, Vector *u) {
    double *u_ptr = calloc(n, sizeof(double));
    memcpy(u_ptr, v, n * sizeof(double));
    int ret = LAPACKE_dgtsv(LAPACK_ROW_MAJOR, n, 1, KL, K, KU, v, n);
    *u = u_ptr;
    return ret;
}

int main() {
    /* int N = 3; */
    /* int M = 4; */
    /*  */
    /* Matrix A = zero_matrix(M, N); */
    /*  */
    /* for (int i = 0; i < 3; i++) { */
    /*     set_element(A, N, M, i, i, 1); */
    /* } */
    /*  */
    /* double x[4] = {5.0, 3.0, 1.0, -1.0}; */
    /* double *x_ptr = x; */
    /*  */
    /* matmul(A, M, N, &x_ptr); */
    /*  */
    /* printf("A = ["); */
    /* for (int i = 0; i < N * M; i++) { */
    /*     printf("%f,", A[i]); */
    /* } */
    /* printf("]\n"); */
    /*  */
    /* printf("x = ["); */
    /* for (int i = 0; i < M; i++) { */
    /*     printf("%f,", x_ptr[i]); */
    /* } */
    /* printf("]\n"); */

    double KL[2] = {0.0, 0.0};
    double K[3] = {1.0, 1.0, 1.0};
    double KU[2] = {0.0, 0.0};
    double v[3] = {10.0, 20.0, 30.0};
    Vector u = NULL;
    solve_tridiag_system(KL, K, KU, 3, v, &u);
    for (int i = 0; i < 3; i++) {
        printf("%f,", u[i]);
    }
    printf("\n");
}
