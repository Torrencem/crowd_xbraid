
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

Vector zero_vector(const int size) {
    double *result;

    result = (double *)calloc(size, sizeof(double));

    return result;
}

inline void matrix_destroy(Matrix a) { free(a); }

inline void vec_destroy(Vector a) { free(a); }

void matrix_copy(const int m, const int n, const Matrix a, Matrix b) {
    memcpy(b, a, m * n * sizeof(double));
}

void vec_copy(const int size, const double *a, double *b) {
    memcpy(b, a, size * sizeof(double));
}

inline void set_element(Matrix a, const int m, const int n, const int i,
                        const int j, const double val) {
    a[i * m + j] = val;
}

void vec_axpy(const int size, const double alpha, const Vector x,
              const double beta, Vector y) {
    for (int i = 0; i < size; i++) {
        y[i] = alpha * x[i] + beta * y[i];
    }
}

void matrix_axpy(const int m, const int n, const double alpha, const Matrix a,
                 const double beta, Matrix b) {
    for (int i = 0; i < m * n; i++) {
        b[i] = alpha * a[i] + beta * b[i];
    }
}

void vec_scale(const int size, const double alpha, Vector x) {
    int i;
    for (i = 0; i < size; i++) {
        x[i] = alpha * x[i];
    }
}

// x should be of size n
inline void matmul(const Matrix a, const int m, const int n, Vector *x) {
    double *y = calloc(m, sizeof(double));
    cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, a, n, *x, 1, 0.0, y, 1);
    *x = y;
}

// Ku = v, solve for u, K is tridiagonal
int solve_tridiag_system(Vector KL, Vector K, Vector KU, const int n,
                         Vector v) {
    int ret = LAPACKE_dgtsv(LAPACK_ROW_MAJOR, n, 1, KL, K, KU, v, 1);
    return ret;
}

#if 0
int main() {
    int N = 3;
    int M = 4;

    Matrix A = zero_matrix(M, N);

    for (int i = 0; i < 3; i++) {
        set_element(A, N, M, i, i, 1);
    }

    double x[4] = {5.0, 3.0, 1.0, -1.0};
    double *x_ptr = x;

    matmul(A, M, N, &x_ptr);

    printf("A = [");
    for (int i = 0; i < N * M; i++) {
        printf("%f,", A[i]);
    }
    printf("]\n");

    printf("x = [");
    for (int i = 0; i < M; i++) {
        printf("%f,", x_ptr[i]);
    }
    printf("]\n");

    double KL[3] = {2, 2, 3};
    double K[4] = {1, 1, 5, 1};
    double KU[3] = {2, 2, 2};
    double v[4] = {14, 28, 27, 13};
    int info = solve_tridiag_system(KL, K, KU, 4, v);
    printf("Returned: %i\n", info);
    for (int i = 0; i < 4; i++) {
        printf("%f,", v[i]);
    }
    printf("\n");
}
#endif
