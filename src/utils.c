
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

Matrix tridiag_to_matrix(const Vector al, const Vector a, const Vector au, const int n) {
    Matrix a_as_mat = zero_matrix(n, n);
    for (int row = 0; row < n; row++) {
        if (row == 0) {
            set_element(a_as_mat, n, n, 0, 0, a[0]);
            set_element(a_as_mat, n, n, 0, 1, au[0]);
        } else if (row == n - 1) {
            set_element(a_as_mat, n, n, n - 1, n - 1, a[n - 1]);
            set_element(a_as_mat, n, n, n - 1, n - 2, al[n - 2]);
        } else {
            set_element(a_as_mat, n, n, row, row - 1, al[row - 1]);
            set_element(a_as_mat, n, n, row, row, a[row]);
            set_element(a_as_mat, n, n, row, row + 1, au[row]);
        }
    }

    return a_as_mat;
}

// x should be of size n
inline void matmul(const Matrix a, const int m, const int n, Vector *x) {
    double *y = calloc(m, sizeof(double));
    cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, a, n, *x, 1, 0.0, y, 1);
    *x = y;
}

// x should be of size x. al, a, and au should be of length n - 1, n, and n - 1 respectively
void matmul_tridiag(const Vector al, const Vector a, const Vector au, const int n, Vector *x) {
    Matrix a_as_mat = tridiag_to_matrix(al, a, au, n);
    matmul(a_as_mat, n, n, x);
}

// Ku = v, solve for u, K is tridiagonal
int solve_tridiag_system(Vector KL, Vector K, Vector KU, const int n,
                         Vector v) {
    int ret = LAPACKE_dgtsv(LAPACK_ROW_MAJOR, n, 1, KL, K, KU, v, 1);
    return ret;
}

// TODO: Turn these into tests, include wolfram alpha link (like https://www.wolframalpha.com/input/?i=%7B%7B1%2C+2%2C+0%2C+0%7D%2C+%7B2%2C+1%2C+2%2C+0%7D%2C+%7B0%2C+2%2C+5%2C+2%7D%2C+%7B0%2C+0%2C+3%2C+1%7D%7D+*+%7B1%2C+2%2C+3%2C+4%7D)
#if 0
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
    /*  */
    double KL[3] = {2, 2, 3};
    double K[4] = {1, 1, 5, 1};
    double KU[3] = {2, 2, 2};
    double v[4] = {14, 28, 27, 13};
    /* int info = solve_tridiag_system(KL, K, KU, 4, v); */
    /* printf("Returned: %i\n", info); */
    /* for (int i = 0; i < 4; i++) { */
    /*     printf("%f,", v[i]); */
    /* } */
    /* printf("\n"); */
    double x[4] = {1, 2, 3, 4};
    double *x_ptr = x;
    matmul_tridiag(KL, K, KU, 4, &x_ptr);
    for (int i = 0; i < 4; i++) {
        printf("%f, ", x_ptr[i]);
    }
    printf("\n");
}
#endif
