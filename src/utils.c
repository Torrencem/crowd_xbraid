
#include "cblas.h"
#include "lapacke.h"

#include "string.h"

#ifdef NOINLINE
#define utils_inline
#else
#define utils_inline inline
#endif

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

void matrix_destroy(Matrix a) { free(a); }

void vec_destroy(Vector a) { free(a); }

void matrix_copy(const int m, const int n, const Matrix a, Matrix b) {
    memcpy(b, a, m * n * sizeof(double));
}

void vec_copy(const int size, const double *a, double *b) {
    memcpy(b, a, size * sizeof(double));
}

// i is the row, j is the column
utils_inline void set_element(Matrix a, const int m, const int n, const int i,
                              const int j, const double val) {
    a[i * n + j] = val;
}

utils_inline double get_element(const Matrix a, const int m, const int n,
                                const int i, const int j) {
    return a[i * n + j];
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

Matrix tridiag_to_matrix(const Vector al, const Vector a, const Vector au,
                         const int n) {
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
utils_inline void matmul(const Matrix a, const int m, const int n, Vector *x) {
    double *y = calloc(m, sizeof(double));
    cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, a, n, *x, 1, 0.0, y, 1);
    *x = y;
}

// x should be of size x. al, a, and au should be of length n - 1, n, and n - 1
// respectively
utils_inline void matmul_tridiag(const Vector al, const Vector a,
                                 const Vector au, const int n, Vector *x) {
    Matrix a_as_mat = tridiag_to_matrix(al, a, au, n);
    matmul(a_as_mat, n, n, x);
}

// Ku = v, solve for u, K is tridiagonal
utils_inline int solve_tridiag_system(Vector KL, Vector K, Vector KU,
                                      const int n, Vector v) {
    int ret = LAPACKE_dgtsv(LAPACK_ROW_MAJOR, n, 1, KL, K, KU, v, 1);
    return ret;
}

#ifdef TESTS
#include <assert.h>
#include <math.h>

void test_row_major_set_get() {
    int N = 3;
    int M = 4;
    Matrix A = zero_matrix(M, N);

    set_element(A, M, N, 0, 1, 1.0);
    set_element(A, M, N, 1, 0, 2.0);
    set_element(A, M, N, 1, 1, 3.0);
    assert(get_element(A, M, N, 0, 1) == 1.0);
    assert(get_element(A, M, N, 1, 0) == 2.0);

    assert(A[1] == 1.0);
    assert(A[3] == 2.0);
    assert(A[4] == 3.0);
}

void test_simple_matmul() {
    int N = 3;
    int M = 4;
    Matrix A = zero_matrix(M, N);

    set_element(A, M, N, 0, 0, 2.0);
    set_element(A, M, N, 0, 1, 1.0);
    set_element(A, M, N, 1, 1, -1.0);

    double x_arr[4] = {5.0, 3.0, 1.0, -1.0};
    double *x = x_arr;

    matmul(A, M, N, &x);

    assert(x[0] == 13.0);
    assert(x[1] == -3.0);
    assert(x[2] == 0.0);
    assert(x[3] == 0.0);
}

double EPSILON = 1e-5;

void test_tridiag_solve() {
    // https://www.wolframalpha.com/input/?i=%7B%7B1%2C+2%2C+0%2C+0%7D%2C+%7B2%2C+1%2C+2%2C+0%7D%2C+%7B0%2C+2%2C+5%2C+2%7D%2C+%7B0%2C+0%2C+3%2C+1%7D%7D+*+%7B3%2C+14%2C+-15%2C+9%7D
    double KL[3] = {2, 2, 3};
    double K[4] = {1, 1, 5, 1};
    double KU[3] = {2, 2, 2};
    double v[4] = {31, -10, -29, -36};
    int info = solve_tridiag_system(KL, K, KU, 4, v);
    assert(info == 0);
    assert(fabs(v[0] - 3.0) <= EPSILON);
    assert(fabs(v[1] - 14.0) <= EPSILON);
    assert(fabs(v[2] - -15.0) <= EPSILON);
    assert(fabs(v[3] - 9.0) <= EPSILON);
}

void test_tridiag_mul() {
    double KL[3] = {2, 2, 3};
    double K[4] = {1, 1, 5, 1};
    double KU[3] = {2, 2, 2};
    double x_arr[4] = {3.0, 14.0, -15.0, 9.0};
    double *x = x_arr;
    matmul_tridiag(KL, K, KU, 4, &x);
    assert(fabs(x[0] - 31.0) <= EPSILON);
    assert(fabs(x[1] - -10.0) <= EPSILON);
    assert(fabs(x[2] - -29.0) <= EPSILON);
    assert(fabs(x[3] - -36.0) <= EPSILON);
}

int main() {
    test_row_major_set_get();
    test_simple_matmul();
    test_tridiag_solve();
    test_tridiag_mul();
}

#endif
