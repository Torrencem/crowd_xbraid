
#include "cblas.h"
#include "lapacke.h"

#include "string.h"

#include <assert.h>

#include "utils.hpp"

#ifdef NOINLINE
#define utils_inline
#else
#define utils_inline extern inline
#endif

Vector_::Vector_(const int len_) noexcept : len(len_) {
    this->data = new double[len_] { 0.0 };
}

// Copy constructor
Vector_::Vector_(const Vector_ &v2) noexcept {
    this->data = new double[v2.len] { 0.0 };
    this->len = v2.len;
    for (int i = 0; i < this->len; i++) {
        this->data[i] = v2.data[i];
    }
}

Vector_::~Vector_() noexcept{
    delete []this->data;
}

double &Vector_::operator[](int i) noexcept {
    assert(0 <= i && i <= this->len);
    return this->data[i];
}

Vector Vector_::operator*(double val) noexcept {
    Vector res = new Vector_(this->len);
    for (int i = 0; i < len; i++) {
        (*res)[i] = this->data[i] * val;
    }
    return res;
}

Vector Vector_::operator+(double val) noexcept {
    Vector res = new Vector_(this->len);
    for (int i = 0; i < len; i++) {
        (*res)[i] = this->data[i] + val;
    }
    return res;
}

Vector Vector_::operator+(const Vector &m2) noexcept {
    assert(m2->len == this->len);
    Vector res = new Vector_(this->len);
    for (int i = 0; i < this->len; i++) {
        (*res)[i] = this->data[i] + m2->data[i];
    }
    return res;
}

Vector Vector_::operator-(double val) noexcept {
    Vector res = new Vector_(this->len);
    for (int i = 0; i < len; i++) {
        (*res)[i] = this->data[i] - val;
    }
    return res;
}

Vector Vector_::operator-(const Vector &m2) noexcept {
    assert(m2->len == this->len);
    Vector res = new Vector_(this->len);
    for (int i = 0; i < this->len; i++) {
        (*res)[i] = this->data[i] - m2->data[i];
    }
    return res;
}

Tridiag_Matrix_::Tridiag_Matrix_(int n_) noexcept : n(n_) {
    assert(n == 16);
    this->al = new double[n - 1] { 0.0 };
    this->a = new double[n] { 0.0 };
    this->au = new double[n - 1] { 0.0 };
}

Tridiag_Matrix_::Tridiag_Matrix_(const Tridiag_Matrix &m2) noexcept {
    this->n = m2->n;
    this->al = (double *) malloc((n - 1) * sizeof(double));
    this->a = (double *) malloc(n * sizeof(double));
    this->au = (double *) malloc((n - 1) * sizeof(double));
    for (int i = 0; i < n; i++) {
        if (i != n - 1) {
            this->al[i] = m2->al[i];
            this->au[i] = m2->au[i];
        }
        this->a[i] = m2->a[i];
    }
}

Tridiag_Matrix_::~Tridiag_Matrix_() noexcept {
    delete []al;
    delete []a;
    delete []au;
}

// NOTE: this function invalidates the tridiag_matrix due
// to the nature of how lapacke does the solve in-place.
// A copy can be used to prevent this
Vector Tridiag_Matrix_::operator*(const Vector &x) noexcept {
    assert(x->len == n);
    Matrix as_mat = new Matrix_(this);
    return (*as_mat) * x;
}

// ---------- Matrix ----------

// Create a new matrix of M rows and N columns
Matrix_::Matrix_(const int m_, const int n_) noexcept : m(m_), n(n_) {
    this->data = new double[m_ * n_] { 0.0 };
}

// Copy constructor
Matrix_::Matrix_(const Matrix &m2) noexcept {
    this->data = new double[m2->m * m2->n];
    this->m = m2->m;
    this->n = m2->n;
    for (int i = 0; i < m * n; i++) {
        this->data[i] = m2->data[i];
    }
}

Matrix_::Matrix_(const Tridiag_Matrix &m2) noexcept {
    this->m = m2->n;
    this->n = m2->n;
    this->data = new double[m2->n * m2->n] { 0.0 };
    for (int row = 0; row < n; row++) {
        if (row == 0) {
            element(0, 0) = m2->a[0];
            element(0, 1) = m2->au[0];
        } else if (row == n - 1) {
            double val = m2->a[n - 1];
            element(n - 1, n - 1) = val;
            element(n - 1, n - 2) = m2->al[n - 2];
        } else {
            element(row, row - 1) = m2->al[row - 1];
            element(row, row) = m2->a[row];
            element(row, row + 1) =  m2->au[row];
        }
    }
}

Matrix_::~Matrix_() noexcept { delete[] this->data; }

double &Matrix_::element(const int row, const int col) noexcept {
    return this->data[row * this->n + col];
}

Matrix Matrix_::operator*(double val) noexcept {
    Matrix res = new Matrix_(m, n);
    for (int i = 0; i < m * n; i++) {
        res->data[i] = this->data[i] * val;
    }
    return res;
}

Matrix Matrix_::operator+(double val) noexcept {
    Matrix res = new Matrix_(m, n);
    for (int i = 0; i < m * n; i++) {
        res->data[i] = this->data[i] + val;
    }
    return res;
}

Matrix Matrix_::operator+(const Matrix &m2) noexcept {
    assert(this->m == m2->m);
    assert(this->n == m2->n);
    Matrix res = new Matrix_(this->m, this->n);
    for (int i = 0; i < m * n; i++) {
        res->data[i] = this->data[i] + m2->data[i];
    }
    return res;
}

Matrix Matrix_::operator-(double val) noexcept {
    Matrix res = new Matrix_(m, n);
    for (int i = 0; i < m * n; i++) {
        res->data[i] = this->data[i] - val;
    }
    return res;
}

Matrix Matrix_::operator-(const Matrix &m2) noexcept {
    assert(this->m == m2->m);
    assert(this->n == m2->n);
    Matrix res = new Matrix_(this->m, this->n);
    for (int i = 0; i < m * n; i++) {
        res->data[i] = this->data[i] - m2->data[i];
    }
    return res;
}

Vector Matrix_::operator*(const Vector &v2) noexcept {
    assert(v2->len == this->n);
    Vector res(v2);
    matmul(this->data, m, n, res->data);
    return res;
}

// x should be of size n
utils_inline void matmul(const double *a, const int m, const int n, double *&x) noexcept {
    double *y = new double[m] { 0.0 };
    cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, a, n, x, 1, 0.0, y, 1);
    delete [] x;
    // free(*x);
    x = y;
}

// Ku = v, solve for u, K is tridiagonal
utils_inline int solve_tridiag_system(double *KL, double *K, double *KU,
                                      const int n, double *v) noexcept {
    int ret = LAPACKE_dgtsv(LAPACK_ROW_MAJOR, n, 1, KL, K, KU, v, 1);
    return ret;
}

#ifdef TESTS

// Tests for functions in utils.c. These make sure that little details are
// implemented correctly!

#include <assert.h>
#include <math.h>

void test_row_major_set_get() {
    int N = 3;
    int M = 4;
    Matrix A (M, N);

    A.element(0, 1) = 1.0;
    A.element(1, 0) = 2.0;
    A.element(1, 1) = 3.0;
    assert(A.element(0, 1) == 1.0);
    assert(A.element(1, 0) == 2.0);

    assert(A.data[1] == 1.0);
    assert(A.data[3] == 2.0);
    assert(A.data[4] == 3.0);
}

void test_simple_matmul() {
    int N = 3;
    int M = 4;
    Matrix A (M, N);

    A.element(0, 0) = 2.0;
    A.element(0, 1) = 1.0;
    A.element(1, 1) = -1.0;

    Vector x_vec (3);
    delete [] x_vec.data;
    x_vec.data = new double[3] {5.0, 3.0, 1.0};

    Vector y = A * x_vec;

    assert(y[0] == 13.0);
    assert(y[1] == -3.0);
    assert(y[2] == 0.0);
    assert(y[3] == 0.0);
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
    double *KL = new double[3] {2, 2, 3};
    double *K = new double[4] {1, 1, 5, 1};
    double *KU = new double[3] {2, 2, 2};
    Vector x_vec (4);
    delete [] x_vec.data;
    x_vec.data = new double[4] {3.0, 14.0, -15.0, 9.0};

    Tridiag_Matrix K_mat (4);
    delete [] K_mat.al;
    delete [] K_mat.a;
    delete [] K_mat.au;
    K_mat.al = KL;
    K_mat.a = K;
    K_mat.au = KU;

    Vector y = K_mat * x_vec;

    assert(fabs(y[0] - 31.0) <= EPSILON);
    assert(fabs(y[1] - -10.0) <= EPSILON);
    assert(fabs(y[2] - -29.0) <= EPSILON);
    assert(fabs(y[3] - -36.0) <= EPSILON);
}

int main() {
    test_row_major_set_get();
    test_simple_matmul();
    test_tridiag_solve();
    test_tridiag_mul();
}

#endif
