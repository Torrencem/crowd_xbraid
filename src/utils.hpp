
#ifndef UTILS_HPP
#define UTILS_HPP

#ifdef NOINLINE
#define utils_inline
#else
#define utils_inline extern inline
#endif

class Vector {
public:
    double *data;
    int len;
    Vector (const int len_);
    Vector (const Vector &v2);
    ~Vector();
    double &operator[](int i);
    const Vector operator*(double val);
    const Vector operator+(double val);
    const Vector operator+(const Vector &m2);
};

class Tridiag_Matrix {
public:
    double *al, *a, *au;
    int n;

    Tridiag_Matrix(int n_);
    Tridiag_Matrix(const Tridiag_Matrix &m2);
    ~Tridiag_Matrix();

    Vector operator*(const Vector &x);
};

class Matrix {
public:
    // The raw data of the matrix, stored in row major order
    double *data;
    int m, n;

    Matrix(const int m_, const int n_);
    Matrix(const Matrix &m2);
    Matrix(const Tridiag_Matrix &m2);
    ~Matrix();

    double &element(const int row, const int col);
    
    Matrix operator*(double val);

    Matrix operator+(double val);

    Matrix operator+(const Matrix &m2);

    Vector operator*(const Vector &v2);
};

utils_inline void matmul(const double *a, const int m, const int n, double *&x);
utils_inline int solve_tridiag_system(double *KL, double *K, double *KU, const int n, double *v);

#endif
