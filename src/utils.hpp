
#ifndef UTILS_HPP
#define UTILS_HPP

#ifdef NOINLINE
#define utils_inline
#else
#define utils_inline extern inline
#endif

class Vector_;

typedef Vector_ *Vector;

class Vector_ {
public:
    double *data;
    int len;
    Vector_ (const int len_) noexcept;
    Vector_ (const Vector_ &v2) noexcept;
    ~Vector_() noexcept;
    double &operator[](int i) noexcept;
    Vector operator*(double val) noexcept;
    Vector operator+(double val) noexcept;
    Vector operator+(const Vector &m2) noexcept;
    Vector operator-(double val) noexcept;
    Vector operator-(const Vector &m2) noexcept;
};

class Tridiag_Matrix_;

typedef Tridiag_Matrix_ *Tridiag_Matrix;

class Tridiag_Matrix_ {
public:
    double *al, *a, *au;
    int n;

    Tridiag_Matrix_(int n_) noexcept;
    Tridiag_Matrix_(const Tridiag_Matrix &m2) noexcept;
    ~Tridiag_Matrix_() noexcept;

    Vector operator*(const Vector &x) noexcept;
};

class Matrix_;

typedef Matrix_ *Matrix;

class Matrix_ {
public:
    // The raw data of the matrix, stored in row major order
    double *data;
    int m, n;

    Matrix_(const int m_, const int n_) noexcept;
    Matrix_(const Matrix &m2) noexcept;
    Matrix_(const Tridiag_Matrix &m2) noexcept;
    ~Matrix_() noexcept;

    double &element(const int row, const int col) noexcept;
    
    Matrix operator*(double val) noexcept;

    Matrix operator+(double val) noexcept;

    Matrix operator+(const Matrix &m2) noexcept;
    
    Matrix operator-(double val) noexcept;

    Matrix operator-(const Matrix &m2) noexcept;

    Vector operator*(const Vector &v2) noexcept;
};

utils_inline void matmul(const double *a, const int m, const int n, double *&x) noexcept;
utils_inline int solve_tridiag_system(double *KL, double *K, double *KU, const int n, double *v) noexcept;

#endif
