#include <iostream>
#include <Eigen/Sparse>

using Eigen::SparseMatrix;
using Eigen::MatrixXd;
using Eigen::VectorXd;

typedef SparseMatrix<double> Sparse;
typedef VectorXd Vector;

Sparse block_diag(const Sparse &diag, const int ncopies) {
    auto m2 = SparseMatrix<double> (diag.rows() * ncopies, diag.cols() * ncopies);

    for (int l = 0; l < ncopies; l++) {
        // For each non-zero entry in diag
        for (int k = 0; k < diag.outerSize(); k++)
            for (SparseMatrix<double>::InnerIterator it(diag, k); it; ++it) {
                m2.insert(l * diag.rows() + it.row(), l * diag.cols() + it.col())
                    = it.value();
            }
    }

    return m2;
}

/// Compute the equivelent of [a; b] from matlab
Sparse jointb(const Sparse &a,
                            const Sparse &b) {
    assert(a.cols() == b.cols());
    auto m2 = Sparse (a.rows() + b.rows(), a.cols());
    
    // For each non-zero entry in a
    for (int k = 0; k < a.outerSize(); k++)
        for (Sparse::InnerIterator it(a, k); it; ++it) {
            m2.insert(it.row(), it.col())
                = it.value();
        }
    
    // For each non-zero entry in b
    for (int k = 0; k < b.outerSize(); k++)
        for (Sparse::InnerIterator it(b, k); it; ++it) {
            m2.insert(it.row() + a.rows(), it.col())
                = it.value();
        }

    return m2;
}

/// Compute the equivelent of [a, b] from matlab
Sparse joinlr(const Sparse &a,
                            const Sparse &b) {
    assert(a.rows() == b.rows());
    auto m2 = Sparse (a.rows(), a.cols() + b.cols());
    
    // For each non-zero entry in a
    for (int k = 0; k < a.outerSize(); k++)
        for (Sparse::InnerIterator it(a, k); it; ++it) {
            m2.insert(it.row(), it.col())
                = it.value();
        }
    
    // For each non-zero entry in b
    for (int k = 0; k < b.outerSize(); k++)
        for (Sparse::InnerIterator it(b, k); it; ++it) {
            m2.insert(it.row(), it.col() + a.cols())
                = it.value();
        }

    return m2;
}

void test_block_diag() {
    Sparse m(3, 3);

    m.insert(0, 0) = 1;
    m.insert(1, 2) = 2;
    m.insert(1, 0) = 3;
    m.insert(2, 0) = 4;
    std::cout << "M is:" << std::endl;
    std::cout << m << std::endl;
    
    auto m2 = block_diag(m, 3);
    
    std::cout << std::endl << "block_diag(M, 3) is:" << std::endl;
    std::cout << m2 << std::endl;
}

Sparse get_As(int mspace, int ntime) {
    Sparse As(mspace, mspace);
    As.insert(0, 0) = 1.0;

    for (int i = 1; i < mspace; i++) {
        As.insert(i, i) = 0.5;
        As.insert(i, i - 1) = 0.5;
    }
    return block_diag(As, ntime);
}

Sparse get_At(int mspace, int ntime) {
    Sparse At(mspace * ntime, mspace * ntime);

    for (int i = 0; i < mspace; i++) {
        At.insert(i, i) = 0.5;
    }

    for (int i = 1; i < ntime; i++) {
        for (int j = 0; j < mspace; j++) {
            At.insert(mspace * (i - 1) + j, mspace * (i - 1) + j) = 0.5;
            At.insert(mspace * (i - 1) + j, mspace * (i - 2) + j) = 0.5;
        }
    }

    return At;
}

void calc_fixed_matrices(int mspace,
        int ntime,
        Sparse &D,
        Sparse &D1,
        Sparse &D2,
        Sparse &As,
        Sparse &At) {
    // TODO: D1 and D2 matrices
    D = joinlr(D1, D2);
    As = get_As(mspace, ntime);
    At = get_At(mspace, ntime);
}

void test_join() {
    Sparse m(3, 5);

    m.insert(0, 0) = 1;
    m.insert(1, 2) = 2;
    m.insert(1, 3) = 3;
    m.insert(2, 4) = 4;
    
    Sparse m2(3, 5);

    m2.insert(0, 0) = 1;
    m2.insert(1, 0) = 2;
    m2.insert(1, 1) = 3;
    m2.insert(2, 2) = 4;

    std::cout << "M is:" << std::endl;
    std::cout << m << std::endl;

    std::cout << "M2 is:" << std::endl;
    std::cout << m2 << std::endl;

    std::cout << "[M; M2] is:" << std::endl;
    std::cout << jointb(m, m2) << std::endl;
    
    std::cout << "[M, M2] is:" << std::endl;
    std::cout << joinlr(m, m2) << std::endl;
}

double GR = 1.61803398875;

template <typename F>
double minarg(F f, double a, double b, double tol = 1e-5) {
    double c = b - (b - a) / GR;
    double d = a + (b - a) / GR;
    while (abs(c - d) > tol) {
        if (f(c) < f(d)) {
            b = d;
        } else {
            a = c;
        }

        c = b - (b - a) / GR;
        d = a + (b - a) / GR;
    }

    return (b + a) / 2;
}

void test_min() {
    auto square = [](double x) { return (x - 1.5) * (x - 1.5); };
    double min = minarg(square, -10.0, 10.0, 1e-6);
    std::cout << "Min is around: " << min << std::endl;
}

int main() {
    // test_block_diag();
    // test_join();
    test_min();
    exit(0);
    int mspace = 20;
    int ntime = 12;
    
    int rho_size = mspace * ntime;
    int lambda_size = mspace * (ntime + 1);

    Vector m(rho_size);
    Vector rho(rho_size);
    Vector lambda(lambda_size);
    m.setConstant(1.1);
    rho.setConstant(1.1);
    lambda.setConstant(1.1);
    Sparse q(lambda_size, 1);

    int time = 1;
    int iters = 20;

    double h = time / ntime;

    for (int i = 0; i < mspace / 2; i++) {
        q.insert(i, 0) = 1.0;
    }
    for (int i = mspace / 2; i < mspace; i++) {
        q.insert(i, 0) = 0.1;
    }
    for (int i = mspace * (ntime - 1); i < mspace * (ntime - 1) + mspace / 2; i++) {
        q.insert(i, 0) = 0.1;
    }
    for (int i = mspace * (ntime - 1) + mspace / 2; i < mspace * ntime; i++) {
        q.insert(i, 0) = 1.0;
    }

    q /= h;
    
    Sparse D, D1, D2, As, At;

    calc_fixed_matrices(mspace, ntime, D, D1, D2, As, At);
}
