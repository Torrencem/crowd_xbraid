#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <iostream>
#include <math.h>

using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;

// typedef SparseMatrix<double> Sparse;
typedef SparseMatrix<double> Sparse;
typedef VectorXd Vector;

Vector get_GwL(Vector &m, Vector &rho, Vector &lambda, Sparse &As, Sparse &At,
               Sparse &D1, Sparse &D2) {
    Vector GmL1 =
        2.0 * m.asDiagonal() * As.transpose() * At * rho.cwiseInverse();
    Vector GmL2 = D1.transpose() * lambda;

    Vector GmL = GmL1 + GmL2;
    Vector GrhoL =
        rho.cwiseProduct(rho).eval().cwiseInverse().eval().asDiagonal() *
            (-At.transpose()) * As * m.cwiseProduct(m) +
        D2.transpose() * lambda;
    Vector result(GmL.size() + GrhoL.size());
    for (int i = 0; i < GmL.size() + GrhoL.size(); i++) {
        if (i >= GmL.size()) {
            result[i] = GrhoL[i - GmL.size()];
        } else {
            result[i] = GmL[i];
        }
    }
    return result;
}

Vector get_GlambdaL(Vector &m, Vector &rho, Vector &q, Sparse &D) {
    Vector x(m.size() + rho.size());
    for (int i = 0; i < m.size() + rho.size(); i++) {
        if (i >= m.size()) {
            x[i] = rho[i - m.size()];
        } else {
            x[i] = m[i];
        }
    }
    return D * x - q;
}

/// In MATLAB, this would be diag(1./diag(a))
Sparse invertDiagonal(const Sparse &a) {
    Sparse m2(a.rows(), a.cols());
    for (int i = 0; i < a.rows(); i++) {
        if (a.coeff(i, i) != 0) {
            m2.insert(i, i) = 1.0 / a.coeff(i, i);
        }
    }
    return m2;
}

Sparse block_diag(const Sparse &diag, const int ncopies) {
    auto m2 =
        SparseMatrix<double>(diag.rows() * ncopies, diag.cols() * ncopies);

    for (int l = 0; l < ncopies; l++) {
        // For each non-zero entry in diag
        for (int k = 0; k < diag.outerSize(); k++)
            for (SparseMatrix<double>::InnerIterator it(diag, k); it; ++it) {
                m2.insert(l * diag.rows() + it.row(),
                          l * diag.cols() + it.col()) = it.value();
            }
    }

    return m2;
}

/// Compute the equivelent of [a; b] from matlab
Sparse jointb(const Sparse &a, const Sparse &b) {
    assert(a.cols() == b.cols());
    auto m2 = Sparse(a.rows() + b.rows(), a.cols());

    // For each non-zero entry in a
    for (int k = 0; k < a.outerSize(); k++)
        for (Sparse::InnerIterator it(a, k); it; ++it) {
            m2.insert(it.row(), it.col()) = it.value();
        }

    // For each non-zero entry in b
    for (int k = 0; k < b.outerSize(); k++)
        for (Sparse::InnerIterator it(b, k); it; ++it) {
            m2.insert(it.row() + a.rows(), it.col()) = it.value();
        }

    return m2;
}

/// Compute the equivelent of [a, b] from matlab
Sparse joinlr(const Sparse &a, const Sparse &b) {
    assert(a.rows() == b.rows());
    auto m2 = Sparse(a.rows(), a.cols() + b.cols());

    // For each non-zero entry in a
    for (int k = 0; k < a.outerSize(); k++)
        for (Sparse::InnerIterator it(a, k); it; ++it) {
            m2.insert(it.row(), it.col()) = it.value();
        }

    // For each non-zero entry in b
    for (int k = 0; k < b.outerSize(); k++)
        for (Sparse::InnerIterator it(b, k); it; ++it) {
            m2.insert(it.row(), it.col() + a.cols()) = it.value();
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
    Sparse As(mspace, mspace + 1);

    for (int i = 0; i < mspace; i++) {
        As.insert(i, i) = 0.5;
        As.insert(i, i + 1) = 0.5;
    }
    return block_diag(As, ntime);
}

Sparse get_At(int mspace, int ntime) {
    Sparse At(mspace * ntime, mspace * (ntime + 1));

    for (int i = 0; i < ntime; i++) {
        for (int j = 0; j < mspace; j++) {
            At.insert(mspace * i + j, mspace * i + j) = 0.5;
            At.insert(mspace * i + j, mspace * (i + 1) + j) = 0.5;
        }
    }

    return At;
}

Sparse get_derivative_matrix_space(int mspace, int ntime, double h) {
    Sparse top_bottom(mspace, (mspace + 1) * ntime);
    // TODO
    Sparse top_bottom2(mspace, (mspace + 1) * ntime);
    Sparse interior_block(mspace, mspace + 1);
    for (int i = 0; i < mspace; i++) {
        interior_block.insert(i, i) = -1.0;
        interior_block.insert(i, i + 1) = 1.0;
    }
    Sparse interior_cell = block_diag(interior_block, ntime);
    Sparse D = jointb(top_bottom, jointb(interior_cell, top_bottom2));
    D /= h;
    return D;
}

Sparse get_derivative_matrix_time(int mspace, int ntime, double h) {
    Sparse top(mspace, mspace * (ntime + 1));
    Sparse bottom(mspace, mspace * (ntime + 1));
    for (int i = 0; i < mspace; i++) {
        top.insert(i, i) = 1.0;
        bottom.insert(i, mspace * ntime + i) = -1.0;
    }

    Sparse center(mspace * ntime, mspace * (ntime + 1));
    for (int i = 0; i < ntime; i++) {
        for (int j = 0; j < mspace; j++) {
            center.insert(i * mspace + j, i * mspace + j) = -1.0;
            center.insert(i * mspace + j, (i + 1) * mspace + j) = 1.0;
        }
    }

    Sparse D = jointb(top, jointb(center, bottom));
    D /= h;

    return D;
}

Sparse get_A_hat(Vector &rho, Vector &m, Sparse &As, Sparse &At) {
    Vector vec_1 = 2 * As.transpose() * At * rho.cwiseInverse();
    Vector vec_2 = 2 * At.transpose() * As * (m.cwiseProduct(m));
    Vector vec_3 = (rho.cwiseProduct(rho.cwiseProduct(rho))).cwiseInverse();
    int dim = vec_1.rows() + vec_2.rows();
    Sparse A(dim, dim);
    for (int i = 0; i < vec_1.rows(); i++) {
        A.insert(i, i) = vec_1[i];
    }
    for (int i = 0; i < vec_2.rows(); i++) {
        A.insert(i + vec_1.rows(), i + vec_1.rows()) = vec_2[i] * vec_3[i];
    }
    return A;
}

Sparse get_zero_matrix(Sparse &S, Sparse &A_hat) {
    Sparse result(S.rows(), A_hat.cols());

    result.setZero();

    return result;
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

template <typename F> double minarg(F f, double a, double b, double tol) {
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

    return (a + b) / 2;
}

double reward(double x, Vector &dm, Vector &drho, Vector &dlambda, Vector &m,
              Vector &rho, Vector &lambda, Sparse &As, Sparse &At, Sparse &D1,
              Sparse &D2, Vector &q, Sparse &D) {
    Vector m2 = m + x * dm;
    Vector rho2 = rho + x * drho;
    Vector lambda2 = lambda + x * dlambda;
    Vector GwL = get_GwL(m2, rho2, lambda2, As, At, D1, D2);
    Vector GlambdaL = get_GlambdaL(m2, rho2, q, D);
    return GwL.squaredNorm() + GlambdaL.squaredNorm();
}

double line_search(Vector &dm, Vector &drho, Vector &dlambda, Vector &m,
                   Vector &rho, Vector &lambda, Sparse &As, Sparse &At,
                   Sparse &D1, Sparse &D2, Vector &q, Sparse &D) {
    auto f = [&](double x) {
        return reward(x, dm, drho, dlambda, m, rho, lambda, As, At, D1, D2, q,
                      D);
    };
    double tmp = minarg(f, -0.1, 0.1, 1e-8);
    return tmp;
}

void test_min() {
    auto square = [](double x) { return (x - 1.5) * (x - 1.5); };
    double min = minarg(square, -100, 100, 1e-8);
    std::cout << "Min is around: " << min << std::endl;
}
