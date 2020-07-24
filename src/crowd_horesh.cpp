#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

using Eigen::SparseMatrix;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// typedef SparseMatrix<double> Sparse;
typedef SparseMatrix<double> Sparse;
typedef VectorXd Vector;

Sparse inverse(const Sparse &A) {
    Eigen::SparseLU<Sparse> solver;
    solver.compute(A);
    Sparse I (A.rows(), A.cols());
    I.setIdentity();
    return solver.solve(I);
}

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

// /// Lay a vector out on the main diagonal of a square matrix
// Sparse diag(const Vector &v) {
//     Sparse result (v.size(), v.size());
//
//     for (int i = 0; i < v.size(); i++) {
//         result.insert(i, i) = v[i];
//     }
//
//     return result;
// }

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
    Sparse top_bottom (mspace, (mspace + 1) * ntime);
    // TODO
    Sparse top_bottom2 (mspace, (mspace + 1) * ntime);
    Sparse interior_block (mspace, mspace + 1);
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
    Sparse top (mspace, mspace * (ntime + 1));
    Sparse bottom (mspace, mspace * (ntime + 1));
    for (int i = 0; i < mspace; i++) {
        top.insert(i, i) = 1.0;
        bottom.insert(i, mspace * ntime + i) = -1.0;
    }

    Sparse center (mspace * ntime, mspace * (ntime + 1));
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

void calc_fixed_matrices(int mspace,
        int ntime,
        double h,
        Sparse &D,
        Sparse &D1,
        Sparse &D2,
        Sparse &As,
        Sparse &At) {
    D1 = get_derivative_matrix_space(mspace, ntime, h);
    D2 = get_derivative_matrix_time(mspace, ntime, h);
    D = joinlr(D1, D2);
    As = get_As(mspace, ntime);
    At = get_At(mspace, ntime);
}

Sparse get_A_hat(Vector &rho, Vector &m, Sparse &As, Sparse &At) {
    Vector vec_1 = 2 * As.transpose() * At * rho.cwiseInverse();
    Vector vec_2 = 2 * At.transpose() * As * (m.cwiseProduct(m));
    Vector vec_3 = (rho.cwiseProduct(rho.cwiseProduct(rho))).cwiseInverse();
    int dim = vec_1.rows() + vec_2.rows();
    Sparse A (dim, dim);
    for (int i = 0; i < vec_1.rows(); i++) {
        A.insert(i, i) = vec_1[i];
    }
    for (int i = 0; i < vec_2.rows(); i++) {
        A.insert(i + vec_1.rows(), i + vec_1.rows()) = vec_2[i] * vec_3[i];
    }
    return A;
}

Vector get_GwL(Vector &m,
        Vector &rho,
        Vector &lambda,
        Sparse &As,
        Sparse &At,
        Sparse &D1,
        Sparse &D2) {
    Vector GmL1 = 2.0 * m.asDiagonal() * As.transpose() * At * rho.cwiseInverse();
    Vector GmL2 = D1.transpose() * lambda;

    // std::cout << "As: " << As << std::endl;
    // exit(0);
    
    // 
    Vector GmL = GmL1 + GmL2;
    // Vector GmRho =
    //     rho.cwiseInverse()
    //        .array()
    //        .square()
    //        .eval()
    //        .matrix()
    //        .asDiagonal()
    //        * (- At.transpose()) * As * m.cwiseProduct(m) + D2.transpose() * lambda;
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
    Vector x (m.size() + rho.size());
    for (int i = 0; i < m.size() + rho.size(); i++) {
        if (i >= m.size()) {
            x[i] = rho[i - m.size()];
        } else {
            x[i] = m[i];
        }
    }
    return D * x - q;
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

template <typename F>
double minarg(F f, double a, double b, double tol) {
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

double reward(double x, Vector &dm, Vector &drho, Vector &dlambda, Vector &m, Vector &rho, Vector &lambda, Sparse &As, Sparse &At, Sparse &D1, Sparse &D2, Vector &q, Sparse &D) {
    Vector m2 = m + x * dm;
    Vector rho2 = rho + x * drho;
    Vector lambda2 = lambda + x * dlambda;
    Vector GwL = get_GwL(m2, rho2, lambda2, As, At, D1, D2);
    Vector GlambdaL = get_GlambdaL(m2, rho2, q, D);
    return GwL.squaredNorm() + GlambdaL.squaredNorm();
}

double line_search(Vector &dm, Vector &drho, Vector &dlambda, Vector &m, Vector &rho, Vector &lambda, Sparse &As, Sparse &At, Sparse &D1, Sparse &D2, Vector &q, Sparse &D) {
    auto f = [&](double x) {
        return reward(x, dm, drho, dlambda, m, rho, lambda, As, At, D1, D2, q, D);
    };
    double tmp = minarg(f, -0.1, 0.1, 1e-8);
    // printf("READ THIS: %f\n", tmp);
    return tmp;
}

void test_min() {
    auto square = [](double x) { return (x - 1.5) * (x - 1.5); };
    double min = minarg(square, -100, 100, 1e-8);
    std::cout << "Min is around: " << min << std::endl;
}

int main() {
    // test_block_diag();
    // test_join();
    // test_min();
    // exit(0);
    int mspace = 45;

    // assert(mspace % 2 == 0);

    int ntime = 45;

    Vector m((mspace + 1) * ntime);
    Vector rho(mspace * (ntime + 1));
    Vector lambda(mspace * (ntime + 2));
    m.setConstant(0.1);
    rho.setConstant(0.5);
    lambda.setConstant(0.1);
    Vector q(mspace * (ntime + 2), 1);

    double time = 1.0;
    int iters = 15;

    double h = time / ntime;

    // for (int i = 0; i < mspace / 2; i++) {
    //     q.insert(i, 0) = 0.1;
    // }
    // for (int i = mspace / 2; i < mspace; i++) {
    //     q.insert(i, 0) = 1.0;
    // }
    // for (int i = mspace * (ntime + 1); i < mspace / 2 * (ntime * 2 + 3); i++) {
    //     q.insert(i, 0) = 1.0;
    // }
    // for (int i = mspace / 2 * (ntime * 2 + 3); i < mspace * (ntime + 2); i++) {
    //     q.insert(i, 0) = 0.1;
    // }

    for (int i = 0; i < mspace; i++) {
        q[i] = 0.5;
    }
    double acc = 0.0;
    for (int i = mspace * (ntime + 1); i < mspace * (ntime + 2); i++) {
        q[i] = -0.5 + -sin(acc * 2.0 * 3.141593) * 0.5;
        acc += 1.0 / ((double) mspace - 1.0);
    }

    q /= h;
    
    Sparse D, D1, D2, As, At;

    calc_fixed_matrices(mspace, ntime, h, D, D1, D2, As, At);

    for (int i = 0; i < iters; i++) {
        Sparse A_hat = get_A_hat(rho, m, As, At);
        // A_hat is diagonal so this is the same as inverse(A_hat)
        Sparse S = -D * A_hat.cwiseInverse() * D.transpose();

        Sparse A = jointb(joinlr(A_hat, D.transpose()), joinlr(get_zero_matrix(S, A_hat), S));
        int dim_D_1 = D.rows();
        int dim_D_2 = D.cols();
        Sparse tl(dim_D_2, dim_D_2);
        tl.setIdentity();
        Sparse br(dim_D_1, dim_D_1);
        br.setIdentity();
        Sparse tr(dim_D_2, dim_D_1);
        tr.setZero();
        Sparse bl = -D * A_hat.cwiseInverse();

        Sparse preconditioner = jointb(joinlr(tl, tr), joinlr(bl, br));
        
        Vector GwL = get_GwL(m, rho, lambda, As, At, D1, D2);
        Vector GlambdaL = get_GlambdaL(m, rho, q, D);
        Vector b(GwL.size() + GlambdaL.size());
        for (int i = 0; i < GwL.size() + GlambdaL.size(); i++) {
            if (i >= GwL.size()) {
                b[i] = GlambdaL[i - GwL.size()];
            } else {
                b[i] = GwL[i];
            }
        }
        b = -preconditioner * b;
        // std::cout << "b is: " << b << std::endl;
        // exit(0);
        // Check convergence
        std::cout << "norm(GwL): " << GwL.norm() / (mspace * ntime) << std::endl;
        std::cout << "norm(GlambdaL): " << GlambdaL.norm() / (mspace * ntime) << std::endl;
        std::cout << "Norm: " << b.norm() / (mspace * ntime) << std::endl;
        
        // Solve A(solution) = b
        Eigen::BiCGSTAB<Sparse, Eigen::IncompleteLUT<double>> solver;
        // Eigen::BiCGSTAB<Sparse> solver;
        // solver.preconditioner().setDroptol(.001);
        // solver.setMaxIterations(100000);
        solver.compute(A);
        if (solver.info() != Eigen::Success) {
            printf("Error decomposing A!\n");
            exit(1);
        }
        Vector solution = solver.solve(b);
        if (solver.info() != Eigen::Success) {
            printf("Error solving Ax = b!\n");
            exit(1);
        }
        // Vector solution = inverse(A) * b;
        
        Eigen::Map<Vector> dm(solution.data(), (mspace + 1) * ntime);
        int start = (mspace + 1) * ntime;
        int len = mspace * (ntime + 1);
        Eigen::Map<Vector> drho(solution.data() + start, len);
        start += len;
        len = solution.size() - start;
        Eigen::Map<Vector> dlambda(solution.data() + start, len);

        Vector dm_(dm);
        Vector drho_(drho);
        Vector dlambda_(dlambda);

        double alpha = line_search(dm_, drho_, dlambda_, m, rho, lambda, As, At, D1, D2, q, D);
        
        std::cout << "drho norm: " << drho.norm() << std::endl;
        std::cout << "alpha: " << alpha << std::endl;

        m += alpha * dm;
        rho += alpha * drho;
        lambda += alpha * dlambda;
    }

    for (int x = 0; x < ntime; x++) {
        for (int y = 0; y < mspace; y++) {
            std::cout 
                << x
                << ", " 
                << y
                << ": " 
                << rho[x * mspace + y]
                << std::endl;
        }
    }
}
