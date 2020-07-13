#include <iostream>
#include <Eigen/Sparse>

using Eigen::SparseMatrix;
using Eigen::MatrixXd;

SparseMatrix<double> block_diag(const SparseMatrix<double> &diag, const int ncopies) {
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
SparseMatrix<double> join(const SparseMatrix<double> &a, const SparseMatrix<double> &b) {
    assert(a.cols() == b.cols());
    auto m2 = SparseMatrix<double> (a.rows() + b.rows(), a.cols());
    
    // For each non-zero entry in a
    for (int k = 0; k < a.outerSize(); k++)
        for (SparseMatrix<double>::InnerIterator it(a, k); it; ++it) {
            m2.insert(it.row(), it.col())
                = it.value();
        }
    
    // For each non-zero entry in b
    for (int k = 0; k < b.outerSize(); k++)
        for (SparseMatrix<double>::InnerIterator it(b, k); it; ++it) {
            m2.insert(it.row() + a.rows(), it.col())
                = it.value();
        }

    return m2;
}

void test_block_diag() {
    SparseMatrix<double> m(3, 3);

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

void test_join() {
    SparseMatrix<double> m(3, 5);

    m.insert(0, 0) = 1;
    m.insert(1, 2) = 2;
    m.insert(1, 3) = 3;
    m.insert(2, 4) = 4;
    
    SparseMatrix<double> m2(3, 5);

    m2.insert(0, 0) = 1;
    m2.insert(1, 0) = 2;
    m2.insert(1, 1) = 3;
    m2.insert(2, 2) = 4;

    std::cout << "M is:" << std::endl;
    std::cout << m << std::endl;

    std::cout << "M2 is:" << std::endl;
    std::cout << m2 << std::endl;

    std::cout << "[M; M2] is:" << std::endl;
    std::cout << join(m, m2) << std::endl;
}

int main() {
    test_block_diag();
    test_join();
}
