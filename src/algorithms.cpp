#include "algorithms.h"
#include "error_utils.h"
#include "hilbert.h"
#include "matrix.h"

typedef std::complex<double> cd;

std::vector<Matrix> qrDecompose(Matrix A) {
        size_t k = 0;
        size_t pivot = 0;

        HilbertSpace space;
        StateVector state({cd{0, 0}});
        Matrix H({{cd{0, 0}}});
        Matrix Q(A.identity());

        for (size_t i = 0; i < A.getRows(); ++i) {
                if (isZeroColumn(A, i)) {
                        continue;
                }

                k = A.getRows() - pivot;
                state.resize(k);
                H.resize(k, k);

                for (size_t j = 0; j < k; ++j) {
                        state[j] = A[j + pivot][i];
                }

                if (std::abs(std::real(state[0])) < 1e-8 && std::abs(std::imag(state[0])) < 1e-8) {
                        state[0] += cd{space.norm(state), 0};
                } else {
                        state[0] += (state[0] / std::abs(state[0])) * cd{space.norm(state), 0};
                }

                H = H.identity() + ((space.outerProduct(state, state.conjugate()) * cd{-2, 0}) * (cd{1, 0} / space.innerProduct(state, state)));
                H = identityPad(H, A.getRows());

                A = H * A;
                Q = Q * H;
        }

        for (size_t i = 0; i < A.getColumns(); ++i) {
                size_t pivot_row = 0;
                double max_val = 0.0;

                if (isZeroColumn(A, i)) {
                        for (size_t j = 0; j < Q.getRows(); ++j) {
                                if (std::abs(Q[j][i]) > max_val) {
                                        max_val = std::abs(Q[j][i]);
                                        pivot_row = j;
                                }
                        }

                        if (max_val > 0 && std::real(Q[pivot_row][i]) < 0) {
                                for (size_t j = 0; j < A.getColumns(); ++j) {
                                        A[pivot_row][j] *= cd{-1, 0};
                                }

                                for (size_t j = 0; j < Q.getRows(); ++j) {
                                        Q[j][pivot_row] *= cd{-1, 0};
                                }
                        }
                } else {
                        for (size_t j = 0; j < A.getRows(); ++j) {
                                if (std::abs(A[j][i]) > max_val) {
                                        max_val = std::abs(A[j][i]);
                                        pivot_row = j;
                                }
                        }

                        if (max_val > 0 && std::real(A[pivot_row][i]) < 0) {
                                for (size_t j = 0; j < A.getColumns(); ++j) {
                                        A[pivot_row][j] *= cd{-1, 0};
                                }

                                for (size_t j = 0; j < Q.getRows(); ++j) {
                                        Q[j][pivot_row] *= cd{-1, 0};
                                }
                        }
                }
        }

        return std::vector<Matrix>{Q, A};
}

bool offDiagonalValues(Matrix A, std::complex<double> tolerance) {
        std::complex<double> sum = {0, 0};

        for (size_t i = 0; i < A.getRows(); ++i) {
                for (size_t j = 0; j < A.getColumns(); ++j) {
                        if (i != j) {
                                sum += A[i][j];
                        }
                }
        }

        if (std::abs(std::real(sum)) < std::real(tolerance) && std::abs(std::imag(sum)) < std::imag(tolerance)) {
                return true;
        }

        return false;
}

Matrix algorithmQR(Matrix A) {
        size_t maxIterations = 500;
        std::complex<double> tolerance = {1e-10, 1e-10};
        Matrix accumulatedQ(A.identity());
        Matrix copyA(A);

        for (size_t i = 0; i < maxIterations; ++i) {
                if (offDiagonalValues(copyA, tolerance)) {
                        if (i == 0) {
                                accumulatedQ = accumulatedQ * qrDecompose(copyA)[0];
                        }

                        break;
                }


                std::vector<Matrix> QR;
                QR = qrDecompose(copyA);

                Matrix Q = QR[0];
                Matrix R = QR[1];

                copyA = R * Q;

                accumulatedQ = accumulatedQ * Q;
        }

        return accumulatedQ;
}

std::vector<Matrix> algorithmSVD(StateVector state) {
        HilbertSpace space;
        Matrix A(split(state));
        Matrix B(A * A.conjugate().transpose());

        std::complex<double> B_trace = B[0][0] + B[1][1];
        std::complex<double> B_determinant = B[0][0] * B[1][1] - B[0][1] * std::conj(B[0][1]);
        std::complex<double> eigenvalue1 = (B_trace + std::sqrt(pow(B_trace, 2) - cd{4, 0} * B_determinant)) / cd{2, 0};
        std::complex<double> eigenvalue2 = (B_trace - std::sqrt(pow(B_trace, 2) - cd{4, 0} * B_determinant)) / cd{2, 0};

        StateVector eigenvector1(2);
        StateVector eigenvector2(2);

        Matrix U(2,2);

        if (std::abs(B[0][1]) < 1e-10) {
                U = U.identity();
        } else if (std::abs(B[0][1]) > std::abs(eigenvalue1 - B[0][0])) {
                eigenvector1[0] = B[0][1];
                eigenvector1[1] = eigenvalue1 - B[0][0];

                eigenvector1 = space.normalize(eigenvector1);

                eigenvector2[0] = cd{-1, 0} * std::conj(eigenvalue1 - B[0][0]);
                eigenvector2[1] = std::conj(B[0][1]);

                eigenvector2 = space.normalize(eigenvector2);

                U[0][0] = eigenvector1[0];
                U[0][1] = eigenvector2[0];
                U[1][0] = eigenvector1[1];
                U[1][1] = eigenvector2[1];
        } else {
                eigenvector1[0] = eigenvalue1 - B[1][1];
                eigenvector1[1] = std::conj(B[0][1]);

                eigenvector1 = space.normalize(eigenvector1);

                eigenvector2[0] = cd{-1, 0} * std::conj(B[0][1]);
                eigenvector2[1] = std::conj(eigenvalue1 - B[1][1]);

                eigenvector2 = space.normalize(eigenvector2);

                U[0][0] = eigenvector1[0];
                U[0][1] = eigenvector2[0];
                U[1][0] = eigenvector1[1];
                U[1][1] = eigenvector2[1];
        }

        if (std::abs(eigenvalue1) < std::abs(eigenvalue2)) {
                std::swap(eigenvector1, eigenvector2);
        }

        std::complex<double> singularValue1 = std::sqrt(eigenvalue1);
        std::complex<double> singularValue2 = std::sqrt(eigenvalue2);

        Matrix sigma(A);
        size_t rowCheck = sigma.getRows();
        size_t columnCheck = sigma.getColumns();

        for (size_t i = 0; i < rowCheck; ++i) {
                for (size_t j = 0; j < columnCheck; ++j) {
                        sigma[i][j] = {0, 0};
                }
        }

        if (columnCheck > 0) {
                sigma[0][0] = singularValue1;
        }

        if (rowCheck > 1 && columnCheck > 1) {
                sigma[1][1] = singularValue2;
        }

        Matrix V(algorithmQR(A.conjugate().transpose() * A).conjugate().transpose());

        return std::vector<Matrix>{U, sigma, V};
}
