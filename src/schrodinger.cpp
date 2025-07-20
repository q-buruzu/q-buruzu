#include "hilbert.h"
#include "matrix.h"
#include "quantumstate.h"
#include "schrodinger.h"

#include <complex>
#include <cmath>
#include <iostream>

typedef std::complex<double> cd;

double hBar = 1.054571817e-34;

Matrix split(StateVector state) {
	std::vector<std::complex<double>> vector = state.get();
	Matrix resultMatrix(2, 0.5 * vector.size());

	for (size_t i = 0; i < 0.5 * vector.size(); ++i) {
		resultMatrix[0][i] = vector[i];
	}

	for (size_t i = 0.5 * vector.size(); i < vector.size(); ++i) {
		resultMatrix[1][i] = vector[i];
	}

	return resultMatrix;
}

Matrix outerProduct(StateVector state) {
        StateVector conjugateVector(state.size());

        for (size_t i = 0; i < state.size(); ++i) {
                conjugateVector[i] = std::conj(state[i]);
        }

        Matrix outerProduct(state.size(), state.size());

        for (size_t i = 0; i < state.size(); ++i) {
                for (size_t j = 0; j < state.size(); ++j) {
                        outerProduct[i][j] = state[i] * conjugateVector[j];
                }
        }

        return outerProduct;
}

Matrix identityPad(Matrix A, size_t size) {
	if (A.getRows() == size) {
		return A;
	}

	Matrix B(size, size);
	B = B.identity();

	for (size_t i = 0; i < A.getRows(); ++i) {
		for (size_t j = 0; j < A.getColumns(); ++j) {
			B[i + (size - A.getRows())][j + (size - A.getColumns())] = A[i][j];
		}
	}

	return B;
}

std::vector<Matrix> qrDecompose(Matrix A) {
	size_t k;
	StateVector state({cd{0, 0}});
	HilbertSpace space;

	Matrix H({{cd{0, 0}}});
	Matrix Q(A.getRows(), A.getRows());
	Q = Q.identity();

	for (size_t i = 0; i < A.getRows() - 1; ++i) {
		k = A.getRows() - i;
		state.resize(k);
		H.resize(k, k);

		for (size_t j = 0; j < k; ++j) {
			state[j] = A[j + i][i];
		}

		state[0] += (state[0] / std::abs(state[0])) * space.norm(state);

		H = H.identity() + ((outerProduct(state) * cd{-2, 0}) * (cd{1, 0} / space.innerProduct(state, state)));

		H = identityPad(H, A.getRows());

		A = H * A;
		Q = Q * H;
	}

	Q.roundValues();
	A.roundValues();

	return std::vector<Matrix>{Q, A};
}

void SVD(StateVector state) {
        Matrix A(split(state));
        Matrix B(A * A.conjugate().transpose());

        std::complex<double> B_trace = B[0][0] + B[1][1];
        std::complex<double> B_determinant = B[0][0] * B[1][1] - B[0][1] * std::conj(B[0][1]);
        std::complex<double> eigenvalue1 = (B_trace + std::sqrt(pow(B_trace, 2) - cd{4, 0} * B_determinant)) / cd{2, 0};
        std::complex<double> eigenvalue2 = (B_trace - std::sqrt(pow(B_trace, 2) - cd{4, 0} * B_determinant)) / cd{2, 0};

        StateVector eigenvector1(2);
        StateVector eigenvector2(2);
	HilbertSpace space;

        if (std::abs(B[0][1]) > std::abs(eigenvalue1 - B[0][0])) {
                eigenvector1[0] = B[0][1];
                eigenvector1[1] = eigenvalue1 - B[0][0];

                space.normalize(eigenvector1);

                eigenvector2[0] = cd{-1, 0} * std::conj(eigenvalue1 - B[0][0]);
                eigenvector2[1] = std::conj(B[0][1]);
        } else {
                eigenvector1[0] = eigenvalue1 - B[1][1];
                eigenvector1[1] = std::conj(B[0][1]);

		space.normalize(eigenvector1);

		eigenvector2[0] = cd{-1, 0} * std::conj(B[0][1]);
                eigenvector2[1] = std::conj(eigenvalue1 - B[1][1]);
        }

        std::vector<std::vector<std::complex<double>>> U(2, std::vector<std::complex<double>>(2, cd{0, 0}));

 	U[0][0] = eigenvector1.get()[0];
        U[0][1] = eigenvector2.get()[0];
        U[1][0] = eigenvector1.get()[1];
        U[1][1] = eigenvector2.get()[1];

        std::complex<double> singularValue1 = std::sqrt(eigenvalue1);
        std::complex<double> singularValue2 = std::sqrt(eigenvalue2);

        Matrix sigma(A);

        for (size_t i = 0; i < 2; ++i) {
                for (size_t j = 0; j < sigma.getColumns(); ++j) {
                        sigma[i][j] = {0, 0};
                }
        }

        sigma[0][0] = singularValue1;
        sigma[1][1] = singularValue2;

        Matrix C(A.conjugate().transpose() * A);

        // Matrix V();
}

Matrix findEvolution(double step) {
	Matrix O({{cd{1, 0}, cd{0, 1}}, {cd{0, 1}, cd{1, 0}}});

	O = (((O * cd{step, 0}) * cd{0, 1}) * cd{1 / hBar, 0}) * cd{-1, 0};

	O = O.identity() + O;

	return O;
}

void evolve(QuantumState state, Matrix evolution) {
	state.set(evolution * state.get());
}

