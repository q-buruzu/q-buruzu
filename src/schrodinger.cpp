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
	Matrix resultMatrix(2, vector.size() / 2);

	for (size_t i = 0; i < vector.size() / 2; ++i) {
		resultMatrix[0][i] = vector[i];
		resultMatrix[1][i] = vector[i + (vector.size() / 2)];
	}

	return resultMatrix;
}

Matrix outerProduct(StateVector vector) {
        Matrix resultMatrix(vector.size(), vector.size());

        for (size_t i = 0; i < vector.size(); ++i) {
                for (size_t j = 0; j < vector.size(); ++j) {
                        resultMatrix[i][j] = vector[i] * std::conj(vector[j]);
                }
        }

        return resultMatrix;
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

	return std::vector<Matrix>{Q, A};
}

bool offDiagonalValues(Matrix A, std::complex<double> tolerance) {
        std::complex<double> sum{0, 0};

        for (size_t i = 0; i < A.getRows(); ++i) {
                for (size_t j = 0; j < A.getColumns(); ++j) {
                        if (i != j) {
                                sum += A[i][j];
                        }
                }
        }

	if (std::real(sum) < std::real(tolerance) && std::imag(sum) < std::imag(tolerance)) {
                return true;
        }

        return false;
}

Matrix algorithmQR(Matrix A) {
        size_t maxIterations{500};
	std::complex<double> tolerance = {1e-10, 1e-10};
        Matrix accumulatedQ(A.identity());

        Matrix copyA(A);

        for (size_t i = 0; i < maxIterations; ++i) {
		if (offDiagonalValues(copyA, tolerance)) {
			break;
		}

                std::vector<Matrix> QR_i;
		QR_i = qrDecompose(copyA);

                Matrix Q_i = QR_i[0];
                Matrix R_i = QR_i[1];

                copyA = R_i * Q_i;

		accumulatedQ = accumulatedQ * Q_i;
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

	Matrix U(2, 2);
	U[0][0] = eigenvector1[0];
        U[0][1] = eigenvector2[0];
        U[1][0] = eigenvector1[1];
        U[1][1] = eigenvector2[1];

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

std::vector<StateVector> extractQubitStates(StateVector state) {
	size_t nQubits = static_cast<size_t>(std::log2(state.size()));
	std::vector<StateVector> individualQubits;
	individualQubits.reserve(nQubits);

	StateVector copyState(state);

	for (size_t i = 0; i < nQubits; ++i) {
		std::vector<Matrix> USV;
		USV = algorithmSVD(copyState);
		Matrix& U(USV[0]);
		Matrix& V(USV[2]);

		StateVector firstQubitState(2);
		firstQubitState[0] = U[0][0];
		firstQubitState[1] = U[1][0];
		individualQubits.push_back(firstQubitState);

		size_t width = V.getRows();
		StateVector remainingQubitStates(width);

		for (size_t j = 0; j < width; ++j) {
			remainingQubitStates[j] = V[j][0];
		}

		copyState.set(remainingQubitStates.get());
	}

	return individualQubits;
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
