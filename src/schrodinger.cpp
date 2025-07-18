#include "matrix.h"
#include "quantumstate.h"
#include "schrodinger.h"

#include <complex>
#include <cmath>
#include <iostream>

typedef std::complex<double> cd;

double hBar = 1.054571817e-34;
Matrix I({{cd{1, 0}, cd{0, 0}}, {cd{0, 0}, cd{1, 0}}});

Matrix resize(StateVector state) {
	std::vector<std::complex<double>> vector = state.get();
	Matrix resultMatrix(2, 0.5 * vector.size());

	for (int i = 0; i < 0.5 * vector.size(); ++i) {
		resultMatrix[0][i] = vector[i];
	}

	for (int i = 0.5 * vector.size(); i < vector.size(); ++i) {
		resultMatrix[1][i] = vector[i];
	}

	return resultMatrix;
}

void tempNormalize(StateVector vector) {
	HilbertSpace space;

        double factor = 1 / space.norm(vector);

        for (size_t i = 0; i < vector.get().size(); ++i) {
                vector[i] *= {factor, 0};
        }
}

double tempNorm(StateVector vector) {
	HilbertSpace space;

	return space.norm(vector);
}

void qrDecompose(Matrix A) {
	int index;
	StateVector state({{0, 0}});
	double norm;

	for (int i = 0; i < A.getColumns; ++i) {
		index = A.getRows - i;
		state.resize(index);

		for (int j = 0; j < index; ++j) {
			state[i] = A[j][i];
		}

		norm = tempNorm(state);

		state[0] += (state[0] / std::abs(state[0])) * norm;

		tempNormalize(state);

		
	}
}


Matrix SVD(StateVector state) {
        Matrix A(resize(state));
        Matrix B(A * A.conjugate().transpose());

        std::complex<double> B_trace = B[0][0] + B[1][1];
        std::complex<double> B_determinant = B[0][0] * B[1][1] - B[0][1] * std::conj(B[0][1]);
        std::complex<double> eigenvalue1 = (B_trace + std::sqrt(pow(B_trace, 2) - 4 * B_determinant)) / 2;
        std::complex<double> eigenvalue2 = (B_trace - std::sqrt(pow(B_trace, 2) - 4 * B_determinant)) / 2;

        StateVector eigenvector1(2);
        StateVector eigenvector2(2);

        if (std::abs(B[0][1]) > std::abs(eigenvalue1 - B[0][0])) {
                eigenvector1[0] = B[0][1];
                eigenvector1[1] = eigenvalue1 - B[0][0];

                normalize(eigenvector1);

                eigenvector2[0] = {-1, 0} * std::conj(eigenvalue1 - B[0][0]);
                eigenvector2[1] = std::conj(B[0][1]);
        } else {
                eigenvector1[0] = eigenvalue1 - B[1][1];
                eigenvector1[1] = std::conj(B.matrix[0][1]);

		normalize(eigenvector1);

		eigenvector2[0] = {-1, 0} * std::conj(B[0][1]);
                eigenvector2[1] = std::conj(eigenvalue1 - B[1][1]);
        }

        std::vector<std::vector<std::complex<double>>> U(2, 2);

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

        Matrix V();
}

Matrix findEvolution(double step) {
	Matrix O({{cd{1, 0}, cd{0, 1}}, {cd{0, 1}, cd{1, 0}}});

	O = (((O * cd{step, 0}) * cd{0, 1}) * cd{1 / hBar, 0}) * cd{-1, 0};

	O = I + O;

	return O;
}

void evolve(QuantumState state, Matrix evolution) {
	state.set(evolution * state.get());
}
