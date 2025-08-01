#include "error_utils.h"
#include "hilbert.h"
#include "matrix.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <vector>

Matrix::Matrix(std::vector<std::vector<std::complex<double>>> otherMatrix) {
	matrix = otherMatrix;

	update();
}

Matrix::Matrix(const Matrix& otherMatrix) {
	matrix = otherMatrix.matrix;

	update();
}

Matrix::Matrix(size_t rows, size_t columns) {
	matrix.resize(rows, std::vector<std::complex<double>>(columns, {0, 0}));

	update();
}

const std::vector<std::vector<std::complex<double>>>& Matrix::get() const {
	return matrix;
}

void Matrix::set(const Matrix& otherMatrix) {
	matrix = otherMatrix.matrix;

	update();
}

void Matrix::resize(size_t rows, size_t columns) {
	matrix.resize(rows);

	for (size_t i = 0; i < rows; ++i) {
		matrix[i].resize(columns, {0, 0});
	}

	update();
}

size_t Matrix::getRows() const {
	return rows;
}

size_t Matrix::getColumns() const {
	return columns;
}

void Matrix::update() {
	rows = matrix.size();
	columns = matrix[0].size();
}

void Matrix::swapRows(size_t i, size_t j) {
	if (i >= 0 && j >= 0 && i < getRows() && j < getRows()) {
                std::swap(matrix[i], matrix[j]);
        } else {
                std::cout << "ERROR CANNOT SWAP NONEXISTENT ROWS";
        }
}

void Matrix::swapColumns(size_t i, size_t j) {
        if (i < getColumns() && j < getColumns()) {
                std::cout << "ERROR CANNOT SWAP NONEXISTENT COLUMNS";
        } else {
                for (size_t k = 0; k < getRows(); ++k) {
                        std::swap(matrix[k][i], matrix[k][j]);
                }
        }
}

Matrix Matrix::operator+(const Matrix& otherMatrix) const {
	throwError(sameDimensions(*this, otherMatrix), "MATRICES MUST BE SAME SIZE (ADDITION)");

	Matrix resultMatrix(rows, columns);

	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < columns; ++j) {
			resultMatrix.matrix[i][j] = matrix[i][j] + otherMatrix.matrix[i][j];
		}
	}

	return resultMatrix;
}

Matrix Matrix::operator*(std::complex<double> scalar) const {
	Matrix resultMatrix(rows, columns);

	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < columns; ++j) {
			resultMatrix.matrix[i][j] = scalar * matrix[i][j];
		}
	}

	return resultMatrix;
}

Matrix Matrix::operator*(const Matrix& otherMatrix) const {
	throwError(multiplyApplicable(*this, otherMatrix), "MATRICES MUST HAVE CORRECT SIZES (MULTIPLICATION)");

	Matrix resultMatrix(rows, otherMatrix.columns);

	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < otherMatrix.columns; ++j) {
			for (size_t k = 0; k < columns; ++k) {
				resultMatrix.matrix[i][j] += matrix[i][k] * otherMatrix.matrix[k][j];
			}
		}
	}

	return resultMatrix;
}

StateVector Matrix::operator*(const StateVector& vector) const {
	StateVector resultVector(rows);

	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < columns; ++j) {
			resultVector[i] += matrix[i][j] * vector[j];
		}
	}

	return resultVector;
}

std::vector<std::complex<double>>& Matrix::operator[](size_t row) {
	return matrix.at(row);
}

const std::vector<std::complex<double>>& Matrix::operator[](size_t row) const {
	return matrix.at(row);
}

Matrix Matrix::transpose() const {
	Matrix resultMatrix(columns, rows);

	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < columns; ++j) {
			resultMatrix.matrix[j][i] = matrix[i][j];
		}
	}

	return resultMatrix;
}

Matrix Matrix::invert() const {
	throwError(isSquare(*this), "MATRIX MUST BE SQUARE (INVERSION)");

	Matrix augment(*this);

	for (size_t i = 0; i < rows; ++i) {
		augment.matrix[i].resize(2 * rows, {0, 0});
		augment.matrix[i][i + rows] = {1, 0};
	}

	augment.update();
	augment.orderRows();

	for (size_t i = 0; i < rows; ++i) {
		std::complex<double> pivot = augment.matrix[i][i];

		for (size_t j = 0; j < 2 * rows; ++j) {
			augment.matrix[i][j] /= pivot;
		}

		for (size_t k = 0; k < rows; ++k) {
			if (k == i) continue;

			std::complex<double> factor = augment.matrix[k][i];

			for (size_t j = 0; j < 2 * rows; ++j) {
				augment.matrix[k][j] -= factor * augment.matrix[i][j];
			}
		}
	}

	Matrix inverseMatrix(rows, rows);

	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < rows; ++j) {
			inverseMatrix.matrix[i][j] = augment.matrix[i][j + rows];
		}
	}

	inverseMatrix.roundValues();

	return inverseMatrix;
}

Matrix Matrix::conjugate() const {
	Matrix resultMatrix(rows, columns);

	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < columns; ++j) {
			resultMatrix.matrix[i][j] = std::conj(matrix[i][j]);
		}
	}

	return resultMatrix;
}

Matrix Matrix::identity() const {
	Matrix resultMatrix(rows, columns);

	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < columns; ++j) {
			resultMatrix.matrix[i][j] = {0, 0};
		}

		resultMatrix.matrix[i][i] = {1, 0};
	}

	return resultMatrix;
}

void Matrix::roundValues() {
	double power = 1e8;

	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < columns; ++j) {
			double real = std::round(matrix[i][j].real() * power) / power;
			double imag = std::round(matrix[i][j].imag() * power) / power;

			matrix[i][j] = {real, imag};
		}
	}
}

void Matrix::print() const {
	for (size_t i = 0; i < rows; ++i) {
        	for (size_t j = 0; j < columns; ++j) {
			std::cout << matrix[i][j] << "\t";
		}

		std::cout << "\n";
	}

	std::cout << "\n";
}

void Matrix::orderRows() {
	size_t maxRow;
	double maxValue;

	for (size_t i = 0; i < rows; ++i) {
		throwError(isZeroColumn(*this, i), "MATRIX CAN NOT HAVE ZERO COLUMN (INVERSION)");

		maxRow = i;
		maxValue = std::abs(matrix[i][i]);

		for (size_t j = i + 1; j < rows; ++j) {
			if (std::abs(matrix[j][i]) > maxValue) {
				maxRow = j;
				maxValue = std::abs(matrix[j][i]);
			}
		}

		if (maxRow != i) {
			std::swap(matrix[i], matrix[maxRow]);
		}
	}
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

Matrix kroneckerProduct(Matrix matrix1, Matrix matrix2) {
	size_t row1 = matrix1.getRows();
	size_t column1 = matrix1.getColumns();
	size_t row2 = matrix2.getRows();
	size_t column2 = matrix2.getColumns();

	Matrix resultMatrix(row1 * row2, column1 * column2);

	for (size_t i = 0; i < row1; ++i) {
		for (size_t j = 0; j < column1; ++j) {
			for (size_t k = 0; k < row2; ++k) {
				for (size_t l = 0; l < column2; ++l) {
					resultMatrix[i * row2 + k][j * column2 + l] = matrix1[i][j] * matrix2[k][l];
				}
			}
		}
	}

	return resultMatrix;
}
