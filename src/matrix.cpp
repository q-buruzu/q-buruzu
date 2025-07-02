#include "matrix.h"
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>

Matrix::Matrix(std::vector<std::vector<std::complex<double>>> initMatrix) {
	matrix = initMatrix;
	update();
}

Matrix::Matrix(const Matrix& otherMatrix) {
	matrix = otherMatrix.matrix;
	update();
}

Matrix::Matrix(int rows, int columns) {
	matrix.resize(rows, std::vector<std::complex<double>>(columns));
	update();
}

const std::vector<std::vector<std::complex<double>>>& Matrix::get() const {
	return matrix;
}

void Matrix::set(const Matrix& otherMatrix) {
	matrix = otherMatrix.matrix;
	update();
}

void Matrix::update() {
	rows = matrix.size();
	columns = matrix[0].size();
}

Matrix Matrix::operator+(const Matrix& otherMatrix) const {
	if (!sameDimensions(otherMatrix)) {
		std::cout << "\nMATRIX ADDITION FAILED\n";
		return *this;
	}

	Matrix resultMatrix(rows, columns);

	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			resultMatrix.matrix[i][j] = matrix[i][j] + otherMatrix.matrix[i][j];
		}
	}

	resultMatrix.roundValues();

	return resultMatrix;
}

Matrix Matrix::operator*(std::complex<double> scalar) const {
	Matrix resultMatrix(rows, columns);

	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			resultMatrix.matrix[i][j] = scalar * matrix[i][j];
		}
	}

	resultMatrix.roundValues();

	return resultMatrix;
}

Matrix Matrix::operator*(const Matrix& otherMatrix) const {
	if (!multiplyApplicable(otherMatrix)) {
		std::cout << "\nMATRIX MULTIPLICATION FAILED\n";
		return *this;
	}

	Matrix resultMatrix(rows, otherMatrix.columns);

	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < otherMatrix.columns; ++j) {
			for (int k = 0; k < columns; ++k) {
				resultMatrix.matrix[i][j] += matrix[i][k] * otherMatrix.matrix[k][j];
			}
		}
	}

	resultMatrix.roundValues();

	return resultMatrix;
}

Matrix Matrix::transpose() const {
	Matrix resultMatrix(columns, rows);

	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			resultMatrix.matrix[j][i] = matrix[i][j];
		}
	}

	return resultMatrix;
}

Matrix Matrix::invert() const {
	if (!isSquare()) {
		std::cout << "\nONLY SQUARE MATRICES CAN BE INVERTED\n";
		return *this;
	}

	Matrix augment(*this);

	for (int i = 0; i < rows; ++i) {
		augment.matrix[i].resize(2 * rows, {0, 0});
		augment.matrix[i][i + rows] = {1, 0};
	}

	augment.update();
	augment.orderRows();

	for (int i = 0; i < rows; ++i) {
		std::complex<double> pivot = augment.matrix[i][i];

		for (int j = 0; j < 2 * rows; ++j) {
			augment.matrix[i][j] /= pivot;
		}

		for (int k = 0; k < rows; ++k) {
			if (k == i) continue;

			std::complex<double> factor = augment.matrix[k][i];

			for (int j = 0; j < 2 * rows; ++j) {
				augment.matrix[k][j] -= factor * augment.matrix[i][j];
			}
		}
	}

	Matrix inverseMatrix(rows, rows);

	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < rows; ++j) {
			inverseMatrix.matrix[i][j] = augment.matrix[i][j + rows];
		}
	}

	inverseMatrix.roundValues();

	return inverseMatrix;
}

void Matrix::print() const {
	for (int i = 0; i < rows; ++i) {
        	for (int j = 0; j < columns; ++j) {
			std::cout << matrix[i][j] << "\t";
		}

		std::cout << "\n";
	}

	std::cout << "\n";
}

bool Matrix::sameDimensions(const Matrix& otherMatrix) const {
	return columns == otherMatrix.columns && rows == otherMatrix.rows;
}

bool Matrix::multiplyApplicable(const Matrix& otherMatrix) const {
	return columns == otherMatrix.rows;
}

bool Matrix::isSquare() const {
	return rows == columns;
}

bool Matrix::isZeroColumn(int column) const {
	for (int i = 0; i < rows; ++i) {
		if (matrix[i][column] != std::complex<double>{0.0, 0.0}) {
			return false;
		}
	}

	return true;
}

void Matrix::roundValues() {
	double power = 1e8;

	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j) {
			double real = std::round(matrix[i][j].real() * power) / power;
			double imag = std::round(matrix[i][j].imag() * power) / power;

			matrix[i][j] = {real, imag};
		}
	}
}

void Matrix::orderRows() {
	int maxRow;
	double maxValue;

	for (int i = 0; i < rows; ++i) {
		if (isZeroColumn(i)) {
			std::cout << "\nMATRIX INVERSION FAILED\n";
			return;
		}

		maxRow = i;
		maxValue = std::abs(matrix[i][i]);

		for (int j = i + 1; j < rows; ++j) {
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

