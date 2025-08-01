#ifndef MATRIX_H
#define MATRIX_H

#include "hilbert.h"

#include <complex>
#include <vector>

class StateVector;

class Matrix {
	public:
		Matrix(std::vector<std::vector<std::complex<double>>> initMatrix);
		Matrix(const Matrix& otherMatrix);
		Matrix(size_t rows, size_t columns);

		const std::vector<std::vector<std::complex<double>>>& get() const;
		void set(const Matrix& otherMatrix);
		void resize(size_t rows, size_t columns);
		size_t getRows() const;
		size_t getColumns() const;
		void swapRows(size_t i, size_t j);
		void swapColumns(size_t i, size_t j);
		Matrix identity() const;
		void roundValues();
		void print() const;

		Matrix operator+(const Matrix& otherMatrix) const;
		Matrix operator*(std::complex<double> scalar) const;
		Matrix operator*(const Matrix& otherMatrix) const;
		StateVector operator*(const StateVector& vector) const;
		std::vector<std::complex<double>>& operator[](size_t row);
		const std::vector<std::complex<double>>& operator[](size_t row) const;

		Matrix transpose() const;
		Matrix invert() const;
		Matrix conjugate() const;

	private:
		std::vector<std::vector<std::complex<double>>> matrix;
		size_t rows, columns;

		void update();
		void orderRows();
};

Matrix identityPad(Matrix A, size_t size);
Matrix kroneckerProduct(Matrix matrix1, Matrix matrix2);

#endif
