#ifndef MATRIX_H
#define MATRIX_H

#include "hilbert.h"

#include <complex>
#include <vector>

class Matrix {
	public:
		Matrix(std::vector<std::vector<std::complex<double>>> initMatrix);
		Matrix(const Matrix& otherMatrix);
		Matrix(size_t rows, size_t columns);

		const std::vector<std::vector<std::complex<double>>>& get() const;
		void set(const Matrix& otherMatrix);
		size_t getRows() const;
		size_t getColumns() const;
		void print() const;

		Matrix operator+(const Matrix& otherMatrix) const;
		Matrix operator*(std::complex<double> scalar) const;
		Matrix operator*(const Matrix& otherMatrix) const;
		StateVector operator*(const StateVector& vector) const;
		Matrix transpose() const;
		Matrix invert() const;
		Matrix conjugate() const;

	private:
		std::vector<std::vector<std::complex<double>>> matrix;
		size_t rows, columns;

		void update();
		void orderRows();
		void roundValues();
};

#endif
