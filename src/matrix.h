#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <cmath>
#include <vector>
#include <complex>

class Matrix {
public:
    Matrix(std::vector<std::vector<std::complex<double>>> initMatrix);
    Matrix(const Matrix& otherMatrix);
    Matrix(int rows, int columns);

    const std::vector<std::vector<std::complex<double>>>& get() const;
    void set(const Matrix& otherMatrix);
    void print() const;

    Matrix operator+(const Matrix& otherMatrix) const;
    Matrix operator*(std::complex<double> scalar) const;
    Matrix operator*(const Matrix& otherMatrix) const;
    Matrix transpose() const;
    Matrix invert() const;

private:
    std::vector<std::vector<std::complex<double>>> matrix;
    int rows, columns;

    void update();
    bool sameDimensions(const Matrix& otherMatrix) const;
    bool multiplyApplicable(const Matrix& otherMatrix) const;
    bool isSquare() const;
    bool isZeroColumn(int column) const;
    void orderRows();
    void roundValues();
};

#endif
