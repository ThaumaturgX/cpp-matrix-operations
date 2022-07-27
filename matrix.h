#ifndef MATRIX_H
#define MATRIX_H
#include <vector>
#include <cmath>

using std::vector;

// Returns a transposed matrix
vector<vector<double>> transpose(vector<vector<double>> &mat);

// Minor of an element (i, j) in a matrix
vector<vector<double>> minor(vector<vector<double>> &mat, int i, int j);


// The functions below only work with square matrices

// Determinant
double det(vector<vector<double>> &mat);

// Returns an inversed matrix
vector<vector<double>> inverse(vector<vector<double>> &mat);

// Solves a linear system Ax = b by using
// an inversed matrix and returns a vector
// of solutions x
vector<double> solveSystem(vector<vector<double>> &mat, vector<double> &b);

#endif // MATRIX_H