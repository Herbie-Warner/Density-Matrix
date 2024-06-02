//Utilities.h
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <string>

#ifndef UTILITIES_H
#define UTILITIES_H

namespace Utilities {
  typedef Eigen::MatrixXcd Matrix;
  typedef Eigen::MatrixXcd Operator;

  Matrix partialTrace(const Matrix& rho, int qubit);
  void normalise_matrix(Matrix& rho);
  Matrix tensorProduct(const Matrix& A, const Matrix& B);
  std::ostream& operator<<(std::ostream& os, const Matrix& matrix);
  bool isPure(const Matrix& matrix);
  bool validateDensityMatrix(Matrix& rho);
  std::vector<std::string> getPurificationString(const std::vector<std::string>& original);
}

#endif