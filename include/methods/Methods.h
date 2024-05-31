//Methods.h
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <vector>

#ifndef METHODS_H
#define METHODS_H

namespace Methods {
  typedef Eigen::MatrixXcd DensityMatrix;
  typedef Eigen::MatrixXcd Operator;


  double computeEntropy(const DensityMatrix& rho);
  void normalise_matrix(DensityMatrix& rho);
  DensityMatrix tensorProduct(const DensityMatrix& A, const DensityMatrix& B);
  std::ostream& operator<<(std::ostream& os, const DensityMatrix& matrix);
  bool isPure(const DensityMatrix& matrix);
  DensityMatrix canonicalPurification(const DensityMatrix& rho);
  DensityMatrix partialTrace(const DensityMatrix& rho, int qubit);
  double ReflectedEntropy(const DensityMatrix& rho);
}

#endif