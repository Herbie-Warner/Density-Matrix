//g++ -I "C:\Users\herbi\eigen-3.4.0" -o den '.\Density Matrix.cpp'
//Herbie Warner 30/05/2024

#include"methods/Utilities.h"
#include<stdexcept>

using namespace Eigen;
using namespace std;

namespace Utilities
{
  double computeEntropy(const Matrix& rho) {
    SelfAdjointEigenSolver<MatrixXcd> es(rho);
    VectorXd eigenvalues = es.eigenvalues().real();
    double entropy = 0.0;
    for (int i = 0; i < eigenvalues.size(); ++i) {
      double lambda = eigenvalues(i);
      if (lambda > 0) {
        entropy -= lambda * log2(lambda);
      }
    }
    return entropy;
  }

  void normalise_matrix(Matrix& rho) {
    double trace = std::abs(rho.trace().real());
    if (trace != 1.0) {
      rho /= trace;
    }
  }


  Matrix tensorProduct(const Matrix& A, const Matrix& B) {
    Matrix C(A.rows() * B.rows(), A.cols() * B.cols());
    for (int i = 0; i < A.rows(); ++i) {
      for (int j = 0; j < A.cols(); ++j) {
        C.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) = A(i, j) * B;
      }
    }
    return C;
  }


  bool isPure(const Matrix& matrix) {
    auto square = matrix * matrix;
    return std::abs(square.trace().real() - 1.0) < 1e-10;
  }

  ostream& operator<<(std::ostream& os, const Matrix& matrix) {
    for (int i = 0; i < matrix.rows(); ++i) {
      for (int j = 0; j < matrix.cols(); ++j) {
        os << matrix(i, j) << "\t";
      }
      os << endl;
    }
    return os;
  }

  Matrix partialTrace(const Matrix& rho, int qubit) {

    int totalQubits = std::log2(rho.rows());
    qubit = totalQubits - qubit;
   

    int dim = 1 << totalQubits;
    int subDim = 1 << (totalQubits - 1);
    Matrix reducedDensityMatrix = Matrix::Zero(subDim, subDim);

    for (int i = 0; i < dim; ++i) {
      for (int j = 0; j < dim; ++j) {
        int i0 = (i & ((1 << qubit) - 1)) | ((i >> (qubit + 1)) << qubit);
        int j0 = (j & ((1 << qubit) - 1)) | ((j >> (qubit + 1)) << qubit);
        int i1 = (i >> qubit) & 1;
        int j1 = (j >> qubit) & 1;

        if (i1 == j1) {
          reducedDensityMatrix(i0, j0) += rho(i, j);
        }
      }
    }

    return reducedDensityMatrix;
  }

  bool validateDensityMatrix(Matrix& rho)
  {
    constexpr double tolerance = - 0.0000001;
    normalise_matrix(rho);
    SelfAdjointEigenSolver<MatrixXcd> es(rho);
    VectorXd eigenvalues = es.eigenvalues().real();
    for (const auto& eigenvalue : eigenvalues)
    {
      if (eigenvalue <= tolerance)
      {
        throw std::invalid_argument("-VE eigenvalues");
        return false;
      }
    }
    if(rho.determinant().real() <= tolerance) { throw std::invalid_argument("-VE det"); }

    if(rho.adjoint() != rho) {throw std::invalid_argument("non hermitian"); }
    return true;
  }

  std::vector<std::string> getPurificationString(const std::vector<std::string>& original)
  {
    std::vector<std::string> purified;
    for (const auto& string : original) {
      purified.push_back(string + "'");
    }
    return purified;
  }
}