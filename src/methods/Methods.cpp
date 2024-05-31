//g++ -I "C:\Users\herbi\eigen-3.4.0" -o den '.\Density Matrix.cpp'
//Herbie Warner 30/05/2024

#include"methods/Methods.h"

using namespace Eigen;
using namespace std;

namespace Methods
{
  double computeEntropy(const DensityMatrix& rho) {
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


  void normalise_matrix(DensityMatrix& rho) {
    double trace = std::abs(rho.trace().real());
    if (trace != 1.0) {
      rho /= trace;
    }
  }


  DensityMatrix tensorProduct(const DensityMatrix& A, const DensityMatrix& B) {
    DensityMatrix C(A.rows() * B.rows(), A.cols() * B.cols());
    for (int i = 0; i < A.rows(); ++i) {
      for (int j = 0; j < A.cols(); ++j) {
        C.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) = A(i, j) * B;
      }
    }
    return C;
  }

  DensityMatrix partialTraceSecondSubsystem(const DensityMatrix& rho, int dimA, int dimB) {
    DensityMatrix tracedRho = DensityMatrix::Zero(dimA, dimA);

    for (int i = 0; i < dimA; ++i) {
      for (int j = 0; j < dimA; ++j) {
        for (int k = 0; k < dimB; ++k) {
          tracedRho(i, j) += rho(i * dimB + k, j * dimB + k);
        }
      }
    }

    return tracedRho;
  }

  DensityMatrix partialTraceFirstSubsystem(const DensityMatrix& rho, int dimA, int dimB) {
    DensityMatrix tracedRho = DensityMatrix::Zero(dimB, dimB);

    for (int i = 0; i < dimB; ++i) {
      for (int j = 0; j < dimB; ++j) {
        for (int k = 0; k < dimA; ++k) {
          tracedRho(i, j) += rho(k * dimB + i, k * dimB + j);
        }
      }
    }

    return tracedRho;
  }


  bool isPure(const DensityMatrix& matrix) {
    auto square = matrix * matrix;
    return std::abs(square.trace().real() - 1.0) < 1e-10;
  }

  DensityMatrix canonicalPurification(const DensityMatrix& rho) {
    if (rho.determinant().real() < 0) { throw std::invalid_argument("-ve det cannot purify"); }

    SelfAdjointEigenSolver<MatrixXcd> es(rho);
    VectorXd eigenvalues = es.eigenvalues().real();
    MatrixXcd eigenvectors = es.eigenvectors();
    int dim = static_cast<int>(eigenvalues.size());

    DensityMatrix psi = DensityMatrix::Zero(dim * dim, 1);
    for (int i = 0; i < dim; ++i) {
      if (eigenvalues(i) > 0) {
        for (int j = 0; j < dim; ++j) {
          psi(i * dim + j, 0) = sqrt(eigenvalues(i)) * eigenvectors(j, i);
        }
      }
    }



    DensityMatrix purifiedState = psi * psi.adjoint();
    return purifiedState;
  }

  ostream& operator<<(std::ostream& os, const DensityMatrix& matrix) {
    for (int i = 0; i < matrix.rows(); ++i) {
      for (int j = 0; j < matrix.cols(); ++j) {
        os << real(matrix(i, j)) << "\t";
      }
      os << endl;
    }
    return os;
  }

  DensityMatrix partialTrace(const DensityMatrix& rho, const vector<int>& dims, int subsystem) {
    int dimA = 1, dimB = dims[subsystem], dimC = 1;
    for (int i = 0; i < subsystem; ++i) dimA *= dims[i];
    for (size_t i = subsystem + 1; i < dims.size(); ++i) dimC *= dims[i];

    int newDim = dimA * dimC;
    DensityMatrix tracedRho = DensityMatrix::Zero(newDim, newDim);

    for (int i = 0; i < dimA; ++i) {
      for (int j = 0; j < dimA; ++j) {
        for (int k = 0; k < dimC; ++k) {
          for (int l = 0; l < dimC; ++l) {
            for (int m = 0; m < dimB; ++m) {
              int row = (i * dimB + m) * dimC + k;
              int col = (j * dimB + m) * dimC + l;
              tracedRho(i * dimC + k, j * dimC + l) += rho(row, col);
            }
          }
        }
      }
    }

    return tracedRho;
  }


}