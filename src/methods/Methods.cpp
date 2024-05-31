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
        os << matrix(i, j) << "\t";
      }
      os << endl;
    }
    return os;
  }

  DensityMatrix partialTrace(const DensityMatrix& rho, int qubit) {

    int totalQubits = std::log2(rho.rows());
    qubit = totalQubits - qubit;
   

    int dim = 1 << totalQubits;
    int subDim = 1 << (totalQubits - 1);
    DensityMatrix reducedDensityMatrix = DensityMatrix::Zero(subDim, subDim);

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

  double ReflectedEntropy(const DensityMatrix& rho) {  // Only for rho 2 qubit
    DensityMatrix ABC = canonicalPurification(rho);
    DensityMatrix AC = partialTrace(ABC, 2);
    DensityMatrix ACApCp = tensorProduct(AC, AC);
    DensityMatrix ACAp = partialTrace(partialTrace(ACApCp,6),5);
    DensityMatrix AAp = partialTrace(partialTrace(ACAp,2),2);
    return computeEntropy(AAp);  
  }
}