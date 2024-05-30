//g++ -I "C:\Users\herbi\eigen-3.4.0" -o den '.\Density Matrix.cpp'
//Herbie Warner 30/05/2024
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>

using namespace Eigen;
using namespace std;

typedef MatrixXcd DensityMatrix;

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

void print()
{

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

DensityMatrix partialTraceFirstSubsystem(const DensityMatrix& rho, int dimA, int dimB) {
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

DensityMatrix partialTraceSecondSubsystem(const DensityMatrix& rho, int dimA, int dimB) {
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

ostream& operator<<(ostream& os, const DensityMatrix& matrix) {
  for (int i = 0; i < matrix.rows(); ++i) {
    for (int j = 0; j < matrix.cols(); ++j) {
      os << real(matrix(i, j)) << "\t";
    }
    os << endl;
  }
  return os;
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

int main() {

  /*
  DensityMatrix rho(4, 4);
  rho <<  5, 2, 0.0, 3,
          2, 5, 0.0, 0.0,
          0.0, 0.0, 5, 0,
          3, 0.0, 0, 5;

  normalise_matrix(rho);

  cout << "Pure: " << std::boolalpha << isPure(rho) << std::endl;
  cout << "rho:" << endl << rho << endl;

  DensityMatrix psi = canonicalPurification(rho);
  cout << "Canonical Purification psi:" << endl << psi << endl;
  cout << "Pure: " << std::boolalpha << isPure(psi) << std::endl;

  DensityMatrix rhoA = partialTraceSecondSubsystem(psi, 4, 4);
  cout << "RhoA" << endl << rhoA << endl;
  double entropyA = computeEntropy(rhoA);
  cout << "Entropy A: " << entropyA << endl;

  DensityMatrix rho_a = partialTraceSecondSubsystem(rhoA, 2, 2);
  cout << "RhoA" << endl << rho_a << endl;
  entropyA = computeEntropy(rho_a);
  cout << "Entropy A: " << entropyA << endl;

  */


  /*
  DensityMatrix rhoA(2, 2);
  rhoA << 1,0,0,1;


  DensityMatrix rhoB(2, 2);
  rhoB << 2,1,1,2;

  normalise_matrix(rhoA);
  normalise_matrix(rhoB);

  auto prod = tensorProduct(rhoA, rhoB);

  double entropy_AB = computeEntropy(prod);
  double entropy_A = computeEntropy(rhoA);
  double entropy_B = computeEntropy(rhoB);

  cout<<entropy_B+entropy_AB<<endl;
  cout<<entropy_A<<endl;
  cout<<entropy_B<<endl;
  */

  DensityMatrix rhoA(2, 2);
  rhoA << 1, 0, 0, 1;
  normalise_matrix(rhoA);


  DensityMatrix psi = canonicalPurification(rhoA);

  DensityMatrix rhoB = partialTraceFirstSubsystem(psi, 2, 2);

  cout << psi << endl;
  cout << computeEntropy(psi) << endl;
  cout << computeEntropy(rhoA) << endl;
  cout << computeEntropy(rhoB) << endl;

  cout<<partialTraceSecondSubsystem(psi,2,2)<<endl;


  return 0;
}
