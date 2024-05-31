//DiscordCalc.cpp

#include "methods/DiscordCalc.h"

using namespace Methods;
using namespace Eigen;
using std::cout, std::endl;

namespace DiscordCalc
{
  double mutualInformationCalc(const Methods::DensityMatrix& rhoAB, int dimA, int dimB)
  {
    double S_A = computeEntropy(partialTraceSecondSubsystem(rhoAB, dimA, dimB));
    double S_B = computeEntropy(partialTraceFirstSubsystem(rhoAB, dimA, dimB));
    double S_AB = computeEntropy(rhoAB);

    return S_A + S_B - S_AB;
  }

  double J_PI_AB_PauliX(const Methods::DensityMatrix& rhoAB, int dimA, int dimB)
  {
    
    Operator identityMatrix = MatrixXd::Identity(dimA, dimA);

    Operator X_d = MatrixXd::Zero(dimB, dimB);
    for (int i = 0; i < dimB; ++i) {
      int j = (i + 1) % dimB;
      X_d(i, j) = 1;
    }

    Operator operate = tensorProduct(identityMatrix, X_d);
    DensityMatrix rhoAB_p = operate * rhoAB * operate.adjoint();

    DensityMatrix rhoA = partialTraceSecondSubsystem(rhoAB, dimA, dimB);
    DensityMatrix rhoA_x = partialTraceSecondSubsystem(rhoAB_p, dimA, dimB);
    normalise_matrix(rhoA_x);

    cout<<rhoAB<<endl;
    cout<<rhoAB_p<<endl;


    double S_A = computeEntropy(rhoA);
    double S_A_p = computeEntropy(rhoA_x);

    cout<<rhoA<<endl;
    cout<<S_A<<endl;

    cout<<S_A_p<<endl;
    cout<<rhoA_x<<endl;

    return S_A - S_A_p; //For only pauliX measurement e.g

  }

  double J_PI_AB(const Methods::DensityMatrix& rhoAB, int dimA, int dimB)
  {
    Operator identityMatrix_A = MatrixXd::Identity(dimA, dimA);
    Operator identityMatrix_B = MatrixXd::Identity(dimB, dimB);




    Operator sigmaX = MatrixXd::Zero(dimB, dimB);
    for (int i = 0; i < dimB; ++i) {
      int j = (i + 1) % dimB;
      sigmaX(i, j) = 1;
    }

    Operator sigmaY = MatrixXcd::Zero(dimB, dimB);
    for (int i = 0; i < dimB - 1; ++i) {
      sigmaY(i, i + 1) = std::complex<double>(0, -1);
      sigmaY(i + 1, i) = std::complex<double>(0, 1);
    }


    Operator sigmaZ = MatrixXcd::Zero(dimB, dimB);
    for (int i = 0; i < dimB; ++i) {
      sigmaZ(i, i) = std::complex<double>(1.0 - 2.0 * i / (dimB - 1), 0);
    }

    Operator operate_identity = tensorProduct(identityMatrix_A, identityMatrix_B);
    Operator operate_sigmaX = tensorProduct(identityMatrix_A, sigmaX);
    Operator operate_sigmaY = tensorProduct(identityMatrix_A, sigmaY);
    Operator operate_sigmaZ = tensorProduct(identityMatrix_A, sigmaZ);

    
    DensityMatrix rhoAB_p_I = operate_identity * rhoAB * operate_identity.adjoint();
    DensityMatrix rhoAB_p_X = operate_sigmaX * rhoAB * operate_sigmaX.adjoint();
    DensityMatrix rhoAB_p_Y = operate_sigmaY * rhoAB * operate_sigmaY.adjoint();
    DensityMatrix rhoAB_p_Z = operate_sigmaZ * rhoAB * operate_sigmaZ.adjoint();


    normalise_matrix(rhoAB_p_I);
    normalise_matrix(rhoAB_p_X);
    normalise_matrix(rhoAB_p_Y);
    normalise_matrix(rhoAB_p_Z);

    cout<<rhoAB_p_X<<endl;
    cout<<rhoAB_p_Y<<endl;
    cout<<rhoAB_p_Z<<endl;

    DensityMatrix rhoA = partialTraceSecondSubsystem(rhoAB, dimA, dimB);


    DensityMatrix rhoA_I = partialTraceSecondSubsystem(rhoAB_p_I, dimA, dimB);
    DensityMatrix rhoA_X = partialTraceSecondSubsystem(rhoAB_p_X, dimA, dimB);
    DensityMatrix rhoA_Y = partialTraceSecondSubsystem(rhoAB_p_Y, dimA, dimB);
    DensityMatrix rhoA_Z = partialTraceSecondSubsystem(rhoAB_p_Z, dimA, dimB);

    cout<<"---"<<endl;
    cout<<partialTraceFirstSubsystem(rhoAB_p_I,dimA,dimB)<<endl;

    cout << partialTraceFirstSubsystem(rhoAB_p_X, dimA, dimB) << endl;

    cout << partialTraceFirstSubsystem(rhoAB_p_Y, dimA, dimB) << endl;

    cout << partialTraceFirstSubsystem(rhoAB_p_Z, dimA, dimB) << endl;

    cout<<"---"<<endl;


    double SA_I = computeEntropy(rhoA_I);
    double SA_X = computeEntropy(rhoA_X);
    double SA_Y = computeEntropy(rhoA_Y);
    double SA_Z = computeEntropy(rhoA_Z);


    cout<<"000000"<<"\n"<<endl;
    cout<<SA_I<<endl;
    cout<<rhoA_I<<endl;
    cout<<SA_X<<endl;
    cout<<rhoA_X<<endl;
    cout<<SA_Y<<endl;
    cout<<rhoA_Y<<endl;
    cout<<SA_Z<<endl;
    cout<<rhoA_Z<<endl;

    double sum = 0.25 * SA_I + 0.25 * SA_X + 0.25 * SA_Y + 0.25 * SA_Z;
    double S_A = computeEntropy(rhoA);

    return S_A - sum;
  }
}