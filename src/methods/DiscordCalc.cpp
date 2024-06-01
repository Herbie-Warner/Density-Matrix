//DiscordCalc.cpp

#include "methods/DiscordCalc.h"

using namespace Methods;
using namespace Eigen;
using std::cout, std::endl, std::abs;

namespace DiscordCalc
{
  double mutualInformationCalc(const Methods::DensityMatrix& rhoAB)
  {
    double S_A = computeEntropy(partialTrace(rhoAB,2));
    double S_B = computeEntropy(partialTrace(rhoAB, 1));
    double S_AB = computeEntropy(rhoAB);

    return S_A + S_B - S_AB;
  }

  

  double J_AB_2Qubits(const Methods::DensityMatrix& rhoAB)
  {

    Eigen::Vector2cd ket0, ket1, basis0, basis1;
    ket0 << 1, 0;
    ket1 << 0, 1;
    const std::complex<double> I(0, 1);


    double theta = 0;
    double phi = 0;
    double maxTheta = 2*PI;
    double maxPhi = 2 * PI;
    int size = 100;

    double best_entropy = 10;
    double best_theta = 0;
    double best_phi = 0;


    for (int i = 0; i < size; ++i)
    {
     
      for (int j = 0; j < size; ++j)
      {

        basis0 = (std::cos(phi) + I * std::sin(phi)) * std::cos(theta) * ket0 + std::sin(theta) * ket1;
        basis1 = (std::cos(phi) - I * std::sin(phi)) * std::sin(theta) * ket0 - std::cos(theta) * ket1;

        Operator A = basis0 * basis0.adjoint();
        Operator B = basis1 * basis1.adjoint();
        std::vector<Operator> operators;

        Operator identityMatrix_A = MatrixXd::Identity(2, 2);
        operators.push_back(tensorProduct(identityMatrix_A, A));
        operators.push_back(tensorProduct(identityMatrix_A, B));


        double entropy_sum = 0;


        for (const auto& oper : operators) {

          double prob = abs((oper * rhoAB).trace());
          DensityMatrix rhoAB_post_measurement = oper * (rhoAB * oper.adjoint()) / prob;
          entropy_sum += prob * computeEntropy(partialTrace(rhoAB_post_measurement, 2));
        }
        if (entropy_sum < best_entropy)
        {
          best_entropy = entropy_sum;
          best_theta = theta;
          best_phi = phi;
        }
        phi += maxPhi/size;
      }
      theta += maxTheta / size;
    }

    DensityMatrix rhoA = partialTrace(rhoAB,2);
    double S_A = computeEntropy(rhoA);
    return S_A - best_entropy;
    /*
    std::complex<double> I(0, 1);

    // Define qubit basis states
    Eigen::Vector2cd ket0, ket1, ketPlus, ketMinus, ketPlusI, ketMinusI;

    ket0 << 1, 0;
    ket1 << 0, 1;
    ketPlus << 1 / std::sqrt(2), 1 / std::sqrt(2);
    ketMinus << 1 / std::sqrt(2), -1 / std::sqrt(2);
    ketPlusI << 1 / std::sqrt(2), I / std::sqrt(2);
    ketMinusI << 1 / std::sqrt(2), -I / std::sqrt(2);

    // Define density matrices
    Operator A = (ket0 * ket0.adjoint());
    Operator B = ket1 * ket1.adjoint();
    Operator C = (ketPlus * ketPlus.adjoint());
    Operator D = (ketMinus * ketMinus.adjoint());
    Operator E = (ketPlusI * ketPlusI.adjoint());
    Operator F = (ketMinusI * ketMinusI.adjoint());

    std::vector<Operator> operators;
    operators.push_back(A);
    operators.push_back(B);
    operators.push_back(C);
    operators.push_back(D);
    operators.push_back(E);
    operators.push_back(F);

    Operator identityMatrix_A = MatrixXd::Identity(dimA, dimA);
    double entropy_sum = 0;


    for (const auto& oper : operators){
      Operator full_space_oper = tensorProduct(identityMatrix_A, oper);
      double prob = abs((full_space_oper * rhoAB).trace());
      cout<<prob<<endl;
      DensityMatrix rhoAB_post_measurement = full_space_oper * (rhoAB * full_space_oper.adjoint()) / prob;
      cout<<rhoAB_post_measurement<<endl;
      entropy_sum += prob* computeEntropy(partialTraceSecondSubsystem(rhoAB_post_measurement, dimA, dimB));
    }

    
 
    
    Operator spin_up_B(2,2);
    spin_up_B<<1,0,0,0;


    Operator spin_down_B(2, 2);
    spin_down_B << 0.5, 0.5, 0.5, 0.5;

    Operator spin_upB = tensorProduct(identityMatrix_A, spin_up_B);
    Operator spin_downB = tensorProduct(identityMatrix_A, spin_down_B);

    double prob_up = abs((spin_upB*rhoAB).trace());
    double prob_down = abs((spin_downB * rhoAB).trace());

    DensityMatrix rhoAB_post_up = spin_upB*(rhoAB*spin_upB.adjoint())/prob_up;
    DensityMatrix rhoAB_post_down = spin_downB * (rhoAB * spin_downB.adjoint()) / prob_down;

   

    DensityMatrix rhoA = partialTraceSecondSubsystem(rhoAB, dimA, dimB);
    double S_A = computeEntropy(rhoA);


    double sum = prob_up*computeEntropy(partialTraceSecondSubsystem(rhoAB_post_up, dimA, dimB)) + prob_down * computeEntropy(partialTraceSecondSubsystem(rhoAB_post_down, dimA, dimB));
    

    DensityMatrix rhoA = partialTraceSecondSubsystem(rhoAB, dimA, dimB);
    double S_A = computeEntropy(rhoA);

    return S_A - entropy_sum;
    */
    
  }

  double compute_DW(const DensityMatrix& rhoAB)
  {
    DensityMatrix rhoB = partialTrace(rhoAB, 1);
    double dW = computeEntropy(rhoB) - computeEntropy(rhoAB) + ReflectedEntropy(rhoAB) / 2;
    return dW;
  }

  double compute_D(const DensityMatrix& rhoAB)
  {
    double res = J_AB_2Qubits(rhoAB);
    double mut_info = mutualInformationCalc(rhoAB);
    double D = mut_info - res;
    return D;
  }

  double compute_DT(const DensityMatrix& rhoAB)
  {
    DensityMatrix rhoCAB = canonicalPurification(rhoAB);
    DensityMatrix rhoCB = partialTrace(rhoCAB, 3);
    DensityMatrix rhoCpBpCB = canonicalPurification(rhoCB);
    DensityMatrix rhoCpBpB = partialTrace(partialTrace(rhoCpBpCB, 5), 4);
    DensityMatrix BpB = partialTrace(partialTrace(rhoCpBpB, 1), 1);
    double S_BBp = computeEntropy(BpB);
    double mut_info = mutualInformationCalc(rhoAB);
    return mut_info - S_BBp;
  }

  double compute_MarkovGap(const Methods::DensityMatrix& rhoAB)
  {
    double mut = mutualInformationCalc(rhoAB);
    double reflected_entropy = ReflectedEntropy(rhoAB);
    return reflected_entropy - mut;
  }
}