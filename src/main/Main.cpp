/*

3.4 - done
3.5 - done


*/




#include "methods/Methods.h"
#include "methods/DiscordCalc.h"

using namespace Methods;
using namespace std;
using namespace DiscordCalc;

void prove_D_less_than_DW()
{
  DensityMatrix rhoAB(4, 4);

  rhoAB << 1, 0, 0, 0,
    0, 2, 0.1, 0,
    0, 0.1, 1, 0.5,
    0, 0, 0.5, 1;

  normalise_matrix(rhoAB);
  validateDensityMatrix(rhoAB);

  /*
  DensityMatrix rhoA(2, 2);
  rhoA << 1, 0.5, 0.5, 1;

  DensityMatrix rhoB(2, 2);
  rhoB << 1, 0.5, 0.5, 1;




  normalise_matrix(rhoA);
  normalise_matrix(rhoB);




  rhoAB = tensorProduct(rhoA, rhoB);

  */


  cout << rhoAB.determinant() << endl;

  double res = J_AB_2Qubits(rhoAB);

  double mut_info = mutualInformationCalc(rhoAB);
  double D = mut_info - res;

  DensityMatrix rhoB = partialTrace(rhoAB, 1);
  double dW = computeEntropy(rhoB) - computeEntropy(rhoAB) + ReflectedEntropy(rhoAB) / 2;



  cout << computeEntropy(rhoB) << endl;
  cout << computeEntropy(rhoAB) << endl;
  cout << ReflectedEntropy(rhoAB) / 2 << endl;

  std::cout << "Discord: " << D << endl;
  std::cout << "DW: " << dW << endl;
 // cout << "Pure: " << std::boolalpha << isPure(rho) << std::endl;
}

void verify_3_5()
{
  DensityMatrix rhoAB(4, 4);

  rhoAB << 1, 0, 0, 0,
    0, 2, 0.1, 0,
    0, 0.1, 1, 0.5,
    0, 0, 0.5, 1;

  DensityMatrix rhoA(2, 2);
  rhoA << 1, 0, 0, 1;

  DensityMatrix rhoB(2, 2);
  rhoB << 1, 0, 0, 1;


  normalise_matrix(rhoA);
  normalise_matrix(rhoB);




  rhoAB = tensorProduct(rhoA, rhoB);

  normalise_matrix(rhoAB);
  validateDensityMatrix(rhoAB);

  double D = compute_D(rhoAB);
  double DW = compute_DW(rhoAB);
  double DT = compute_DT(rhoAB);

  std::cout << "Discord: " << D << endl;
  std::cout << "DW: " << DW << endl;
  std::cout << "DT: " << DT << endl;
}

void verify_7_2()
{
  DensityMatrix rhoAB(4, 4);

  rhoAB << 1, 0, 0, 0,
    0, 2, 0.1, 0,
    0, 0.1, 1, 0.5,
    0, 0, 0.5, 1;

  validateDensityMatrix(rhoAB);

  double D = compute_D(rhoAB);


  double DW = compute_DW(rhoAB);
  double Q = DW - mutualInformationCalc(rhoAB)/2 ;

  std::cout << "Discord: " << D << endl;
  std::cout << "DW: " << DW << endl;
  std::cout << "DeltaQ: " << Q << endl;
}

void verify_7_3()
{
  DensityMatrix rhoAB(4, 4);

  rhoAB << 1, 0, 0, 0,
    0, 2, 0.1, 0,
    0, 0.1, 1, 0.5,
    0, 0, 0.5, 1;

  validateDensityMatrix(rhoAB);

  double D = compute_D(rhoAB);


  double DW = compute_DW(rhoAB);
  double Q = DW - mutualInformationCalc(rhoAB) / 2;

  std::cout << "Discord: " << D << endl;
  std::cout << "HAC/2: " << Q << endl;
  std::cout << "HAC: " << 2*Q << endl;
}

int main() {
  verify_7_3();
  return 0;

}

