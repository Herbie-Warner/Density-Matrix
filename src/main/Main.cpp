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

  double mut_info = mutualInformationCalc(rhoAB, 2, 2);
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

double compute_DW(const DensityMatrix& rhoAB)
{
  DensityMatrix rhoB = partialTrace(rhoAB, 1);
  double dW = computeEntropy(rhoB) - computeEntropy(rhoAB) + ReflectedEntropy(rhoAB) / 2;
  return dW;
}

double compute_D(const DensityMatrix& rhoAB)
{
  double res = J_AB_2Qubits(rhoAB);
  double mut_info = mutualInformationCalc(rhoAB, 2, 2);
  double D = mut_info - res;
  return D;
}

double compute_DT(const DensityMatrix& rhoAB)
{
  DensityMatrix rhoCAB = canonicalPurification(rhoAB);
  DensityMatrix rhoCB = partialTrace(rhoCAB,3);
  DensityMatrix rhoCpBpCB = canonicalPurification(rhoCB);
  DensityMatrix rhoCpBpB = partialTrace(partialTrace(rhoCpBpCB, 5), 4); 
  DensityMatrix BpB = partialTrace(partialTrace(rhoCpBpB, 1), 1);
  double S_BBp = computeEntropy(BpB);
  double mut_info = mutualInformationCalc(rhoAB, 2, 2);
  return mut_info - S_BBp;
}

int main() {
  DensityMatrix rhoAB(4, 4);

  rhoAB << 1, 0, 0, 0,
    0, 2, 0.1, 0,
    0, 0.1, 1, 0.5,
    0, 0, 0.5, 1;

  normalise_matrix(rhoAB);
  validateDensityMatrix(rhoAB);

  double D = compute_D(rhoAB);
  double DW = compute_DW(rhoAB);
  double DT = compute_DT(rhoAB);

  std::cout << "Discord: " << D << endl;
  std::cout << "DW: " << DW << endl;
  std::cout << "DT: " << DT << endl;

  return 0;

}

