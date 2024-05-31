#include "methods/Methods.h"
#include "methods/DiscordCalc.h"

using namespace Methods;
using namespace std;
using namespace DiscordCalc;

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

  /*
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
  */

  DensityMatrix rhoAB(4, 4);
  rhoAB << 5, 2, 1, 3,
           2, 2, 1, 0,
           1, 1, 4, 2,
           3, 0, 2, 8;

  

  normalise_matrix(rhoAB);
  cout<<rhoAB<<endl;

  /*
  DensityMatrix rhoABC = canonicalPurification(rhoAB);

  DensityMatrix rhoBC = partialTraceFirstSubsystem(rhoABC, 4, 4);
  
  DensityMatrix rhoBC_Bp_Cp = canonicalPurification(rhoBC);

  cout << "Pure: " << std::boolalpha << isPure(rhoBC_Bp_Cp) << endl;
  cout<< rhoBC_Bp_Cp<<endl;

  cout<<computeEntropy(rhoBC_Bp_Cp);
 

  DensityMatrix rhoAB_retraced = Methods::partialTrace(rhoABC, { 4, 4, 4 }, 2);

  cout << "rhoAB after tracing out C from rhoABC:" << endl;
  cout << rhoAB_retraced << endl;
  */

  double mut_info = mutualInformationCalc(rhoAB, 2, 2);




  double res = J_PI_AB(rhoAB, 2, 2);
  cout<<mut_info<<endl;
  cout<<res<<endl;
  cout<<mut_info- res<<endl;

  return 0;
}

