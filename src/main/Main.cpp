#include<iostream>
#include"methods/DensityMatrix.h"
#include"methods/Utilities.h"
#include"methods/EntropicQuantities.h"

using Utilities::Matrix;

int main() {
  //verify_7_3();
  
  
  Matrix rhoA(2, 2);
  rhoA << 1, 0, 0, 1;

  Matrix rhoB(2, 2);
  rhoB << 1, 0.1, 0.1, 1;

  /*
  Matrix rhoC(2,2);
  rhoC << 1,0.5,0.5,1;
  */

  auto rhoAB = Utilities::tensorProduct(rhoA, rhoB);
  //auto rhoABC = Utilities::tensorProduct(rhoAB, rhoC);


  DensityMatrix dm(rhoAB);


  dm.defineSubsystem("A", {1}); 
  dm.defineSubsystem("B", {2});
  //dm.defineSubsystem("C", {3});
  
  double reflectedEntropy = EntropicQuantities::ReflectedEntropy(dm, { "A" }, { "B" });
  std::cout<<reflectedEntropy<<std::endl;

 
  return 0;

}

