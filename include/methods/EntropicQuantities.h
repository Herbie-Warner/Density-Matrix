//DiscordCalc.h
#include "methods/Utilities.h"
#include "methods/DensityMatrix.h"
#include <Eigen/Dense>

#ifndef ENTROPICQUANTITIES_H
#define ENTROPICQUANTITIES_H

namespace EntropicQuantities {
  const double PI = 3.14159265358979323846;
  /*
  double mutualInformationCalc(const Utilities::Matrix& rhoAB);
  double J_AB_2Qubits(const Utilities::Matrix& rhoAB);
  double compute_DW(const Utilities::Matrix& rhoAB);
  double compute_D(const Utilities::Matrix& rhoAB);
  double compute_DT(const Utilities::Matrix& rhoAB);
  double compute_MarkovGap(const Utilities::Matrix& rhoAB);
  */
  double ReflectedEntropy(const DensityMatrix& rho, std::vector<std::string> systemA, std::vector<std::string> systemB);
}


#endif
