//DiscordCalc.h
#include "methods/Methods.h"
#include <Eigen/Dense>

#ifndef DISCORDCALC_H
#define DISCORDCALC_H

namespace DiscordCalc {
  const double PI = 3.14159265358979323846;
  double mutualInformationCalc(const Methods::DensityMatrix& rhoAB);
  double J_AB_2Qubits(const Methods::DensityMatrix& rhoAB);
  double compute_DW(const Methods::DensityMatrix& rhoAB);
  double compute_D(const Methods::DensityMatrix& rhoAB);
  double compute_DT(const Methods::DensityMatrix& rhoAB);
  double compute_MarkovGap(const Methods::DensityMatrix& rhoAB);
}

#endif