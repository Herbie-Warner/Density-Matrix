//DiscordCalc.h
#include "methods/Methods.h"
#include <Eigen/Dense>

#ifndef DISCORDCALC_H
#define DISCORDCALC_H

namespace DiscordCalc {
  const double PI = 3.14159265358979323846;
  double mutualInformationCalc(const Methods::DensityMatrix& rhoAB, int dimA, int dimB);
  double J_AB_2Qubits(const Methods::DensityMatrix& rhoAB);
}

#endif