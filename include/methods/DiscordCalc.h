//DiscordCalc.h
#include "methods/Methods.h"
#include <Eigen/Dense>

#ifndef DISCORDCALC_H
#define DISCORDCALC_H

namespace DiscordCalc {
  double mutualInformationCalc(const Methods::DensityMatrix& rhoAB, int dimA, int dimB);
  double J_PI_AB_PauliX(const Methods::DensityMatrix& rhoAB, int dimA, int dimB);
  double J_PI_AB(const Methods::DensityMatrix& rhoAB, int dimA, int dimB);
}

#endif