#include<Eigen/Dense>
#include<complex>
#include<vector>
#include<map>
#include "methods/Utilities.h"

#ifndef DENSITYMATRIX_H
#define DENSITYMATRIX_H

class DensityMatrix {
public:
  DensityMatrix(const Utilities::Matrix& mat);
  ~DensityMatrix() {}
  DensityMatrix(const DensityMatrix& other);
  DensityMatrix& operator=(const DensityMatrix& other); 
  void normalise();
  void defineSubsystem(const std::string& name, const std::vector<int>& qubits); //Error
  DensityMatrix partialTrace(std::string subsystem) const;
  DensityMatrix partialTrace(std::vector<std::string> subsystem) const;
  DensityMatrix getSubsystem(std::vector<std::string> subsystem) const;
  Utilities::Matrix getDensityMatrix() const {return matrix;}
  void print() const;
  void printSubSystems() const;
  double computeEntropy() const;
  double computeEntropy(std::string system) const;
  DensityMatrix purify() const;
  bool isPure() const;
private:
  int numQubits;
  Utilities::Matrix matrix;
  std::map<std::string, std::vector<int>> subsystems;
  std::map<std::string, std::vector<int>> edit_subsystem_definitions(int removed_qubit) const;
  DensityMatrix performPartialTrace(const std::string& subsystem) const; 
};

#endif