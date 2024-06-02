#include"methods/Utilities.h"
#include"methods/DensityMatrix.h"
#include<iostream>
#include<sstream>
#include<string>
#include<algorithm>
#include<iterator>
#include<set>

using Utilities::Matrix;


DensityMatrix::DensityMatrix(const Matrix& mat) : matrix(mat), numQubits(std::log2(mat.rows())) {
  Utilities::validateDensityMatrix(matrix);
}

DensityMatrix::DensityMatrix(const DensityMatrix& other) : numQubits(other.numQubits), matrix(other.matrix), subsystems(other.subsystems) {}

DensityMatrix& DensityMatrix::operator=(const DensityMatrix& other) {
    if (this == &other) {
      return *this;
    }
    numQubits = other.numQubits;
    matrix = other.matrix;
    subsystems = other.subsystems;
    return *this;
}

void DensityMatrix::normalise()
{
  Utilities::normalise_matrix(matrix);
}

void DensityMatrix::defineSubsystem(const std::string& name, const std::vector<int>& qubits) {
  subsystems[name] = qubits;
}

DensityMatrix DensityMatrix::partialTrace(std::string subsystem) const {
  if (subsystems.find(subsystem) == subsystems.end()) {throw std::invalid_argument("Cannot find subsystem");}
  return performPartialTrace(subsystem);
}

DensityMatrix DensityMatrix::partialTrace(std::vector<std::string> subsystem) const {
  DensityMatrix traced = *this;  

  for (const auto& subsys : subsystem) {
    if (subsystems.find(subsys) == subsystems.end()) {throw std::invalid_argument("Cannot find subsystem");}
    traced = traced.performPartialTrace(subsys); 
  }
  return traced;
}

DensityMatrix DensityMatrix::getSubsystem(std::vector<std::string> subsystem) const {
  std::set<std::string> vecSet(subsystem.begin(), subsystem.end());

  std::vector<std::string> difference;

  for (const auto& pair : subsystems) {
    if (vecSet.find(pair.first) == vecSet.end()) {
      difference.push_back(pair.first);
    }
  };
  return partialTrace(difference);
}


void DensityMatrix::print() const
{
  std::cout << matrix << std::endl;
}


std::map<std::string, std::vector<int>> DensityMatrix::edit_subsystem_definitions(int removed_qubit) const
{
  std::map<std::string, std::vector<int>> new_subSystems;
  for (const auto& pair : subsystems)
  {
    std::vector<int> new_qubits;
    for (const auto& qubit : pair.second)
    {
      if (qubit < removed_qubit)
      {
        new_qubits.push_back(qubit);
      }
      else if (qubit > removed_qubit)
      {
        new_qubits.push_back(qubit - 1);
      }

    }
    if (!new_qubits.empty()) { new_subSystems[pair.first] = new_qubits; }
  }
  return new_subSystems;
}

DensityMatrix DensityMatrix::performPartialTrace(const std::string& subsystem) const {
  DensityMatrix subsys = *this;
  while (subsys.subsystems.find(subsystem) != subsys.subsystems.end())
  {
    int qubit_to_remove = subsys.subsystems.at(subsystem)[0];
    subsys.matrix = Utilities::partialTrace(subsys.matrix, qubit_to_remove);
    subsys.subsystems = subsys.edit_subsystem_definitions(qubit_to_remove);
  }
  return subsys;
}

double DensityMatrix::computeEntropy() const {
  Eigen::SelfAdjointEigenSolver<Matrix> es(matrix);
  Eigen::VectorXd eigenvalues = es.eigenvalues().real();
  double entropy = 0.0;
  for (int i = 0; i < eigenvalues.size(); ++i) {
    double lambda = eigenvalues(i);
    if (lambda > 0) {
      entropy -= lambda * log2(lambda);
    }
  }
  return entropy;
}

void DensityMatrix::printSubSystems() const {
  for (const auto& pair : subsystems) {
    std::stringstream qubits;
    qubits<<pair.first<<": ";
    for (const auto& qubit : pair.second)
    {
      qubits<<qubit<<" ";
    }
    std::cout<<qubits.str()<<std::endl;
  }
}

double DensityMatrix::computeEntropy(std::string qubit) const {
  std::vector<std::string> qubits_to_remove;
  for (const auto& pair : subsystems)
  {
    if (pair.first != qubit) {
      qubits_to_remove.push_back(pair.first);
    }
  }
  DensityMatrix subSys = *this;
  for (const auto& qubit : qubits_to_remove)
  {
    subSys.partialTrace(qubit);
  }
  Eigen::SelfAdjointEigenSolver<Matrix> es(subSys.matrix);
  Eigen::VectorXd eigenvalues = es.eigenvalues().real();
  double entropy = 0.0;
  for (int i = 0; i < eigenvalues.size(); ++i) {
    double lambda = eigenvalues(i);
    if (lambda > 0) {
      entropy -= lambda * log2(lambda);
    }
  }
  return entropy;
}

DensityMatrix DensityMatrix::purify() const {
  Eigen::SelfAdjointEigenSolver<Matrix> es(matrix);
  Eigen::VectorXd eigenvalues = es.eigenvalues().real();
  Eigen::MatrixXcd eigenvectors = es.eigenvectors();
  int dim = static_cast<int>(eigenvalues.size());

  Matrix psi = Matrix::Zero(dim * dim, 1);
  for (int i = 0; i < dim; ++i) {
    if (eigenvalues(i) > 0) {
      for (int j = 0; j < dim; ++j) {
        psi(i * dim + j, 0) = sqrt(eigenvalues(i)) * eigenvectors(j, i);
      }
    }
  }
  /*
  
  Make key complex numbers
  Each prime gives an imaginary component
  If purifying over system already with primes demote to real numbers in the sequence
  Problem for multiple puficiations where one wishes to search for specific purification system?
  Meh

  */


  Matrix purifiedState = psi * psi.adjoint();
  DensityMatrix purifiedMatrix(purifiedState);

  std::map<std::string, std::vector<int>> subsystemsPrime = subsystems;
  for (const auto& pair : subsystems) {
    std::cout<<pair.first<<", "<<pair.second[0]<<std::endl;
    std::string sysPrime = pair.first + "'";
    subsystemsPrime[sysPrime] = pair.second;
    for(auto& elem : subsystemsPrime[pair.first]) {elem += numQubits;}
  }
  std::cout<<"\n"<<std::endl;
  purifiedMatrix.subsystems.insert(subsystemsPrime.begin(), subsystemsPrime.end());
  return purifiedMatrix;
}

bool DensityMatrix::isPure() const
{
  return Utilities::isPure(matrix);
}