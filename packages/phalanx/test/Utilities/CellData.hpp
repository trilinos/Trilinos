#ifndef PHX_EXAMPLE_CELL_DATA_HPP
#define PHX_EXAMPLE_CELL_DATA_HPP

#include "Phalanx_ConfigDefs.hpp" // for std::vector
#include "AlgebraicTypes.hpp"

class CellData {
  
public:

  CellData();
  
  virtual ~CellData() {}
  
  std::vector<double>& getBasisFunctions();
  
  std::vector< MyVector<double> >& getBasisFunctionGradients();
  
  std::vector< MyVector<double> >& getNodeCoordinates();
  
private:
  
  std::vector< MyVector<double> > coords_;
  
  std::vector<double> phi_;
  
  std::vector< MyVector<double> > grad_phi_;
};

#endif
