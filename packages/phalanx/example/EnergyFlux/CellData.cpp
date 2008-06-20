#include "CellData.hpp"

//**********************************************************************
CellData::CellData() :
  phi_(4, 0.25),
  grad_phi_(4, 0.25)
{ }

//**********************************************************************
std::vector< MyVector<double> >& CellData::getNodeCoordinates()
{
  return coords_;
}

//**********************************************************************
std::vector<double>& CellData::getBasisFunctions()
{
  return phi_;
}

//**********************************************************************
std::vector< MyVector<double> >& CellData::getBasisFunctionGradients()
{
  return grad_phi_;
}

//**********************************************************************
