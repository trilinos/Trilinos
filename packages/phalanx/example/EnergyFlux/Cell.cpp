// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#include "Cell.hpp"

//**********************************************************************
MyCell::MyCell() :
  phi_(4),
  grad_phi_(4)
{ 
  for (std::size_t i=0; i < phi_.size(); ++i)
    phi_[i].resize(4,0.25);

  for (std::size_t i=0; i < grad_phi_.size(); ++i)
    grad_phi_[i].resize(4,MyVector<double>(0.25,0.25,0.25));
}

//**********************************************************************
std::vector< MyVector<double> >& MyCell::getNodeCoordinates()
{
  return coords_;
}

//**********************************************************************
std::vector< std::vector<double> >& MyCell::getBasisFunctions()
{
  return phi_;
}

//**********************************************************************
std::vector< std::vector< MyVector<double> > >& 
MyCell::getBasisFunctionGradients()
{
  return grad_phi_;
}

//**********************************************************************
std::size_t MyCell::localIndex()
{
  return local_index_;
}

//**********************************************************************
void MyCell::setLocalIndex(std::size_t index)
{
  local_index_ = index;
}

//**********************************************************************
