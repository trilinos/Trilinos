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
  m_phi_mem(Teuchos::arcp<double>(4*4)),
  m_grad_phi_mem(Teuchos::arcp<double>(4*4*3)),
  m_phi(&(m_phi_mem[0]),4,4),
  m_grad_phi(&(m_grad_phi_mem[0]),4,4,3)
{ 
  using namespace Teuchos;

  for (ArrayRCP<double>::Ordinal i = 0; i < m_phi_mem.size(); ++i)
    m_phi_mem[i] = 0.25;

  for (ArrayRCP<double>::Ordinal i = 0; i < m_grad_phi_mem.size(); ++i)
    m_grad_phi_mem[i] = 0.25;
}

//**********************************************************************
PHX::Array<double,PHX::NaturalOrder,Node,Dim>& 
MyCell::getNodeCoordinates()
{
  return m_coords;
}

//**********************************************************************
PHX::Array<double,PHX::NaturalOrder,QuadPoint,Node>& 
MyCell::getBasisFunctions()
{
  return m_phi;
}

//**********************************************************************
PHX::Array<double,PHX::NaturalOrder,QuadPoint,Node,Dim>& 
MyCell::getBasisFunctionGradients()
{
  return m_grad_phi;
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
