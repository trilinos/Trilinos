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

#ifndef PHX_EXAMPLE_MYCELL_HPP
#define PHX_EXAMPLE_MYCELL_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Dimension.hpp"
#include "Phalanx_Array.hpp"

class MyCell {
  
public:

  MyCell();
  
  virtual ~MyCell() {}
  
  PHX::Array<double,PHX::NaturalOrder,Node,Dim>& getNodeCoordinates();
  
  PHX::Array<double,PHX::NaturalOrder,QuadPoint,Node>& 
  getBasisFunctions();
  
  PHX::Array<double,PHX::NaturalOrder,QuadPoint,Node,Dim>& 
  getBasisFunctionGradients();
  
  std::size_t localIndex();

  void setLocalIndex(std::size_t index);

private:
  
  std::size_t local_index_;

  Teuchos::ArrayRCP<double> m_coords_mem;
  
  Teuchos::ArrayRCP<double> m_phi_mem;
  
  Teuchos::ArrayRCP<double> m_grad_phi_mem;

  PHX::Array<double,PHX::NaturalOrder,Node,Dim> m_coords;
  
  PHX::Array<double,PHX::NaturalOrder,QuadPoint,Node> m_phi;

  PHX::Array<double,PHX::NaturalOrder,QuadPoint,Node,Dim> m_grad_phi;

};

#endif
