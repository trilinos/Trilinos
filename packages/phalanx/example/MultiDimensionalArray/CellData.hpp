// @HEADER
// ************************************************************************
// 
//            Phalanx: A Partial Differential Equation Assembly 
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

#ifndef PHX_EXAMPLE_CELL_DATA_HPP
#define PHX_EXAMPLE_CELL_DATA_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Dimension.hpp"
#include "Phalanx_Array.hpp"

class CellData {
  
public:

  CellData();
  
  virtual ~CellData() {}
  
  phdmesh::ArrayNatural<double,Node,Dim>& getNodeCoordinates();
  
  phdmesh::ArrayNatural<double,QuadPoint,Node>& getBasisFunctions();
  
  phdmesh::ArrayNatural<double,QuadPoint,Node,Dim>& 
  getBasisFunctionGradients();
  
private:
  
  std::vector<double> m_coords_mem;
  
  std::vector<double> m_phi_mem;
  
  std::vector<double> m_grad_phi_mem;

  phdmesh::ArrayNatural<double,Node,Dim> m_coords;
  
  phdmesh::ArrayNatural<double,QuadPoint,Node> m_phi;

  phdmesh::ArrayNatural<double,QuadPoint,Node,Dim> m_grad_phi;

};

#endif
