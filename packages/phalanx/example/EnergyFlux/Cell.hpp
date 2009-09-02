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

#ifndef PHX_EXAMPLE_CELL_HPP
#define PHX_EXAMPLE_CELL_HPP

#include <vector>
#include "Phalanx_ConfigDefs.hpp"
#include "AlgebraicTypes.hpp"

//! A representation of a finite element cell.  This is not a realistic element, but is meant to represent what an element would act like (the actual basis functions values and node coordinates are fake).
class MyCell {
  
public:

  MyCell();
  
  virtual ~MyCell() {}
  
  std::vector< MyVector<double> >& getNodeCoordinates();
  
  std::vector< std::vector<double> >& getBasisFunctions();
  
  std::vector< std::vector< MyVector<double> > >& getBasisFunctionGradients();

  std::size_t localIndex();

  void setLocalIndex(std::size_t index);

private:
  
  std::size_t local_index_;

  std::vector< MyVector<double> > coords_;
  
  std::vector< std::vector<double> > phi_;
  
  std::vector< std::vector< MyVector<double> > > grad_phi_;
};

#endif
