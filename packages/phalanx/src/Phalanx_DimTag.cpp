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

#include <sstream>
#include "Phalanx_DimTag.hpp"

//**********************************************************************
PHX::DimTag::~DimTag() {}

//**********************************************************************
std::string PHX::DimTag::to_string( size_t n , int i ) const
{
  //array_bounds_checking( n , i );
  std::ostringstream s ;
  s << i ;
  return s.str();
}

//**********************************************************************
int PHX::DimTag::to_index( size_t n , const std::string & s ) const
{
  const int i = atoi( s.c_str() );
  //array_bounds_checking( n , i );
  return i ;
}

//**********************************************************************
