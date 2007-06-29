// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

// taylor_example
//
//  usage: 
//     taylor_example
//
//  output:  
//     prints a summary line and one line "Hello" for each process to standard out

#include "Sacado.hpp"

template <typename ScalarT>
void func(const ScalarT& a, const ScalarT& b) {
  //ScalarT t1 = a*(b-1.);
  ScalarT t1 = asin(a);

  std::cout << t1 << std::endl;
}

int main(int argc, char **argv)
{

  Sacado::Tay::CacheTaylor<double> a_dtaylor(3, 0.5);
  a_dtaylor.fastAccessCoeff(1) = 0.5;
  a_dtaylor.fastAccessCoeff(2) = 0.5;
  a_dtaylor.fastAccessCoeff(3) = 0.5;
  Sacado::Tay::CacheTaylor<double> b_dtaylor(3, 3.0);
  func(a_dtaylor, b_dtaylor);

  Sacado::Tay::Taylor<double> a_staylor(3, 0.5);
  a_staylor.fastAccessCoeff(1) = 0.5;
  a_staylor.fastAccessCoeff(2) = 0.5;
  a_staylor.fastAccessCoeff(3) = 0.5;
  Sacado::Tay::Taylor<double> b_staylor(3, 3.0);
  func(a_staylor, b_staylor);

  return 0;
}

  

