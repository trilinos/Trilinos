// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

// sacado_test
//
//  usage: 
//     sacado_test
//
//  output:  
//     prints a summary line and one line "Hello" for each process to standard out

#include "Sacado.hpp"

template <typename ScalarT>
void func(const ScalarT& a, const ScalarT& b) {
  ScalarT t1 = a*(b-1.);

  std::cout << t1 << std::endl;
}

int main(int argc, char **argv)
{

  Sacado::Taylor::DTaylor<double> a_dtaylor(3, 2.0);
  Sacado::Taylor::DTaylor<double> b_dtaylor(3, 3.0);
  func(a_dtaylor, b_dtaylor);

  return 0;
}

  

