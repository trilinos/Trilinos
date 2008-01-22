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

// hermite_example
//
//  usage: 
//     hermite_example
//
//  output:  
//     prints the Hermite Polynomial Chaos Expansion of a simple function

#include <iostream>
#include <iomanip>

#include "Stokhos_HermiteExpansion.hpp"
#include "Stokhos_HermiteEBasis.hpp"

int main(int argc, char **argv)
{
  try {
    const unsigned int d = 7;
    Stokhos::HermiteExpansion< double,Stokhos::HermiteEBasis<double> > he(d);
    Stokhos::HermitePoly<double> u(d), v(d), w(d);
    u[0] = 1.0;
    u[1] = 0.4;
    u[2] = 0.06;
    u[3] = 0.002;

    he.log(v,u);
    he.times(w,v,v);
    he.plusEqual(w,1.0);
    he.divide(v,1.0,w);

    std::cout << "u (hermite basis) = " << u << std::endl;
    std::cout << "v (hermite basis) = " << v << std::endl;
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
