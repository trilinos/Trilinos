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

// pce_example
//
//  usage: 
//     pce_example
//
//  output:  
//     prints the Hermite Polynomial Chaos Expansion of a simple function
//     using the class Sacado::PCE::Hermite

#include <iostream>

#include "Sacado.hpp"
#include "Sacado_PCE_OrthogPoly.hpp"
#include "Stokhos.hpp"

int main(int argc, char **argv)
{
  try {
    const unsigned int d = 7;
    typedef Stokhos::HermiteBasis<int,double> basis_type;
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(1); 
    bases[0] = Teuchos::rcp(new basis_type(d));
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
    Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expansion = 
      Teuchos::rcp(new Stokhos::DerivOrthogPolyExpansion<int,double>(basis));
    Sacado::PCE::OrthogPoly<double> u(expansion);
    u.fastAccessCoeff(0) = 1.0;
    u.fastAccessCoeff(1) = 0.4;
    u.fastAccessCoeff(2) = 0.06;
    u.fastAccessCoeff(3) = 0.002;

    Sacado::PCE::OrthogPoly<double> w = std::log(u);
    Sacado::PCE::OrthogPoly<double> v = std::sinh(1.0/(w*w + 1.0));

    std::cout << "u (hermite basis) = " << u;
    std::cout << "v (hermite basis) = " << v;

    std::cout << "\nExample passed!" << std::endl;
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
