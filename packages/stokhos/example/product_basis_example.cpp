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

#include "Stokhos_OrthogPolyExpansion.hpp"
#include "Stokhos_TayOrthogPolyExpansion.hpp"
#include "Stokhos_HermiteEBasis.hpp"
#include "Stokhos_HermiteEBasis2.hpp"
#include "Stokhos_HermiteBasis.hpp"
#include "Stokhos_UnitHermiteBasis.hpp"
#include "Stokhos_LegendreBasis.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"

typedef Stokhos::HermiteEBasis2<double> basis_type;
//typedef Stokhos::LegendreBasis<double> basis_type;

int main(int argc, char **argv)
{
  try {
    const unsigned int d = 1;
    const unsigned int p = 4;
    std::vector< Teuchos::RCP<const Stokhos::OrthogPolyBasis<double> > > bases(d); 
    std::vector<double> deriv_coeffs(d);
    for (unsigned int i=0; i<d; i++) {
      bases[i] = Teuchos::rcp(new basis_type(p));
      deriv_coeffs[i] = 1.0;
    }
    Teuchos::RCP< Stokhos::CompletePolynomialBasis<double> > basis = 
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<double>(bases,
								deriv_coeffs));
//     Teuchos::RCP<basis_type> basis = Teuchos::rcp(new basis_type(p));
    std::cout << *basis << std::endl;
    Stokhos::TayOrthogPolyExpansion<double> he(basis);
    unsigned int sz = basis->size();
    Stokhos::OrthogPolyApprox<double> u(sz),v(sz),w(sz),x(sz);
    u.term(*basis, 0,0) = 1.0;

    u.term(*basis, 1,0) = 0.4;
    u.term(*basis, 2,0) = 0.06;
    u.term(*basis, 3,0) = 0.002;

//     u.term(*basis, 0,1) = 0.4;
//     u.term(*basis, 0,2) = 0.06;
//     u.term(*basis, 0,3) = 0.002;

//     u.term(*basis, 0,0) = 0.5;
//     u.term(*basis, 1,0) = 0.05;
//     u.term(*basis, 0,1) = 0.05;

    std::cout << "u = " << std::endl;
    u.print(*basis, std::cout);

//     he.exp(v,u);
//     he.exp(w,v);
//     he.log(x,v);
    he.log(v,u);

//     he.times(w,v,v);
//     he.plusEqual(w,1.0);
//     he.divide(v,1.0,w);
//     he.sinh(w,v);

    Stokhos::OrthogPolyApprox<double> du(sz),dv(sz);

//     he.derivative(du,u);
//     he.derivative(dv,v);
//     he.times(w,v,du);

    std::cout.precision(12);
    std::cout << "v = " << std::endl;
    v.print(*basis, std::cout);

//     std::cout << "w = " << std::endl;
//     w.print(*basis, std::cout);

//     std::cout << "x = " << std::endl;
//     x.print(*basis, std::cout);

//     std::cout.precision(12);
//     std::cout << "du = " << std::endl;
//     du.print(*basis, std::cout);

//     std::cout.precision(12);
//     std::cout << "dv = " << std::endl;
//     dv.print(*basis, std::cout);

//     std::cout.precision(12);
//     std::cout << "v*du = " << std::endl;
//     w.print(*basis, std::cout);
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
