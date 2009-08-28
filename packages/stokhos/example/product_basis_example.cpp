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

#include <iostream>
#include <iomanip>

#include "Stokhos.hpp"

// This example demonstrates computing a multi-variate polynomial chaos
// expansion using a product basis.

int main(int argc, char **argv)
{
  try {
    const int d = 3;
    const int p = 5;
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d); 
    for (int i=0; i<d; i++) {
      bases[i] = Teuchos::rcp(new Stokhos::HermiteBasis<int,double>(p));
    }
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis = 
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad = 
        Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));
    Stokhos::QuadOrthogPolyExpansion<int,double> expn(basis, quad);
    Stokhos::OrthogPolyApprox<int,double> u(basis), v(basis), w(basis);
    u.term(0,0) = 1.0;

    for (int i=0; i<d; i++) {
      u.term(i,1) = 0.4 / d;
      u.term(i,2) = 0.06 / d;
      u.term(i,3) = 0.002 / d;
    }
    u.print(std::cout);

    expn.log(v,u);
    expn.times(w,v,v);
    expn.plusEqual(w,1.0);
    expn.divide(v,1.0,w);

    std::cout.precision(16);
    v.print(std::cout);

    double mean = v[0];
    double std_dev = 0.0;
    const Teuchos::Array<double> nrm2 = basis->norm_squared();
    for (int i=1; i<basis->size(); i++)
      std_dev += v[i]*v[i]*nrm2[i];
    std_dev = std::sqrt(std_dev);

    std::cout << "Mean =      " << mean << std::endl;
    std::cout << "Std. Dev. = " << std_dev << std::endl;

  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
