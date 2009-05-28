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
//     prints the Hermite Polynomial Chaos Expansion of the simple function
//
//     w = sinh(1/(log(u)^2+1)
//
//     where u = 1 + 0.4*H_1(x) + 0.06*H_2(x) + 0.002*H_3(x), x is a zero-mean
//     and unit-variance Gaussian random variable, and H_i(x) is the i-th
//     Hermite polynomial.
//
//     The expansion for w is computed from an order of 0 to 10 with the 
//     corresponding mean and standard deviation.

#include <iostream>
#include <iomanip>

#include "Stokhos.hpp"

int main(int argc, char **argv)
{
  try {
    const int pmax = 10;
    std::vector<double> sd(pmax+1);
    std::vector<double> me(pmax+1);
    for (int p=0; p<=pmax; p++) {
      std::vector< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(1); 
      bases[0] = Teuchos::rcp(new Stokhos::HermiteBasis<int,double>(p));
      Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis = 
        Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
      Stokhos::DerivOrthogPolyExpansion<int,double> expn(basis);
      Stokhos::OrthogPolyApprox<int,double> u(p+1),v(p+1),w(p+1);
      u[0] = 1.0;
      if (p >= 1)
        u[1] = 0.4;
      if (p >= 2)
        u[2] = 0.06;
      if (p >= 3)
        u[3] = 0.002;
      
      expn.log(v,u);
      expn.times(w,v,v);
      expn.plusEqual(w,1.0);
      expn.divide(v,1.0,w);
      expn.sinh(w,v);
      
      std::cout.precision(16);
      std::cout << "w = " << w << std::endl;

      double mean = w[0];
      double std_dev = 0.0;
      const std::vector<double> nrm2 = basis->norm_squared();
      for (int i=1; i<basis->size(); i++)
        std_dev += w[i]*w[i]*nrm2[i];
      std_dev = std::sqrt(std_dev);

      std::cout << "Mean =      " << mean << std::endl;
      std::cout << "Std. Dev. = " << std_dev << std::endl;

      me[p] = mean;
      sd[p] = std_dev;
    }

    for (int p=0; p<=pmax; p++)
      std::cout << me[p] << "  " << sd[p] << std::endl;
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
