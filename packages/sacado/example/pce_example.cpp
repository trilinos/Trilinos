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
#include <iomanip>
#include <math.h>

#include "Sacado.hpp"
#include "Sacado_PCE_Hermite.hpp"

#ifdef HAVE_SACADO_STOKHOS
#include "Sacado_PCE_OrthogPoly.hpp"
#include "Stokhos_HermiteEBasis2.hpp"
#endif

// The function to differentiate
template <typename ScalarT>
ScalarT func(const ScalarT& a) {
  ScalarT r = std::exp(a);

  return r;
}

int main(int argc, char **argv)
{
  try {
    const unsigned int d = 7;
    Sacado::PCE::Hermite<double>::initWorkspace(d);
    Sacado::PCE::Hermite<double> u(d);
    u.fastAccessCoeff(0) = 1.0;
    u.fastAccessCoeff(1) = 0.4;
    u.fastAccessCoeff(2) = 0.06;
    u.fastAccessCoeff(3) = 0.002;

    Sacado::PCE::Hermite<double> v = std::log(u);
    //Sacado::PCE::Hermite<double> v = std::sinh(1.0/(w*w + 1.0));

    std::cout.precision(12);
    std::cout << "u (hermite basis) = " << u << std::endl;
    std::cout << "u (standard basis) = " << u.toStandardBasis() << std::endl;

    std::cout << "v (hermite basis) = " << v << std::endl;
    std::cout << "v (standard basis) = " << v.toStandardBasis() << std::endl;

    Sacado::PCE::StandardPoly<double> us = u.toStandardBasis();
    Sacado::Tay::Taylor<double> ut(d);
    for (unsigned int i=0; i<=d; i++)
      ut.fastAccessCoeff(i) = us.coeff(i);

    Sacado::Tay::Taylor<double> vt = std::log(ut);
    //Sacado::Tay::Taylor<double> vt = std::sinh(1.0/(wt*wt + 1.0));

    std::cout.precision(12);
    std::cout << "u (taylor basis) = " << ut << std::endl;
    std::cout << "v (taylor basis) = " << vt << std::endl;

#ifdef HAVE_SACADO_STOKHOS
    typedef Stokhos::HermiteEBasis2<double> basis_type;
    typedef Sacado::PCE::OrthogPoly<double>::expansion_type expansion_type;
    Teuchos::RCP<basis_type> basis = Teuchos::rcp(new basis_type(d));
    Teuchos::RCP<expansion_type> expansion = 
      Teuchos::rcp(new expansion_type(basis));
    Sacado::PCE::OrthogPoly<double>::initExpansion(expansion);
    Sacado::PCE::OrthogPoly<double> ue(d+1);
    for (unsigned int i=0; i<=d; i++)
      ue.fastAccessCoeff(i) = u.coeff(i);

    Sacado::PCE::OrthogPoly<double> ve = std::log(ue);
    //Sacado::PCE::OrthogPoly<double> ve = std::sinh(1.0/(we*we + 1.0));

    std::cout << "ue (hermite basis) = " << ue << std::endl;
    std::cout << "ue (standard basis) = " << ue.toStandardBasis() << std::endl;
    std::cout << "ve (hermite basis) = " << ve << std::endl;
    std::cout << "ve (standard basis) = " << ve.toStandardBasis() << std::endl;
#endif
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
