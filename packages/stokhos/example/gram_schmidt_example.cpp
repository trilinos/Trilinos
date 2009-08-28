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
#include "Stokhos_StieltjesPCEBasis.hpp"
#include "Stokhos_GramSchmidtBasis.hpp"
#include "Stokhos_UserDefinedQuadrature.hpp"

typedef Stokhos::LegendreBasis<int,double> basis_type;

// This example demonstrates computing an orthogonal basis via Gram-Schmidt

struct stieltjes_pce_quad_func {
  stieltjes_pce_quad_func(
   const Stokhos::OrthogPolyApprox<int,double>& pce_,
   const Stokhos::OrthogPolyBasis<int,double>& basis_) :
    pce(pce_), basis(basis_), vec(2) {}
  
  double operator() (const double& a, double& b) const {
    vec[0] = a;
    vec[1] = b;
    return pce.evaluate(vec);
  }
  const Stokhos::OrthogPolyApprox<int,double>& pce;
  const Stokhos::OrthogPolyBasis<int,double>& basis;
  mutable Teuchos::Array<double> vec;
};

int main(int argc, char **argv)
{
  try {

    const unsigned int d = 1;
    const unsigned int p = 7;

    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d); 

    // Create product basis
    for (unsigned int i=0; i<d; i++)
      bases[i] = Teuchos::rcp(new basis_type(p));
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis = 
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));

    // Tensor product quadrature
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad = 
      Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));
    // Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad = 
    //   Teuchos::rcp(new Stokhos::SparseGridQuadrature<int,double>(basis,
    // 							     p+1));

    // Create approximation
    Stokhos::OrthogPolyApprox<int,double> x(basis), u(basis), v(basis), 
      w(basis), w2(basis);
    for (unsigned int i=0; i<d; i++) {
      x.term(i, 1) = 1.0;
    }
	  
    Teuchos::RCP< Stokhos::GramSchmidtBasis<int,double> > gs_basis = 
      Teuchos::rcp(new Stokhos::GramSchmidtBasis<int,double>(basis,
							     quad->getQuadPoints(),
							     quad->getQuadWeights()));

    std::cout << *gs_basis << std::endl;
	
    // User-defined quadrature
    Teuchos::RCP< const Teuchos::Array< Teuchos::Array<double> > > points =
      Teuchos::rcp(&(quad->getQuadPoints()),false);
    Teuchos::RCP< const Teuchos::Array<double> > weights =
      Teuchos::rcp(&(quad->getQuadWeights()),false);
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > gs_quad = 
      Teuchos::rcp(new Stokhos::UserDefinedQuadrature<int,double>(gs_basis,
								  points,
								  weights));
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
