// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Stokhos_SparseGridQuadrature.hpp"

#ifdef HAVE_STOKHOS_DAKOTA

void 
Stokhos::getMyPoints( int order, int np, double p[], double x[] ){

  long int pointer = p[0];
  const OneDOrthogPolyBasis<int,double>* basis = (const OneDOrthogPolyBasis<int,double>*) pointer;
  Teuchos::Array<double> quad_points;
  Teuchos::Array<double> quad_weights;
  Teuchos::Array< Teuchos::Array<double> > quad_values;
  basis->getQuadPoints(2*order-1, quad_points, quad_weights, quad_values);
  for(std::size_t i = 0; i<quad_points.size(); i++){
    x[i] = quad_points[i];
  }
  
}


void 
Stokhos::getMyWeights( int order, int np, double p[], double w[] ){

  long int pointer = p[0];
  const OneDOrthogPolyBasis<int,double>* basis = (const OneDOrthogPolyBasis<int,double>*) pointer;
  Teuchos::Array<double> quad_points;
  Teuchos::Array<double> quad_weights;
  Teuchos::Array< Teuchos::Array<double> > quad_values;
  basis->getQuadPoints(2*order-1, quad_points, quad_weights, quad_values);

  for(std::size_t i = 0; i<quad_points.size(); i++){
    w[i] = quad_weights[i];
  }

}

#endif // HAVE_STOKHOS_DAKOTA
