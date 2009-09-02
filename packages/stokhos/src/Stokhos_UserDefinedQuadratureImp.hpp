// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2008) Sandia Corporation
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

template <typename ordinal_type, typename value_type>
Stokhos::UserDefinedQuadrature<ordinal_type, value_type>::
UserDefinedQuadrature(
const Teuchos::RCP<const OrthogPolyBasis<ordinal_type,value_type> >& basis,
const Teuchos::RCP<const Teuchos::Array< Teuchos::Array<value_type> > >& points,
     const Teuchos::RCP<const Teuchos::Array<value_type> >& weights)
  : quad_points(points),
    quad_weights(weights)
{
  ordinal_type nqp = points->size();
  Teuchos::RCP<Teuchos::Array< Teuchos::Array<value_type> > > qv = 
    Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(nqp));
  for (ordinal_type i=0; i<nqp; i++) {
    (*qv)[i].resize(basis->size());
    basis->evaluateBases((*points)[i], (*qv)[i]);
  }
  quad_values = qv;
}

template <typename ordinal_type, typename value_type>
Stokhos::UserDefinedQuadrature<ordinal_type, value_type>::
UserDefinedQuadrature(
const Teuchos::RCP<const Teuchos::Array< Teuchos::Array<value_type> > >& points,
const Teuchos::RCP<const Teuchos::Array<value_type> >& weights,
const Teuchos::RCP<const Teuchos::Array< Teuchos::Array<value_type> > >& values)
  : quad_points(points),
    quad_weights(weights),
    quad_values(values)
{
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array< Teuchos::Array<value_type> >&
Stokhos::UserDefinedQuadrature<ordinal_type, value_type>::
getQuadPoints() const
{
  return *quad_points;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array<value_type>&
Stokhos::UserDefinedQuadrature<ordinal_type, value_type>::
getQuadWeights() const
{
  return *quad_weights;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array< Teuchos::Array<value_type> >&
Stokhos::UserDefinedQuadrature<ordinal_type, value_type>::
getBasisAtQuadPoints() const
{
  return *quad_values;
}
