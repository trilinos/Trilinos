// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

template <typename ordinal_type, typename value_type>
std::ostream& 
Stokhos::UserDefinedQuadrature<ordinal_type,value_type>::
print(std::ostream& os) const
{
  ordinal_type nqp = quad_weights->size();
  os << "Sparse Grid Quadrature with " << nqp << " points:"
     << std::endl << "Weight : Points" << std::endl;
  for (ordinal_type i=0; i<nqp; i++) {
    os << i << ": " << (*quad_weights)[i] << " : ";
    for (ordinal_type j=0; 
	 j<static_cast<ordinal_type>((*quad_points)[i].size()); 
	 j++)
      os << (*quad_points)[i][j] << " ";
    os << std::endl;
  }
  os << "Basis values at quadrature points:" << std::endl;
  for (ordinal_type i=0; i<nqp; i++) {
    os << i << " " << ": ";
    for (ordinal_type j=0; 
	 j<static_cast<ordinal_type>((*quad_values)[i].size()); 
	 j++)
      os << (*quad_values)[i][j] << " ";
    os << std::endl;
  }

  return os;
}
