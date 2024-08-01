// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_TimeMonitor.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::TensorProductQuadrature<ordinal_type, value_type>::
TensorProductQuadrature(const Teuchos::RCP<const ProductBasis<ordinal_type,value_type> >& product_basis) 
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::TensorProductQuadrature -- Quad Grid Generation");
#endif
  ordinal_type d = product_basis->dimension();
  ordinal_type sz = product_basis->size();

  Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type,value_type> > > coordinate_bases = product_basis->getCoordinateBases();

  // Compute quad points, weights, values
  Teuchos::Array< Teuchos::Array<value_type> > gp(d);
  Teuchos::Array< Teuchos::Array<value_type> > gw(d);
  Teuchos::Array< Teuchos::Array< Teuchos::Array<value_type> > > gv(d);
  Teuchos::Array<ordinal_type> n(d);
  ordinal_type ntot = 1;
  for (ordinal_type i=0; i<d; i++) {
    coordinate_bases[i]->getQuadPoints(2*(coordinate_bases[i]->order()), 
				       gp[i], gw[i], gv[i]);
    n[i] = gp[i].size();
    ntot *= n[i];
  }
  quad_points.resize(ntot);
  quad_weights.resize(ntot);
  quad_values.resize(ntot);
  Teuchos::Array<ordinal_type> index(d);
  for (ordinal_type i=0; i<d; i++)
    index[i] = 0;
  ordinal_type cnt = 0;
  while (cnt < ntot) {
    quad_points[cnt].resize(d);
    quad_weights[cnt] = value_type(1.0);
    quad_values[cnt].resize(sz);
    for (ordinal_type j=0; j<d; j++) {
      quad_points[cnt][j] = gp[j][index[j]];
      quad_weights[cnt] *= gw[j][index[j]];
    }
    for (ordinal_type k=0; k<sz; k++) {
      quad_values[cnt][k] = value_type(1.0);
      const MultiIndex<ordinal_type>& term = product_basis->term(k);
      for (ordinal_type j=0; j<d; j++) 
        quad_values[cnt][k] *= gv[j][index[j]][term[j]];
    }
    ++index[0];
    ordinal_type i = 0;
    while (i < d-1 && index[i] == n[i]) {
      index[i] = 0;
      ++i;
      ++index[i];
    }
    ++cnt;
  }

  //std::cout << "Number of quadrature points = " << ntot << std::endl;
}

template <typename ordinal_type, typename value_type>
Stokhos::TensorProductQuadrature<ordinal_type, value_type>::
TensorProductQuadrature(const Teuchos::RCP<const ProductBasis<ordinal_type,value_type> >& product_basis,const ordinal_type& quad_order) 
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::TensorProductQuadrature -- Quad Grid Generation");
#endif
  ordinal_type d = product_basis->dimension();
  ordinal_type sz = product_basis->size();

  Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type,value_type> > > coordinate_bases = product_basis->getCoordinateBases();

  // Compute quad points, weights, values
  Teuchos::Array< Teuchos::Array<value_type> > gp(d);
  Teuchos::Array< Teuchos::Array<value_type> > gw(d);
  Teuchos::Array< Teuchos::Array< Teuchos::Array<value_type> > > gv(d);
  Teuchos::Array<ordinal_type> n(d);
  ordinal_type ntot = 1;
  for (ordinal_type i=0; i<d; i++) {
    coordinate_bases[i]->getQuadPoints(quad_order, 
				       gp[i], gw[i], gv[i]);
    n[i] = gp[i].size();
    ntot *= n[i];
  }
  quad_points.resize(ntot);
  quad_weights.resize(ntot);
  quad_values.resize(ntot);
  Teuchos::Array<ordinal_type> index(d);
  for (ordinal_type i=0; i<d; i++)
    index[i] = 0;
  ordinal_type cnt = 0;
  while (cnt < ntot) {
    quad_points[cnt].resize(d);
    quad_weights[cnt] = value_type(1.0);
    quad_values[cnt].resize(sz);
    for (ordinal_type j=0; j<d; j++) {
      quad_points[cnt][j] = gp[j][index[j]];
      quad_weights[cnt] *= gw[j][index[j]];
    }
    for (ordinal_type k=0; k<sz; k++) {
      quad_values[cnt][k] = value_type(1.0);
      MultiIndex<ordinal_type> term = product_basis->term(k);
      for (ordinal_type j=0; j<d; j++) 
        quad_values[cnt][k] *= gv[j][index[j]][term[j]];
    }
    ++index[0];
    ordinal_type i = 0;
    while (i < d-1 && index[i] == n[i]) {
      index[i] = 0;
      ++i;
      ++index[i];
    }
    ++cnt;
  }

  //std::cout << "Number of quadrature points = " << ntot << std::endl;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array< Teuchos::Array<value_type> >&
Stokhos::TensorProductQuadrature<ordinal_type, value_type>::
getQuadPoints() const
{
  return quad_points;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array<value_type>&
Stokhos::TensorProductQuadrature<ordinal_type, value_type>::
getQuadWeights() const
{
  return quad_weights;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array< Teuchos::Array<value_type> >&
Stokhos::TensorProductQuadrature<ordinal_type, value_type>::
getBasisAtQuadPoints() const
{
  return quad_values;
}

template <typename ordinal_type, typename value_type>
std::ostream& 
Stokhos::TensorProductQuadrature<ordinal_type,value_type>::
print(std::ostream& os) const
{
  ordinal_type nqp = quad_weights.size();
  os << "Tensor Product Quadrature with " << nqp << " points:"
     << std::endl << "Weight : Points" << std::endl;
  for (ordinal_type i=0; i<nqp; i++) {
    os << i << ": " << quad_weights[i] << " : ";
    for (ordinal_type j=0; j<static_cast<ordinal_type>(quad_points[i].size()); 
	 j++)
      os << quad_points[i][j] << " ";
    os << std::endl;
  }
  os << "Basis values at quadrature points:" << std::endl;
  for (ordinal_type i=0; i<nqp; i++) {
    os << i << " " << ": ";
    for (ordinal_type j=0; j<static_cast<ordinal_type>(quad_values[i].size()); 
	 j++)
      os << quad_values[i][j] << " ";
    os << std::endl;
  }

  return os;
}
