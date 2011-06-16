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

#include "TriKota_ConfigDefs.hpp"
#include "sandia_sgmg.H"
#include "sandia_rules.H"
#include "pecos_global_defs.hpp"
#include "Teuchos_TimeMonitor.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::SparseGridQuadrature<ordinal_type, value_type>::
SparseGridQuadrature(
  const Teuchos::RCP<const ProductBasis<ordinal_type,value_type> >& product_basis,
  ordinal_type sparse_grid_level,
  value_type duplicate_tol) :
  coordinate_bases(product_basis->getCoordinateBases())
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Sparse Grid Generation");
#endif

  ordinal_type d = product_basis->dimension();
  ordinal_type p = product_basis->order();
  ordinal_type sz = product_basis->size();
  ordinal_type level = sparse_grid_level;

  // Mike's heuristic formula for computing the level
  if (level == 0) {
    level = static_cast<ordinal_type>(std::ceil(0.5*(p+d-1)));
    if (level < d)
      level = p;
  }

  std::cout << "Sparse grid level = " << level << std::endl;

  // Compute quad points, weights, values
  Teuchos::Array<int> rules(d), growthRules(d);
  
  Teuchos::Array< void (*) ( int order, int dim, double x[] ) > compute1DPoints(d);
  Teuchos::Array< void (*) ( int order, int dim, double w[] ) > compute1DWeights(d);
  for (ordinal_type i=0; i<d; i++) {
    compute1DPoints[i] = &(getMyPoints);
    compute1DWeights[i] = &(getMyWeights);
    rules[i] = coordinate_bases[i]->getSparseGridRule();
    growthRules[i] = coordinate_bases[i]->getSparseGridGrowthRule();
  }

  // Set the static sparse grid quadrature pointer to this
  // (this will cause a conflict if another sparse grid quadrature object
  // is trying to access the VPISparseGrid library, but that's all we can
  // do at this point).
  sgq = this;

  int num_total_pts  =
    webbur::sgmg_size_total(d, level, &rules[0], &growthRules[0]);
  ordinal_type ntot =
    webbur::sgmg_size(d, level, &rules[0], &compute1DPoints[0], 1e-15,
		      &growthRules[0]);
  Teuchos::Array<int> sparse_order(ntot*d);
  Teuchos::Array<int> sparse_index(ntot*d);
  Teuchos::Array<int> sparse_unique_index(num_total_pts);
  quad_points.resize(ntot);
  quad_weights.resize(ntot);
  quad_values.resize(ntot);
  Teuchos::Array<value_type> gp(ntot*d);

  webbur::sgmg_unique_index(d, level, &rules[0], &compute1DPoints[0],
			    duplicate_tol, ntot, num_total_pts, 
			    &growthRules[0],
			    &sparse_unique_index[0]);
  webbur::sgmg_index(d, level, &rules[0], ntot, num_total_pts, 
		     &sparse_unique_index[0],
		     &growthRules[0], 
		     &sparse_order[0], &sparse_index[0]);
  webbur::sgmg_weight(d, level, &rules[0], &compute1DWeights[0],
		      ntot, num_total_pts, 
		      &sparse_unique_index[0], 
		      &growthRules[0],
		      &quad_weights[0]);
  webbur::sgmg_point(d, level, &rules[0], &compute1DPoints[0],
		     ntot, &sparse_order[0], 
		     &sparse_index[0], 
		     &growthRules[0], 
		     &gp[0]);
  
  for (ordinal_type i=0; i<ntot; i++) {
    quad_values[i].resize(sz);
    quad_points[i].resize(d);
    for (ordinal_type j=0; j<d; j++)
      quad_points[i][j] = gp[i*d+j];
    product_basis->evaluateBases(quad_points[i], quad_values[i]);
  }

  std::cout << "Number of quadrature points = " << ntot << std::endl;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array< Teuchos::Array<value_type> >&
Stokhos::SparseGridQuadrature<ordinal_type, value_type>::
getQuadPoints() const
{
  return quad_points;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array<value_type>&
Stokhos::SparseGridQuadrature<ordinal_type, value_type>::
getQuadWeights() const
{
  return quad_weights;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array< Teuchos::Array<value_type> >&
Stokhos::SparseGridQuadrature<ordinal_type, value_type>::
getBasisAtQuadPoints() const
{
  return quad_values;
}

template <typename ordinal_type, typename value_type>
void 
Stokhos::SparseGridQuadrature<ordinal_type,value_type>::
getMyPoints( int order, int dim, double x[] )
{
  if (sgq->coordinate_bases[dim]->getSparseGridRule() == 1)
    webbur::clenshaw_curtis_compute_points(order, x);
  else {
    Teuchos::Array<double> quad_points;
    Teuchos::Array<double> quad_weights;
    Teuchos::Array< Teuchos::Array<double> > quad_values;
    sgq->coordinate_bases[dim]->getQuadPoints(2*order-1, quad_points, 
					      quad_weights, quad_values);
    for (int i = 0; i<quad_points.size(); i++) {
      x[i] = quad_points[i];
    }
  }
}

template <typename ordinal_type, typename value_type>
void 
Stokhos::SparseGridQuadrature<ordinal_type,value_type>::
getMyWeights( int order, int dim, double w[] )
{
  if (sgq->coordinate_bases[dim]->getSparseGridRule() == 1) {
    webbur::clenshaw_curtis_compute_weights(order, w);
    for (int i=0; i<order; i++)
      w[i] *= 0.5;
  }
  else {
    Teuchos::Array<double> quad_points;
    Teuchos::Array<double> quad_weights;
    Teuchos::Array< Teuchos::Array<double> > quad_values;
    sgq->coordinate_bases[dim]->getQuadPoints(2*order-1, quad_points, 
					      quad_weights, quad_values);
    for (int i = 0; i<quad_points.size(); i++) {
      w[i] = quad_weights[i];
    }
  }
}

template <typename ordinal_type, typename value_type>
std::ostream& 
Stokhos::SparseGridQuadrature<ordinal_type,value_type>::
print(std::ostream& os) const
{
  ordinal_type nqp = quad_weights.size();
  os << "Sparse Grid Quadrature with " << nqp << " points:"
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
