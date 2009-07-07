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

#include "sparse_grid_mixed_growth.H"

// namespace webbur {
//   int sparse_grid_mixed_size ( int dim_num, int level_max, int rule[] );
//   void sparse_grid_mixed_point ( int dim_num, int level_max, int rule[], 
//                                  int point_num, double sparse_point[] );
//   void sparse_grid_mixed_weight ( int dim_num, int level_max, int rule[], 
//                                   int point_num, double sparse_weight[] );
// }

template <typename ordinal_type, typename value_type>
Stokhos::SparseGridQuadrature<ordinal_type, value_type>::
SparseGridQuadrature(Teuchos::RCP<const OrthogPolyBasis<ordinal_type,value_type> >& product_basis,
		     ordinal_type sparse_grid_level) 
{
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

  const std::vector< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type,value_type> > >& coordinate_bases = product_basis->getCoordinateBases();

  // Compute quad points, weights, values
  std::vector<int> rules(d);
  std::vector< void (*) ( int order, int np, double p[], double x[] ) > compute1DPoints(d);
  std::vector< void (*) ( int order, int np, double p[], double w[] ) > compute1DWeights(d);
  for (ordinal_type i=0; i<d; i++) {
    rules[i] = coordinate_bases[i]->getRule();
    if (rules[i] == 5) {
      compute1DPoints[i] = webbur::hermite_compute_points_np;
      compute1DWeights[i] = webbur::hermite_compute_weights_np;
    }
    else if (rules[i] == 1) {
      compute1DPoints[i] = webbur::clenshaw_curtis_compute_points_np;
      compute1DWeights[i] = webbur::clenshaw_curtis_compute_weights_np;
    }
  }
  std::vector<int> nparams(d);
  std::vector<double> params(d);
  int num_total_pts  =
    webbur::sparse_grid_mixed_growth_size_total(d, level, &rules[0],
						webbur::level_to_order_default);
  ordinal_type ntot =
    webbur::sparse_grid_mixed_growth_size(d, level, &rules[0], 
					  &nparams[0], &params[0], 
					  &compute1DPoints[0],
					  1e-15,
					  webbur::level_to_order_default);
  std::vector<int> sparse_order(ntot*d);
  std::vector<int> sparse_index(ntot*d);
  std::vector<int> sparse_unique_index(num_total_pts);
  quad_points.resize(ntot);
  quad_weights.resize(ntot);
  quad_values.resize(ntot);
  std::vector<value_type> gp(ntot*d);

  webbur::sparse_grid_mixed_growth_unique_index(d, level, &rules[0],
						&nparams[0], 
						&params[0],
						&compute1DPoints[0],
						1e-15, ntot, num_total_pts, 
						webbur::level_to_order_default,
						&sparse_unique_index[0]);
  webbur::sparse_grid_mixed_growth_index(d, level, &rules[0], ntot, 
					 num_total_pts, &sparse_unique_index[0],
					 webbur::level_to_order_default, 
					 &sparse_order[0], &sparse_index[0]);
  webbur::sparse_grid_mixed_growth_weight(d, level, &rules[0], 
					  &nparams[0], &params[0], 
					  &compute1DWeights[0],
					  ntot, num_total_pts, 
					  &sparse_unique_index[0], 
					  webbur::level_to_order_default,
					  &quad_weights[0]);
  webbur::sparse_grid_mixed_growth_point(d, level, &rules[0], 
					 &nparams[0], &params[0], 
					 &compute1DPoints[0],
					 ntot, &sparse_order[0], 
					 &sparse_index[0], 
					 webbur::level_to_order_default, 
					 &gp[0]);
  
  value_type weight_factor = 1.0;
  std::vector<value_type> point_factor(d);
  for (ordinal_type i=0; i<d; i++) {
    weight_factor *= coordinate_bases[i]->getQuadWeightFactor();
    point_factor[i] = coordinate_bases[i]->getQuadPointFactor();
  }
  for (ordinal_type i=0; i<ntot; i++) {
    quad_values[i].resize(sz);
    quad_points[i].resize(d);
    quad_weights[i] *= weight_factor;
    for (ordinal_type j=0; j<d; j++)
      quad_points[i][j] = gp[i*d+j]*point_factor[j];
    quad_values[i] = product_basis->evaluateBases(quad_points[i]);
  }

  //std::cout << "ntot = " << ntot << std::endl;

//   std::cout << "Sparse grid quadrature points, weights, values = " << std::endl;
//   for (int i=0; i<n; i++) {
//     std::cout << "\t" << this->quad_points[i][0] 
//               << "\t" << this->quad_weights[i];
//     for (ordinal_type j=0; j<sz; j++)
//       std::cout << "\t" << this->quad_values[i][j];
//     cout << std::endl;
//   }

  // // Check monomial exactness
  // std::cout << "\n---------- Monomial exactness check ----------\n";
  // double max_error = 0.0;
  // for (ordinal_type k=0; k<sz; k++) {
  //   double sgi = 0.0;
  //   for (ordinal_type gp=0; gp<quad_points.size(); gp++) {
  //     double v = 1.0;
  //     for (ordinal_type j=0; j<d; j++)
  // 	v *= std::pow(quad_points[gp][j],static_cast<int>(terms[k][j]));
  //     sgi += v*quad_weights[gp];
  //   }
  //   double exact = 1.0;
  //   for (ordinal_type j=0; j<d;j ++) {
  //     exact *= 0.5*(1.0  - std::pow(-1.0, static_cast<int>(terms[k][j]+1))) / 
  // 	(terms[k][j] + 1.0);
  //   }
  //   double error = std::abs(exact-sgi);
  //   if (error > max_error)
  //     max_error = error;
  //   std::cout << "(";
  //   for (ordinal_type j=0; j<d; j++) {
  //     if (j != 0)
  // 	std::cout << ", ";
  //     std::cout << terms[k][j];
  //   }
  //   std::cout << ") -- " << sgi << "\t" << exact << "\t" << error << "\n";
  // }
  // std::cout << "max error = " << max_error << std::endl;
}

template <typename ordinal_type, typename value_type>
const std::vector< std::vector<value_type> >&
Stokhos::SparseGridQuadrature<ordinal_type, value_type>::
getQuadPoints() const
{
  return quad_points;
}

template <typename ordinal_type, typename value_type>
const std::vector<value_type>&
Stokhos::SparseGridQuadrature<ordinal_type, value_type>::
getQuadWeights() const
{
  return quad_weights;
}

template <typename ordinal_type, typename value_type>
const std::vector< std::vector<value_type> >&
Stokhos::SparseGridQuadrature<ordinal_type, value_type>::
getBasisAtQuadPoints() const
{
  return quad_values;
}
