// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_TimeMonitor.hpp"

template <typename ordinal_type, typename value_type, typename point_compare_type>
template <typename index_set_type>
Stokhos::SmolyakSparseGridQuadrature<ordinal_type, value_type, point_compare_type>::
SmolyakSparseGridQuadrature(
  const Teuchos::RCP<const ProductBasis<ordinal_type,value_type> >& product_basis,
  const index_set_type& index_set,
  const value_type duplicate_tol,
  const point_compare_type& point_compare)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos: Sparse Grid Generation");
#endif

  typedef SmolyakBasis<ordinal_type,value_type> smolyak_basis_type;
  smolyak_basis_type smolyak_basis(
    product_basis->getCoordinateBases(), index_set, duplicate_tol);

  typedef SmolyakPseudoSpectralOperator<ordinal_type,value_type,point_compare_type> smolyak_operator_type;
  smolyak_operator_type smolyak_operator(smolyak_basis, true, true,
                                         point_compare);
  ordinal_type nqp = smolyak_operator.point_size();
  ordinal_type npc = product_basis->size();

  // Compute quad points, weights, values
  quad_points.resize(nqp);
  quad_weights.resize(nqp);
  quad_values.resize(nqp);
  typedef typename smolyak_operator_type::const_set_iterator const_iterator;
  ordinal_type i = 0;
  for (const_iterator it = smolyak_operator.set_begin();
       it != smolyak_operator.set_end(); ++it, ++i) {
    quad_points[i] = it->first.getTerm();
    quad_weights[i] = it->second.first;
    quad_values[i].resize(npc);
    product_basis->evaluateBases(quad_points[i], quad_values[i]);
  }
}

template <typename ordinal_type, typename value_type, typename point_compare_type>
const Teuchos::Array< Teuchos::Array<value_type> >&
Stokhos::SmolyakSparseGridQuadrature<ordinal_type, value_type, point_compare_type>::
getQuadPoints() const
{
  return quad_points;
}

template <typename ordinal_type, typename value_type, typename point_compare_type>
const Teuchos::Array<value_type>&
Stokhos::SmolyakSparseGridQuadrature<ordinal_type, value_type, point_compare_type>::
getQuadWeights() const
{
  return quad_weights;
}

template <typename ordinal_type, typename value_type, typename point_compare_type>
const Teuchos::Array< Teuchos::Array<value_type> >&
Stokhos::SmolyakSparseGridQuadrature<ordinal_type, value_type, point_compare_type>::
getBasisAtQuadPoints() const
{
  return quad_values;
}

template <typename ordinal_type, typename value_type, typename point_compare_type>
std::ostream&
Stokhos::SmolyakSparseGridQuadrature<ordinal_type, value_type, point_compare_type>::
print(std::ostream& os) const
{
  ordinal_type nqp = quad_weights.size();
  os << "Smolyak Sparse Grid Quadrature with " << nqp << " points:"
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
