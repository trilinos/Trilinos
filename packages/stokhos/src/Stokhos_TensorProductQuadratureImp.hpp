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
Stokhos::TensorProductQuadrature<ordinal_type, value_type>::
TensorProductQuadrature(Teuchos::RCP<const OrthogPolyBasis<ordinal_type,value_type> >& product_basis) 
{
  ordinal_type d = product_basis->dimension();
  ordinal_type sz = product_basis->size();

  const std::vector< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type,value_type> > >& coordinate_bases = product_basis->getCoordinateBases();

  

  // Compute quad points, weights, values
  std::vector< std::vector<value_type> > gp(d);
  std::vector< std::vector<value_type> > gw(d);
  std::vector< std::vector< std::vector<value_type> > > gv(d);
  std::vector<ordinal_type> n(d);
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
  std::vector<ordinal_type> index(d);
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
      std::vector<ordinal_type> term = product_basis->getTerm(k);
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

  //std::cout << "ntot = " << ntot << std::endl;

//   std::cout << "Tensor product quadrature points, weights, values = " << std::endl;
//   for (int i=0; i<n; i++) {
//     std::cout << "\t" << this->quad_points[i][0] 
//               << "\t" << this->quad_weights[i];
//     for (ordinal_type j=0; j<sz; j++)
//       std::cout << "\t" << this->quad_values[i][j];
//     cout << std::endl;
//   }
}

template <typename ordinal_type, typename value_type>
Stokhos::TensorProductQuadrature<ordinal_type, value_type>::
TensorProductQuadrature(Teuchos::RCP<const OrthogPolyBasis<ordinal_type,value_type> >& product_basis,const ordinal_type& quad_order) 
{
  ordinal_type d = product_basis->dimension();
  ordinal_type sz = product_basis->size();

  const std::vector< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type,value_type> > >& coordinate_bases = product_basis->getCoordinateBases();

  

  // Compute quad points, weights, values
  std::vector< std::vector<value_type> > gp(d);
  std::vector< std::vector<value_type> > gw(d);
  std::vector< std::vector< std::vector<value_type> > > gv(d);
  std::vector<ordinal_type> n(d);
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
  std::vector<ordinal_type> index(d);
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
      std::vector<ordinal_type> term = product_basis->getTerm(k);
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

  //std::cout << "ntot = " << ntot << std::endl;

//   std::cout << "Tensor product quadrature points, weights, values = " << std::endl;
//   for (int i=0; i<n; i++) {
//     std::cout << "\t" << this->quad_points[i][0] 
//               << "\t" << this->quad_weights[i];
//     for (ordinal_type j=0; j<sz; j++)
//       std::cout << "\t" << this->quad_values[i][j];
//     cout << std::endl;
//   }
}

template <typename ordinal_type, typename value_type>
const std::vector< std::vector<value_type> >&
Stokhos::TensorProductQuadrature<ordinal_type, value_type>::
getQuadPoints() const
{
  return quad_points;
}

template <typename ordinal_type, typename value_type>
const std::vector<value_type>&
Stokhos::TensorProductQuadrature<ordinal_type, value_type>::
getQuadWeights() const
{
  return quad_weights;
}

template <typename ordinal_type, typename value_type>
const std::vector< std::vector<value_type> >&
Stokhos::TensorProductQuadrature<ordinal_type, value_type>::
getBasisAtQuadPoints() const
{
  return quad_values;
}
