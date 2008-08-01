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

#include "Stokhos_Polynomial.hpp"

template <typename BasisT>
Stokhos::TripleProduct<BasisT>::
TripleProduct(const Teuchos::RCP<const BasisT>& basis_) :
  l(basis_->size()),
  basis(basis_),
  Bij(l*l),
  Cijk(l*l*l)
{
  compute();
}

template <typename BasisT>
Stokhos::TripleProduct<BasisT>::
~TripleProduct()
{
}

template <typename BasisT>
const typename Stokhos::TripleProduct<BasisT>::value_type&
Stokhos::TripleProduct<BasisT>::
double_deriv(unsigned int i, unsigned int j) const
{
  return Bij[ l*j + i ];
}

template <typename BasisT>
const typename Stokhos::TripleProduct<BasisT>::value_type&
Stokhos::TripleProduct<BasisT>::
triple_value(unsigned int i, unsigned int j, unsigned int k) const
{
  return Cijk[ l*(l*k + j) + i ];
}

template <typename BasisT>
const typename Stokhos::TripleProduct<BasisT>::value_type&
Stokhos::TripleProduct<BasisT>::
triple_deriv(unsigned int i, unsigned int j, unsigned int k) const
{
  return Dijk[ l*(l*k + j) + i ];
}

template <typename BasisT>
const typename Stokhos::TripleProduct<BasisT>::value_type&
Stokhos::TripleProduct<BasisT>::
norm_squared(unsigned int i) const
{
  return basis->norm_squared()[i];
}

template <typename BasisT>
void
Stokhos::TripleProduct<BasisT>::
compute()  
{
  // Compute Cijk = < \Psi_i \Psi_j \Psi_k >
  std::vector<value_type> a(2*l);
  for (unsigned int i=0; i<l; i++) {
    for (unsigned int j=0; j<l; j++) {
      basis->projectProduct(i, j, a);
      for (unsigned int k=0; k<l; k++)
	Cijk[ l*(l*k+j) + i ] = a[k];
    }
  }

  // Compute Dijk = < \Psi_i \Psi_j \Psi_k' >
  const std::vector<value_type>& nrm_sq = basis->norm_squared();
  for (unsigned int k=0; k<l; k++) {
    std::vector<value_type> b(l);
    basis->projectDerivative(k, b);
    for (unsigned int m=0; m<l; m++)
      Bij[ l*k+m ] = b[m]*nrm_sq[m];
    for (unsigned int i=0; i<l; i++)
      for (unsigned int j=0; j<l; j++) {
	Dijk[ l*(l*k+j) + i ] = value_type(0.0);
	for (unsigned int m=0; m<l; m++)
	  Dijk[ l*(l*k+j) + i ] += b[m]*Cijk[ l*(l*m+j) + i ];
      }
  }
}

template <typename BasisT>
unsigned int
Stokhos::TripleProduct<BasisT>::
num_values(unsigned int k) const
{
  return l*l;
}

template <typename BasisT>
void
Stokhos::TripleProduct<BasisT>::
triple_value(unsigned int k, unsigned int ll, 
	     unsigned int& i, unsigned int& j, value_type& c) const
{
  j = ll/l;
  i = ll-j*l;
  c = triple_value(i,j,k);
}
