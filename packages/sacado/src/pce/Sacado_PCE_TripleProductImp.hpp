// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Sacado_PCE_StandardPoly.hpp"

template <typename BasisT>
Sacado::PCE::TripleProduct<BasisT>::
TripleProduct(unsigned int degree) :
  l(degree+1),
  basis(2*degree),
  Cijk(l*l*l)
{
  compute();
}

template <typename BasisT>
Sacado::PCE::TripleProduct<BasisT>::
TripleProduct(const Sacado::PCE::TripleProduct<BasisT>& tp) :
  l(tp.l),
  basis(tp.basis),
  Cijk(tp.Cijk)
{
}

template <typename BasisT>
Sacado::PCE::TripleProduct<BasisT>::
~TripleProduct()
{
}

template <typename BasisT>
Sacado::PCE::TripleProduct<BasisT>&
Sacado::PCE::TripleProduct<BasisT>::
operator=(const Sacado::PCE::TripleProduct<BasisT>& tp)
{
  if (this != &tp) {
    l = tp.l;
    basis = tp.basis;
    Cijk = tp.Cijk;
  }
  return *this;
}

template <typename BasisT>
const typename Sacado::PCE::TripleProduct<BasisT>::value_type&
Sacado::PCE::TripleProduct<BasisT>::
value(unsigned int i, unsigned int j, unsigned int k) const
{
  return Cijk[ l*(l*k + j) + i ];
}

template <typename BasisT>
const typename Sacado::PCE::TripleProduct<BasisT>::value_type&
Sacado::PCE::TripleProduct<BasisT>::
norm_squared(unsigned int i) const
{
  return basis.norm_squared()[i];
}

template <typename BasisT>
void
Sacado::PCE::TripleProduct<BasisT>::
resize(unsigned int degree)
{
  if (degree+1 != l) {
    l = degree+1;
    basis = BasisT(2*degree);
    Cijk.resize(l*l*l);
    compute();
  }
}

// Important note:  To get the correct value for <\Psi_i \Psi_j \Psi_k>,
// we have to expand \Psi_i*\Psi_j in the full degree_i+degree_j basis, not
// just the degree d basis.  There for the basis needs to be of size 2*d
template <typename BasisT>
void
Sacado::PCE::TripleProduct<BasisT>::
compute()  
{
  const std::vector<value_type>& nrm_sq = basis.norm_squared();
  StandardPoly<value_type> pij(2*(l-1));
  std::vector<value_type> a(2*l);
  for (unsigned int i=0; i<l; i++) {
    const StandardPoly<value_type>& pi = basis.getBasisPoly(i);
    for (unsigned int j=0; j<l; j++) {
      const StandardPoly<value_type>& pj = basis.getBasisPoly(j);
      pij.multiply(1.0, pi, pj, 0.0);
      basis.project(pij, a);
      for (unsigned int k=0; k<l; k++)
	Cijk[ l*(l*k+j) + i ] = a[k]*nrm_sq[k];
    }
  }
}
