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

extern "C" {
  void ncijk_(int*, int*);
  void triple_(int*, int*, int*, int*, double*);
  void psinorm_(int*, double*);
}

template <typename BasisT>
Stokhos::TayTripleProduct<BasisT>::
TayTripleProduct()
{
}

template <typename BasisT>
Stokhos::TayTripleProduct<BasisT>::
~TayTripleProduct()
{
}

template <typename BasisT>
unsigned int
Stokhos::TayTripleProduct<BasisT>::
num_values(unsigned int k) const
{
  int kk = k;
  int n;
  ncijk_(&kk, &n);
  return n;
}

template <typename BasisT>
void
Stokhos::TayTripleProduct<BasisT>::
triple_value(unsigned int k, unsigned int l, 
	     unsigned int& i, unsigned int& j, 
	     typename Stokhos::TayTripleProduct<BasisT>::value_type& c) const
{
  int kk = k;
  int ll = l+1;
  int ii, jj;
  triple_(&kk, &ll, &ii, &jj, &c);
  i = ii;
  j = jj;
}

template <typename BasisT>
typename Stokhos::TayTripleProduct<BasisT>::value_type
Stokhos::TayTripleProduct<BasisT>::
norm_squared(unsigned int i) const
{
  int ii = i;
  double nrm;
  psinorm_(&ii, &nrm);
  return nrm;
}

