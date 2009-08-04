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

template <typename BasisT>
Sacado::PCE::Workspace<BasisT>::
Workspace(unsigned int sz_) :
  sz(sz_),
  A(2*sz,2*sz),
  b(2*sz,2),
  piv(2*sz),
  Cijk(sz-1),
  lapack()
{
}

template <typename BasisT>
void
Sacado::PCE::Workspace<BasisT>::
resize(unsigned int sz_)
{
  if (sz_ != sz) {
    sz = sz_;
    A.shape(2*sz,2*sz);
    b.shape(2*sz,2);
    piv.resize(2*sz);
    Cijk.resize(sz-1);
  }
}

extern "C" {
  double F77_BLAS_MANGLE(dlange,DLANGE)(char*, int*, int*, double*, int*, double*);
}

template <typename BasisT>
typename Sacado::PCE::Workspace<BasisT>::ordinal_type
Sacado::PCE::Workspace<BasisT>::
solve(typename Sacado::PCE::Workspace<BasisT>::ordinal_type s,
      typename Sacado::PCE::Workspace<BasisT>::ordinal_type nrhs)
{
  ordinal_type info;
//   lapack.GESV(s, nrhs, A.values(), A.numRows(), &(piv[0]), b.values(), 
// 	      b.numRows(), &info);
  lapack.GETRF(s, s, A.values(), A.numRows(), &(piv[0]), &info);
  value_type norm, rcond;
  std::vector<ordinal_type> iwork(4*s);
  std::vector<value_type> work(4*s);
  norm = 1.0;
  ordinal_type n = A.numRows();
  char t = '1';
  norm = F77_BLAS_MANGLE(dlange,DLANGE)(&t, &s, &s, A.values(), &n, &work[0]);
  lapack.GECON('1', s, A.values(), A.numRows(), norm, &rcond, &work[0], 
	       &iwork[0], &info);
  std::cout << "condition number = " << 1.0/rcond << std::endl;
  lapack.GETRS('N', s, nrhs, A.values(), A.numRows(), &(piv[0]), b.values(), 
	       b.numRows(), &info);
  return info;
}
