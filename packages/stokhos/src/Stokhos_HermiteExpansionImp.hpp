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

#include "Teuchos_TestForException.hpp"
#include "Stokhos_DynamicArrayTraits.hpp"

template <typename T, typename BasisT> 
Stokhos::HermiteExpansion<T,BasisT>::
HermiteExpansion(unsigned int d) :
  sz(d+1),
  A(2*sz,2*sz),
  B(2*sz,2),
  piv(2*sz),
  Cijk(sz-1),
  lapack()
{
}

template <typename T, typename BasisT> 
void
Stokhos::HermiteExpansion<T,BasisT>::
resize(unsigned int d)
{
  if (sz != d+1) {
    sz = d+1;
    A.shape(2*sz,2*sz);
    B.shape(2*sz,2);
    piv.resize(2*sz);
    Cijk.resize(sz-1);
  }
}

extern "C" {
  double dlange_(char*, int*, int*, double*, int*, double*);
}

template <typename T, typename BasisT>
typename Stokhos::HermiteExpansion<T,BasisT>::ordinal_type
Stokhos::HermiteExpansion<T,BasisT>::
solve(typename Stokhos::HermiteExpansion<T,BasisT>::ordinal_type s,
      typename Stokhos::HermiteExpansion<T,BasisT>::ordinal_type nrhs)
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
  norm = dlange_(&t, &s, &s, A.values(), &n, &work[0]);
  lapack.GECON('1', s, A.values(), A.numRows(), norm, &rcond, &work[0], 
	       &iwork[0], &info);
  std::cout << "condition number = " << 1.0/rcond << std::endl;
  lapack.GETRS('N', s, nrhs, A.values(), A.numRows(), &(piv[0]), B.values(), 
	       B.numRows(), &info);
  return info;
}

template <typename T, typename BasisT> 
void
Stokhos::HermiteExpansion<T,BasisT>::
unaryMinus(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a)
{
  T* cc = c.coeff();
  const T* ca = a.coeff();
  unsigned int da = a.degree();

  for (unsigned int i=0; i<=da; i++)
    cc[i] = -ca[i];
}

template <typename T, typename BasisT> 
void
Stokhos::HermiteExpansion<T,BasisT>::
plusEqual(Stokhos::HermitePoly<T>& c, const T& val)
{
  c[0] += val;
}

template <typename T, typename BasisT> 
void
Stokhos::HermiteExpansion<T,BasisT>::
minusEqual(Stokhos::HermitePoly<T>& c, const T& val)
{
  c[0] -= val;
}

template <typename T, typename BasisT> 
void
Stokhos::HermiteExpansion<T,BasisT>::
timesEqual(Stokhos::HermitePoly<T>& c, const T& val)
{
  unsigned int dc = c.degree();
  T* cc = c.coeff();
  for (unsigned int i=0; i<dc; i++)
    cc[i] *= val;
}

template <typename T, typename BasisT> 
void
Stokhos::HermiteExpansion<T,BasisT>::
divideEqual(Stokhos::HermitePoly<T>& c, const T& val)
{
  unsigned int dc = c.degree();
  T* cc = c.coeff();
  for (unsigned int i=0; i<=dc; i++)
    cc[i] /= val;
}

template <typename T, typename BasisT> 
void
Stokhos::HermiteExpansion<T,BasisT>::
plusEqual(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& x)
{
  unsigned int d = c.degree();
  unsigned int xd = x.degree();
  unsigned int dmin = xd < d ? xd : d;
  if (xd > d)
    c.resize(xd,true);

  T* cc = c.coeff();
  T* xc = x.coeff();
  for (unsigned int i=0; i<=dmin; i++)
    cc[i] += xc[i];
  if (d < xd)
    for (unsigned int i=d+1; i<=xd; i++)
      cc[i] = xc[i];
}

template <typename T, typename BasisT> 
void
Stokhos::HermiteExpansion<T,BasisT>::
minusEqual(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& x)
{
  unsigned int d = c.degree();
  unsigned int xd = x.degree();
  unsigned int dmin = xd < d ? xd : d;
  if (xd > d)
    c.resize(xd,true);

  T* cc = c.coeff();
  T* xc = x.coeff();
  for (unsigned int i=0; i<=dmin; i++)
    cc[i] -= xc[i];
  if (d < xd)
    for (unsigned int i=d+1; i<=xd; i++)
      cc[i] = -xc[i];
}

template <typename T, typename BasisT> 
void
Stokhos::HermiteExpansion<T,BasisT>::
timesEqual(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& x)
{
#ifdef STOKHOS_DEBUG
  const char* func = "Stokhos::HermiteExpansion::timesEqual()";
  TEST_FOR_EXCEPTION((x.degree() != c.degree()) && (x.degree() != 0) && 
		     (c.degree() != 0), std::logic_error,
		     func << ":  Attempt to assign with incompatible degrees");
#endif
  unsigned int d = c.degree();
  unsigned int xd = x.degree();

  T* cc = c.coeff();
  T* xc = x.coeff();
  
  if (d > 0 && xd > 0) {
#ifdef STOKHOS_DEBUG
    TEST_FOR_EXCEPTION(degree() < d, std::logic_error,
		       func << ":  Expansion degree (" << degree() 
		       << ") is too small for computation (" << d 
		       << " needed).");
#endif

    // Copy c coefficients into temporary array
    T* tc = Stokhos::ds_array<T>::get_and_fill(cc,d+1);
    T tmp;
    for (unsigned int k=0; k<=d; k++) {
      tmp = T(0.0);
      for (unsigned int i=0; i<=d; i++) {
	for (unsigned int j=0; j<=xd; j++)
	  tmp += Cijk.value(i,j,k)*tc[i]*xc[j];
      }
      cc[k] = tmp / Cijk.norm_squared(k);
    }
  }
  else if (d > 0) {
    for (unsigned int i=0; i<=d; i++)
      cc[i] *= xc[0];
  }
  else  {
    if (xd > d)
      c.resize(xd,true);
    for (int i=xd; i>=0; i--)
	cc[i] = cc[0]*xc[i];
  }
}

template <typename T, typename BasisT> 
void
Stokhos::HermiteExpansion<T,BasisT>::
divideEqual(Stokhos::HermitePoly<T>& c, const HermitePoly<T>& x)
{
  const char* func = "Stokhos::HermiteExpansion::divideEquals()";

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION((x.degree() != c.degree()) && (x.degree() != 0) && 
		     (c.degree() != 0), std::logic_error,
		     func << ":  Attempt to assign with incompatible degrees");
#endif

  unsigned int d = c.degree();
  unsigned int xd = x.degree();

  T* cc = c.coeff();
  T* xc = x.coeff();
  
  if (xd > 0) {
    if (xd > d)
      c.resize(xd,true);

#ifdef STOKHOS_DEBUG
    TEST_FOR_EXCEPTION(degree() < xd, std::logic_error,
		       func << ":  Expansion degree (" << degree() 
		       << ") is too small for computation (" << xd
		       << " needed).");
#endif
    
    // Fill A
    for (unsigned int i=0; i<=xd; i++) {
      for (unsigned int j=0; j<=xd; j++) {
	A(i,j) = 0.0;
	for (unsigned int k=0; k<=xd; k++) {
	  A(i,j) += Cijk.value(i,j,k)*xc[k];
	}
	A(i,j) /= Cijk.norm_squared(i);
      }
    }
    
    // Fill B
    for (unsigned int i=0; i<=xd; i++)
      B(i,0) = cc[i];

    // Solve system
    int info = solve(xd+1, 1);

    TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		       func << ":  Argument " << info 
		            << " for solve had illegal value");
    TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		       func << ":  Diagonal entry " << info 
		            << " in LU factorization is exactly zero");

    // Get coefficients
    for (unsigned int i=0; i<=xd; i++)
      cc[i] = B(i,0);
  }
  else {
    for (unsigned int i=0; i<=d; i++)
      cc[i] /= xc[0];
  }
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
plus(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a, 
     const Stokhos::HermitePoly<T>& b)
{
#ifdef STOKHOS_DEBUG
  const char* func = "Stokhos::HermiteExpansion::plus()";
  TEST_FOR_EXCEPTION((a.degree() != b.degree()) && (a.degree() != 0) && 
		     (b.degree() != 0), 
		     std::logic_error,
		     func << ":  Arguments have incompatible degrees!");
#endif

  unsigned int da = a.degree();
  unsigned int db = b.degree();
  unsigned int dc = da > db ? da : db;
  if (dc != c.degree())
    c.resize(dc);

  const T* ca = a.coeff();
  const T* cb = b.coeff();
  T* cc = c.coeff();

  if (da > 0 && db > 0) {
    for (unsigned int i=0; i<=dc; i++)
      cc[i] = ca[i] + cb[i];
  }
  else if (da > 0) {
    cc[0] = ca[0] + cb[0];
    for (unsigned int i=1; i<=dc; i++)
      cc[i] = ca[i];
  }
  else if (db >= 0) {
    cc[0] = ca[0] + cb[0];
    for (unsigned int i=1; i<=dc; i++)
      cc[i] = cb[i];
  }
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
plus(Stokhos::HermitePoly<T>& c, const T& a, const Stokhos::HermitePoly<T>& b)
{
  unsigned int dc = b.degree();
  if (dc != c.degree())
    c.resize(dc);

  const T* cb = b.coeff();
  T* cc = c.coeff();

  cc[0] = a + cb[0];
  for (unsigned int i=1; i<=dc; i++)
    cc[i] = cb[i];
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
plus(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a, const T& b)
{
  unsigned int dc = a.degree();
  if (dc != c.degree())
    c.resize(dc);

  const T* ca = a.coeff();
  T* cc = c.coeff();

  cc[0] = ca[0] + b;
  for (unsigned int i=1; i<=dc; i++)
    cc[i] = ca[i];
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
minus(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a, 
      const Stokhos::HermitePoly<T>& b)
{
#ifdef STOKHOS_DEBUG
  const char* func = "Stokhos::HermiteExpansion::minus()";
  TEST_FOR_EXCEPTION((a.degree() != b.degree()) && (a.degree() != 0) && 
		     (b.degree() != 0), 
		     std::logic_error,
		     func << ":  Arguments have incompatible degrees!");
#endif

  unsigned int da = a.degree();
  unsigned int db = b.degree();
  unsigned int dc = da > db ? da : db;
  if (dc != c.degree())
    c.resize(dc);

  const T* ca = a.coeff();
  const T* cb = b.coeff();
  T* cc = c.coeff();

  if (da > 0 && db > 0) {
    for (unsigned int i=0; i<=dc; i++)
      cc[i] = ca[i] - cb[i];
  }
  else if (da > 0) {
    cc[0] = ca[0] - cb[0];
    for (unsigned int i=1; i<=dc; i++)
      cc[i] = ca[i];
  }
  else if (db >= 0) {
    cc[0] = ca[0] - cb[0];
    for (unsigned int i=1; i<=dc; i++)
      cc[i] = -cb[i];
  }
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
minus(Stokhos::HermitePoly<T>& c, const T& a, const Stokhos::HermitePoly<T>& b)
{
  unsigned int dc = b.degree();
  if (dc != c.degree())
    c.resize(dc);

  const T* cb = b.coeff();
  T* cc = c.coeff();

  cc[0] = a - cb[0];
  for (unsigned int i=1; i<=dc; i++)
    cc[i] = -cb[i];
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
minus(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a, const T& b)
{
  unsigned int dc = a.degree();
  if (dc != c.degree())
    c.resize(dc);

  const T* ca = a.coeff();
  T* cc = c.coeff();

  cc[0] = ca[0] - b;
  for (unsigned int i=1; i<=dc; i++)
    cc[i] = ca[i];

  return c;
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
times(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a, 
      const Stokhos::HermitePoly<T>& b)
{
#ifdef STOKHOS_DEBUG
  const char* func = "Stokhos::HermiteExpansion::times()";
  TEST_FOR_EXCEPTION((a.degree() != b.degree()) && (a.degree() != 0) && 
		     (b.degree() != 0), 
		     std::logic_error,
		     func << ":  Arguments have incompatible degrees!");
#endif

  unsigned int da = a.degree();
  unsigned int db = b.degree();
  unsigned int dc = da > db ? da : db;
  if (dc != c.degree())
    c.resize(dc);

  const T* ca = a.coeff();
  const T* cb = b.coeff();
  T* cc = c.coeff();

  if (da > 0 && db > 0) {
#ifdef STOKHOS_DEBUG
    TEST_FOR_EXCEPTION(degree() < dc, std::logic_error,
		       func << ":  Expansion degree (" << degree()
		       << ") is too small for computation (" << dc
		       << " needed).");
#endif

    T tmp;
    for (unsigned int k=0; k<=dc; k++) {
      tmp = T(0.0);
      for (unsigned int i=0; i<=da; i++) {
	for (unsigned int j=0; j<=db; j++)
	  tmp += Cijk.value(i,j,k)*ca[i]*cb[j];
      }
      cc[k] = tmp / Cijk.norm_squared(k);
    }
  }
  else if (da > 0) {
    for (unsigned int i=0; i<=dc; i++)
      cc[i] = ca[i]*cb[0];
  }
  else if (db >= 0) {
    for (unsigned int i=0; i<=dc; i++)
      cc[i] = ca[0]*cb[i];
  }
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
times(Stokhos::HermitePoly<T>& c, const T& a, const Stokhos::HermitePoly<T>& b)
{
  unsigned int dc = b.degree();
  if (dc != c.degree())
    c.resize(dc);

  const T* cb = b.coeff();
  T* cc = c.coeff();

  for (unsigned int i=0; i<=dc; i++)
    cc[i] = a*cb[i];
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
times(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a, const T& b)
{
  unsigned int dc = a.degree();
  if (dc != c.degree())
    c.resize(dc);

  const T* ca = a.coeff();
  T* cc = c.coeff();

  for (unsigned int i=0; i<=dc; i++)
    cc[i] = ca[i]*b;

  return c;
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
divide(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a, 
       const Stokhos::HermitePoly<T>& b)
{
  const char* func = "Stokhos::HermiteExpansion::divide()";

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION((a.degree() != b.degree()) && (a.degree() != 0) && 
		     (b.degree() != 0), 
		     std::logic_error,
		     func << ":  Arguments have incompatible degrees!");
#endif

  unsigned int da = a.degree();
  unsigned int db = b.degree();
  unsigned int dc = da > db ? da : db;
  if (dc != c.degree())
    c.resize(dc);

  const T* ca = a.coeff();
  const T* cb = b.coeff();
  T* cc = c.coeff();

  if (db > 0) {

#ifdef STOKHOS_DEBUG
    TEST_FOR_EXCEPTION(degree() < dc, std::logic_error,
		       func << ":  Expansion degree (" << degree() 
		       << ") is too small for computation (" << dc
		       << " needed).");
#endif

    // Fill A
    for (unsigned int i=0; i<=dc; i++) {
      for (unsigned int j=0; j<=dc; j++) {
	A(i,j) = 0.0;
	for (unsigned int k=0; k<=db; k++) {
	  A(i,j) += Cijk.value(i,j,k)*cb[k];
	}
	A(i,j) /= Cijk.norm_squared(i);
      }
    }
    
    // Fill B
    for (unsigned int i=0; i<=da; i++)
      B(i,0) = ca[i];
    for (unsigned int i=da+1; i<=dc; i++)
      B(i,0) = T(0.0);

    // Solve system
    int info = solve(dc+1, 1);

    TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		       func << ":  Argument " << info 
		            << " for solve had illegal value");
    TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		       func << ":  Diagonal entry " << info 
		            << " in LU factorization is exactly zero");

    // Get coefficients
    for (unsigned int i=0; i<=dc; i++)
      cc[i] = B(i,0);
  }
  else {
    for (unsigned int i=0; i<=da; i++)
      cc[i] = ca[i]/cb[0];
  }
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
divide(Stokhos::HermitePoly<T>& c, const T& a, 
       const Stokhos::HermitePoly<T>& b)
{
  const char* func = "Stokhos::HermiteExpansion::divide()";

  unsigned int dc = b.degree();
  if (dc != c.degree())
    c.resize(dc);

  const T* cb = b.coeff();
  T* cc = c.coeff();

  if (dc > 0) {
#ifdef STOKHOS_DEBUG
    TEST_FOR_EXCEPTION(degree() < dc, std::logic_error,
		       func << ":  Expansion degree (" << degree()
		       << ") is too small for computation (" << dc
		       << " needed).");
#endif
   
    // Fill A
    for (unsigned int i=0; i<=dc; i++) {
      for (unsigned int j=0; j<=dc; j++) {
	A(i,j) = 0.0;
	for (unsigned int k=0; k<=dc; k++) {
	  A(i,j) += Cijk.value(i,j,k)*cb[k];
	}
	A(i,j) /= Cijk.norm_squared(i);
      }
    }

    // Fill B
    B(0,0) = a;
    for (unsigned int i=1; i<=dc; i++)
      B(i,0) = T(0.0);

    // Solve system
    int info = solve(dc+1, 1);

    TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		       func << ":  Argument " << info 
		            << " for solve had illegal value");
    TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		       func << ":  Diagonal entry " << info 
		            << " in LU factorization is exactly zero");

    // Get coefficients
    for (unsigned int i=0; i<=dc; i++)
      cc[i] = B(i,0);
  }
  else 
    cc[0] = a / cb[0];
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
divide(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a, 
       const T& b)
{
  unsigned int dc = a.degree();
  if (dc != c.degree())
    c.resize(dc);

  const T* ca = a.coeff();
  T* cc = c.coeff();

  for (unsigned int i=0; i<=dc; i++)
    cc[i] = ca[i]/b;
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
exp(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a)
{
  const char* func = "Stokhos::HermiteExpansion::exp()";

  unsigned int dc = a.degree();
  if (dc != c.degree())
    c.resize(dc);

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(degree() < dc, std::logic_error,
		     func << ":  Expansion degree (" << degree() 
		     << ") is too small for computation (" << dc
		     << " needed).");
#endif

  const T* ca = a.coeff();
  T* cc = c.coeff();

  const typename tp_type::basis_type& basis = Cijk.getBasis();

  // Fill A and B
  for (unsigned int i=1; i<=dc; i++) {
    for (unsigned int j=1; j<=dc; j++) {
      A(i-1,j-1) = 0.0;
      for (unsigned int k=1; k<=dc; k++) {
	if (i+j >= k) {
	  A(i-1,j-1) += Cijk.value(i-1,j,k-1)*basis.derivCoeff(k)*ca[k];
	}
      }
      A(i-1,j-1) /= Cijk.norm_squared(i-1)*basis.derivCoeff(i);
    }
    A(i-1,i-1) -= T(1.0);
    B(i-1,0) = -ca[i];
  }

  // Solve system
  int info = solve(dc, 1);

  TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		     func << ":  Argument " << info 
		     << " for solve had illegal value");
  TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		     func << ":  Diagonal entry " << info 
		     << " in LU factorization is exactly zero");

  // Compute degree-0 coefficient
  T s = basis.getBasisPoly(0).coeff(0) * ca[0];
  T t = T(1.0);
  for (unsigned int i=1; i<=dc; i++) {
    s += basis.getBasisPoly(i).coeff(0) * ca[i];
    t += basis.getBasisPoly(i).coeff(0) * B(i-1,0);
  }
  s = std::exp(s);
  cc[0] = (s/t);

  // Compute remaining coefficients
  for (unsigned int i=1; i<=dc; i++)
    cc[i] = B(i-1,0) * cc[0];

  cc[0] /= basis.getBasisPoly(0).coeff(0);
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
log(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a)
{
   const char* func = "Stokhos::HermiteExpansion::log()";

  unsigned int dc = a.degree();
  if (dc != c.degree())
    c.resize(dc);

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(degree() < dc, std::logic_error,
		     func << ":  Expansion degree (" << degree() 
		     << ") is too small for computation (" << dc
		     << " needed).");
#endif

  const T* ca = a.coeff();
  T* cc = c.coeff();

  const typename tp_type::basis_type& basis = Cijk.getBasis();

  // Fill A and B
  for (unsigned int i=1; i<=dc; i++) {
    for (unsigned int j=1; j<=dc; j++) {
      A(i-1,j-1) = 0.0;
      for (unsigned int k=0; k<=dc; k++) {
	if (i+j-2 >= k) {
	  A(i-1,j-1) += Cijk.value(i-1,j-1,k)*basis.derivCoeff(j)*ca[k];
	}
      }
      A(i-1,j-1) /= Cijk.norm_squared(i-1)*basis.derivCoeff(i);
    }
    B(i-1,0) = ca[i];
  }

  // Solve system
  int info = solve(dc, 1);

  TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		     func << ":  Argument " << info 
		     << " for solve had illegal value");
  TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		     func << ":  Diagonal entry " << info 
		     << " in LU factorization is exactly zero");

  // Compute degree-0 coefficient
  T s = basis.getBasisPoly(0).coeff(0) * ca[0];
  T t = T(0.0);
  for (unsigned int i=1; i<=dc; i++) {
    s += basis.getBasisPoly(i).coeff(0) * ca[i];
    t += basis.getBasisPoly(i).coeff(0) * B(i-1,0);
  }
  cc[0] = (std::log(s) - t) / basis.getBasisPoly(0).coeff(0);

  // Compute remaining coefficients
  for (unsigned int i=1; i<=dc; i++)
    cc[i] = B(i-1,0);
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
log10(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a)
{
  log(c,a);
  divide(c,a,std::log(10.0));
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
sqrt(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a)
{
  log(c,a);
  timesEqual(c,T(0.5));
  exp(c,c);
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
pow(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a,
    const Stokhos::HermitePoly<T>& b)
{
  log(c,a);
  timesEqual(c,b);
  exp(c,c);
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
pow(Stokhos::HermitePoly<T>& c, const T& a, const Stokhos::HermitePoly<T>& b)
{
  times(c,std::log(a),b);
  exp(c,c);
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
pow(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a, const T& b)
{
  log(c,a);
  timesEqual(c,b);
  exp(c,c);
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
sincos(Stokhos::HermitePoly<T>& s, Stokhos::HermitePoly<T>& c, 
       const Stokhos::HermitePoly<T>& a)
{
  const char* func = "Stokhos::HermiteExpansion::sincos()";
  unsigned int dc = a.degree();

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(degree() < dc, std::logic_error,
		     func << ":  Expansion degree (" << degree() 
		     << ") is too small for computation (" << dc 
		     << " needed).");
#endif

  if (s.degree() != dc)
    s.resize(dc);
  if (c.degree() != dc)
    c.resize(dc);

  const T* ca = a.coeff();
  T* cs = s.coeff();
  T* cc = c.coeff();

  const typename tp_type::basis_type& basis = Cijk.getBasis();

  // Fill A and b
  A.putScalar(T(0.0));
  B.putScalar(T(0.0));
  for (unsigned int i=1; i<=dc; i++) {
    for (unsigned int j=1; j<=dc; j++) {
      T tmp = T(0.0);
      for (unsigned int k=1; k<=dc; k++) {
	if (i+j >= k) {
	  tmp += Cijk.value(i-1,j,k-1)*basis.derivCoeff(k)*ca[k];
	}
      }
      tmp /= Cijk.norm_squared(i-1)*basis.derivCoeff(i);
      A(i-1,j-1+dc) = -tmp;
      A(i-1+dc,j-1) = tmp;
    }
    A(i-1,i-1) = T(1.0);
    A(i-1+dc,i-1+dc) = T(1.0);
    B(i-1,0) = ca[i];
    B(i-1+dc,1) = ca[i];
  }

  // Solve system
  int info = solve(2*dc, 2);

  TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		     func << ":  Argument " << info 
		     << " for solve had illegal value");
  TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		     func << ":  Diagonal entry " << info 
		     << " in LU factorization is exactly zero");

  // Compute degree-0 coefficients
  T t = basis.getBasisPoly(0).coeff(0) * ca[0];
  T a00 = T(0.0);
  T a01 = T(0.0);
  T a10 = T(0.0);
  T a11 = T(0.0);
  T b0 = B(0,0);
  T b1 = B(1,0);
  for (unsigned int i=1; i<=dc; i++) {
    t += basis.getBasisPoly(i).coeff(0) * ca[i];
    a00 += basis.getBasisPoly(i).coeff(0) * B(i-1,1);
    a01 += basis.getBasisPoly(i).coeff(0) * B(i-1,0);
    a10 += basis.getBasisPoly(i).coeff(0) * B(i-1+dc,1);
    a11 += basis.getBasisPoly(i).coeff(0) * B(i-1+dc,0);
  }
  A(0,0) = T(1.0) - a00;
  A(0,1) = a01;
  A(1,0) = -a10;
  A(1,1) = T(1.0) + a11;
  B(0,0) = std::sin(t);
  B(1,0) = std::cos(t);
  info = solve(2, 1);
  TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		     func << ":  Argument " << info 
		     << " for (2x2) solve had illegal value");
  TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		     func << ":  Diagonal entry " << info 
		     << " in (2x2) LU factorization is exactly zero");
  cs[0] = B(0,0);
  cc[0] = B(1,0);

  // Compute remaining coefficients
  B(0,0) = b0;
  B(1,0) = b1;
  for (unsigned int i=1; i<=dc; i++) {
    cs[i] = cc[0]*B(i-1,0) - cs[0]*B(i-1,1);
    cc[i] = cc[0]*B(i-1+dc,0) - cs[0]*B(i-1+dc,1);
  }

  cs[0] /= basis.getBasisPoly(0).coeff(0);
  cc[0] /= basis.getBasisPoly(0).coeff(0);
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
sin(Stokhos::HermitePoly<T>& s, const Stokhos::HermitePoly<T>& a)
{
  unsigned int dc = a.degree();
  HermitePoly<T> c(dc);
  sincos(s, c, a);
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
cos(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a)
{
  unsigned int dc = a.degree();
  HermitePoly<T> s(dc);
  sincos(s, c, a);
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
tan(Stokhos::HermitePoly<T>& t, const Stokhos::HermitePoly<T>& a)
{
  unsigned int dc = a.degree();
  HermitePoly<T> c(dc);
  
  sincos(t, c, a);
  divideEqual(t,c);
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
sinhcosh(Stokhos::HermitePoly<T>& s, Stokhos::HermitePoly<T>& c, 
	 const Stokhos::HermitePoly<T>& a)
{
  const char* func = "Stokhos::HermiteExpansion::sinhcosh()";
  unsigned int dc = a.degree();

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(degree() < dc, std::logic_error,
		     func << ":  Expansion degree (" << degree() 
		     << ") is too small for computation (" << dc
		     << " needed).");
#endif

  if (s.degree() != dc)
    s.resize(dc);
  if (c.degree() != dc)
    c.resize(dc);

  const T* ca = a.coeff();
  T* cs = s.coeff();
  T* cc = c.coeff();

  const typename tp_type::basis_type& basis = Cijk.getBasis();

  // Fill A and B
  A.putScalar(T(0.0));
  B.putScalar(T(0.0));
  for (unsigned int i=1; i<=dc; i++) {
    for (unsigned int j=1; j<=dc; j++) {
      T tmp = T(0.0);
      for (unsigned int k=1; k<=dc; k++) {
	if (i+j >= k) {
	  tmp += Cijk.value(i-1,j,k-1)*basis.derivCoeff(k)*ca[k];
	}
      }
      tmp /= Cijk.norm_squared(i-1)*basis.derivCoeff(i);
      A(i-1,j-1+dc) = -tmp;
      A(i-1+dc,j-1) = -tmp;
    }
    A(i-1,i-1) = T(1.0);
    A(i-1+dc,i-1+dc) = T(1.0);
    B(i-1,0) = ca[i];
    B(i-1+dc,1) = ca[i];
  }

  // Solve system
  int info = solve(2*dc, 2);

  TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		     func << ":  Argument " << info 
		     << " for solve had illegal value");
  TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		     func << ":  Diagonal entry " << info 
		     << " in LU factorization is exactly zero");

  // Compute degree-0 coefficients
  T t = basis.getBasisPoly(0).coeff(0) * ca[0];
  T a00 = T(1.0);
  T a01 = T(0.0);
  T a10 = T(0.0);
  T a11 = T(1.0);
  T b0 = B(0,0);
  T b1 = B(1,0);
  for (unsigned int i=1; i<=dc; i++) {
    t += basis.getBasisPoly(i).coeff(0) * ca[i];
    a00 += basis.getBasisPoly(i).coeff(0) * B(i-1,1);
    a01 += basis.getBasisPoly(i).coeff(0) * B(i-1,0);
    a10 += basis.getBasisPoly(i).coeff(0) * B(i-1+dc,1);
    a11 += basis.getBasisPoly(i).coeff(0) * B(i-1+dc,0);
  }
  A(0,0) = a00;
  A(0,1) = a01;
  A(1,0) = a10;
  A(1,1) = a11;
  B(0,0) = std::sinh(t);
  B(1,0) = std::cosh(t);
  info = solve(2, 1);
  TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		     func << ":  Argument " << info 
		     << " for (2x2) solve had illegal value");
  TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		     func << ":  Diagonal entry " << info 
		     << " in (2x2) LU factorization is exactly zero");
  cs[0] = B(0,0);
  cc[0] = B(1,0);

  // Compute remaining coefficients
  B(0,0) = b0;
  B(1,0) = b1;
  for (unsigned int i=1; i<=dc; i++) {
    cs[i] = cc[0]*B(i-1,0) + cs[0]*B(i-1,1);
    cc[i] = cc[0]*B(i-1+dc,0) + cs[0]*B(i-1+dc,1);
  }

  cs[0] /= basis.getBasisPoly(0).coeff(0);
  cc[0] /= basis.getBasisPoly(0).coeff(0);
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
sinh(Stokhos::HermitePoly<T>& s, const Stokhos::HermitePoly<T>& a)
{
  unsigned int dc = a.degree();
  HermitePoly<T> c(dc);
  sinhcosh(s, c, a);
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
cosh(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a)
{
  unsigned int dc = a.degree();
  HermitePoly<T> s(dc);
  sinhcosh(s, c, a);
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
tanh(Stokhos::HermitePoly<T>& t, const Stokhos::HermitePoly<T>& a)
{
  unsigned int dc = a.degree();
  HermitePoly<T> c(dc);
  
  sinhcosh(t, c, a);
  divideEqual(t,c);
}

template <typename T, typename BasisT>
template <typename OpT>
void
Stokhos::HermiteExpansion<T,BasisT>::
quad(const OpT& quad_func,
     Stokhos::HermitePoly<T>& c, 
     const Stokhos::HermitePoly<T>& a,
     const Stokhos::HermitePoly<T>& b)
{
  const char* func = "Stokhos::HermiteExpansion::quad()";
  unsigned int dc = a.degree();
  if (dc != c.degree())
    c.resize(dc);

  const T* ca = a.coeff();
  const T* cb = b.coeff();
  T* cc = c.coeff();

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(degree() < dc, std::logic_error,
		     func << ":  Expansion degree (" << degree()
		     << ") is too small for computation (" << dc
		     << " needed).");
#endif

  const typename tp_type::basis_type& basis = Cijk.getBasis();
    
  // Fill A
  for (unsigned int i=1; i<=dc; i++) {
    for (unsigned int j=1; j<=dc; j++) {
      A(i-1,j-1) = 0.0;
      for (unsigned int k=0; k<=dc; k++) {
	A(i-1,j-1) += Cijk.value(i-1,j-1,k)*cb[k];
      }
      A(i-1,j-1) = A(i-1,j-1)*basis.derivCoeff(j) / 
	(Cijk.norm_squared(i-1)*basis.derivCoeff(i));
    }
  }
    
  // Fill B
  for (unsigned int i=1; i<=dc; i++)
    B(i-1,0) = ca[i];
  
  // Solve system
  int info = solve(dc, 1);

  TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		     func << ":  Argument " << info 
		     << " for solve had illegal value");
  TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		     func << ":  Diagonal entry " << info 
		     << " in LU factorization is exactly zero");

  // Compute degree-0 coefficient
  T s = basis.getBasisPoly(0).coeff(0) * ca[0];
  T t = T(0.0);
  for (unsigned int i=1; i<=dc; i++) {
    s += basis.getBasisPoly(i).coeff(0) * ca[i];
    t += basis.getBasisPoly(i).coeff(0) * B(i-1,0);
  }
  cc[0] = (quad_func(s) - t) / basis.getBasisPoly(0).coeff(0);

  // Get coefficients
  for (unsigned int i=1; i<=dc; i++)
    cc[i] = B(i-1,0);
}

template <typename T>
struct acos_quad_func { 
  T operator() (const T& a) const { return std::acos(a); } 
};

template <typename T>
struct asin_quad_func { 
  T operator() (const T& a) const { return std::asin(a); } 
};

template <typename T>
struct atan_quad_func { 
  T operator() (const T& a) const { return std::atan(a); } 
};

template <typename T>
struct acosh_quad_func { 
  T operator() (const T& a) const { return std::log(a+std::sqrt(a*a-T(1.0))); }
};

template <typename T>
struct asinh_quad_func { 
  T operator() (const T& a) const { return std::log(a+std::sqrt(a*a+T(1.0))); }
};

template <typename T>
struct atanh_quad_func { 
  T operator() (const T& a) const { return 0.5*std::log((T(1.0)+a)/(T(1.0)-a)); } 
};

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
acos(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a)
{
  times(c,a,a);
  minus(c,T(1.0),c);
  sqrt(c,c);
  timesEqual(c,T(-1.0));
  quad(acos_quad_func<T>(), c, a, c);
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
asin(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a)
{
  times(c,a,a);
  minus(c,T(1.0),c);
  sqrt(c,c);
  quad(asin_quad_func<T>(), c, a, c);
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
atan(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a)
{
  times(c,a,a);
  plusEqual(c,T(1.0));
  quad(atan_quad_func<T>(), c, a, c);
}

// template <typename T, typename BasisT>
// Hermite<T>
// atan2(const Hermite<T>& a,
//       const Hermite<T>& b)
// {
//   Hermite<T> c = atan(a/b);
//   c.fastAccessCoeff(0) = atan2(a.coeff(0),b.coeff(0));
// }

// template <typename T, typename BasisT>
// Hermite<T>
// atan2(const T& a,
//       const Hermite<T>& b)
// {
//   Hermite<T> c = atan(a/b);
//   c.fastAccessCoeff(0) = atan2(a,b.coeff(0));
// }

// template <typename T, typename BasisT>
// Hermite<T>
// atan2(const Hermite<T>& a,
//       const T& b)
// {
//   Hermite<T> c = atan(a/b);
//   c.fastAccessCoeff(0) = atan2(a.coeff(0),b);
// }

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
acosh(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a)
{
  times(c,a,a);
  minusEqual(c,T(1.0));
  sqrt(c,c);
  quad(acosh_quad_func<T>(), c, a, c);
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
asinh(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a)
{
  times(c,a,a);
  plusEqual(c,T(1.0));
  sqrt(c,c);
  quad(asinh_quad_func<T>(), c, a, c);
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
atanh(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a)
{
  times(c,a,a);
  minus(c,T(1.0),c);
  quad(atanh_quad_func<T>(), c, a, c);
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
fabs(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a)
{
  if (a[0] >= 0)
    c = a;
  else
    unaryMinus(c,a);
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
abs(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a)
{
  if (a[0] >= 0)
    c = a;
  else
    unaryMinus(c,a);
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
max(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a,
    const Stokhos::HermitePoly<T>& b)
{
  if (a[0] >= b[0])
    c = a;
  else
    c = b;
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
max(Stokhos::HermitePoly<T>& c, const T& a, const Stokhos::HermitePoly<T>& b)
{
  if (a >= b[0])
    c = HermitePoly<T>(b.degree(), a);
  else
    c = b;
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
max(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a, const T& b)
{
  if (a[0] >= b)
    c = a;
  else
    c = HermitePoly<T>(a.degree(), b);
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
min(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a,
    const Stokhos::HermitePoly<T>& b)
{
  if (a[0] <= b[0])
    return c = a;
  else
    return c = b;
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
min(Stokhos::HermitePoly<T>& c, const T& a, const Stokhos::HermitePoly<T>& b)
{
  if (a <= b[0])
    c = HermitePoly<T>(b.degree(), a);
  else
    c = b;
}

template <typename T, typename BasisT>
void
Stokhos::HermiteExpansion<T,BasisT>::
min(Stokhos::HermitePoly<T>& c, const Stokhos::HermitePoly<T>& a, const T& b)
{
  if (a[0] <= b)
    c = a;
  else
    c = HermitePoly<T>(a.degree(), b);
}
