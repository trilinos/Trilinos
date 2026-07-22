// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Assert.hpp"
#include "Stokhos_DynamicArrayTraits.hpp"

template <typename ordinal_type, typename value_type> 
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
DerivOrthogPolyExpansion(
  const Teuchos::RCP<const Stokhos::DerivBasis<ordinal_type, value_type> >& basis_,
  const Teuchos::RCP<const Teuchos::SerialDenseMatrix<ordinal_type, value_type> >& Bij_,
      const Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk_,
      const Teuchos::RCP<const Stokhos::Dense3Tensor<ordinal_type, value_type> >& Dijk_) :
  basis(basis_),
  Bij(Bij_),
  Cijk(Cijk_),
  Dijk(Dijk_),
  sz(basis->size()),
  A(2*sz,2*sz),
  B(2*sz,2),
  piv(2*sz),
  lapack()
{
}

extern "C" {
  double dlange_(char*, int*, int*, double*, int*, double*);
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
solve(ordinal_type s, ordinal_type nrhs)
{
  if (s == 0 || nrhs == 0)
    return 0;

  ordinal_type info;
//   lapack.GESV(s, nrhs, A.values(), A.numRows(), &(piv[0]), b.values(), 
// 	      b.numRows(), &info);
  lapack.GETRF(s, s, A.values(), A.numRows(), &(piv[0]), &info);
  value_type norm, rcond;
  Teuchos::Array<ordinal_type> iwork(4*s);
  Teuchos::Array<value_type> work(4*s);
  norm = 1.0;
  ordinal_type n = A.numRows();
  char t = '1';
  norm = dlange_(&t, &s, &s, A.values(), &n, &work[0]);
  lapack.GECON('1', s, A.values(), A.numRows(), norm, &rcond, &work[0], 
	       &iwork[0], &info);
  //std::cout << "condition number = " << 1.0/rcond << std::endl;
  lapack.GETRS('N', s, nrhs, A.values(), A.numRows(), &(piv[0]), B.values(), 
	       B.numRows(), &info);
  return info;
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
unaryMinus(
  Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
  const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  ordinal_type pc = a.size();
  if (c.size() != pc)
    c.resize(pc);

  value_type* cc = c.coeff();
  const value_type* ca = a.coeff();

  for (ordinal_type i=0; i<pc; i++)
    cc[i] = -ca[i];
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
plusEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	  const value_type& val)
{
  c[0] += val;
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
minusEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	   const value_type& val)
{
  c[0] -= val;
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
timesEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	   const value_type& val)
{
  ordinal_type pc = c.size();
  value_type* cc = c.coeff();
  for (ordinal_type i=0; i<pc; i++)
    cc[i] *= val;
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
divideEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	    const value_type& val)
{
  ordinal_type pc = c.size();
  value_type* cc = c.coeff();
  for (ordinal_type i=0; i<pc; i++)
    cc[i] /= val;
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
plusEqual(
  Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
  const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& x)
{
  ordinal_type xp = x.size();
  if (c.size() < xp)
    c.resize(xp);

  value_type* cc = c.coeff();
  const value_type* xc = x.coeff();
  for (ordinal_type i=0; i<xp; i++)
    cc[i] += xc[i];
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
minusEqual(
  Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
  const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& x)
{
  ordinal_type xp = x.size();
  if (c.size() < xp)
    c.resize(xp);

  value_type* cc = c.coeff();
  const value_type* xc = x.coeff();
  for (ordinal_type i=0; i<xp; i++)
    cc[i] -= xc[i];
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
timesEqual(
  Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
  const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& x)
{
  ordinal_type p = c.size();
  ordinal_type xp = x.size();
  ordinal_type pc;
  if (p > 1 && xp > 1)
    pc = sz;
  else
    pc = p*xp;
  TEUCHOS_TEST_FOR_EXCEPTION(sz < pc, std::logic_error,
		     "Stokhos::DerivOrthogPolyExpansion::timesEqual()" <<
		     ":  Expansion size (" << sz << 
		     ") is too small for computation.");
  if (c.size() != pc)
    c.resize(pc);

  value_type* cc = c.coeff();
  const value_type* xc = x.coeff();
  
  if (p > 1 && xp > 1) {
    // Copy c coefficients into temporary array
    value_type* tc = Stokhos::ds_array<value_type>::get_and_fill(cc,p);
    value_type tmp, cijk;
    ordinal_type i,j;
    for (ordinal_type k=0; k<pc; k++) {
      tmp = value_type(0.0);
      ordinal_type n = Cijk->num_values(k);
      for (ordinal_type l=0; l<n; l++) {
	Cijk->value(k,l,i,j,cijk);
	if (i < p && j < xp)
	  tmp += cijk*tc[i]*xc[j];
      }
      cc[k] = tmp / basis->norm_squared(k);
    }
  }
  else if (p > 1) {
    for (ordinal_type i=0; i<p; i++)
      cc[i] *= xc[0];
  }
  else if (xp > 1) {
    for (ordinal_type i=1; i<xp; i++)
      cc[i] = cc[0]*xc[i];
    cc[0] *= xc[0];
  }
  else {
    cc[0] *= xc[0];
  }
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
divideEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	    const OrthogPolyApprox<ordinal_type, value_type, node_type>& x)
{
  const char* func = "Stokhos::DerivOrthogPolyExpansion::divide()";
  ordinal_type p = c.size();
  ordinal_type xp = x.size();
  ordinal_type pc;
  if (xp > 1)
    pc = sz;
  else
    pc = p;
  TEUCHOS_TEST_FOR_EXCEPTION(sz < pc, std::logic_error,
		     "Stokhos::DerivOrthogPolyExpansion::divideEqual()" <<
		     ":  Expansion size (" << sz << 
		     ") is too small for computation.");
  if (c.size() != pc)
    c.resize(pc);

  value_type* cc = c.coeff();
  const value_type* xc = x.coeff();
  
  if (xp > 1) {
    
    // Fill A
    value_type cijk;
    ordinal_type i,j;
    for (ordinal_type i=0; i<pc; i++)
      for (ordinal_type j=0; j<pc; j++)
	A(i,j) = 0.0;
    for (ordinal_type k=0; k<xp; k++) {
      ordinal_type n = Cijk->num_values(k);
      for (ordinal_type l=0; l<n; l++) {
	Cijk->value(k,l,i,j,cijk);
	A(i,j) += cijk*xc[k]/basis->norm_squared(i);
      }
    }
    
    // Fill B
    for (ordinal_type i=0; i<p; i++)
      B(i,0) = cc[i];
    for (ordinal_type i=p; i<pc; i++)
      B(i,0) = value_type(0.0);

    // Solve system
    int info = solve(pc, 1);

    TEUCHOS_TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		       func << ":  Argument " << info 
		            << " for solve had illegal value");
    TEUCHOS_TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		       func << ":  Diagonal entry " << info 
		            << " in LU factorization is exactly zero");

    // Get coefficients
    for (ordinal_type i=0; i<pc; i++)
      cc[i] = B(i,0);
  }
  else {
    for (ordinal_type i=0; i<p; i++)
      cc[i] /= xc[0];
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
plus(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  ordinal_type pa = a.size();
  ordinal_type pb = b.size();
  ordinal_type pc = pa > pb ? pa : pb;
  if (c.size() != pc)
    c.resize(pc);

  const value_type* ca = a.coeff();
  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  if (pa > pb) {
    for (ordinal_type i=0; i<pb; i++)
      cc[i] = ca[i] + cb[i];
    for (ordinal_type i=pb; i<pc; i++)
      cc[i] = ca[i];
  }
  else {
    for (ordinal_type i=0; i<pa; i++)
      cc[i] = ca[i] + cb[i];
    for (ordinal_type i=pa; i<pc; i++)
      cc[i] = cb[i];
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
plus(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const value_type& a, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  ordinal_type pc = b.size();
  if (c.size() != pc)
    c.resize(pc);

  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  cc[0] = a + cb[0];
  for (ordinal_type i=1; i<pc; i++)
    cc[i] = cb[i];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
plus(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
     const value_type& b)
{
  ordinal_type pc = a.size();
  if (c.size() != pc)
    c.resize(pc);

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();

  cc[0] = ca[0] + b;
  for (ordinal_type i=1; i<pc; i++)
    cc[i] = ca[i];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
minus(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  ordinal_type pa = a.size();
  ordinal_type pb = b.size();
  ordinal_type pc = pa > pb ? pa : pb;
  if (c.size() != pc)
    c.resize(pc);

  const value_type* ca = a.coeff();
  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  if (pa > pb) {
    for (ordinal_type i=0; i<pb; i++)
      cc[i] = ca[i] - cb[i];
    for (ordinal_type i=pb; i<pc; i++)
      cc[i] = ca[i];
  }
  else {
    for (ordinal_type i=0; i<pa; i++)
      cc[i] = ca[i] - cb[i];
    for (ordinal_type i=pa; i<pc; i++)
      cc[i] = -cb[i];
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
minus(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const value_type& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  ordinal_type pc = b.size();
  if (c.size() != pc)
    c.resize(pc);

  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  cc[0] = a - cb[0];
  for (ordinal_type i=1; i<pc; i++)
    cc[i] = -cb[i];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
minus(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const value_type& b)
{
  ordinal_type pc = a.size();
  if (c.size() != pc)
    c.resize(pc);

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();

  cc[0] = ca[0] - b;
  for (ordinal_type i=1; i<pc; i++)
    cc[i] = ca[i];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
times(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  ordinal_type pa = a.size();
  ordinal_type pb = b.size();
  ordinal_type pc;
  if (pa > 1 && pb > 1)
    pc = sz;
  else
    pc = pa*pb;
  TEUCHOS_TEST_FOR_EXCEPTION(sz < pc, std::logic_error,
		     "Stokhos::DerivOrthogPolyExpansion::times()" <<
		     ":  Expansion size (" << sz << 
		     ") is too small for computation.");
  if (c.size() != pc)
    c.resize(pc);

  const value_type* ca = a.coeff();
  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  if (pa > 1 && pb > 1) {
    value_type tmp, cijk;
    ordinal_type i,j;
    for (ordinal_type k=0; k<pc; k++) {
      tmp = value_type(0.0);
      ordinal_type n = Cijk->num_values(k);
      for (ordinal_type l=0; l<n; l++) {
    	Cijk->value(k,l,i,j,cijk);
	if (i < pa && j < pb)
	  tmp += cijk*ca[i]*cb[j];
      }
      cc[k] = tmp / basis->norm_squared(k);
    }
  }
  else if (pa > 1) {
    for (ordinal_type i=0; i<pc; i++)
      cc[i] = ca[i]*cb[0];
  }
  else if (pb > 1) {
    for (ordinal_type i=0; i<pc; i++)
      cc[i] = ca[0]*cb[i];
  }
  else {
    cc[0] = ca[0]*cb[0];
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
times(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const value_type& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  ordinal_type pc = b.size();
  if (c.size() != pc)
    c.resize(pc);

  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  for (ordinal_type i=0; i<pc; i++)
    cc[i] = a*cb[i];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
times(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const value_type& b)
{
  ordinal_type pc = a.size();
  if (c.size() != pc)
    c.resize(pc);

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();

  for (ordinal_type i=0; i<pc; i++)
    cc[i] = ca[i]*b;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
divide(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  const char* func = "Stokhos::DerivOrthogPolyExpansion::divide()";
  ordinal_type pa = a.size();
  ordinal_type pb = b.size();
  ordinal_type pc;
  if (pb > 1)
    pc = sz;
  else
    pc = pa;
  TEUCHOS_TEST_FOR_EXCEPTION(sz < pc, std::logic_error,
		     "Stokhos::DerivOrthogPolyExpansion::divide()" <<
		     ":  Expansion size (" << sz << 
		     ") is too small for computation.");
  if (c.size() != pc)
    c.resize(pc);

  const value_type* ca = a.coeff();
  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  if (pb > 1) {

    // Fill A
    value_type cijk;
    ordinal_type i,j;
    for (ordinal_type i=0; i<pc; i++)
      for (ordinal_type j=0; j<pc; j++)
	A(i,j) = 0.0;

    for (ordinal_type k=0; k<pb; k++) {
      ordinal_type n = Cijk->num_values(k);
      for (ordinal_type l=0; l<n; l++) {
	Cijk->value(k,l,i,j,cijk);
	A(i,j) += cijk*cb[k] / basis->norm_squared(i);
      }
    }
    
    // Fill B
    for (ordinal_type i=0; i<pa; i++)
      B(i,0) = ca[i];
    for (ordinal_type i=pa; i<pc; i++)
      B(i,0) = value_type(0.0);

    // Solve system
    int info = solve(pc, 1);

    TEUCHOS_TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		       func << ":  Argument " << info 
		            << " for solve had illegal value");
    TEUCHOS_TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		       func << ":  Diagonal entry " << info 
		            << " in LU factorization is exactly zero");

    // Get coefficients
    for (ordinal_type i=0; i<pc; i++)
      cc[i] = B(i,0);
  }
  else {
    for (ordinal_type i=0; i<pa; i++)
      cc[i] = ca[i]/cb[0];
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
divide(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
       const value_type& a, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  const char* func = "Stokhos::DerivOrthogPolyExpansion::divide()";
  ordinal_type pb = b.size();
  ordinal_type pc;
  if (pb > 1)
    pc = sz;
  else
    pc = 1;
  if (c.size() != pc)
    c.resize(pc);

  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  if (pb > 1) {
   
    // Fill A
    value_type cijk;
    ordinal_type i,j;
    for (ordinal_type i=0; i<pc; i++)
      for (ordinal_type j=0; j<pc; j++)
	A(i,j) = 0.0;

    for (ordinal_type k=0; k<pb; k++) {
      ordinal_type n = Cijk->num_values(k);
      for (ordinal_type l=0; l<n; l++) {
	Cijk->value(k,l,i,j,cijk);
	A(i,j) += cijk*cb[k] / basis->norm_squared(i);
      }
    }

    // Fill B
    B(0,0) = a;
    for (ordinal_type i=1; i<pc; i++)
      B(i,0) = value_type(0.0);

    // Solve system
    int info = solve(pc, 1);

    TEUCHOS_TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		       func << ":  Argument " << info 
		            << " for solve had illegal value");
    TEUCHOS_TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		       func << ":  Diagonal entry " << info 
		            << " in LU factorization is exactly zero");

    // Get coefficients
    for (ordinal_type i=0; i<pc; i++)
      cc[i] = B(i,0);
  }
  else 
    cc[0] = a / cb[0];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
divide(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
       const value_type& b)
{
  ordinal_type pc = a.size();
  if (c.size() != pc)
    c.resize(pc);

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();

  for (ordinal_type i=0; i<pc; i++)
    cc[i] = ca[i]/b;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
exp(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  const char* func = "Stokhos::DerivOrthogPolyExpansion::exp()";
  ordinal_type pa = a.size();
  ordinal_type pc;
  if (pa > 1)
    pc = sz;
  else
    pc = 1;
  if (c.size() != pc)
    c.resize(pc);

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();

  if (pa > 1) {
    value_type psi_0 = basis->evaluateZero(0);

    // Fill A and B
    for (ordinal_type i=1; i<pc; i++) {
      B(i-1,0) = 0.0;
      for (ordinal_type j=1; j<pc; j++) {
	A(i-1,j-1) = (*Bij)(i-1,j);
	for (ordinal_type k=1; k<pa; k++)
	  A(i-1,j-1) -= ca[k]*(*Dijk)(i-1,j,k);
	B(i-1,0) += ca[j]*(*Bij)(i-1,j);
      }
      B(i-1,0) *= psi_0;
    }

    // Solve system
    int info = solve(pc-1, 1);

    TEUCHOS_TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		       func << ":  Argument " << info 
		       << " for solve had illegal value");
    TEUCHOS_TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		       func << ":  Diagonal entry " << info 
		       << " in LU factorization is exactly zero");

    // Compute order-0 coefficient
    value_type s = psi_0 * ca[0];
    value_type t = psi_0;
    for (ordinal_type i=1; i<pc; i++) {
      s += basis->evaluateZero(i) * ca[i];
      t += basis->evaluateZero(i) * B(i-1,0);
    }
    s = std::exp(s);
    cc[0] = (s/t);

    // Compute remaining coefficients
    for (ordinal_type i=1; i<pc; i++)
      cc[i] = B(i-1,0) * cc[0];
  }
  else
    cc[0] = std::exp(ca[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
log(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  const char* func = "Stokhos::DerivOrthogPolyExpansion::log()";
  ordinal_type pa = a.size();
  ordinal_type pc;
  if (pa > 1)
    pc = sz;
  else
    pc = 1;
  if (c.size() != pc)
    c.resize(pc);

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();

  if (pa > 1) {
    value_type psi_0 = basis->evaluateZero(0);

    // Fill A and B
    for (ordinal_type i=1; i<pc; i++) {
      B(i-1,0) = 0.0;
      for (ordinal_type j=1; j<pc; j++) {
	A(i-1,j-1) = 0.0;
	for (ordinal_type k=0; k<pa; k++)
	  A(i-1,j-1) += ca[k]*(*Dijk)(i-1,k,j);
	B(i-1,0) += ca[j]*(*Bij)(i-1,j); 
      }
    }

    // Solve system
    int info = solve(pc-1, 1);
    
    TEUCHOS_TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		       func << ":  Argument " << info 
		       << " for solve had illegal value");
    TEUCHOS_TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		       func << ":  Diagonal entry " << info 
		       << " in LU factorization is exactly zero");
    
    // Compute order-0 coefficient
    value_type s = psi_0 * ca[0];
    value_type t = value_type(0.0);
    for (ordinal_type i=1; i<pc; i++) {
      s += basis->evaluateZero(i) * ca[i];
      t += basis->evaluateZero(i) * B(i-1,0);
    }
    cc[0] = (std::log(s) - t) / psi_0;

    // Compute remaining coefficients
    for (ordinal_type i=1; i<pc; i++)
      cc[i] = B(i-1,0);
  }
  else
    cc[0] = std::log(ca[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
log10(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() > 1) {
    log(c,a);
    divide(c,c,std::log(10.0));
  }
  else {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::log10(a[0]);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
sqrt(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() > 1) {
    log(c,a);
    timesEqual(c,value_type(0.5));
    exp(c,c);
  }
  else {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::sqrt(a[0]);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
cbrt(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() > 1) {
    log(c,a);
    timesEqual(c,value_type(1.0/3.0));
    exp(c,c);
  }
  else {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::cbrt(a[0]);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
pow(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (a.size() > 1 || b.size() > 1) {
    log(c,a);
    timesEqual(c,b);
    exp(c,c);
  }
  else {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::pow(a[0], b[0]);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
pow(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const value_type& a, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (b.size() > 1) {
    times(c,std::log(a),b);
    exp(c,c);
  }
  else {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::pow(a, b[0]);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
pow(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
    const value_type& b)
{
  if (a.size() > 1) {
    log(c,a);
    timesEqual(c,b);
    exp(c,c);
  }
  else {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::pow(a[0], b);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
sincos(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& s, 
       Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  const char* func = "Stokhos::DerivOrthogPolyExpansion::sincos()";
  ordinal_type pa = a.size();
  ordinal_type pc;
  if (pa > 1)
    pc = sz;
  else
    pc = 1;
  if (c.size() != pc)
    c.resize(pc);
  if (s.size() != pc)
    s.resize(pc);

  const value_type* ca = a.coeff();
  value_type* cs = s.coeff();
  value_type* cc = c.coeff();

  if (pa > 1) {
    value_type psi_0 = basis->evaluateZero(0);
    ordinal_type offset = pc-1;

    // Fill A and b
    B.putScalar(value_type(0.0));
    value_type tmp, tmp2;
    for (ordinal_type i=1; i<pc; i++) {
      tmp2 = value_type(0.0);
      for (ordinal_type j=1; j<pc; j++) {
	tmp = (*Bij)(i-1,j);
	A(i-1,j-1) = tmp;
	A(i-1+offset,j-1+offset) = tmp;
	tmp = value_type(0.0);
	for (ordinal_type k=1; k<pa; k++)
	  tmp += ca[k]*(*Dijk)(i-1,j,k);
	A(i-1+offset,j-1) = tmp;
	A(i-1,j-1+offset) = -tmp;
	tmp2 += ca[j]*(*Bij)(i-1,j);
      }
      B(i-1,0) = tmp2*psi_0;
      B(i-1+offset,1) = tmp2*psi_0;
    }

    // Solve system
    int info = solve(2*pc-2, 2);

    TEUCHOS_TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		       func << ":  Argument " << info 
		       << " for solve had illegal value");
    TEUCHOS_TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		       func << ":  Diagonal entry " << info 
		       << " in LU factorization is exactly zero");
    
    // Compute order-0 coefficients
    value_type t = psi_0 * ca[0];
    value_type a00 = psi_0;
    value_type a01 = value_type(0.0);
    value_type a10 = value_type(0.0);
    value_type a11 = psi_0;
    value_type b0 = B(0,0);
    value_type b1 = B(1,0);
    for (ordinal_type i=1; i<pc; i++) {
      t += basis->evaluateZero(i) * ca[i];
      a00 -= basis->evaluateZero(i) * B(i-1,1);
      a01 += basis->evaluateZero(i) * B(i-1,0);
      a10 -= basis->evaluateZero(i) * B(i-1+offset,1);
      a11 += basis->evaluateZero(i) * B(i-1+offset,0);
    }
    A(0,0) = a00;
    A(0,1) = a01;
    A(1,0) = a10;
    A(1,1) = a11;
    B(0,0) = std::sin(t);
    B(1,0) = std::cos(t);
    
    info = solve(2, 1);

    TEUCHOS_TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		       func << ":  Argument " << info 
		       << " for (2x2) solve had illegal value");
    TEUCHOS_TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		       func << ":  Diagonal entry " << info 
		       << " in (2x2) LU factorization is exactly zero");
    cs[0] = B(0,0);
    cc[0] = B(1,0);

    // Compute remaining coefficients
    B(0,0) = b0;
    B(1,0) = b1;
    for (ordinal_type i=1; i<pc; i++) {
      cs[i] = cc[0]*B(i-1,0) - cs[0]*B(i-1,1);
      cc[i] = cc[0]*B(i-1+offset,0) - cs[0]*B(i-1+offset,1);
    }
  }
  else {
    cs[0] = std::sin(ca[0]);
    cc[0] = std::cos(ca[0]);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
sin(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& s, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() > 1) {
    OrthogPolyApprox<ordinal_type, value_type, node_type> c(s);
    sincos(s, c, a);
  }
  else {
    if (s.size() != 1)
      s.resize(1);
    s[0] = std::sin(a[0]);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
cos(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() > 1) {
    OrthogPolyApprox<ordinal_type, value_type, node_type> s(c);
    sincos(s, c, a);
  }
  else {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::cos(a[0]);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
tan(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& t, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() > 1) {
    OrthogPolyApprox<ordinal_type, value_type, node_type> c(t);
    sincos(t, c, a);
    divideEqual(t,c);
  }
  else {
    if (t.size() != 1)
      t.resize(1);
    t[0] = std::tan(a[0]);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
sinhcosh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& s, 
	 Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	 const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  const char* func = "Stokhos::DerivOrthogPolyExpansion::sinhcosh()";
  ordinal_type pa = a.size();
  ordinal_type pc;
  if (pa > 1)
    pc = sz;
  else
    pc = 1;
  if (c.size() != pc)
    c.resize(pc);
  if (s.size() != pc)
    s.resize(pc);

  const value_type* ca = a.coeff();
  value_type* cs = s.coeff();
  value_type* cc = c.coeff();

  if (pa > 1) {
    value_type psi_0 = basis->evaluateZero(0);
    ordinal_type offset = pc-1;

    // Fill A and b
    B.putScalar(value_type(0.0));
    value_type tmp, tmp2;
    for (ordinal_type i=1; i<pc; i++) {
      tmp2 = value_type(0.0);
      for (ordinal_type j=1; j<pc; j++) {
	tmp = (*Bij)(i-1,j);
	A(i-1,j-1) = tmp;
	A(i-1+offset,j-1+offset) = tmp;
	tmp = value_type(0.0);
	for (ordinal_type k=1; k<pa; k++)
	  tmp += ca[k]*(*Dijk)(i-1,j,k);
	A(i-1+offset,j-1) = -tmp;
	A(i-1,j-1+offset) = -tmp;
	tmp2 += ca[j]*(*Bij)(i-1,j);
      }
      B(i-1,0) = tmp2*psi_0;
      B(i-1+offset,1) = tmp2*psi_0;
    }

    // Solve system
    int info = solve(2*pc-2, 2);

    TEUCHOS_TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		       func << ":  Argument " << info 
		       << " for solve had illegal value");
    TEUCHOS_TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		       func << ":  Diagonal entry " << info 
		       << " in LU factorization is exactly zero");
    
    // Compute order-0 coefficients
    value_type t = psi_0 * ca[0];
    value_type a00 = psi_0;
    value_type a01 = value_type(0.0);
    value_type a10 = value_type(0.0);
    value_type a11 = psi_0;
    value_type b0 = B(0,0);
    value_type b1 = B(1,0);
    for (ordinal_type i=1; i<pc; i++) {
      t += basis->evaluateZero(i) * ca[i];
      a00 += basis->evaluateZero(i) * B(i-1,1);
      a01 += basis->evaluateZero(i) * B(i-1,0);
      a10 += basis->evaluateZero(i) * B(i-1+offset,1);
      a11 += basis->evaluateZero(i) * B(i-1+offset,0);
    }
    A(0,0) = a00;
    A(0,1) = a01;
    A(1,0) = a10;
    A(1,1) = a11;
    B(0,0) = std::sinh(t);
    B(1,0) = std::cosh(t);
    info = solve(2, 1);
    TEUCHOS_TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		       func << ":  Argument " << info 
		       << " for (2x2) solve had illegal value");
    TEUCHOS_TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		       func << ":  Diagonal entry " << info 
		       << " in (2x2) LU factorization is exactly zero");
    cs[0] = B(0,0);
    cc[0] = B(1,0);

    // Compute remaining coefficients
    B(0,0) = b0;
    B(1,0) = b1;
    for (ordinal_type i=1; i<pc; i++) {
      cs[i] = cc[0]*B(i-1,0) + cs[0]*B(i-1,1);
      cc[i] = cc[0]*B(i-1+offset,0) + cs[0]*B(i-1+offset,1);
    }
  }
  else {
    cs[0] = std::sinh(ca[0]);
    cc[0] = std::cosh(ca[0]);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
sinh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& s, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() > 1) {
    OrthogPolyApprox<ordinal_type, value_type, node_type> c(s);
    sinhcosh(s, c, a);
  }
  else {
    if (s.size() != 1)
      s.resize(1);
    s[0] = std::sinh(a[0]);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
cosh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() > 1) {
    OrthogPolyApprox<ordinal_type, value_type, node_type> s(c);
    sinhcosh(s, c, a);
  }
  else {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::cosh(a[0]);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
tanh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& t, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() > 1) {
    OrthogPolyApprox<ordinal_type, value_type, node_type> c(t);
    sinhcosh(t, c, a);
    divideEqual(t,c);
  }
  else {
    if (t.size() != 1)
      t.resize(1);
    t[0] = std::tanh(a[0]);
  }
}

template <typename ordinal_type, typename value_type>
template <typename OpT>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
quad(const OpT& quad_func,
     Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  const char* func = "Stokhos::DerivOrthogPolyExpansion::quad()";
  ordinal_type pa = a.size();
  ordinal_type pb = b.size();
  ordinal_type pc = sz;
  if (c.size() != pc)
    c.resize(pc);

  const value_type* ca = a.coeff();
  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  value_type psi_0 = basis->evaluateZero(0);
    
  // Fill A and B
  for (ordinal_type i=1; i<pc; i++) {
    B(i-1,0) = 0.0;
    for (ordinal_type j=1; j<pc; j++) {
      A(i-1,j-1) = 0.0;
      for (ordinal_type k=0; k<pb; k++)
	A(i-1,j-1) += cb[k]*(*Dijk)(i-1,k,j);
    }
    for (ordinal_type j=1; j<pa; j++) {
      B(i-1,0) += ca[j]*(*Bij)(i-1,j);
    }
  }
    
  // Solve system
  int info = solve(pc-1, 1);
  
  TEUCHOS_TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		     func << ":  Argument " << info 
		     << " for solve had illegal value");
  TEUCHOS_TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		     func << ":  Diagonal entry " << info 
		     << " in LU factorization is exactly zero");
  
  // Compute order-0 coefficient
  value_type s = psi_0 * ca[0];
  value_type t = value_type(0.0);
  for (ordinal_type i=1; i<pc; i++) {
    s += basis->evaluateZero(i) * ca[i];
    t += basis->evaluateZero(i) * B(i-1,0);
  }
  cc[0] = (quad_func(s) - t) / psi_0;
  
  // Get coefficients
  for (ordinal_type i=1; i<pc; i++)
    cc[i] = B(i-1,0);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
acos(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() > 1) {
    times(c,a,a);
    minus(c,value_type(1.0),c);
    sqrt(c,c);
    timesEqual(c,value_type(-1.0));
    quad(acos_quad_func(), c, a, c);
  }
  else {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::acos(a[0]);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
asin(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() > 1) {
    times(c,a,a);
    minus(c,value_type(1.0),c);
    sqrt(c,c);
    quad(asin_quad_func(), c, a, c);
  }
  else {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::asin(a[0]);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
atan(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() > 1) {
    times(c,a,a);
    plusEqual(c,value_type(1.0));
    quad(atan_quad_func(), c, a, c);
  }
  else {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::atan(a[0]);
  }
}

// template <typename ordinal_type, typename value_type>
// Hermite<ordinal_type, value_type>
// atan2(const Hermite<ordinal_type, value_type>& a,
//       const Hermite<ordinal_type, value_type>& b)
// {
//   Hermite<ordinal_type, value_type> c = atan(a/b);
//   c.fastAccessCoeff(0) = atan2(a.coeff(0),b.coeff(0));
// }

// template <typename ordinal_type, typename value_type>
// Hermite<ordinal_type, value_type>
// atan2(const T& a,
//       const Hermite<ordinal_type, value_type>& b)
// {
//   Hermite<ordinal_type, value_type> c = atan(a/b);
//   c.fastAccessCoeff(0) = atan2(a,b.coeff(0));
// }

// template <typename ordinal_type, typename value_type>
// Hermite<ordinal_type, value_type>
// atan2(const Hermite<ordinal_type, value_type>& a,
//       const T& b)
// {
//   Hermite<ordinal_type, value_type> c = atan(a/b);
//   c.fastAccessCoeff(0) = atan2(a.coeff(0),b);
// }

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
acosh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() > 1) {
    times(c,a,a);
    minusEqual(c,value_type(1.0));
    sqrt(c,c);
    quad(acosh_quad_func(), c, a, c);
  }
  else {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::log(a[0]+std::sqrt(a[0]*a[0]-value_type(1.0)));
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
asinh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() > 1) {
    times(c,a,a);
    plusEqual(c,value_type(1.0));
    sqrt(c,c);
    quad(asinh_quad_func(), c, a, c);
  }
  else {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::log(a[0]+std::sqrt(a[0]*a[0]+value_type(1.0)));
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
atanh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() > 1) {
    times(c,a,a);
    minus(c,value_type(1.0),c);
    quad(atanh_quad_func(), c, a, c);
  }
  else {
    if (c.size() != 1)
      c.resize(1);
    c[0] = 0.5*std::log((value_type(1.0)+a[0])/(value_type(1.0)-a[0]));
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
fabs(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a[0] >= 0)
    c = a;
  else
    unaryMinus(c,a);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
abs(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a[0] >= 0)
    c = a;
  else
    unaryMinus(c,a);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
max(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (a[0] >= b[0])
    c = a;
  else
    c = b;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
max(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const value_type& a, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (a >= b[0]) {
    c = OrthogPolyApprox<ordinal_type, value_type, node_type>(basis);
    c[0] = a;
  }
  else
    c = b;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
max(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
    const value_type& b)
{
  if (a[0] >= b)
    c = a;
  else {
    c = OrthogPolyApprox<ordinal_type, value_type, node_type>(basis);
    c[0] = b;
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
min(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (a[0] <= b[0])
    c = a;
  else
    c = b;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
min(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const value_type& a, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (a <= b[0]) {
    c = OrthogPolyApprox<ordinal_type, value_type, node_type>(basis);
    c[0] = a;
  }
  else
    c = b;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
min(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
    const value_type& b)
{
  if (a[0] <= b)
    c = a;
  else {
    c = OrthogPolyApprox<ordinal_type, value_type, node_type>(basis);
    c[0] = b;
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::DerivOrthogPolyExpansion<ordinal_type, value_type>::
derivative(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  ordinal_type pc = a.size();

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();

  for (ordinal_type i=0; i<pc; i++) {
    cc[i] = 0.0;
    for (ordinal_type j=0; j<pc; j++)
      cc[i] += ca[j]*(*Bij)(i,j);
    cc[i] /= basis->norm_squared(i);
  }
}
