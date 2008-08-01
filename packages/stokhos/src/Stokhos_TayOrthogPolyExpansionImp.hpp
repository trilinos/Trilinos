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

extern "C" {
  void uq_prep_(int*, int*, int*);
  void uq_prod2_(const double*, const double*, double*, int*);
  void uq_div_(const double*, const double*, double*, int*);
  void uq_exp_(const double*, double*, int*, int*, double*, int*);
  void uq_log_(const double*, double*, int*, int*, double*, int*);
  void uq_sqrt_(const double*, double*, int*, int*);
  void uq_exp_int_(const double*, double*, int*);
  void uq_log_int_(const double*, double*, int*);
}

template <typename T> 
Stokhos::TayOrthogPolyExpansion<T>::
TayOrthogPolyExpansion(
	     const Teuchos::RCP<const Stokhos::OrthogPolyBasis<T> >& basis_) :
  basis(basis_),
  Cijk(),
  rtol(1.0e-12)
{
  order = basis->order();
  dim = basis->dimension();
  int nup;

  std::cout << "starting uq_prep..." << std::endl;
  uq_prep_(&order, &dim, &nup);
  std::cout << "finished uq_prep..." << std::endl;
  sz = nup+1;
}

template <typename T> 
void
Stokhos::TayOrthogPolyExpansion<T>::
unaryMinus(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
	   const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a)
{
  value_type* cc = c.coeff();
  const value_type* ca = a.coeff();
  unsigned int pa = a.size();

  for (unsigned int i=0; i<pa; i++)
    cc[i] = -ca[i];
}

template <typename T> 
void
Stokhos::TayOrthogPolyExpansion<T>::
plusEqual(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, const typename Stokhos::TayOrthogPolyExpansion<T>::value_type& val)
{
  c[0] += val;
}

template <typename T> 
void
Stokhos::TayOrthogPolyExpansion<T>::
minusEqual(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, const typename Stokhos::TayOrthogPolyExpansion<T>::value_type& val)
{
  c[0] -= val;
}

template <typename T> 
void
Stokhos::TayOrthogPolyExpansion<T>::
timesEqual(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, const typename Stokhos::TayOrthogPolyExpansion<T>::value_type& val)
{
  unsigned int pc = c.size();
  value_type* cc = c.coeff();
  for (unsigned int i=0; i<pc; i++)
    cc[i] *= val;
}

template <typename T> 
void
Stokhos::TayOrthogPolyExpansion<T>::
divideEqual(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, const typename Stokhos::TayOrthogPolyExpansion<T>::value_type& val)
{
  unsigned int pc = c.size();
  value_type* cc = c.coeff();
  for (unsigned int i=0; i<pc; i++)
    cc[i] /= val;
}

template <typename T> 
void
Stokhos::TayOrthogPolyExpansion<T>::
plusEqual(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
	  const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& x)
{
  unsigned int p = c.size();
  unsigned int xp = x.size();
  unsigned int pmin = xp < p ? xp : p;
  if (xp > p)
    c.resize(xp);

  value_type* cc = c.coeff();
  const value_type* xc = x.coeff();
  for (unsigned int i=0; i<pmin; i++)
    cc[i] += xc[i];
  if (p < xp)
    for (unsigned int i=p; i<xp; i++)
      cc[i] = xc[i];
}

template <typename T> 
void
Stokhos::TayOrthogPolyExpansion<T>::
minusEqual(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
	   const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& x)
{
  unsigned int p = c.size();
  unsigned int xp = x.size();
  unsigned int pmin = xp < p ? xp : p;
  if (xp > p)
    c.resize(xp);

  value_type* cc = c.coeff();
  const value_type* xc = x.coeff();
  for (unsigned int i=0; i<pmin; i++)
    cc[i] -= xc[i];
  if (p < xp)
    for (unsigned int i=p; i<xp; i++)
      cc[i] = -xc[i];
}

template <typename T> 
void
Stokhos::TayOrthogPolyExpansion<T>::
timesEqual(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
	   const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& x)
{
#ifdef STOKHOS_DEBUG
  const char* func = "Stokhos::TayOrthogPolyExpansion::timesEqual()";
  TEST_FOR_EXCEPTION((x.size() != c.size()) && (x.size() != 1) && 
		     (c.size() != 1), std::logic_error,
		     func << ":  Attempt to assign with incompatible orders");
#endif
  unsigned int p = c.size();
  unsigned int xp = x.size();

  value_type* cc = c.coeff();
  const value_type* xc = x.coeff();
  
  if (p > 1 && xp > 1) {
#ifdef STOKHOS_DEBUG
    TEST_FOR_EXCEPTION(size() < p, std::logic_error,
		       func << ":  Expansion size (" << size() 
		       << ") is too small for computation (" << p 
		       << " needed).");
#endif

    // Copy c coefficients into temporary array
    value_type* tc = Stokhos::ds_array<value_type>::get_and_fill(cc,p);
   
    int nup = p-1;
    uq_prod2_(cc, xc, tc, &nup);
    Stokhos::ds_array<value_type>::copy(tc, cc, p);
  }
  else if (p > 1) {
    for (unsigned int i=0; i<p; i++)
      cc[i] *= xc[0];
  }
  else  {
    if (xp > p)
      c.resize(xp);
    for (int i=xp-1; i>=0; i--)
	cc[i] = cc[0]*xc[i];
  }
}

template <typename T> 
void
Stokhos::TayOrthogPolyExpansion<T>::
divideEqual(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
	    const OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& x)
{
  const char* func = "Stokhos::TayOrthogPolyExpansion::divideEquals()";

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION((x.size() != c.size()) && (x.size() != 1) && 
		     (c.size() != 1), std::logic_error,
		     func << ":  Attempt to assign with incompatible sizes");
#endif

  unsigned int p = c.size();
  unsigned int xp = x.size();

  value_type* cc = c.coeff();
  const value_type* xc = x.coeff();
  
  if (xp > 1) {
    if (xp > p)
      c.resize(xp);

#ifdef STOKHOS_DEBUG
    TEST_FOR_EXCEPTION(size() < xp, std::logic_error,
		       func << ":  Expansion size (" << size() 
		       << ") is too small for computation (" << xp
		       << " needed).");
#endif
    
    // Copy c coefficients into temporary array
    value_type* tc = Stokhos::ds_array<value_type>::get_and_fill(cc,p);
   
    int nup = xp-1;
    uq_div_(cc, xc, tc, &nup);
    Stokhos::ds_array<value_type>::copy(tc, cc, xp);
    
  }
  else {
    for (unsigned int i=0; i<p; i++)
      cc[i] /= xc[0];
  }
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
plus(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
     const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a, 
     const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& b)
{
#ifdef STOKHOS_DEBUG
  const char* func = "Stokhos::TayOrthogPolyExpansion::plus()";
  TEST_FOR_EXCEPTION((a.size() != b.size()) && (a.size() != 1) && 
		     (b.size() != 1), 
		     std::logic_error,
		     func << ":  Arguments have incompatible sizes!");
#endif

  unsigned int pa = a.size();
  unsigned int pb = b.size();
  unsigned int pc = pa > pb ? pa : pb;
  if (pc != c.size())
    c.resize(pc);

  const value_type* ca = a.coeff();
  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  if (pa > 1 && pb > 1) {
    for (unsigned int i=0; i<pc; i++)
      cc[i] = ca[i] + cb[i];
  }
  else if (pa > 1) {
    cc[0] = ca[0] + cb[0];
    for (unsigned int i=1; i<pc; i++)
      cc[i] = ca[i];
  }
  else if (pb >= 1) {
    cc[0] = ca[0] + cb[0];
    for (unsigned int i=1; i<pc; i++)
      cc[i] = cb[i];
  }
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
plus(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
     const typename Stokhos::TayOrthogPolyExpansion<T>::value_type& a, 
     const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& b)
{
  unsigned int pc = b.size();
  if (pc != c.size())
    c.resize(pc);

  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  cc[0] = a + cb[0];
  for (unsigned int i=1; i<pc; i++)
    cc[i] = cb[i];
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
plus(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
     const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a, 
     const typename Stokhos::TayOrthogPolyExpansion<T>::value_type& b)
{
  unsigned int pc = a.size();
  if (pc != c.size())
    c.resize(pc);

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();

  cc[0] = ca[0] + b;
  for (unsigned int i=1; i<pc; i++)
    cc[i] = ca[i];
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
minus(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
      const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a, 
      const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& b)
{
#ifdef STOKHOS_DEBUG
  const char* func = "Stokhos::TayOrthogPolyExpansion::minus()";
  TEST_FOR_EXCEPTION((a.size() != b.size()) && (a.size() != 1) && 
		     (b.size() != 1), 
		     std::logic_error,
		     func << ":  Arguments have incompatible sizes!");
#endif

  unsigned int pa = a.size();
  unsigned int pb = b.size();
  unsigned int pc = pa > pb ? pa : pb;
  if (pc != c.size())
    c.resize(pc);

  const value_type* ca = a.coeff();
  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  if (pa > 1 && pb > 1) {
    for (unsigned int i=0; i<pc; i++)
      cc[i] = ca[i] - cb[i];
  }
  else if (pa > 1) {
    cc[0] = ca[0] - cb[0];
    for (unsigned int i=1; i<pc; i++)
      cc[i] = ca[i];
  }
  else if (pb >= 1) {
    cc[0] = ca[0] - cb[0];
    for (unsigned int i=1; i<pc; i++)
      cc[i] = -cb[i];
  }
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
minus(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
      const typename Stokhos::TayOrthogPolyExpansion<T>::value_type& a, 
      const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& b)
{
  unsigned int pc = b.size();
  if (pc != c.size())
    c.resize(pc);

  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  cc[0] = a - cb[0];
  for (unsigned int i=1; i<pc; i++)
    cc[i] = -cb[i];
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
minus(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
      const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a, const typename Stokhos::TayOrthogPolyExpansion<T>::value_type& b)
{
  unsigned int pc = a.size();
  if (pc != c.size())
    c.resize(pc);

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();

  cc[0] = ca[0] - b;
  for (unsigned int i=1; i<pc; i++)
    cc[i] = ca[i];
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
times(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
      const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a, 
      const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& b)
{
#ifdef STOKHOS_DEBUG
  const char* func = "Stokhos::TayOrthogPolyExpansion::times()";
  TEST_FOR_EXCEPTION((a.size() != b.size()) && (a.size() != 1) && 
		     (b.size() != 1), 
		     std::logic_error,
		     func << ":  Arguments have incompatible sizes!");
#endif

  unsigned int pa = a.size();
  unsigned int pb = b.size();
  unsigned int pc = pa > pb ? pa : pb;
  if (pc != c.size())
    c.resize(pc);

  const value_type* ca = a.coeff();
  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  if (pa > 1 && pb > 1) {
#ifdef STOKHOS_DEBUG
    TEST_FOR_EXCEPTION(size() < pc, std::logic_error,
		       func << ":  Expansion size (" << size()
		       << ") is too small for computation (" << pc
		       << " needed).");
#endif

    int nup = pc-1;
    uq_prod2_(ca, cb, cc, &nup);
  }
  else if (pa > 1) {
    for (unsigned int i=0; i<pc; i++)
      cc[i] = ca[i]*cb[0];
  }
  else if (pb > 1) {
    for (unsigned int i=0; i<pc; i++)
      cc[i] = ca[0]*cb[i];
  }
  else {
    cc[0] = ca[0]*cb[0];
  }
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
times(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
      const typename Stokhos::TayOrthogPolyExpansion<T>::value_type& a, 
      const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& b)
{
  unsigned int pc = b.size();
  if (pc != c.size())
    c.resize(pc);

  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  for (unsigned int i=0; i<pc; i++)
    cc[i] = a*cb[i];
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
times(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
      const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a, 
      const typename Stokhos::TayOrthogPolyExpansion<T>::value_type& b)
{
  unsigned int pc = a.size();
  if (pc != c.size())
    c.resize(pc);

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();

  for (unsigned int i=0; i<pc; i++)
    cc[i] = ca[i]*b;
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
divide(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
       const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a, 
       const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& b)
{
  const char* func = "Stokhos::TayOrthogPolyExpansion::divide()";

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION((a.size() != b.size()) && (a.size() != 1) && 
		     (b.size() != 1), 
		     std::logic_error,
		     func << ":  Arguments have incompatible sizes!");
#endif

  unsigned int pa = a.size();
  unsigned int pb = b.size();
  unsigned int pc = pa > pb ? pa : pb;
  if (pc != c.size())
    c.resize(pc);

  const value_type* ca = a.coeff();
  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  if (pb > 1) {

#ifdef STOKHOS_DEBUG
    TEST_FOR_EXCEPTION(size() < pc, std::logic_error,
		       func << ":  Expansion size (" << size() 
		       << ") is too small for computation (" << pc
		       << " needed).");
#endif

    int nup = pc-1;
    uq_div_(ca, cb, cc, &nup);
  }
  else {
    for (unsigned int i=0; i<pa; i++)
      cc[i] = ca[i]/cb[0];
  }
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
divide(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
       const typename Stokhos::TayOrthogPolyExpansion<T>::value_type& a, 
       const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& b)
{
  const char* func = "Stokhos::TayOrthogPolyExpansion::divide()";

  unsigned int pc = b.size();
  if (pc != c.size())
    c.resize(pc);

  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  if (pc > 1) {
#ifdef STOKHOS_DEBUG
    TEST_FOR_EXCEPTION(size() < pc, std::logic_error,
		       func << ":  Expansion size (" << size()
		       << ") is too small for computation (" << pc
		       << " needed).");
#endif
   
    value_type* ca = Stokhos::ds_array<value_type>::get_and_fill(pc);
    cc[0] = a;
    int nup = pc-1;
    uq_div_(ca, cb, cc, &nup);
  }
  else 
    cc[0] = a / cb[0];
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
divide(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
       const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a, 
       const typename Stokhos::TayOrthogPolyExpansion<T>::value_type& b)
{
  unsigned int pc = a.size();
  if (pc != c.size())
    c.resize(pc);

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();

  for (unsigned int i=0; i<pc; i++)
    cc[i] = ca[i]/b;
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
exp(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
    const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a)
{
  unsigned int pc = a.size();
  if (pc != c.size())
    c.resize(pc);

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();

  int nrm = 1;
  int nup = pc-1;
  uq_exp_(ca, cc, &dim, &nup, &rtol, &nrm);
  //uq_exp_int_(ca, cc, &nup);
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
log(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
    const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a)
{
  unsigned int pc = a.size();
  if (pc != c.size())
    c.resize(pc);

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();

  int nrm = 1;
  int nup = pc-1;
  //uq_log_(ca, cc, &dim, &nup, &rtol, &nrm);
  uq_log_int_(ca, cc, &nup);
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
log10(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
      const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a)
{
  log(c,a);
  divide(c,c,std::log(10.0));
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
sqrt(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
     const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a)
{
  unsigned int pc = a.size();
  if (pc != c.size())
    c.resize(pc);

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();
  int iguess = 0;

  int nup = pc-1;
  uq_sqrt_(ca, cc, &nup, &iguess);
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
pow(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
    const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a,
    const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& b)
{
  log(c,a);
  timesEqual(c,b);
  exp(c,c);
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
pow(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
    const typename Stokhos::TayOrthogPolyExpansion<T>::value_type& a, 
    const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& b)
{
  times(c,std::log(a),b);
  exp(c,c);
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
pow(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
    const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a, 
    const typename Stokhos::TayOrthogPolyExpansion<T>::value_type& b)
{
  log(c,a);
  timesEqual(c,b);
  exp(c,c);
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
sin(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& s, 
    const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a)
{
  throw "Not Defined";
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
cos(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
    const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a)
{
  throw "Not Defined";
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
tan(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& t, 
    const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a)
{
  throw "Not Defined";
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
sinh(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& s, 
     const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a)
{
  throw "Not Defined";
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
cosh(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
     const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a)
{
  throw "Not Defined";
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
tanh(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& t, 
     const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a)
{
  throw "Not Defined";
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
acos(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
     const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a)
{
  throw "Not Defined";
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
asin(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
     const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a)
{
  throw "Not Defined";
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
atan(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
     const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a)
{
  throw "Not Defined";
}

// template <typename T>
// Hermite<value_type>
// atan2(const Hermite<value_type>& a,
//       const Hermite<value_type>& b)
// {
//   Hermite<value_type> c = atan(a/b);
//   c.fastAccessCoeff(0) = atan2(a.coeff(0),b.coeff(0));
// }

// template <typename T>
// Hermite<value_type>
// atan2(const T& a,
//       const Hermite<value_type>& b)
// {
//   Hermite<value_type> c = atan(a/b);
//   c.fastAccessCoeff(0) = atan2(a,b.coeff(0));
// }

// template <typename T>
// Hermite<value_type>
// atan2(const Hermite<value_type>& a,
//       const T& b)
// {
//   Hermite<value_type> c = atan(a/b);
//   c.fastAccessCoeff(0) = atan2(a.coeff(0),b);
// }

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
acosh(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
      const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a)
{
  throw "Not Defined";
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
asinh(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
      const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a)
{
  throw "Not Defined";
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
atanh(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
      const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a)
{
  throw "Not Defined";
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
fabs(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
     const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a)
{
  if (a[0] >= 0)
    c = a;
  else
    unaryMinus(c,a);
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
abs(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
    const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a)
{
  if (a[0] >= 0)
    c = a;
  else
    unaryMinus(c,a);
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
max(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
    const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a,
    const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& b)
{
  if (a[0] >= b[0])
    c = a;
  else
    c = b;
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
max(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
    const typename Stokhos::TayOrthogPolyExpansion<T>::value_type& a, 
    const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& b)
{
  if (a >= b[0])
    c = OrthogPolyApprox<value_type>(b.size(), a);
  else
    c = b;
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
max(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
    const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a, 
    const typename Stokhos::TayOrthogPolyExpansion<T>::value_type& b)
{
  if (a[0] >= b)
    c = a;
  else
    c = OrthogPolyApprox<value_type>(a.size(), b);
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
min(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
    const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a,
    const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& b)
{
  if (a[0] <= b[0])
    return c = a;
  else
    return c = b;
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
min(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
    const typename Stokhos::TayOrthogPolyExpansion<T>::value_type& a, 
    const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& b)
{
  if (a <= b[0])
    c = OrthogPolyApprox<value_type>(b.size(), a);
  else
    c = b;
}

template <typename T>
void
Stokhos::TayOrthogPolyExpansion<T>::
min(Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& c, 
    const Stokhos::OrthogPolyApprox<typename Stokhos::TayOrthogPolyExpansion<T>::value_type >& a, 
    const typename Stokhos::TayOrthogPolyExpansion<T>::value_type& b)
{
  if (a[0] <= b)
    c = a;
  else
    c = OrthogPolyApprox<value_type>(a.size(), b);
}
