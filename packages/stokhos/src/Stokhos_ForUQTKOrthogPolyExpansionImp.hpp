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
#include "Teuchos_ConfigDefs.hpp"

#define UQ_PREP_F77 F77_FUNC_(uq_prep,UQ_PREP)
#define UQ_PROD2_F77 F77_FUNC_(uq_prod2,UQ_PROD2)
#define UQ_DIV_F77 F77_FUNC_(uq_div,UQ_DIV)
#define UQ_EXP_F77 F77_FUNC_(uq_exp,UQ_EXP)
#define UQ_LOG_F77 F77_FUNC_(uq_log,UQ_LOG)
#define UQ_SQRT_F77 F77_FUNC_(uq_sqrt,UQ_SQRT)
#define UQ_EXP_INT_F77 F77_FUNC_(uq_exp_int,UQ_EXP)
#define UQ_LOG_INT_F77 F77_FUNC_(uq_log_int,UQ_LOG)

extern "C" {
  void UQ_PREP_F77(int*, int*, int*);
  void UQ_PROD2_F77(const double*, const double*, double*, int*);
  void UQ_DIV_F77(const double*, const double*, double*, int*);
  void UQ_EXP_F77(const double*, double*, int*, int*, double*, int*);
  void UQ_LOG_F77(const double*, double*, int*, int*, double*, int*);
  void UQ_SQRT_F77(const double*, double*, int*, int*);
  void UQ_EXP_INT_F77(const double*, double*, int*);
  void UQ_LOG_INT_F77(const double*, double*, int*);
}

template <typename ordinal_type, typename value_type> 
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
ForUQTKOrthogPolyExpansion(
			   const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis_,
			   EXPANSION_METHOD method_,
			   value_type rtol_) :
  basis(basis_),
  rtol(rtol_),
  method(method_)
{
  order = basis->order();
  dim = basis->dimension();
  int nup;
  UQ_PREP_F77(&order, &dim, &nup);
  sz = nup+1;
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
unaryMinus(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
	   const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a)
{
  value_type* cc = c.coeff();
  const value_type* ca = a.coeff();
  ordinal_type pa = a.size();

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(c.size() != pa, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::unaryMinus()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "a.size() = " << pa << ", c.size() = " << c.size()
		     << ".");
#endif

  for (ordinal_type i=0; i<pa; i++)
    cc[i] = -ca[i];
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
plusEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, const value_type& val)
{
  c[0] += val;
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
minusEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, const value_type& val)
{
  c[0] -= val;
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
timesEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, const value_type& val)
{
  ordinal_type pc = c.size();
  value_type* cc = c.coeff();
  for (ordinal_type i=0; i<pc; i++)
    cc[i] *= val;
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
divideEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, const value_type& val)
{
  ordinal_type pc = c.size();
  value_type* cc = c.coeff();
  for (ordinal_type i=0; i<pc; i++)
    cc[i] /= val;
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
plusEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
	  const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& x)
{
  ordinal_type p = c.size();
  ordinal_type xp = x.size();
#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(p < xp, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::plusEqual()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "x.size() = " << xp << ", c.size() = " << p << ".");
#endif

  value_type* cc = c.coeff();
  const value_type* xc = x.coeff();
  for (ordinal_type i=0; i<xp; i++)
    cc[i] += xc[i];
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
minusEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
	   const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& x)
{
  ordinal_type p = c.size();
  ordinal_type xp = x.size();
#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(p < xp, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::minusEqual()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "x.size() = " << xp << ", c.size() = " << p << ".");
#endif

  value_type* cc = c.coeff();
  const value_type* xc = x.coeff();
  for (ordinal_type i=0; i<xp; i++)
    cc[i] -= xc[i];
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
timesEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
	   const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& x)
{
  ordinal_type p = c.size();
  ordinal_type xp = x.size();
#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(p < xp, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::timesEqual()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "x.size() = " << xp << ", c.size() = " << p << ".");
#endif

  value_type* cc = c.coeff();
  const value_type* xc = x.coeff();
  
  if (p > 1 && xp > 1) {
#ifdef STOKHOS_DEBUG
    TEST_FOR_EXCEPTION(p != xp, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::timesEqual()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "x.size() = " << xp << ", c.size() = " << p << ".");
    TEST_FOR_EXCEPTION(size() < p, std::logic_error,
                       "Stokhos::ForUQTKOrthogPolyExpansion::timesEqual()" 
		       << ":  Expansion size (" << size() 
                       << ") is too small for computation (" << p 
                       << " needed).");
#endif

    // Copy c coefficients into temporary array
    value_type* tc = Stokhos::ds_array<value_type>::get_and_fill(cc,p);
   
    int nup = p-1;
    UQ_PROD2_F77(cc, xc, tc, &nup);
    Stokhos::ds_array<value_type>::copy(tc, cc, p);
  }
  else {
    for (ordinal_type i=0; i<p; i++)
      cc[i] *= xc[0];
  }
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
divideEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
	    const OrthogPolyApprox<ordinal_type, value_type>& x)
{
  ordinal_type p = c.size();
  ordinal_type xp = x.size();
#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(p < xp, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::divideEqual()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "x.size() = " << xp << ", c.size() = " << p << ".");
#endif

  value_type* cc = c.coeff();
  const value_type* xc = x.coeff();
  
  if (xp > 1) {

#ifdef STOKHOS_DEBUG
    TEST_FOR_EXCEPTION(p != xp, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::divideEqual()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "x.size() = " << xp << ", c.size() = " << p << ".");
    TEST_FOR_EXCEPTION(size() < xp, std::logic_error,
		       "Stokhos::ForUQTKOrthogPolyExpansion::divideEqual()" 
		       << ":  Expansion size (" << size() 
		       << ") is too small for computation (" << xp
		       << " needed).");
#endif
    
    // Copy c coefficients into temporary array
    value_type* tc = Stokhos::ds_array<value_type>::get_and_fill(cc,p);
   
    int nup = xp-1;
    UQ_DIV_F77(tc, xc, cc, &nup);
    
  }
  else {
    for (ordinal_type i=0; i<p; i++)
      cc[i] /= xc[0];
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
plus(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& b)
{
#ifdef STOKHOS_DEBUG
  const char* func = "Stokhos::ForUQTKOrthogPolyExpansion::plus()";
  TEST_FOR_EXCEPTION((a.size() != b.size()) && (a.size() != 1) && 
		     (b.size() != 1), 
		     std::logic_error,
		     func << ":  Arguments have incompatible sizes!");
#endif

  ordinal_type pa = a.size();
  ordinal_type pb = b.size();
  ordinal_type pc = pa > pb ? pa : pb;

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(c.size() != pc, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::plus()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "a.size() = " << pa << ", b.size() = " << pb
		     << ", c.size() = " << c.size() << ".");
#endif

  const value_type* ca = a.coeff();
  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  if (pa > 1 && pb > 1) {
    for (ordinal_type i=0; i<pc; i++)
      cc[i] = ca[i] + cb[i];
  }
  else if (pa > 1) {
    cc[0] = ca[0] + cb[0];
    for (ordinal_type i=1; i<pc; i++)
      cc[i] = ca[i];
  }
  else if (pb >= 1) {
    cc[0] = ca[0] + cb[0];
    for (ordinal_type i=1; i<pc; i++)
      cc[i] = cb[i];
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
plus(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
     const value_type& a, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& b)
{
  ordinal_type pc = b.size();

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(c.size() != pc, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::plus()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "b.size() = " << pc << ", c.size() = " << c.size()
		     << ".");
#endif

  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  cc[0] = a + cb[0];
  for (ordinal_type i=1; i<pc; i++)
    cc[i] = cb[i];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
plus(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a, 
     const value_type& b)
{
  ordinal_type pc = a.size();

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(c.size() != pc, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::plus()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "a.size() = " << pc << ", c.size() = " << c.size()
		     << ".");
#endif

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();

  cc[0] = ca[0] + b;
  for (ordinal_type i=1; i<pc; i++)
    cc[i] = ca[i];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
minus(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& b)
{
#ifdef STOKHOS_DEBUG
  const char* func = "Stokhos::ForUQTKOrthogPolyExpansion::minus()";
  TEST_FOR_EXCEPTION((a.size() != b.size()) && (a.size() != 1) && 
		     (b.size() != 1), 
		     std::logic_error,
		     func << ":  Arguments have incompatible sizes!");
#endif

  ordinal_type pa = a.size();
  ordinal_type pb = b.size();
  ordinal_type pc = pa > pb ? pa : pb;

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(c.size() != pc, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::minus()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "a.size() = " << pa << ", b.size() = " << pb
		     << ", c.size() = " << c.size() << ".");
#endif

  const value_type* ca = a.coeff();
  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  if (pa > 1 && pb > 1) {
    for (ordinal_type i=0; i<pc; i++)
      cc[i] = ca[i] - cb[i];
  }
  else if (pa > 1) {
    cc[0] = ca[0] - cb[0];
    for (ordinal_type i=1; i<pc; i++)
      cc[i] = ca[i];
  }
  else if (pb >= 1) {
    cc[0] = ca[0] - cb[0];
    for (ordinal_type i=1; i<pc; i++)
      cc[i] = -cb[i];
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
minus(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
      const value_type& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& b)
{
  ordinal_type pc = b.size();

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(c.size() != pc, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::minus()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "b.size() = " << pc << ", c.size() = " << c.size()
		     << ".");
#endif

  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  cc[0] = a - cb[0];
  for (ordinal_type i=1; i<pc; i++)
    cc[i] = -cb[i];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
minus(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a, const value_type& b)
{
  ordinal_type pc = a.size();

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(c.size() != pc, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::minus()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "a.size() = " << pc << ", c.size() = " << c.size()
		     << ".");
#endif

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();

  cc[0] = ca[0] - b;
  for (ordinal_type i=1; i<pc; i++)
    cc[i] = ca[i];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
times(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& b)
{
#ifdef STOKHOS_DEBUG
  const char* func = "Stokhos::ForUQTKOrthogPolyExpansion::times()";
  TEST_FOR_EXCEPTION((a.size() != b.size()) && (a.size() != 1) && 
		     (b.size() != 1), 
		     std::logic_error,
		     func << ":  Arguments have incompatible sizes!");
#endif

  ordinal_type pa = a.size();
  ordinal_type pb = b.size();
  ordinal_type pc = pa > pb ? pa : pb;

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(c.size() != pc, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::times()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "a.size() = " << pa << ", b.size() = " << pb
		     << ", c.size() = " << c.size() << ".");
#endif

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
    UQ_PROD2_F77(ca, cb, cc, &nup);
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
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
times(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
      const value_type& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& b)
{
  ordinal_type pc = b.size();

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(c.size() != pc, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::times()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "b.size() = " << pc << ", c.size() = " << c.size()
		     << ".");
#endif

  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  for (ordinal_type i=0; i<pc; i++)
    cc[i] = a*cb[i];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
times(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a, 
      const value_type& b)
{
  ordinal_type pc = a.size();

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(c.size() != pc, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::times()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "a.size() = " << pc << ", c.size() = " << c.size()
		     << ".");
#endif

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();

  for (ordinal_type i=0; i<pc; i++)
    cc[i] = ca[i]*b;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
divide(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& b)
{
  const char* func = "Stokhos::ForUQTKOrthogPolyExpansion::divide()";

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION((a.size() != b.size()) && (a.size() != 1) && 
		     (b.size() != 1), 
		     std::logic_error,
		     func << ":  Arguments have incompatible sizes!");
#endif

  ordinal_type pa = a.size();
  ordinal_type pb = b.size();
  ordinal_type pc = pa > pb ? pa : pb;

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(c.size() != pc, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::divide()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "a.size() = " << pa << ", b.size() = " << pb
		     << ", c.size() = " << c.size() << ".");
#endif

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
    UQ_DIV_F77(ca, cb, cc, &nup);
  }
  else {
    for (ordinal_type i=0; i<pa; i++)
      cc[i] = ca[i]/cb[0];
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
divide(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
       const value_type& a, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& b)
{
  ordinal_type pc = b.size();

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(c.size() != pc, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::divide()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "b.size() = " << pc << ", c.size() = " << c.size()
		     << ".");
#endif

  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  if (pc > 1) {
#ifdef STOKHOS_DEBUG
    const char* func = "Stokhos::ForUQTKOrthogPolyExpansion::divide()";
    TEST_FOR_EXCEPTION(size() < pc, std::logic_error,
		       func << ":  Expansion size (" << size()
		       << ") is too small for computation (" << pc
		       << " needed).");
#endif
   
    value_type* ca = Stokhos::ds_array<value_type>::get_and_fill(pc);
    ca[0] = a;
    int nup = pc-1;
    UQ_DIV_F77(ca, cb, cc, &nup);
  }
  else 
    cc[0] = a / cb[0];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
divide(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a, 
       const value_type& b)
{
  ordinal_type pc = a.size();

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(c.size() != pc, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::divide()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "a.size() = " << pc << ", c.size() = " << c.size()
		     << ".");
#endif

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();

  for (ordinal_type i=0; i<pc; i++)
    cc[i] = ca[i]/b;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
exp(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a)
{
  ordinal_type pc = a.size();

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(c.size() != pc, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::exp()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "a.size() = " << pc << ", c.size() = " << c.size()
		     << ".");
#endif

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();

  int nup = pc-1;
  if (method == TAYLOR) {
    int nrm = 1;
    UQ_EXP_F77(ca, cc, &dim, &nup, &rtol, &nrm);
  }
  else
    UQ_EXP_INT_F77(ca, cc, &nup);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
log(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a)
{
  ordinal_type pc = a.size();

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(c.size() != pc, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::log()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "a.size() = " << pc << ", c.size() = " << c.size()
		     << ".");
#endif

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();

  int nup = pc-1;
  if (method == TAYLOR) {
    int nrm = 1;
    UQ_LOG_F77(ca, cc, &dim, &nup, &rtol, &nrm);
  }
  else
    UQ_LOG_INT_F77(ca, cc, &nup);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
log10(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a)
{
  log(c,a);
  divide(c,c,std::log(10.0));
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
sqrt(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a)
{
  ordinal_type pc = a.size();

#ifdef STOKHOS_DEBUG
  TEST_FOR_EXCEPTION(c.size() != pc, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::sqrt()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "a.size() = " << pc << ", c.size() = " << c.size()
		     << ".");
#endif

  const value_type* ca = a.coeff();
  value_type* cc = c.coeff();
  int iguess = 0;

  int nup = pc-1;
  UQ_SQRT_F77(ca, cc, &nup, &iguess);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
pow(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a,
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& b)
{
  log(c,a);
  timesEqual(c,b);
  exp(c,c);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
pow(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
    const value_type& a, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& b)
{
  times(c,std::log(a),b);
  exp(c,c);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
pow(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a, 
    const value_type& b)
{
  log(c,a);
  timesEqual(c,b);
  exp(c,c);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
sin(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& s, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a)
{
  throw "Not Defined";
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
cos(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a)
{
  throw "Not Defined";
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
tan(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& t, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a)
{
  throw "Not Defined";
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
sinh(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& s, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a)
{
  // sinh(x) = (exp(x) - exp(-x))/2.0
  Stokhos::OrthogPolyApprox<ordinal_type, value_type> t(a);
  timesEqual(t, -1.0);
  exp(s, a);
  exp(t, t);
  minusEqual(s, t);
  divideEqual(s, 2.0);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
cosh(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a)
{
  // cosh(x) = (exp(x) + exp(-x))/2.0
  Stokhos::OrthogPolyApprox<ordinal_type, value_type> t(a);
  timesEqual(t, -1.0);
  exp(c, a);
  exp(t, t);
  plusEqual(c, t);
  divideEqual(c, 2.0);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
tanh(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& t, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a)
{
  // tanh(x) = (exp(x) - exp(-x))/(exp(x) + exp(-x))
  Stokhos::OrthogPolyApprox<ordinal_type, value_type> s(a);
  Stokhos::OrthogPolyApprox<ordinal_type, value_type> c(a);
  timesEqual(s, -1.0);
  exp(s, s);
  exp(c, a);
  minus(t, c, s);
  plusEqual(c, s);
  divideEqual(t, c);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
acos(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a)
{
  throw "Not Defined";
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
asin(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a)
{
  throw "Not Defined";
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
atan(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a)
{
  throw "Not Defined";
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
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
acosh(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a)
{
  throw "Not Defined";
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
asinh(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a)
{
  throw "Not Defined";
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
atanh(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a)
{
  throw "Not Defined";
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
fabs(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a)
{
  if (a[0] >= 0)
    c = a;
  else
    unaryMinus(c,a);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
abs(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a)
{
  if (a[0] >= 0)
    c = a;
  else
    unaryMinus(c,a);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
max(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a,
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& b)
{
  if (a[0] >= b[0])
    c = a;
  else
    c = b;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
max(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
    const value_type& a, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& b)
{
  if (a >= b[0]) {
    c = OrthogPolyApprox<ordinal_type, value_type>(basis);
    c[0] = a;
  }
  else
    c = b;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
max(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a, 
    const value_type& b)
{
  if (a[0] >= b)
    c = a;
  else {
    c = OrthogPolyApprox<ordinal_type, value_type>(basis);
    c[0] = b;
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
min(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a,
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& b)
{
  if (a[0] <= b[0])
    c = a;
  else
    c = b;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
min(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
    const value_type& a, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& b)
{
  if (a <= b[0]) {
    c = OrthogPolyApprox<ordinal_type, value_type>(basis);
    c[0] = a;
  }
  else
    c = b;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
min(Stokhos::OrthogPolyApprox<ordinal_type, value_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& a, 
    const value_type& b)
{
  if (a[0] <= b)
    c = a;
  else {
    c = OrthogPolyApprox<ordinal_type, value_type>(basis);
    c[0] = b;
  }
}
