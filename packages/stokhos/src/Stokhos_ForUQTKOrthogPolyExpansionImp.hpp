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
  const Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk_,
  EXPANSION_METHOD method_,
  value_type rtol_) :
  OrthogPolyExpansionBase<ordinal_type, value_type, node_type>(basis_, Cijk_),
  rtol(rtol_),
  method(method_)
{
  order = this->basis->order();
  dim = this->basis->dimension();
  int nup;
  UQ_PREP_F77(&order, &dim, &nup);
  sz = nup+1;
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
timesEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	   const value_type& val)
{
  OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::timesEqual(c,val);
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
divideEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	    const value_type& val)
{
  OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::divideEqual(c,val);
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
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
		     "Stokhos::ForUQTKOrthogPolyExpansion::timesEqual()" <<
		     ":  Expansion size (" << sz << 
		     ") is too small for computation.");
  if (c.size() != pc)
    c.resize(pc);

  value_type* cc = c.coeff();
  const value_type* xc = x.coeff();
  
  if (p > 1 && xp > 1) {
    TEUCHOS_TEST_FOR_EXCEPTION(pc != xp, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::timesEqual()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "x.size() = " << xp << ", c.size() = " << pc << ".");

    // Copy c coefficients into temporary array
    value_type* tc = Stokhos::ds_array<value_type>::get_and_fill(cc,pc);
   
    int nup = pc-1;
    UQ_PROD2_F77(cc, xc, tc, &nup);
    Stokhos::ds_array<value_type>::copy(tc, cc, p);
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
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
divideEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	    const OrthogPolyApprox<ordinal_type, value_type, node_type>& x)
{
  ordinal_type p = c.size();
  ordinal_type xp = x.size();
  ordinal_type pc;
  if (xp > 1)
    pc = sz;
  else
    pc = p;
  TEUCHOS_TEST_FOR_EXCEPTION(sz < pc, std::logic_error,
		     "Stokhos::ForUQTKOrthogPolyExpansion::divideEqual()" <<
		     ":  Expansion size (" << sz << 
		     ") is too small for computation.");
  if (c.size() != pc)
    c.resize(pc);

  value_type* cc = c.coeff();
  const value_type* xc = x.coeff();
  
  if (xp > 1) {
    TEUCHOS_TEST_FOR_EXCEPTION(pc != xp, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::divideEqual()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "x.size() = " << xp << ", c.size() = " << pc << ".");
    
    // Copy c coefficients into temporary array
    value_type* tc = Stokhos::ds_array<value_type>::get_and_fill(cc,pc);
   
    int nup = pc-1;
    UQ_DIV_F77(tc, xc, cc, &nup);
    
  }
  else {
    for (ordinal_type i=0; i<pc; i++)
      cc[i] /= xc[0];
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
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
		     "Stokhos::ForUQTKOrthogPolyExpansion::times()" <<
		     ":  Expansion size (" << sz << 
		     ") is too small for computation.");
  if (c.size() != pc)
    c.resize(pc);

  const value_type* ca = a.coeff();
  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  if (pa > 1 && pb > 1) {
    TEUCHOS_TEST_FOR_EXCEPTION(pa != pc || pb != pc, std::logic_error,
		       "Stokhos::ForUQTKOrthogPolyExpansion::times()" 
		       << ":  Arguments have incompatible sizes:  "
		       << "a.size() = " << pa << ", b.size() = " << pb 
		       << ", required size = " << pc << ".");

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
times(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const value_type& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::times(c,a,b);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
times(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const value_type& b)
{
  OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::times(c,a,b);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
divide(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  ordinal_type pa = a.size();
  ordinal_type pb = b.size();
  ordinal_type pc;
  if (pb > 1)
    pc = sz;
  else
    pc = pa;
  TEUCHOS_TEST_FOR_EXCEPTION(sz < pc, std::logic_error,
		     "Stokhos::ForUQTKOrthogPolyExpansion::divide()" <<
		     ":  Expansion size (" << sz << 
		     ") is too small for computation.");
  if (c.size() != pc)
    c.resize(pc);

  const value_type* ca = a.coeff();
  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  if (pb > 1) {
    TEUCHOS_TEST_FOR_EXCEPTION(pa != pc || pb != pc, std::logic_error,
		       "Stokhos::ForUQTKOrthogPolyExpansion::divide()" 
		       << ":  Arguments have incompatible sizes:  "
		       << "a.size() = " << pa << ", b.size() = " << pb 
		       << ", required size = " << pc << ".");

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
divide(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
       const value_type& a, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
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
    TEUCHOS_TEST_FOR_EXCEPTION(pb != pc, std::logic_error,
		       "Stokhos::ForUQTKOrthogPolyExpansion::divide()" 
		       << ":  Arguments have incompatible sizes:  "
		       << "b.size() = " << pb 
		       << ", required size = " << pc << ".");
   
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
divide(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
       const value_type& b)
{
  OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::divide(c,a,b);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
exp(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
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
    TEUCHOS_TEST_FOR_EXCEPTION(pa != pc, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::exp()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "a.size() = " << pa << ", c.size() = " << pc
		     << ".");
    int nup = pc-1;
    if (method == TAYLOR) {
      int nrm = 1;
      UQ_EXP_F77(ca, cc, &dim, &nup, &rtol, &nrm);
    }
    else
      UQ_EXP_INT_F77(ca, cc, &nup);
  }
  else
    cc[0] = std::exp(ca[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
log(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
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
    TEUCHOS_TEST_FOR_EXCEPTION(pa != pc, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::log()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "a.size() = " << pa << ", c.size() = " << pc
		     << ".");
    int nup = pc-1;
    if (method == TAYLOR) {
      int nrm = 1;
      UQ_LOG_F77(ca, cc, &dim, &nup, &rtol, &nrm);
    }
    else
      UQ_LOG_INT_F77(ca, cc, &nup);
  }
  else
    cc[0] = std::log(ca[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
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
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
sqrt(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
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
    TEUCHOS_TEST_FOR_EXCEPTION(pa != pc, std::logic_error,
                     "Stokhos::ForUQTKOrthogPolyExpansion::sqrt()" 
                     << ":  Arguments have incompatible sizes:  "
		     << "a.size() = " << pa << ", c.size() = " << pc
		     << ".");
    int iguess = 0;
    int nup = pc-1;
    UQ_SQRT_F77(ca, cc, &nup, &iguess);
  }
  else
    cc[0] = std::sqrt(ca[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
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
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
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
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
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
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
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
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
sin(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& s, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() == 1) {
    if (s.size() != 1)
      s.resize(1);
    s[0] = std::sin(a[0]);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Stokhos::ForUQTKOrthogPolyExpansion::sin()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
cos(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() == 1) {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::cos(a[0]);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Stokhos::ForUQTKOrthogPolyExpansion::cos()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
tan(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& t, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() == 1) {
    if (t.size() != 1)
      t.resize(1);
    t[0] = std::tan(a[0]);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Stokhos::ForUQTKOrthogPolyExpansion::tan()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
sinh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& s, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() > 1) {
    // sinh(x) = (exp(x) - exp(-x))/2.0
    Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type> t(a);
    timesEqual(t, -1.0);
    exp(s, a);
    exp(t, t);
    this->minusEqual(s, t);
    divideEqual(s, 2.0);
  }
  else {
    if (s.size() != 1)
      s.resize(1);
    s[0] = std::sinh(a[0]);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
cosh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() > 1) {
    // cosh(x) = (exp(x) + exp(-x))/2.0
    Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type> t(a);
    timesEqual(t, -1.0);
    exp(c, a);
    exp(t, t);
    this->plusEqual(c, t);
    divideEqual(c, 2.0);
  }
  else {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::cosh(a[0]);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
tanh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& t, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() > 1) {
    // tanh(x) = (exp(x) - exp(-x))/(exp(x) + exp(-x))
    Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type> s(a);
    Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type> c(a);
    timesEqual(s, -1.0);
    exp(s, s);
    exp(c, a);
    this->minus(t, c, s);
    this->plusEqual(c, s);
    divideEqual(t, c);
  }
  else {
    if (t.size() != 1)
      t.resize(1);
    t[0] = std::tanh(a[0]);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
acos(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() == 1) {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::acos(a[0]);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Stokhos::ForUQTKOrthogPolyExpansion::acos()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
asin(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() == 1) {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::asin(a[0]);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Stokhos::ForUQTKOrthogPolyExpansion::asin()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
atan(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() == 1) {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::atan(a[0]);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Stokhos::ForUQTKOrthogPolyExpansion::atan()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
atan2(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (a.size() == 1 && b.size() == 1) {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::atan2(a[0], b[0]);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Stokhos::ForUQTKOrthogPolyExpansion::atan2()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
atan2(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const value_type& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (b.size() == 1) {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::atan2(a, b[0]);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Stokhos::ForUQTKOrthogPolyExpansion::atan2()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
atan2(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const value_type& b)
{
  if (a.size() == 1) {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::atan2(a[0], b);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Stokhos::ForUQTKOrthogPolyExpansion::atan2()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
acosh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() == 1) {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::log(a[0]+std::sqrt(a[0]*a[0]-value_type(1.0)));
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Stokhos::ForUQTKOrthogPolyExpansion::acosh()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
asinh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
 if (a.size() == 1) {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::log(a[0]+std::sqrt(a[0]*a[0]+value_type(1.0)));
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Stokhos::ForUQTKOrthogPolyExpansion::asinh()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ForUQTKOrthogPolyExpansion<ordinal_type, value_type>::
atanh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() == 1) {
    if (c.size() != 1)
      c.resize(1);
    c[0] = 0.5*std::log((value_type(1.0)+a[0])/(value_type(1.0)-a[0]));
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Stokhos::ForUQTKOrthogPolyExpansion::atanh()" 
		       << ":  Method not implemented!");
}
