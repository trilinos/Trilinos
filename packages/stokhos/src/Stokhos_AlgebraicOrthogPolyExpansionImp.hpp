// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Assert.hpp"

template <typename ordinal_type, typename value_type> 
Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type>::
AlgebraicOrthogPolyExpansion(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis_,
  const Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk_,
  const Teuchos::RCP<Teuchos::ParameterList>& params_) :
  OrthogPolyExpansionBase<ordinal_type, value_type, node_type>(basis_, Cijk_, params_)
{
}

template <typename ordinal_type, typename value_type>
void
Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type>::
exp(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() == 1) {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::exp(a[0]);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Stokhos::AlgebraicOrthogPolyExpansion::exp()" 
		       << ":  Method not implemented!");

}

template <typename ordinal_type, typename value_type>
void
Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type>::
log(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() == 1) {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::log(a[0]);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Stokhos::AlgebraicOrthogPolyExpansion::log()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type>::
log10(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() == 1) {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::log10(a[0]);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Stokhos::AlgebraicOrthogPolyExpansion::log10()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type>::
sqrt(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() == 1) {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::sqrt(a[0]);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Stokhos::AlgebraicOrthogPolyExpansion::sqrt()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type>::
cbrt(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() == 1) {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::cbrt(a[0]);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Stokhos::AlgebraicOrthogPolyExpansion::cbrt()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type>::
pow(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (a.size() == 1 && b.size() == 1) {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::pow(a[0], b[0]);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Stokhos::AlgebraicOrthogPolyExpansion::pow()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type>::
pow(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const value_type& a, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (b.size() == 1) {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::pow(a, b[0]);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Stokhos::AlgebraicOrthogPolyExpansion::pow()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type>::
pow(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
    const value_type& b)
{
  if (a.size() == 1) {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::pow(a[0], b);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Stokhos::AlgebraicOrthogPolyExpansion::pow()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type>::
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
		       "Stokhos::AlgebraicOrthogPolyExpansion::sin()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type>::
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
		       "Stokhos::AlgebraicOrthogPolyExpansion::cos()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type>::
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
		       "Stokhos::AlgebraicOrthogPolyExpansion::tan()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type>::
sinh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& s, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() == 1) {
    if (s.size() != 1)
      s.resize(1);
    s[0] = std::sinh(a[0]);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Stokhos::AlgebraicOrthogPolyExpansion::sinh()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type>::
cosh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() == 1) {
    if (c.size() != 1)
      c.resize(1);
    c[0] = std::cosh(a[0]);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Stokhos::AlgebraicOrthogPolyExpansion::cosh()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type>::
tanh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& t, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a.size() == 1) {
    if (t.size() != 1)
      t.resize(1);
    t[0] = std::tanh(a[0]);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Stokhos::AlgebraicOrthogPolyExpansion::tanh()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type>::
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
		       "Stokhos::AlgebraicOrthogPolyExpansion::acos()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type>::
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
		       "Stokhos::AlgebraicOrthogPolyExpansion::asin()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type>::
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
		       "Stokhos::AlgebraicOrthogPolyExpansion::atan()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type>::
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
		       "Stokhos::AlgebraicOrthogPolyExpansion::atan2()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type>::
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
		       "Stokhos::AlgebraicOrthogPolyExpansion::atan2()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type>::
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
		       "Stokhos::AlgebraicOrthogPolyExpansion::atan2()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type>::
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
		       "Stokhos::AlgebraicOrthogPolyExpansion::acosh()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type>::
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
		       "Stokhos::AlgebraicOrthogPolyExpansion::asinh()" 
		       << ":  Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type>::
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
		       "Stokhos::AlgebraicOrthogPolyExpansion::atanh()" 
		       << ":  Method not implemented!");
}
