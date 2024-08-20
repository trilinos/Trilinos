// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

template <typename ordinal_type, typename value_type> 
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
ConstantOrthogPolyExpansion()
{
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
unaryMinus(
  Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
  const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = -a[0];
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
plusEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	  const value_type& val)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] += val;
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
minusEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	   const value_type& val)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] -= val;
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
timesEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	   const value_type& val)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] *= val;
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
divideEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	    const value_type& val)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] /= val;
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
plusEqual(
  Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
  const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& x)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] += x[0];
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
minusEqual(
  Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
  const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& x)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] -= x[0];
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
timesEqual(
  Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
  const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& x)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] *= x[0];
}

template <typename ordinal_type, typename value_type> 
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
divideEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
            const Stokhos::OrthogPolyApprox<ordinal_type, value_type >& x)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] /= x[0];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
plus(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = a[0] + b[0];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
plus(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const value_type& a, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = a + b[0];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
plus(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
     const value_type& b)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = a[0] + b;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
minus(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = a[0] - b[0];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
minus(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const value_type& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = a - b[0];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
minus(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const value_type& b)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = a[0] - b;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
times(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = a[0] * b[0];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
times(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const value_type& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = a * b[0];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
times(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const value_type& b)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = a[0] * b;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
divide(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = a[0] / b[0];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
divide(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
       const value_type& a, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = a / b[0];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
divide(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
       const value_type& b)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = a[0] / b;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
exp(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::exp(a[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
log(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::log(a[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
log10(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::log10(a[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
sqrt(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::sqrt(a[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
cbrt(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::cbrt(a[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
pow(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::pow(a[0], b[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
pow(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const value_type& a, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::pow(a, b[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
pow(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
    const value_type& b)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::pow(a[0], b);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
sin(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& s, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (s.size() < 1)
    s.resize(1);
  s[0] = std::sin(a[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
cos(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::cos(a[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
tan(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& t, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (t.size() < 1)
    t.resize(1);
  t[0] = std::tan(a[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
sinh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& s, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (s.size() < 1)
    s.resize(1);
  s[0] = std::sinh(a[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
cosh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::cosh(a[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
tanh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& t, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (t.size() < 1)
    t.resize(1);
  t[0] = std::tanh(a[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
acos(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::acos(a[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
asin(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::asin(a[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
atan(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::atan(a[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
atan2(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::atan2(a[0], b[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
atan2(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const value_type& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::atan2(a, b[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
atan2(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const value_type& b)
{
  c[0] = std::atan2(a[0], b);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
acosh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::log(a[0]+std::sqrt(a[0]*a[0]-value_type(1.0))); 
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
asinh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::log(a[0]+std::sqrt(a[0]*a[0]+value_type(1.0))); 
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
atanh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = 0.5*std::log((value_type(1.0)+a[0])/(value_type(1.0)-a[0])); 
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
fabs(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::fabs(a[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
abs(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::fabs(a[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
max(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::max(a[0], b[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
max(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const value_type& a, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::max(a, b[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
max(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
    const value_type& b)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::max(a[0], b);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
min(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::min(a[0], b[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
min(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const value_type& a, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::min(a, b[0]);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ConstantOrthogPolyExpansion<ordinal_type, value_type>::
min(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
    const value_type& b)
{
  if (c.size() < 1)
    c.resize(1);
  c[0] = std::min(a[0], b);
}
