// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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

#include "Teuchos_Assert.hpp"
#include "Stokhos_DynamicArrayTraits.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_Tuple.hpp"

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type> 
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
PseudoSpectralOrthogPolyExpansion(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis_,
  const Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk_,
  const Teuchos::RCP<const PseudoSpectralOperator<ordinal_type, value_type, point_compare_type> >& ps_op_,
  const Teuchos::RCP<Teuchos::ParameterList>& params_) :
  OrthogPolyExpansionBase<ordinal_type, value_type, node_type>(basis_, Cijk_, params_),
  ps_op(ps_op_),
  sz(ps_op->coeff_size()),
  nqp(ps_op->point_size()),
  avals(nqp),
  bvals(nqp),
  fvals(nqp)
{
  Teuchos::RCP<Teuchos::ParameterList> params = params_;
  if (params == Teuchos::null)
    params = Teuchos::rcp(new Teuchos::ParameterList);
  use_quad_for_times = params->get("Use Quadrature for Times", false);
  use_quad_for_division = params->get("Use Quadrature for Division", true);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type> 
template <typename FuncT>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
unary_op(const FuncT& func,
         OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
         const OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  ordinal_type pa = a.size();
  ordinal_type pc;
  if (a.size() == 1)
    pc = 1;
  else
    pc = sz;
  if (c.size() != pc)
    c.resize(pc);

  if (pc == 1) {
    c[0] = func(a[0]);
    return;
  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::PSExp -- Unary Polynomial Evaluation");
#endif

  // Evaluate input
  SDV a_sdv(Teuchos::View, const_cast<value_type*>(a.coeff()), pa);
  ps_op->transformPCE2QP(1.0, a_sdv, avals, 0.0, false);

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::PSExp -- Unary Function Evaluation");
#endif

    // Evaluate function
  for (ordinal_type qp=0; qp<nqp; qp++) {
    fvals[qp] = func(avals[qp]);
  }

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::PSExp -- Unary Polynomial Integration");
#endif

  // Integrate
  SDV c_sdv(Teuchos::View, c.coeff(), pc);
  ps_op->transformQP2PCE(1.0, fvals, c_sdv, 0.0, false);

  }
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type> 
template <typename FuncT>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
binary_op(const FuncT& func,
          OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
          const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
          const OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  ordinal_type pa = a.size();
  ordinal_type pb = b.size();
  ordinal_type pc;
  if (pa == 1 && pb == 1)
    pc = 1;
  else
    pc = sz;
  if (c.size() != pc)
    c.resize(pc);

  if (pc == 1) {
    c[0] = func(a[0], b[0]);
    return;
  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::PSExp -- PP Binary Polynomial Evaluation");
#endif

  // Evaluate input
  SDV a_sdv(Teuchos::View, const_cast<value_type*>(a.coeff()), pa);
  SDV b_sdv(Teuchos::View, const_cast<value_type*>(b.coeff()), pb);
  ps_op->transformPCE2QP(1.0, a_sdv, avals, 0.0, false);
  ps_op->transformPCE2QP(1.0, b_sdv, bvals, 0.0, false);

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::PSExp -- PP Binary Function Evaluation");
#endif

    // Evaluate function
  for (ordinal_type qp=0; qp<nqp; qp++)
    fvals[qp] = func(avals[qp], bvals[qp]);

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::PSExp -- PP Binary Polynomial Integration");
#endif

  // Integrate
  SDV c_sdv(Teuchos::View, c.coeff(), pc);
  ps_op->transformQP2PCE(1.0, fvals, c_sdv, 0.0, false);

  }
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type> 
template <typename FuncT>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
binary_op(const FuncT& func,
          OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
          const value_type& a, 
          const OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  ordinal_type pb = b.size();
  ordinal_type pc;
  if (pb == 1)
    pc = 1;
  else
    pc = sz;
  if (c.size() != pc)
    c.resize(pc);

  if (pc == 1) {
    c[0] = func(a, b[0]);
    return;
  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::PSExp -- CP Binary Polynomial Evaluation");
#endif

  // Evaluate input
  SDV b_sdv(Teuchos::View, const_cast<value_type*>(b.coeff()), pb);
  ps_op->transformPCE2QP(1.0, b_sdv, bvals, 0.0, false);

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::PSExp -- CP Binary Function Evaluation");
#endif

    // Evaluate function
  for (ordinal_type qp=0; qp<nqp; qp++)
    fvals[qp] = func(a, bvals[qp]);

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::PSExp -- CP Binary Polynomial Integration");
#endif

  // Integrate
  SDV c_sdv(Teuchos::View, c.coeff(), pc);
  ps_op->transformQP2PCE(1.0, fvals, c_sdv, 0.0, false);

  }
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type> 
template <typename FuncT>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
binary_op(const FuncT& func,
          OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
          const OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
          const value_type& b)
{
  ordinal_type pa = a.size();
  ordinal_type pc;
  if (pa == 1)
    pc = 1;
  else
    pc = sz;
  if (c.size() != pc)
    c.resize(pc);

  if (pc == 1) {
    c[0] = func(a[0], b);
    return;
  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::PSExp -- PC Binary Polynomial Evaluation");
#endif

  // Evaluate input
  SDV a_sdv(Teuchos::View, const_cast<value_type*>(a.coeff()), pa);
  ps_op->transformPCE2QP(1.0, a_sdv, avals, 0.0, false);

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::PSExp -- PC Binary Function Evaluation");
#endif

    // Evaluate function
  for (ordinal_type qp=0; qp<nqp; qp++)
    fvals[qp] = func(avals[qp], b);

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::PSExp -- PC Binary Polynomial Integration");
#endif

  // Integrate
  SDV c_sdv(Teuchos::View, c.coeff(), pc);
  ps_op->transformQP2PCE(1.0, fvals, c_sdv, 0.0, false);

  }
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type> 
template <typename FuncT>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
nary_op(const FuncT& func,
	OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	const OrthogPolyApprox<ordinal_type, value_type, node_type>** na)
{
  const int N = FuncT::N;
  bool is_constant = true;
  for (int i=0; i<N; i++) {
    if (na[i]->size() > 1) {
      is_constant = false;
      break;
    }
  }
  ordinal_type pc;
  if (is_constant)
    pc = 1;
  else
    pc = sz;
  if (c.size() != pc)
    c.resize(pc);

  if (pc == 1) {
    value_type val[N];
    for (int i=0; i<N; i++)
      val[i] = (*na[i])[0];
    c[0] = func(val);
    return;
  }

  if (N >= navals.size())
    navals.resize(N+1);
  if (navals[N].size() != N) {
    navals[N].resize(N);
    for (int i=0; i<N; i++)
      navals[N][i].resize(nqp);
  }
  
  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::PSExp -- N(" << N << ")-ary Polynomial Evaluation");
#endif

  // Evaluate input
  for (int i=0; i<N; i++) {
    SDV sdv(Teuchos::View, const_cast<value_type*>(na[i]->coeff()), 
	    na[i]->size());
    ps_op->transformPCE2QP(1.0, sdv, navals[N][i], 0.0, false);
  }

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::PSExp -- N(" << N << ")-ary Function Evaluation");
#endif

  // Evaluate function
  value_type val[N];
  for (ordinal_type qp=0; qp<nqp; qp++) {
    for (int i=0; i<N; i++)
      val[i] = navals[N][i][qp];
    fvals[qp] = func(val);
  }

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::PSExp -- N(" << N << ")-ary Polynomial Integration");
#endif

  // Integrate
  SDV c_sdv(Teuchos::View, c.coeff(), pc);
  ps_op->transformQP2PCE(1.0, fvals, c_sdv, 0.0, false);

  }
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type> 
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
timesEqual(
  Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
  const value_type& x)
{
  OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::timesEqual(c,x);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type> 
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
divideEqual(
  Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
  const value_type& x)
{
  OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::divideEqual(c,x);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type> 
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
timesEqual(
  Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
  const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& x)
{
  if (use_quad_for_times)
    binary_op(times_quad_func(), c, c, x);
  else
    OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::timesEqual(c, x);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type> 
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
divideEqual(
  Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
  const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& x)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::divideEqual(OPA)");
#endif
  if (x.size() == 1) {
    ordinal_type p = c.size();
    value_type* cc = c.coeff();
    const value_type* xc = x.coeff();
    for (ordinal_type i=0; i<p; i++)
      cc[i] /= xc[0];
  }
  else {
    if (use_quad_for_division)
      binary_op(div_quad_func(), c, c, x);
    else
      OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::divideEqual(c, x);
  }
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
times(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (use_quad_for_times)
    binary_op(times_quad_func(), c, a, b);
  else
    OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::times(c, a, b);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
times(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const value_type& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::times(c, a, b);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
times(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
      const value_type& b)
{
  OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::times(c, a, b);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
divide(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::OrthogPolyExpansionBase::divide(OPA,OPA)");
#endif
  if (b.size() == 1) {
    ordinal_type pc = a.size();
    if (c.size() != pc)
      c.resize(pc);

    const value_type* ca = a.coeff();
    const value_type* cb = b.coeff();
    value_type* cc = c.coeff();

    for (ordinal_type i=0; i<pc; i++)
      cc[i] = ca[i]/cb[0];
  }
  else {
    if (use_quad_for_division)
      binary_op(div_quad_func(), c, a, b);
    else
      OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::divide(c, a, b);
  }
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
divide(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
       const value_type& a, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (use_quad_for_division)
    binary_op(div_quad_func(), c, a, b);
  else
    OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::divide(c, a, b);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
divide(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
      const value_type& b)
{
  OrthogPolyExpansionBase<ordinal_type, value_type, node_type>::divide(c, a, b);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
exp(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(exp_quad_func(), c, a);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
log(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(log_quad_func(), c, a);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
log10(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(log10_quad_func(), c, a);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
sqrt(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(sqrt_quad_func(), c, a);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
pow(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  binary_op(pow_quad_func(), c, a, b);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
pow(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const value_type& a, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  binary_op(pow_quad_func(), c, a, b);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
pow(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
    const value_type& b)
{
  binary_op(pow_quad_func(), c, a, b);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
sin(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& s, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(sin_quad_func(), s, a);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
cos(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(cos_quad_func(), c, a);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
tan(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& t, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(tan_quad_func(), t, a);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
sinh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& s, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(sinh_quad_func(), s, a);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
cosh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(cosh_quad_func(), c, a);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
tanh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& t, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(tanh_quad_func(), t, a);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
acos(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(acos_quad_func(), c, a);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
asin(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(asin_quad_func(), c, a);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
atan(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(atan_quad_func(), c, a);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
atan2(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  binary_op(atan2_quad_func(), c, a, b);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
atan2(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const value_type& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  binary_op(atan2_quad_func(), c, a, b);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
atan2(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const value_type& b)
{
  binary_op(atan2_quad_func(), c, a, b);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
acosh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(acosh_quad_func(), c, a);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
asinh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(asinh_quad_func(), c, a);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
void
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
atanh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(atanh_quad_func(), c, a);
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
template <typename ExprT1, typename ExprT2>
value_type
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
compute_times_coeff(ordinal_type i, const ExprT1& a, const ExprT2& b) const
{
  ordinal_type pa = a.size();
  ordinal_type pb = b.size();

  if (pa > 1 && pb > 1) {
    ordinal_type k_lim = pa;
    ordinal_type j_lim = pb;
    if (pb < pa) {
      k_lim = pb;
      j_lim = pa;
    }
    typename Cijk_type::i_iterator i_it = this->Cijk->find_i(i);
#ifdef STOKHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(i_it == this->Cijk->i_end(), std::logic_error,
		     "Stokhos::PseudoSpectralOrthogPolyExpansion::compute_times_coeff()" 
		       << ":  Index " << i << " is out of range [0," 
		       << this->Cijk->num_i() << ")!");
#endif
    value_type cc = value_type(0);
    value_type aa, bb, cijk;
    ordinal_type j, k;
    for (typename Cijk_type::ik_iterator k_it = this->Cijk->k_begin(i_it); 
	 k_it != this->Cijk->k_end(i_it); ++k_it) {
      k = index(k_it);
      if (k < k_lim) {
	if (pa < pb) {
	  if (k == 0)
	    aa = a.val();
	  else
	    aa = a.higher_order_coeff(k);
	}
	else {
	  if (k == 0)
	    aa = b.val();
	  else
	    aa = b.higher_order_coeff(k);
	}
	for (typename Cijk_type::ikj_iterator j_it = this->Cijk->j_begin(k_it);
	     j_it != this->Cijk->j_end(k_it); ++j_it) {
	  j = index(j_it);
	  cijk = value(j_it);
	  if (j < j_lim) {
	    if (pa < pb) {
	      if (j == 0)
		bb = b.val();
	      else
		bb = b.higher_order_coeff(j);
	    }
	    else {
	      if (j == 0)
		bb = a.val();
	      else
		bb = a.higher_order_coeff(j);
	    }
	    cc += cijk*aa*bb;
	  }
	}
      }
    }
    return cc / this->basis->norm_squared(i);
  }
  else if (i == 0)
    return a.val() * b.val();
  else if (pa > 1) {
    return a.higher_order_coeff(i)*b.val();
  }
  else {
    return a.val()*b.higher_order_coeff(i);
  }
}

template <typename ordinal_type, typename value_type, typename point_compare_type, typename node_type>
template <typename ExprT1, typename ExprT2>
value_type
Stokhos::PseudoSpectralOrthogPolyExpansion<ordinal_type, value_type, point_compare_type, node_type>::
fast_compute_times_coeff(ordinal_type i, const ExprT1& a, const ExprT2& b) const
{
  typename Cijk_type::i_iterator i_it = this->Cijk->find_i(i);
#ifdef STOKHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(i_it == this->Cijk->i_end(), std::logic_error,
		     "Stokhos::PseudoSpectralOrthogPolyExpansion::fast_ompute_times_coeff()" 
		     << ":  Index " << i << " is out of range [0," 
		     << this->Cijk->num_i() << ")!");
#endif
  value_type cc = value_type(0);
  value_type aa, bb, cijk;
  ordinal_type j, k;
  for (typename Cijk_type::ik_iterator k_it = this->Cijk->k_begin(i_it); 
       k_it != this->Cijk->k_end(i_it); ++k_it) {
    k = index(k_it);
    if (k == 0)
      aa = a.val();
    else
      aa = a.higher_order_coeff(k);
    for (typename Cijk_type::ikj_iterator j_it = this->Cijk->j_begin(k_it);
	 j_it != this->Cijk->j_end(k_it); ++j_it) {
      j = index(j_it);
      cijk = value(j_it);
      if (j == 0)
	bb = b.val();
      else
	bb = b.higher_order_coeff(j);
      cc += cijk*aa*bb;
    }
  }

  return cc / this->basis->norm_squared(i);
}
