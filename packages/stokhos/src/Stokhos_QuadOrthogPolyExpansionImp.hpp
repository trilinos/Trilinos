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

#include "Teuchos_TestForException.hpp"
#include "Stokhos_DynamicArrayTraits.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_Tuple.hpp"

template <typename ordinal_type, typename value_type, typename node_type> 
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
QuadOrthogPolyExpansion(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis_,
  const Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk_,
  const Teuchos::RCP<const Quadrature<ordinal_type, value_type> >& quad_,
  bool use_quad_for_times_) :
  basis(basis_),
  Cijk(Cijk_),
  quad(quad_),
  use_quad_for_times(use_quad_for_times_),
  sz(basis->size()),
  blas(),
  quad_points(quad->getQuadPoints()),
  quad_weights(quad->getQuadWeights()),
  quad_values(quad->getBasisAtQuadPoints()),
  norms(basis->norm_squared()),
  nqp(quad_points.size()),
  avals(nqp),
  bvals(nqp),
  fvals(nqp),
  qv(nqp*sz),
  sqv(nqp*sz)
{
  for (ordinal_type qp=0; qp<nqp; qp++)
    for (ordinal_type i=0; i<sz; i++) {
      qv[qp*sz+i] = quad_values[qp][i];
      sqv[qp*sz+i] = quad_values[qp][i]/norms[i];
    }
}

template <typename ordinal_type, typename value_type, typename node_type> 
template <typename FuncT>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
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
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::QuadExp -- Unary Polynomial Evaluation");
#endif

  // Evaluate input
  blas.GEMV(Teuchos::TRANS, pa, nqp, 1.0, &qv[0], sz, a.coeff(), 1, 0.0,
            &avals[0], 1);

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::QuadExp -- Unary Function Evaluation");
#endif

    // Evaluate function
  for (ordinal_type qp=0; qp<nqp; qp++)
    fvals[qp] = func(avals[qp])*quad_weights[qp];

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::QuadExp -- Unary Polynomial Integration");
#endif

  // Integrate
  blas.GEMV(Teuchos::NO_TRANS, pc, nqp, 1.0, &sqv[0], sz, &fvals[0], 1, 0.0,
            c.coeff(), 1);

  }
}

template <typename ordinal_type, typename value_type, typename node_type> 
template <typename FuncT>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
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
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::QuadExp -- PP Binary Polynomial Evaluation");
#endif

  // Evaluate input
  blas.GEMV(Teuchos::TRANS, pa, nqp, 1.0, &qv[0], sz, a.coeff(), 1, 0.0,
            &avals[0], 1);
  blas.GEMV(Teuchos::TRANS, pb, nqp, 1.0, &qv[0], sz, b.coeff(), 1, 0.0,
            &bvals[0], 1);

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::QuadExp -- PP Binary Function Evaluation");
#endif

    // Evaluate function
  for (ordinal_type qp=0; qp<nqp; qp++)
    fvals[qp] = func(avals[qp], bvals[qp])*quad_weights[qp];

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::QuadExp -- PP Binary Polynomial Integration");
#endif

  // Integrate
  blas.GEMV(Teuchos::NO_TRANS, pc, nqp, 1.0, &sqv[0], sz, &fvals[0], 1, 0.0,
            c.coeff(), 1);

  }
}

template <typename ordinal_type, typename value_type, typename node_type> 
template <typename FuncT>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
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
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::QuadExp -- CP Binary Polynomial Evaluation");
#endif

  // Evaluate input
  blas.GEMV(Teuchos::TRANS, pb, nqp, 1.0, &qv[0], sz, b.coeff(), 1, 0.0,
            &bvals[0], 1);

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::QuadExp -- CP Binary Function Evaluation");
#endif

    // Evaluate function
  for (ordinal_type qp=0; qp<nqp; qp++)
    fvals[qp] = func(a, bvals[qp])*quad_weights[qp];

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::QuadExp -- CP Binary Polynomial Integration");
#endif

  // Integrate
  blas.GEMV(Teuchos::NO_TRANS, pc, nqp, 1.0, &sqv[0], sz, &fvals[0], 1, 0.0,
            c.coeff(), 1);

  }
}

template <typename ordinal_type, typename value_type, typename node_type> 
template <typename FuncT>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
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
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::QuadExp -- PC Binary Polynomial Evaluation");
#endif

  // Evaluate input
  blas.GEMV(Teuchos::TRANS, pa, nqp, 1.0, &qv[0], sz, a.coeff(), 1, 0.0,
            &avals[0], 1);

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::QuadExp -- PC Binary Function Evaluation");
#endif

    // Evaluate function
  for (ordinal_type qp=0; qp<nqp; qp++)
    fvals[qp] = func(avals[qp], b)*quad_weights[qp];

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::QuadExp -- PC Binary Polynomial Integration");
#endif

  // Integrate
  blas.GEMV(Teuchos::NO_TRANS, pc, nqp, 1.0, &sqv[0], sz, &fvals[0], 1, 0.0,
            c.coeff(), 1);

  }
}

template <typename ordinal_type, typename value_type, typename node_type> 
template <typename FuncT>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
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
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::QuadExp -- N(" << N << ")-ary Polynomial Evaluation");
#endif

  // Evaluate input
  for (int i=0; i<N; i++) {
    //navals[i].resize(nqp);
    ordinal_type pa = na[i]->size();
    blas.GEMV(Teuchos::TRANS, pa, nqp, 1.0, &qv[0], sz, na[i]->coeff(), 1, 0.0,
	      &navals[N][i][0], 1);
  }

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::QuadExp -- N(" << N << ")-ary Function Evaluation");
#endif

  // Evaluate function
  value_type val[N];
  for (ordinal_type qp=0; qp<nqp; qp++) {
    for (int i=0; i<N; i++)
      val[i] = navals[N][i][qp];
    fvals[qp] = func(val)*quad_weights[qp];
  }

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::QuadExp -- N(" << N << ")-ary Polynomial Integration");
#endif

  // Integrate
  blas.GEMV(Teuchos::NO_TRANS, pc, nqp, 1.0, &sqv[0], sz, &fvals[0], 1, 0.0,
            c.coeff(), 1);

  }
}

template <typename ordinal_type, typename value_type, typename node_type> 
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
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

template <typename ordinal_type, typename value_type, typename node_type> 
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
plusEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	  const value_type& val)
{
  c[0] += val;
}

template <typename ordinal_type, typename value_type, typename node_type> 
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
minusEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	   const value_type& val)
{
  c[0] -= val;
}

template <typename ordinal_type, typename value_type, typename node_type> 
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
timesEqual(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
	   const value_type& val)
{
  ordinal_type pc = c.size();
  value_type* cc = c.coeff();
  for (ordinal_type i=0; i<pc; i++)
    cc[i] *= val;
}

template <typename ordinal_type, typename value_type, typename node_type> 
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
divideEqual(
  Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
  const value_type& val)
{
  ordinal_type pc = c.size();
  value_type* cc = c.coeff();
  for (ordinal_type i=0; i<pc; i++)
    cc[i] /= val;
}

template <typename ordinal_type, typename value_type, typename node_type> 
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
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

template <typename ordinal_type, typename value_type, typename node_type> 
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
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

template <typename ordinal_type, typename value_type, typename node_type> 
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
timesEqual(
  Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
  const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& x)
{
  if (use_quad_for_times) {
    binary_op(times_quad_func(), c, c, x);
    return;
  }

  ordinal_type p = c.size();
  ordinal_type xp = x.size();
  ordinal_type pc;
  if (p > 1 && xp > 1)
    pc = sz;
  else
    pc = p*xp;
  TEST_FOR_EXCEPTION(sz < pc, std::logic_error,
		     "Stokhos::QuadOrthogPolyExpansion::timesEqual()" <<
		     ":  Expansion size (" << sz << 
		     ") is too small for computation.");
  if (c.size() != pc)
    c.resize(pc);

  value_type* cc = c.coeff();
  const value_type* xc = x.coeff();
  
  typename Cijk_type::k_iterator k_begin = Cijk->k_begin();
  typename Cijk_type::k_iterator k_end = Cijk->k_end();
  typename Cijk_type::k_reverse_iterator k_last = Cijk->k_rbegin();
  if (index(k_last) > pc)
    k_end = ++(Cijk->find_k(pc));
  
  if (p > 1 && xp > 1) {
    // Copy c coefficients into temporary array
    value_type* tc = Stokhos::ds_array<value_type>::get_and_fill(cc,p);
    value_type tmp, cijk;
    ordinal_type i,j,k;
     for (typename Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
      k = index(k_it);
      tmp = value_type(0.0);
      for (typename Cijk_type::kj_iterator j_it = Cijk->j_begin(k_it); 
	   j_it != Cijk->j_end(k_it); ++j_it) {
	j = index(j_it);
	if (j < xp) {
	  for (typename Cijk_type::kji_iterator i_it = Cijk->i_begin(j_it);
	       i_it != Cijk->i_end(j_it); ++i_it) {
	    i = index(i_it);
	    cijk = value(i_it);
	    if (i < p)
	      tmp += cijk*tc[i]*xc[j];
	  }
	}
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

template <typename ordinal_type, typename value_type, typename node_type> 
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
divideEqual(
  Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
  const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& x)
{
  if (x.size() == 1) {
    ordinal_type p = c.size();
    value_type* cc = c.coeff();
    const value_type* xc = x.coeff();
    for (ordinal_type i=0; i<p; i++)
      cc[i] /= xc[0];
  }
  else
    binary_op(div_quad_func(), c, c, x);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
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

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
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

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
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

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
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

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
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

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
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

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
times(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (use_quad_for_times) {
    binary_op(times_quad_func(), c, a, b);
    return;
  }

  ordinal_type pa = a.size();
  ordinal_type pb = b.size();
  ordinal_type pc;
  if (pa > 1 && pb > 1)
    pc = sz;
  else
    pc = pa*pb;
  TEST_FOR_EXCEPTION(sz < pc, std::logic_error,
		     "Stokhos::QuadOrthogPolyExpansion::times()" <<
		     ":  Expansion size (" << sz << 
		     ") is too small for computation.");
  if (c.size() != pc)
    c.resize(pc);

  const value_type* ca = a.coeff();
  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  typename Cijk_type::k_iterator k_begin = Cijk->k_begin();
  typename Cijk_type::k_iterator k_end = Cijk->k_end();
  typename Cijk_type::k_reverse_iterator k_last = Cijk->k_rbegin();
  if (index(k_last) > pc)
    k_end = ++(Cijk->find_k(pc));

  if (pa > 1 && pb > 1) {
    value_type tmp, cijk;
    ordinal_type i,j,k;
    for (typename Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
      k = index(k_it);
      tmp = value_type(0.0);
      for (typename Cijk_type::kj_iterator j_it = Cijk->j_begin(k_it); 
	   j_it != Cijk->j_end(k_it); ++j_it) {
	j = index(j_it);
	if (j < pb) {
	  for (typename Cijk_type::kji_iterator i_it = Cijk->i_begin(j_it);
	       i_it != Cijk->i_end(j_it); ++i_it) {
	    i = index(i_it);
	    cijk = value(i_it);
	    if (i < pa)
	      tmp += cijk*ca[i]*cb[j];
	  }
	}
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

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
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

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
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

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
divide(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
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
  else
    binary_op(div_quad_func(), c, a, b);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
divide(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
       const value_type& a, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  binary_op(div_quad_func(), c, a, b);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
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

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
exp(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(exp_quad_func(), c, a);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
log(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(log_quad_func(), c, a);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
log10(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(log10_quad_func(), c, a);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
sqrt(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(sqrt_quad_func(), c, a);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
pow(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  binary_op(pow_quad_func(), c, a, b);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
pow(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const value_type& a, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  binary_op(pow_quad_func(), c, a, b);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
pow(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
    const value_type& b)
{
  binary_op(pow_quad_func(), c, a, b);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
sin(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& s, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(sin_quad_func(), s, a);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
cos(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(cos_quad_func(), c, a);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
tan(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& t, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(tan_quad_func(), t, a);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
sinh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& s, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(sinh_quad_func(), s, a);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
cosh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(cosh_quad_func(), c, a);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
tanh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& t, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(tanh_quad_func(), t, a);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
acos(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(acos_quad_func(), c, a);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
asin(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(asin_quad_func(), c, a);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
atan(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(atan_quad_func(), c, a);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
atan2(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  binary_op(atan2_quad_func(), c, a, b);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
atan2(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const value_type& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  binary_op(atan2_quad_func(), c, a, b);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
atan2(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const value_type& b)
{
  binary_op(atan2_quad_func(), c, a, b);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
acosh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(acosh_quad_func(), c, a);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
asinh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(asinh_quad_func(), c, a);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
atanh(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  unary_op(atanh_quad_func(), c, a);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
fabs(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
     const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a[0] >= 0)
    c = a;
  else
    unaryMinus(c,a);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
abs(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a)
{
  if (a[0] >= 0)
    c = a;
  else
    unaryMinus(c,a);
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
max(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (a[0] >= b[0])
    c = a;
  else
    c = b;
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
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

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
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

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
min(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
    const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b)
{
  if (a[0] <= b[0])
    c = a;
  else
    c = b;
}

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
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

template <typename ordinal_type, typename value_type, typename node_type>
void
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
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

template <typename ordinal_type, typename value_type, typename node_type>
template <typename ExprT1, typename ExprT2>
value_type
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
compute_times_coeff(ordinal_type k, const ExprT1& a, const ExprT2& b) const
{
  ordinal_type pa = a.size();
  ordinal_type pb = b.size();

  if (pa > 1 && pb > 1) {
    value_type cc = value_type(0);
    value_type aa, bb, cijk;
    ordinal_type i,j;
    cc = value_type(0.0);
    ordinal_type n = Cijk->num_values(k);
    for (ordinal_type l=0; l<n; l++) {
      Cijk->value(k,l,i,j,cijk);
      if (i < pa && j < pb) {
	if (i == 0)
	  aa = a.val();
	else
	  aa = a.higher_order_coeff(i);
	if (j == 0)
	  bb = b.val();
	else
	  bb = b.higher_order_coeff(j);
	cc += cijk*aa*bb;
      }
    }
    return cc / basis->norm_squared(k);
  }
  else if (k == 0)
    return a.val() * b.val();
  else if (pa > 1) {
    return a.higher_order_coeff(k)*b.val();
  }
  else {
    return a.val()*b.higher_order_coeff(k);
  }
}

// template <typename ordinal_type, typename value_type, typename node_type>
// template <typename ExprT1, typename ExprT2>
// value_type
// Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
// fast_compute_times_coeff(ordinal_type k, const ExprT1& a, const ExprT2& b) const
// {
//   value_type cc = value_type(0);
//   value_type cijk;
//   ordinal_type i,j;
//   ordinal_type n = Cijk->num_values(k);
//   for (ordinal_type l=0; l<n; l++) {
//     Cijk->value(k,l,i,j,cijk);
//     cc += cijk*a.fastAccessCoeff(i)*b.fastAccessCoeff(j);
//   }
//   return cc / basis->norm_squared(k);
// }

// template <typename ordinal_type, typename value_type, typename node_type>
// template <typename ExprT1, typename ExprT2>
// value_type
// Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
// fast_compute_times_coeff(ordinal_type k, const ExprT1& a, const ExprT2& b) const
// {
//   value_type cc = value_type(0);
//   value_type aa, bb, cijk;
//   ordinal_type i,j;
//   ordinal_type n = Cijk->num_values(k);
//   for (ordinal_type l=0; l<n; l++) {
//     Cijk->value(k,l,i,j,cijk);
//     aa = a.fastAccessCoeff(i);
//     bb = b.fastAccessCoeff(j);
//     cc += cijk*aa*bb;
//   }
//   return cc / basis->norm_squared(k);
// }

template <typename ordinal_type, typename value_type, typename node_type>
template <typename ExprT1, typename ExprT2>
value_type
Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type, node_type>::
fast_compute_times_coeff(ordinal_type k, const ExprT1& a, const ExprT2& b) const
{
  value_type cc = value_type(0);
  value_type aa, bb, cijk;
  ordinal_type i,j;
  ordinal_type n = Cijk->num_values(k);
  for (ordinal_type l=0; l<n; l++) {
    Cijk->value(k,l,i,j,cijk);
    if (i == 0)
      aa = a.val();
    else
      aa = a.fast_higher_order_coeff(i);
    if (j == 0)
      bb = b.val();
    else
      bb = b.fast_higher_order_coeff(j);
    cc += cijk*aa*bb;
  }
  return cc / basis->norm_squared(k);
}
