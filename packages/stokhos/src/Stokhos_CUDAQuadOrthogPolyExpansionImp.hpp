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

#include <thrust/transform.h>
#include "cublas.h"

namespace Stokhos {

  template <int N, typename array_type> struct make_tuple_N {};
  
  template <typename array_type>
  struct make_tuple_N<1, array_type> {
    typedef typename array_type::value_type::const_iterator T;
    static thrust::tuple<T> 
    begin(const array_type& a) { 
      return thrust::make_tuple(a[0].begin());
    }
    static thrust::tuple<T> 
    end(const array_type& a) { 
      return thrust::make_tuple(a[0].end());
    }
  };

  template <typename array_type>
  struct make_tuple_N<2, array_type> {
    typedef typename array_type::value_type::const_iterator T;
    static thrust::tuple<T,T> 
    begin(const array_type& a) { 
      return thrust::make_tuple(a[0].begin(), a[1].begin());
    }
    static thrust::tuple<T,T> 
    end(const array_type& a) { 
      return thrust::make_tuple(a[0].end(), a[1].end());
    }
  };

  template <typename array_type>
  struct make_tuple_N<3, array_type> {
    typedef typename array_type::value_type::const_iterator T;
    static thrust::tuple<T,T,T> 
    begin(const array_type& a) { 
      return thrust::make_tuple(a[0].begin(), a[1].begin(), a[2].begin());
    }
    static thrust::tuple<T,T,T> 
    end(const array_type& a) { 
      return thrust::make_tuple(a[0].end(), a[1].end(), a[2].end());
    }
  };

  template <typename array_type>
  struct make_tuple_N<4, array_type> {
    typedef typename array_type::value_type::const_iterator T;
    static thrust::tuple<T,T,T,T> 
    begin(const array_type& a) { 
      return thrust::make_tuple(a[0].begin(), a[1].begin(), a[2].begin(),
				a[3].begin());
    }
    static thrust::tuple<T,T,T,T> 
    end(const array_type& a) { 
      return thrust::make_tuple(a[0].end(), a[1].end(), a[2].end(), a[3].end());
    }
  };
}

template <typename FuncT>
void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
unary_op(const FuncT& func,
         OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
         const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a)
{
  int pa = a.size();
  int pc;
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

  float *qv_ptr = thrust::raw_pointer_cast(qv.data());
  float *sqv_ptr = thrust::raw_pointer_cast(sqv.data());
  const float *a_ptr = thrust::raw_pointer_cast(a.coeff());
  float *c_ptr = thrust::raw_pointer_cast(c.coeff());
  float *avals_ptr = thrust::raw_pointer_cast(avals.data());
  float *fvals_ptr = thrust::raw_pointer_cast(fvals.data());

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::QuadExp -- Unary Polynomial Evaluation");
#endif

  // Evaluate input
  cublasSgemv('T', pa, nqp, 1.0, qv_ptr, sz, a_ptr, 1, 0.0, avals_ptr, 1);

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::QuadExp -- Unary Function Evaluation");
#endif

  // Evaluate function
  thrust::transform(avals.begin(), avals.end(), fvals.begin(), func);
  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::QuadExp -- Unary Polynomial Integration");
#endif

  // Integrate
  cublasSgemv('N', pc, nqp, 1.0, sqv_ptr, sz, fvals_ptr, 1, 0.0, c_ptr, 1);

  }
}
 
template <typename FuncT>
void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
binary_op(const FuncT& func,
          OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
          const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
          const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b)
{
  int pa = a.size();
  int pb = b.size();
  int pc;
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

  float *qv_ptr = thrust::raw_pointer_cast(qv.data());
  float *sqv_ptr = thrust::raw_pointer_cast(sqv.data());
  const float *a_ptr = thrust::raw_pointer_cast(a.coeff());
  const float *b_ptr = thrust::raw_pointer_cast(b.coeff());
  float *c_ptr = thrust::raw_pointer_cast(c.coeff());
  float *avals_ptr = thrust::raw_pointer_cast(avals.data());
  float *bvals_ptr = thrust::raw_pointer_cast(bvals.data());
  float *fvals_ptr = thrust::raw_pointer_cast(fvals.data());

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::QuadExp --Binary Polynomial Evaluation");
#endif

  // Evaluate input
  cublasSgemv('T', pa, nqp, 1.0, qv_ptr, sz, a_ptr, 1, 0.0, avals_ptr, 1);
  cublasSgemv('T', pb, nqp, 1.0, qv_ptr, sz, b_ptr, 1, 0.0, bvals_ptr, 1);

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::QuadExp -- Binary Function Evaluation");
#endif

  // Evaluate function
  thrust::transform(avals.begin(), avals.end(), bvals.begin(), fvals.begin(), 
		    func);

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::QuadExp -- Binary Polynomial Integration");
#endif

  // Integrate
  cublasSgemv('N', pc, nqp, 1.0, sqv_ptr, sz, fvals_ptr, 1, 0.0, c_ptr, 1);

  }
}

template <typename FuncT>
void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int,float> >::
nary_op(const FuncT& func,
	OrthogPolyApprox<int, float, Stokhos::CUDAStorage<int,float> >& c, 
	const OrthogPolyApprox<int, float, Stokhos::CUDAStorage<int,float> >** na)
{
  const int N = FuncT::N;
  bool is_constant = true;
  for (int i=0; i<N; i++) {
    if (na[i]->size() > 1) {
      is_constant = false;
      break;
    }
  }
  int pc;
  if (is_constant)
    pc = 1;
  else
    pc = sz;
  if (c.size() != pc)
    c.resize(pc);

  if (pc == 1) {
    float val[N];
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

  float *qv_ptr = thrust::raw_pointer_cast(qv.data());
  float *sqv_ptr = thrust::raw_pointer_cast(sqv.data());
  float *c_ptr = thrust::raw_pointer_cast(c.coeff());
  float *fvals_ptr = thrust::raw_pointer_cast(fvals.data());

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::QuadExp -- N(" << N << ")-ary Polynomial Evaluation");
#endif

  // Evaluate input
  for (int i=0; i<N; i++) {
    int pa = na[i]->size();
    const float *na_ptr = thrust::raw_pointer_cast(na[i]->coeff());
    float *navals_ptr = thrust::raw_pointer_cast(navals[N][i].data());
    cublasSgemv('T', pa, nqp, 1.0, qv_ptr, sz, na_ptr, 1, 0.0, navals_ptr, 1);
  }

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::QuadExp -- N(" << N << ")-ary Function Evaluation");
#endif

  // Evaluate function
  thrust::transform(
    thrust::make_zip_iterator(make_tuple_N<N,navals_type>::begin(navals[N])), 
    thrust::make_zip_iterator(make_tuple_N<N,navals_type>::end(navals[N])), 
    fvals.begin(), 
    func);

  }

  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::QuadExp -- N(" << N << ")-ary Polynomial Integration");
#endif

  // Integrate
  cublasSgemv('N', pc, nqp, 1.0, sqv_ptr, sz, fvals_ptr, 1, 0.0, c_ptr, 1);

  }
}

template <typename ExprT1, typename ExprT2>
float
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int,float> >::
compute_times_coeff(int k, const ExprT1& a, const ExprT2& b) const
{
  int pa = a.size();
  int pb = b.size();

  typename Cijk_type::k_iterator k_it = Cijk->find_k(k);

  if (pa > 1 && pb > 1) {
    float cc = float(0);
    float aa, bb, cijk;
    int i,j;
    cc = float(0.0);
    for (typename Cijk_type::kj_iterator j_it = Cijk->j_begin(k_it); 
	 j_it != Cijk->j_end(k_it); ++j_it) {
      j = index(j_it);
      if (j < pb) {
	if (j == 0)
	  bb = b.val();
	else
	  bb = b.higher_order_coeff(j);
	for (typename Cijk_type::kji_iterator i_it = Cijk->i_begin(j_it);
	     i_it != Cijk->i_end(j_it); ++i_it) {
	  i = index(i_it);
	  cijk = value(i_it);
	  if (i < pa) {
	    if (i == 0)
	      aa = a.val();
	    else
	      aa = a.higher_order_coeff(i);
	  }
	  cc += cijk*aa*bb;
	}
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

template <typename ExprT1, typename ExprT2>
float
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int,float> >::
fast_compute_times_coeff(int k, const ExprT1& a, const ExprT2& b) const
{
  float cc = float(0);
  float aa, bb, cijk;
  int i,j;
  typename Cijk_type::k_iterator k_it = Cijk->find_k(k);
  for (typename Cijk_type::kj_iterator j_it = Cijk->j_begin(k_it); 
       j_it != Cijk->j_end(k_it); ++j_it) {
    j = index(j_it);
    if (j == 0)
      bb = b.val();
    else
      bb = b.fast_higher_order_coeff(j);
    for (typename Cijk_type::kji_iterator i_it = Cijk->i_begin(j_it);
	 i_it != Cijk->i_end(j_it); ++i_it) {
      i = index(i_it);
      cijk = value(i_it);
      if (i == 0)
	aa = a.val();
      else
	aa = a.fast_higher_order_coeff(i);
      cc += cijk*aa*bb;
    }
  }
  return cc / basis->norm_squared(k);
}
