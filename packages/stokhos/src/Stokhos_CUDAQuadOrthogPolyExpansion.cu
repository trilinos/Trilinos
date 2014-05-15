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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Stokhos_CUDAQuadOrthogPolyExpansion.hpp"

#include <thrust/copy.h>
#include <thrust/transform.h>
#include "cublas.h"

#include "Teuchos_Assert.hpp"
#include "Stokhos_DynamicArrayTraits.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace Stokhos {

  struct negate_func { 
    __host__ __device__ float operator() (const float& a) const { 
      return -a; 
    } 
  };

  struct plus_func { 
    __host__ __device__ float operator() (const float& a, const float& b) const { 
      return a + b; 
    } 
  };

  struct plus_lc_func { 
    float a;
    plus_lc_func(const float& a_) : a(a_) {}
    __host__ __device__ float operator() (const float& b) const { 
      return a + b; 
    } 
  };

  struct plus_rc_func { 
    float b;
    plus_rc_func(const float& b_) : b(b_) {}
    __host__ __device__ float operator() (const float& a) const { 
      return a + b; 
    } 
  };

  struct minus_func { 
    __host__ __device__ float operator() (const float& a, const float& b) const { 
      return a - b; 
    } 
  };

  struct minus_lc_func { 
    float a;
    minus_lc_func(const float& a_) : a(a_) {}
    __host__ __device__ float operator() (const float& b) const { 
      return a - b; 
    } 
  };

  struct minus_rc_func { 
    float b;
    minus_rc_func(const float& b_) : b(b_) {}
    __host__ __device__ float operator() (const float& a) const { 
      return a - b; 
    } 
  };

  struct times_func { 
    __host__ __device__ float operator() (const float& a, const float& b) const { 
      return a * b; 
    } 
  };

  struct times_lc_func { 
    float a;
    times_lc_func(const float& a_) : a(a_) {}
    __host__ __device__ float operator() (const float& b) const { 
      return a * b; 
    } 
  };

  struct times_rc_func { 
    float b;
    times_rc_func(const float& b_) : b(b_) {}
    __host__ __device__ float operator() (const float& a) const { 
      return a * b; 
    } 
  };

  struct div_func { 
    __host__ __device__ float operator() (const float& a, const float& b) const { 
      return a / b; 
    } 
  };

  struct div_lc_func { 
    float a;
    div_lc_func(const float& a_) : a(a_) {}
    __host__ __device__ float operator() (const float& b) const { 
      return a / b; 
    } 
  };

  struct div_rc_func { 
    float b;
    div_rc_func(const float& b_) : b(b_) {}
    __host__ __device__ float operator() (const float& a) const { 
      return a / b; 
    } 
  };

  struct div_rc_func_dev { 
    typedef typename thrust::device_vector<float>::const_reference const_reference;
    //const_reference b;
    const float *b;
    div_rc_func_dev(const float* b_) : b(b_) {}
    __host__ __device__ float operator() (const float& a) const { 
      return a / *b; 
    } 
  };

  struct exp_func { 
    __host__ __device__ float operator() (const float& a) const { 
      return std::exp(a); 
    } 
  };

  struct log_func { 
    __host__ __device__ float operator() (const float& a) const { 
      return std::log(a); 
    } 
  };

  struct log10_func { 
    __host__ __device__ float operator() (const float& a) const { 
      return std::log10(a); 
    } 
  };
    
  struct sqrt_func { 
    __host__ __device__ float operator() (const float& a) const { 
      return std::sqrt(a); 
    } 
  };

  struct pow_func { 
    __host__ __device__ float operator() (const float& a, const float& b) const { 
      return std::pow(a,b); 
    } 
  };

  struct pow_lc_func { 
    float a;
    pow_lc_func(const float& a_) : a(a_) {}
    __host__ __device__ float operator() (const float& b) const { 
      return std::pow(a,b); 
    } 
  };

  struct pow_rc_func { 
    float b;
    pow_rc_func(const float& b_) : b(b_) {}
    __host__ __device__ float operator() (const float& a) const { 
      return std::pow(a,b); 
    } 
  };

  struct cos_func { 
    __host__ __device__ float operator() (const float& a) const { 
      return std::cos(a); 
    } 
  };

  struct sin_func { 
    __host__ __device__ float operator() (const float& a) const { 
      return std::sin(a); 
    } 
  };

  struct tan_func { 
    __host__ __device__ float operator() (const float& a) const { 
      return std::tan(a); 
    } 
  };

  struct cosh_func { 
    __host__ __device__ float operator() (const float& a) const { 
      return std::cosh(a); 
    } 
  };

  struct sinh_func { 
    __host__ __device__ float operator() (const float& a) const { 
      return std::sinh(a); 
    } 
  };

  struct tanh_func { 
    __host__ __device__ float operator() (const float& a) const { 
      return std::tanh(a); 
    } 
  };

  struct acos_func { 
    __host__ __device__ float operator() (const float& a) const { 
      return std::acos(a); 
    } 
  };

  struct asin_func { 
    __host__ __device__ float operator() (const float& a) const { 
      return std::asin(a); 
    } 
  };

  struct atan_func { 
    __host__ __device__ float operator() (const float& a) const { 
      return std::atan(a); 
    } 
  };

  struct atan2_func { 
    __host__ __device__ float operator() (const float& a, const float& b) const { 
      return std::atan2(a,b); 
    } 
  };

  struct atan2_lc_func { 
    float a;
    atan2_lc_func(const float& a_) : a(a_) {}
    __host__ __device__ float operator() (const float& b) const { 
      return std::atan2(a,b); 
    } 
  };

  struct atan2_rc_func { 
    float b;
    atan2_rc_func(const float& b_) : b(b_) {}
    __host__ __device__ float operator() (const float& a) const { 
      return std::atan2(a,b); 
    } 
  };

  struct acosh_func { 
    __host__ __device__ float operator() (const float & a) const { 
      return std::log(a+std::sqrt(a*a-float(1.0))); 
    }
  };

  struct asinh_func { 
    __host__ __device__ float operator() (const float& a) const { 
      return std::log(a+std::sqrt(a*a+float(1.0))); 
    }
  };

  struct atanh_func { 
    __host__ __device__ float operator() (const float& a) const { 
      return 0.5*std::log((float(1.0)+a)/(float(1.0)-a)); 
    } 
  };
}

Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
QuadOrthogPolyExpansion(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, float> >& basis_,
  const Teuchos::RCP<const Stokhos::Sparse3Tensor<int, float> >& Cijk_,
  const Teuchos::RCP<const Quadrature<int, float> >& quad_,
  bool use_quad_for_times_) :
  basis(basis_),
  Cijk(Cijk_),
  quad(quad_),
  use_quad_for_times(use_quad_for_times_),
  sz(basis->size()),
  quad_points(quad->getQuadPoints()),
  quad_weights(quad->getQuadWeights()),
  quad_values(quad->getBasisAtQuadPoints()),
  norms(basis->norm_squared()),
  nqp(quad_points.size()),
  avals(nqp),
  bvals(nqp),
  fvals(nqp),
  host_qv(nqp*sz),
  host_sqv(nqp*sz),
  qv(nqp*sz),
  sqv(nqp*sz)
{
  for (int qp=0; qp<nqp; qp++) {
    for (int i=0; i<sz; i++) {
      host_qv[qp*sz+i] = quad_values[qp][i];
      host_sqv[qp*sz+i] = quad_weights[qp]*quad_values[qp][i]/norms[i];
    }
  }
  thrust::copy(host_qv.begin(), host_qv.end(), qv.begin());
  thrust::copy(host_sqv.begin(), host_sqv.end(), sqv.begin());

  cublasInit();
}

Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
~QuadOrthogPolyExpansion()
{
  cublasShutdown();
}
 
void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
unaryMinus(
  Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
  const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a)
{
  int pc = a.size();
  if (c.size() != pc)
    c.resize(pc);

  thrust::transform(a.coeff(), a.coeff()+pc, c.coeff(), negate_func());
}

 
void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
plusEqual(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
	  const float& val)
{
  thrust::transform(c.coeff(), c.coeff()+1, c.coeff(), plus_rc_func(val));
}

 
void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
minusEqual(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
	   const float& val)
{
  thrust::transform(c.coeff(), c.coeff()+1, c.coeff(), minus_rc_func(val));
}

 
void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
timesEqual(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
	   const float& val)
{
  int pc = c.size();
  thrust::transform(c.coeff(), c.coeff()+pc, c.coeff(), times_rc_func(val));
}

 
void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
divideEqual(
  Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
  const float& val)
{
  int pc = c.size();
  thrust::transform(c.coeff(), c.coeff()+pc, c.coeff(), div_rc_func(val));
}

 
void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
plusEqual(
  Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
  const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& x)
{
  int xp = x.size();
  if (c.size() < xp)
    c.resize(xp);

  thrust::transform(c.coeff(), c.coeff()+xp, x.coeff(), c.coeff(), plus_func());
}

 
void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
minusEqual(
  Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
  const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& x)
{
  int xp = x.size();
  if (c.size() < xp)
    c.resize(xp);

  thrust::transform(c.coeff(), c.coeff()+xp, x.coeff(), c.coeff(), 
		    minus_func());
}

 
void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
timesEqual(
  Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
  const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& x)
{
  if (use_quad_for_times) {
    binary_op(times_func(), c, c, x);
    return;
  }

  int p = c.size();
  int xp = x.size();
  int pc;
  if (p > 1 && xp > 1)
    pc = sz;
  else
    pc = p*xp;
  TEUCHOS_TEST_FOR_EXCEPTION(sz < pc, std::logic_error,
		     "Stokhos::QuadOrthogPolyExpansion::timesEqual()" <<
		     ":  Expansion size (" << sz << 
		     ") is too small for computation.");
  if (c.size() != pc)
    c.resize(pc);

  pointer cc = c.coeff();
  const_pointer xc = x.coeff();
  
  if (p > 1 && xp > 1) {
    // Copy c coefficients into temporary array
    thrust::device_vector<float> t(c.coeff(), c.coeff()+p);
    pointer tc = t.data();

    typename Cijk_type::i_iterator i_begin = Cijk->i_begin();
    typename Cijk_type::i_iterator i_end = Cijk->i_end();
    if (pc < Cijk->num_i())
      i_end = Cijk->find_i(pc);
    int k_lim = p;
    int j_lim = xp;
    const_pointer kc = tc;
    const_pointer jc = xc;
    if (xp < p) {
      k_lim = xp;
      j_lim = p;
      kc = xc;
      jc = tc;
    }

    float tmp, cijk;
    int i,j,k;
    for (typename Cijk_type::i_iterator i_it=i_begin; i_it!=i_end; ++i_it) {
      i = index(i_it);
      tmp = float(0.0);
      for (typename Cijk_type::ik_iterator k_it = Cijk->k_begin(i_it); 
	   k_it != Cijk->k_end(i_it); ++k_it) {
	k = index(k_it);
	if (k < k_lim) {
	  for (typename Cijk_type::ikj_iterator j_it = Cijk->j_begin(k_it);
	       j_it != Cijk->j_end(k_it); ++j_it) {
	    j = index(j_it);
	    cijk = value(j_it);
	    if (j < j_lim)
	      tmp += cijk*kc[k]*jc[j];
	  }
	}
      }
      cc[i] = tmp / basis->norm_squared(i);
    }
  }
  else if (p > 1) {
    for (int i=0; i<p; i++)
      cc[i] *= xc[0];
  }
  else if (xp > 1) {
    for (int i=1; i<xp; i++)
      cc[i] = cc[0]*xc[i];
    cc[0] *= xc[0];
  }
  else {
    cc[0] *= xc[0];
  }
}

 
void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
divideEqual(
  Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
  const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& x)
{
  if (x.size() == 1) {
    int p = c.size();
    pointer cc = c.coeff();
    const_pointer xc = x.coeff();
    thrust::transform(cc, cc+p, cc, 
		      div_rc_func_dev(thrust::raw_pointer_cast(xc)));
  }
  else
    binary_op(div_func(), c, c, x);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
plus(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
     const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
     const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b)
{
  int pa = a.size();
  int pb = b.size();
  int pc = pa > pb ? pa : pb;
  if (c.size() != pc)
    c.resize(pc);

  const_pointer ca = a.coeff();
  const_pointer cb = b.coeff();
  pointer cc = c.coeff();

  if (pa > pb) {
    thrust::transform(ca, ca+pb, cb, cc, plus_func());
    thrust::copy(ca+pb, ca+pc, cc+pb);
  }
  else {
    thrust::transform(ca, ca+pa, cb, cc, plus_func());
    thrust::copy(cb+pa, cb+pc, cc+pa);
  }
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
plus(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
     const float& a, 
     const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b)
{
  int pc = b.size();
  if (c.size() != pc)
    c.resize(pc);

  const_pointer cb = b.coeff();
  pointer cc = c.coeff();

  thrust::transform(cb, cb+1, cc, plus_lc_func(a));
  thrust::copy(cb+1, cb+pc, cc+1);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
plus(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
     const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
     const float& b)
{
  int pc = a.size();
  if (c.size() != pc)
    c.resize(pc);

  const_pointer ca = a.coeff();
  pointer cc = c.coeff();

  thrust::transform(ca, ca+1, cc, plus_rc_func(b));
  thrust::copy(ca+1, ca+pc, cc+1);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
minus(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
      const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
      const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b)
{
  int pa = a.size();
  int pb = b.size();
  int pc = pa > pb ? pa : pb;
  if (c.size() != pc)
    c.resize(pc);

  const_pointer ca = a.coeff();
  const_pointer cb = b.coeff();
  pointer cc = c.coeff();

  if (pa > pb) {
    thrust::transform(ca, ca+pb, cb, cc, minus_func());
    thrust::copy(ca+pb, ca+pc, cc+pb);
  }
  else {
    thrust::transform(ca, ca+pa, cb, cc, minus_func());
    thrust::transform(cb+pa, cb+pc, cc+pa, negate_func());
  }
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
minus(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
      const float& a, 
      const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b)
{
  int pc = b.size();
  if (c.size() != pc)
    c.resize(pc);

  const_pointer cb = b.coeff();
  pointer cc = c.coeff();

  thrust::transform(cb, cb+1, cc, minus_lc_func(a));
  thrust::transform(cb+1, cb+pc, cc+1, negate_func());
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
minus(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
      const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
      const float& b)
{
  int pc = a.size();
  if (c.size() != pc)
    c.resize(pc);

  const_pointer ca = a.coeff();
  pointer cc = c.coeff();

  thrust::transform(ca, ca+1, cc, minus_rc_func(b));
  thrust::copy(ca+1, ca+pc, cc+1);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
times(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
      const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
      const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b)
{
  if (use_quad_for_times) {
    binary_op(times_func(), c, a, b);
    return;
  }

  int pa = a.size();
  int pb = b.size();
  int pc;
  if (pa > 1 && pb > 1)
    pc = sz;
  else
    pc = pa*pb;
  TEUCHOS_TEST_FOR_EXCEPTION(sz < pc, std::logic_error,
		     "Stokhos::QuadOrthogPolyExpansion::times()" <<
		     ":  Expansion size (" << sz << 
		     ") is too small for computation.");
  if (c.size() != pc)
    c.resize(pc);

  const_pointer ca = a.coeff();
  const_pointer cb = b.coeff();
  pointer cc = c.coeff();

  if (pa > 1 && pb > 1) {
    typename Cijk_type::i_iterator i_begin = Cijk->i_begin();
    typename Cijk_type::i_iterator i_end = Cijk->i_end();
    if (pc < Cijk->num_i())
      i_end = Cijk->find_i(pc);
    int k_lim = pa;
    int j_lim = pb;
    const_pointer kc = ca;
    const_pointer jc = cb;
    if (pb < pa) {
      k_lim = pb;
      j_lim = pa;
      kc = cb;
      jc = ca;
    }

    float tmp, cijk;
    int i,j,k;
    for (typename Cijk_type::i_iterator i_it=i_begin; i_it!=i_end; ++i_it) {
      i = index(i_it);
      tmp = float(0.0);
      for (typename Cijk_type::ik_iterator k_it = Cijk->k_begin(i_it); 
	   k_it != Cijk->k_end(i_it); ++k_it) {
	k = index(k_it);
	if (k < k_lim) {
	  for (typename Cijk_type::ikj_iterator j_it = Cijk->j_begin(k_it);
	       j_it != Cijk->j_end(k_it); ++j_it) {
	    j = index(j_it);
	    cijk = value(j_it);
	    if (j < j_lim)
	      tmp += cijk*kc[k]*jc[j];
	  }
	}
      }
      cc[i] = tmp / basis->norm_squared(i);
    }
  }
  else if (pa > 1) {
    for (int i=0; i<pc; i++)
      cc[i] = ca[i]*cb[0];
  }
  else if (pb > 1) {
    for (int i=0; i<pc; i++)
      cc[i] = ca[0]*cb[i];
  }
  else {
    cc[0] = ca[0]*cb[0];
  }
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
times(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
      const float& a, 
      const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b)
{
  int pc = b.size();
  if (c.size() != pc)
    c.resize(pc);

  const_pointer cb = b.coeff();
  pointer cc = c.coeff();

  thrust::transform(cb, cb+pc, cc, times_lc_func(a));
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
times(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
      const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
      const float& b)
{
  int pc = a.size();
  if (c.size() != pc)
    c.resize(pc);

  const_pointer ca = a.coeff();
  pointer cc = c.coeff();

  thrust::transform(ca, ca+pc, cc, times_rc_func(b));
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
divide(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
       const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
       const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b)
{
  if (b.size() == 1) {
    int pc = a.size();
    if (c.size() != pc)
      c.resize(pc);

    const_pointer ca = a.coeff();
    const_pointer cb = b.coeff();
    pointer cc = c.coeff();

    thrust::transform(ca, ca+pc, cc, 
		      div_rc_func_dev(thrust::raw_pointer_cast(cb)));
  }
  else
    binary_op(div_func(), c, a, b);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
divide(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
       const float& a, 
       const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b)
{
  unary_op(div_lc_func(a), c, b);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
divide(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
       const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
       const float& b)
{
  int pc = a.size();
  if (c.size() != pc)
    c.resize(pc);

  const_pointer ca = a.coeff();
  pointer cc = c.coeff();

  thrust::transform(ca, ca+pc, cc, div_rc_func(b));
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
exp(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
    const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a)
{
  unary_op(exp_func(), c, a);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
log(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
    const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a)
{
  unary_op(log_func(), c, a);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
log10(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
      const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a)
{
  unary_op(log10_func(), c, a);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
sqrt(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
     const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a)
{
  unary_op(sqrt_func(), c, a);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
pow(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
    const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a,
    const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b)
{
  binary_op(pow_func(), c, a, b);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
pow(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
    const float& a, 
    const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b)
{
  unary_op(pow_lc_func(a), c, b);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
pow(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
    const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
    const float& b)
{
  unary_op(pow_rc_func(b), c, a);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
sin(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& s, 
    const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a)
{
  unary_op(sin_func(), s, a);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
cos(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
    const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a)
{
  unary_op(cos_func(), c, a);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
tan(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& t, 
    const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a)
{
  unary_op(tan_func(), t, a);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
sinh(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& s, 
     const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a)
{
  unary_op(sinh_func(), s, a);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
cosh(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
     const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a)
{
  unary_op(cosh_func(), c, a);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
tanh(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& t, 
     const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a)
{
  unary_op(tanh_func(), t, a);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
acos(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
     const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a)
{
  unary_op(acos_func(), c, a);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
asin(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
     const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a)
{
  unary_op(asin_func(), c, a);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
atan(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
     const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a)
{
  unary_op(atan_func(), c, a);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
atan2(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
      const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a,
      const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b)
{
  binary_op(atan2_func(), c, a, b);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
atan2(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
      const float& a, 
      const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b)
{
  unary_op(atan2_lc_func(a), c, b);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
atan2(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
      const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
      const float& b)
{
  unary_op(atan2_rc_func(b), c, a);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
acosh(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
      const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a)
{
  unary_op(acosh_func(), c, a);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
asinh(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
      const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a)
{
  unary_op(asinh_func(), c, a);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
atanh(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
      const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a)
{
  unary_op(atanh_func(), c, a);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
fabs(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
     const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a)
{
  if (a[0] >= 0)
    c = a;
  else
    unaryMinus(c,a);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
abs(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
    const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a)
{
  if (a[0] >= 0)
    c = a;
  else
    unaryMinus(c,a);
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
max(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
    const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a,
    const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b)
{
  if (a[0] >= b[0])
    c = a;
  else
    c = b;
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
max(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
    const float& a, 
    const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b)
{
  if (a >= b[0]) {
    c = OrthogPolyApprox<int, float, CUDAStorage<int,float> >(basis);
    c[0] = a;
  }
  else
    c = b;
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
max(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
    const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
    const float& b)
{
  if (a[0] >= b)
    c = a;
  else {
    c = OrthogPolyApprox<int, float, CUDAStorage<int,float> >(basis);
    c[0] = b;
  }
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
min(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
    const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a,
    const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b)
{
  if (a[0] <= b[0])
    c = a;
  else
    c = b;
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
min(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
    const float& a, 
    const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b)
{
  if (a <= b[0]) {
    c = OrthogPolyApprox<int, float, CUDAStorage<int,float> >(basis);
    c[0] = a;
  }
  else
    c = b;
}


void
Stokhos::QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >::
min(Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
    const Stokhos::OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
    const float& b)
{
  if (a[0] <= b)
    c = a;
  else {
    c = OrthogPolyApprox<int, float, CUDAStorage<int,float> >(basis);
    c[0] = b;
  }
}
