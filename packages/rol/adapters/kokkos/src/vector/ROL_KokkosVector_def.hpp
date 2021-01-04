// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#pragma once
#ifndef ROL_KOKKOSVECTOR_DEF_HPP
#define ROL_KOKKOSVECTOR_DEF_HPP

namespace ROL {

template<typename Real, typename Device>
void KokkosVector<Real,Device>::applyBinary( const Elementwise::BinaryFunction<Real>& bf, 
                                             const Vector<Real>& x )  {
  data_.sync_device();
  auto binary = Elementwise::KokkosBinaryFunction<Real,Device>::create(bf);
  binary->apply(*this,cast(x));
  data_.modify_device();
}


template<typename Real, typename Device>
void KokkosVector<Real,Device>::applyUnary( const Elementwise::UnaryFunction<Real>& uf ) {
  data_.sync_device();
  auto unary = Elementwise::KokkosUnaryFunction<Real,Device>::create(uf);
  unary->apply(*this);
  data_.modify_device();
}

template<typename Real, typename Device>
void KokkosVector<Real,Device>::axpy( const Real alpha, const Vector<Real>& x )  {
  data_.sync_device();
  device_view_type d  = data_.view_device();
  device_view_type xd = cast(x).view_device();
  Kokkos::parallel_for( dimension(), KOKKOS_LAMBDA( const int i ) {
    d(i) += alpha * xd(i);
  });
  data_.modify_device();
} 

template<typename Real, typename Device>
Ptr<Vector<Real>> KokkosVector<Real,Device>::basis( int i ) const  {
  auto bp = makePtr<KokkosVector>(dimension());
  device_view_type bpd = bp->view_device();
  Kokkos::parallel_for( dimension(), KOKKOS_LAMBDA( const int j ) {
    bpd(j) = static_cast<Real>(i==j);
  });
  return bp;
} 

template<typename Real, typename Device>
Ptr<Vector<Real>> KokkosVector<Real,Device>::clone() const {
  return makePtr<KokkosVector>(dimension());
} 

template<typename Real, typename Device>
int KokkosVector<Real,Device>::dimension() const { 
  return static_cast<int>(data_.extent(0)); 
}

template<typename Real, typename Device>
Real KokkosVector<Real,Device>::dot( const Vector<Real>& x ) const  { 
  data_.sync_device();
  device_view_type d  = data_.view_device();
  device_view_type xd = cast(x).view_device();
  Real result = 0;
  Kokkos::parallel_reduce( dimension(), KOKKOS_LAMBDA( const int i, Real& lresult ) {
    lresult += d(i) * xd(i);
  }, result);
  return result;
} 

template<typename Real, typename Device>
Real KokkosVector<Real,Device>::norm() const  { 
  data_.sync_device();
  device_view_type d = data_.view_device();
  Real result = 0;
  Kokkos::parallel_reduce( dimension(), KOKKOS_LAMBDA( const int i, Real& lresult ) {
    lresult += d(i) * d(i);
  }, result);
  return std::sqrt(result);
} // end norm 


template<typename Real, typename Device>
void KokkosVector<Real,Device>::plus( const Vector<Real>& x )  {
  data_.sync_device();
  device_view_type d  = data_.view_device();
  device_view_type xd = cast(x).view_device();
  Kokkos::parallel_for( dimension(), KOKKOS_LAMBDA( const int i ) {
    d(i) += xd(i);
  });
  data_.modify_device();
} 

template<typename Real, typename Device>
void KokkosVector<Real,Device>::randomize( const Real l, const Real u )  {
  host_view_type h = data_.view_host();
  for( int i=0; i<dimension(); ++i ) {
    auto x = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
    h(i) = u*x + l*(1-x);
  }
  data_.modify_host();
} 
 
template<typename Real, typename Device>
Real KokkosVector<Real,Device>::reduce( const Elementwise::ReductionOp<Real>& r ) const  {
  data_.sync_device();
  auto reducer = Elementwise::KokkosReductionOp<Real,Device>::create(r);
  return reducer->reduce(*this);
} 

template<typename Real, typename Device>
void KokkosVector<Real,Device>::scale( const Real alpha )  {
  device_view_type d  = data_.view_device();
  Kokkos::parallel_for( dimension(), KOKKOS_LAMBDA( const int i ) {
    d(i) *= alpha;
  });
  data_.modify_device();
} 

template<typename Real, typename Device>
void KokkosVector<Real,Device>::set( const Vector<Real>& x )  {
  device_view_type d  = data_.view_device();
  device_view_type xd = cast(x).view_device();
  Kokkos::parallel_for( dimension(), KOKKOS_LAMBDA( const int i ) {
    d(i) = xd(i);
  });
  data_.modify_device();
}

template<typename Real, typename Device>
void KokkosVector<Real,Device>::setScalar( const Real alpha )  {
  device_view_type d  = data_.view_device();
  Kokkos::parallel_for( dimension(), KOKKOS_LAMBDA( const int i ) {
    d(i) = alpha;
  });
  data_.modify_device();
} 

template<typename Real, typename Device>
void KokkosVector<Real,Device>::zero() {
  device_view_type d = data_.view_device();
  Kokkos::parallel_for( data_.extent(0), KOKKOS_LAMBDA( const int i ) {
    d(i) = 0;
  });
  data_.modify_device();
} 


template<typename Real, typename Device>
void KokkosVector<Real,Device>::print( std::ostream& os ) const {
  data_.sync_host();
  host_view_type h = data_.view_host(); 
  for( int i=0; i<dimension(); ++i ) os << data_.h_view(i) << " ";
  os << std::endl;
} 

} // namespace ROL

#endif  // ROL_KOKKOSVECTOR_DEF_HPP
