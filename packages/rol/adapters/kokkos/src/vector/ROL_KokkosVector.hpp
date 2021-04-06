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

#ifndef ROL_KOKKOSVECTOR_HPP
#define ROL_KOKKOSVECTOR_HPP

#include "ROL_Vector.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_DualView.hpp"

/** @ingroup la_group
    \class ROL::KokkosVector
    \brief Defines the ROL interface for Kokkos::View
*/

namespace ROL {

template<typename Real, typename Device>
class KokkosVector 
  : public ROL::Vector<Real> {
public:

  using view_type = Kokkos::DualView<Real*,Device>;
  using host_view_type = typename view_type::t_host;
  using device_view_type = typename view_type::t_dev;

  KokkosVector( int n ) : data_("data",n) {} 
  virtual ~KokkosVector() {}

  void copy_host_to_device() const {
    Kokkos::deep_copy(data_.d_view,data_.h_view);
  }

  void copy_device_to_host() const {
    Kokkos::deep_copy(data_.h_view,data_.d_view);
  }

  host_view_type   get_host_view()   const { return data_.view_host();   }
  device_view_type get_device_view() const { return data_.view_device(); }

  int dimension() const override { return static_cast<int>(data_.extent(0)); }

  void zero() override {
    device_view_type d = get_device_view();
    Kokkos::parallel_for( data_.extent(0), KOKKOS_LAMBDA( const int i ) {
      d(i) = 0;
    });
  } // end zero

  void set( const Vector<Real>& x ) override {
    device_view_type d  = get_device_view();
    device_view_type xd = cast(x).get_device_view();
    
    Kokkos::parallel_for( dimension(), KOKKOS_LAMBDA( const int i ) {
      d(i) = xd(i);
    });
  } // end set

  void plus( const Vector<Real>& x ) override {
    device_view_type d  = get_device_view();
    device_view_type xd = cast(x).get_device_view();
    Kokkos::parallel_for( dimension(), KOKKOS_LAMBDA( const int i ) {
      d(i) += xd(i);
    });
  } // end plus

  void axpy( const Real alpha, const Vector<Real>& x ) override {
    device_view_type d  = get_device_view();
    device_view_type xd = cast(x).get_device_view();
    Kokkos::parallel_for( dimension(), KOKKOS_LAMBDA( const int i ) {
      d(i) += alpha * xd(i);
    });
  } // end axpy

  void setScalar( const Real alpha ) override {
    device_view_type d  = get_device_view();
    Kokkos::parallel_for( dimension(), KOKKOS_LAMBDA( const int i ) {
      d(i) = alpha;
    });
  } // end setScalar

  void scale( const Real alpha ) override {
    device_view_type d  = get_device_view();
    Kokkos::parallel_for( dimension(), KOKKOS_LAMBDA( const int i ) {
      d(i) *= alpha;
    });
  } // end scale

  Real dot( const Vector<Real>& x ) const override { 
    device_view_type d  = get_device_view();
    device_view_type xd = cast(x).get_device_view();
    Real result = 0;
    Kokkos::parallel_reduce( dimension(), KOKKOS_LAMBDA( const int i, Real& lresult ) {
      lresult += d(i) * xd(i);
    }, result);
    return result;
  } // end dot

  Real norm() const override { 
    device_view_type d  = get_device_view();
    Real result = 0;
    Kokkos::parallel_reduce( dimension(), KOKKOS_LAMBDA( const int i, Real& lresult ) {
      lresult += d(i) * d(i);
    }, result);
    return std::sqrt(result);
  } // end norm 

  Ptr<Vector<Real>> clone() const override {
    return makePtr<KokkosVector>(dimension());
  } // end clone

  Ptr<Vector<Real>> basis( int i ) const override {
    auto bp = makePtr<KokkosVector>(dimension());
    device_view_type bpd = bp->get_device_view();
    Kokkos::parallel_for( dimension(), KOKKOS_LAMBDA( const int j ) {
      bpd(j) = static_cast<Real>(i==j);
    });
    return bp;
  } // end basis

  void randomize( const Real l = 0.0, const Real u = 1.0 ) override {
    host_view_type h = get_host_view();
    for( int i=0; i<dimension(); ++i ) {
      auto x = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
      data_.h_view(i) = u*x + l*(1-x);
    }
    copy_host_to_device();
  } // end randomize
 
  /** 
  * Elementwise functions must be executed on the host
  */

  void applyUnary( const Elementwise::UnaryFunction<Real>& f ) override {
    host_view_type h = get_host_view();
    copy_device_to_host();
    for( int i=0; i<dimension(); ++i )  h(i) = f.apply(h(i));
    copy_host_to_device();
  }

  void applyBinary( const Elementwise::BinaryFunction<Real>& f, const Vector<Real>& x ) override {
    copy_device_to_host();
    const auto& xkv = cast(x);
    xkv.copy_device_to_host();
    host_view_type h  = get_host_view();
    host_view_type xh = cast(x).get_host_view();
    for( int i=0; i<dimension(); ++i )  h(i) = f.apply(xh(i),h(i));
    copy_host_to_device();
  }

  Real reduce( const Elementwise::ReductionOp<Real>& r ) const override {
    copy_device_to_host();
    host_view_type h = get_host_view();
    Real result = r.initialValue();
    for( int i=0; i<dimension(); ++i ) r.reduce(h(i),result);
    return result;
  } 

  void print( std::ostream& os ) const override {
    copy_device_to_host();
    for( int i=0; i<dimension(); ++i ) os << data_.h_view(i) << " ";
    os << std::endl;
  } // end print

  static const KokkosVector& cast( const Vector<Real>& x ) {
    return static_cast<const KokkosVector&>(x);
  }

  static KokkosVector& cast( Vector<Real>& x ) {
    return static_cast<KokkosVector&>(x);
  } 

private:

  mutable view_type data_;

}; // KokkosVector

} // namespace ROL



#endif  // ROL_KOKKOSVECTOR_HPP
