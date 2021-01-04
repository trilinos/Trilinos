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
#ifndef ROL_KOKKOSVECTOR_DECL_HPP
#define ROL_KOKKOSVECTOR_DECL_HPP

/** @ingroup la_group
    \class ROL::KokkosVector
    \brief Defines the ROL interface for Kokkos::View
*/

namespace ROL {

template<typename Real, typename Device>
class KokkosVector : public ROL::Vector<Real> {
public:

  using view_type = Kokkos::DualView<Real*,Device>;
  using host_view_type = typename view_type::t_host;
  using device_view_type = typename view_type::t_dev;

  static const KokkosVector& cast( const Vector<Real>& x ) {
    return static_cast<const KokkosVector&>(x);
  }

  static KokkosVector& cast( Vector<Real>& x ) {
    return static_cast<KokkosVector&>(x);
  } 

  KokkosVector( int n ) : data_("data",n) {} 
  virtual ~KokkosVector() {}

  //--------------------------------------------------------------------------------
  // Kokkos::DualView Interface

  inline void modify_device() { data_.modify_device(); }
  inline void modify_host() { data_.modify_host(); } 
  inline bool need_sync_device() { return data_.need_sync_device(); }
  inline bool need_sync_host() { return data_.need_sync_host(); } 
  inline void sync_device() { data_.sync_device(); }
  inline void sync_host() { data_.sync_host(); } 
  KOKKOS_INLINE_FUNCTION host_view_type   view_host()   { return data_.view_host();   }
  KOKKOS_INLINE_FUNCTION device_view_type view_device() { return data_.view_device(); }
  KOKKOS_INLINE_FUNCTION host_view_type   view_host()   const { return data_.view_host();   }
  KOKKOS_INLINE_FUNCTION device_view_type view_device() const { return data_.view_device(); }

  //--------------------------------------------------------------------------------
  // Set values via host access given an aritrary index function f:int->Real
  void set_from_index_function( std::function<Real(int)> f ) {
    host_view_type h = data_.view_host();
    for( int i=0; i<dimension(); ++i ) h(i) = f(i);
    data_.modify_host();
    data_.sync_device();
  }

  //--------------------------------------------------------------------------------
  // ROL::Vector Interface

  void applyBinary( const Elementwise::BinaryFunction<Real>&,  const Vector<Real>& ) override;
  void applyUnary( const Elementwise::UnaryFunction<Real>& ) override;
  void axpy( const Real, const Vector<Real>& ) override;
  Ptr<Vector<Real>> basis( int ) const override;
  Ptr<Vector<Real>> clone() const override;
  int dimension() const override;
  Real dot( const Vector<Real>& ) const override;
  Real norm() const override;
  void plus( const Vector<Real>& ) override;
  void print( std::ostream& ) const override;
  void randomize( const Real = 0.0, const Real = 1.0 ) override;
  Real reduce( const Elementwise::ReductionOp<Real>& ) const override;
  void scale( const Real alpha ) override;
  void set( const Vector<Real>& ) override;
  void setScalar( const Real ) override;
  void zero() override;
 
  //--------------------------------------------------------------------------------

private:

  mutable view_type data_;

}; // KokkosVector

} // namespace ROL

#endif  // ROL_KOKKOSVECTOR_DECL_HPP
