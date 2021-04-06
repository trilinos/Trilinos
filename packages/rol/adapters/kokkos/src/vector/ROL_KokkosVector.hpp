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
#include "Kokkos_Vector.hpp"

/** @ingroup la_group
    \class ROL::KokkosVector
    \brief Defines the ROL interface for Kokkos Vector
*/

namespace ROL {

template<typename Real, typename Device>
class KokkosVector : public ROL::Vector<Real> {
public:

  using container_type = Kokkos::vector<Real,Device>;
  using size_type      = typename container_type::size_type;

  KokkosVector( int n, Real fillValue = Real() ) : data_(n,fillValue) {} 

  Real& operator[] ( int i ) { return data_[i]; }
  const Real& operator[] ( int i ) const { return data_[i]; }

  int dimension() const override { return data_.size(); }

  void zero() override {
    Kokkos::parallel_for( data_.size(), KOKKOS_LAMBDA( const int i ) {
      data_[i] = 0;
    });
  } // end zero

  void set( const Vector<Real>& x ) override {
    const auto& xk = cast(x);
    Kokkos::parallel_for( data_.size(), KOKKOS_LAMBDA( const int i ) {
      data_[i] = xk.data_[i];
    });
  } // end set

  void plus( const Vector<Real>& x ) override {
    const auto& xk = cast(x);
    Kokkos::parallel_for( data_.size(), KOKKOS_LAMBDA( const int i ) {
      data_[i] += xk.data_[i];
    });
  } // end plus

  void axpy( const Real alpha, const Vector<Real>& x ) override {
    const auto& xk = cast(x);
    Kokkos::parallel_for( data_.size(), KOKKOS_LAMBDA( const int i ) {
      data_[i] += alpha * xk.data_[i];
    });
  } // end axpy

  void setScalar( const Real alpha ) override {
    Kokkos::parallel_for( data_.size(), KOKKOS_LAMBDA( const int i ) {
      data_[i] = alpha;
    });
  } // end setScalar

  void scale( const Real alpha ) override {
    Kokkos::parallel_for( data_.size(), KOKKOS_LAMBDA( const int i ) {
      data_[i] *= alpha;
    });
  } // end scale

  Real dot( const Vector<Real>& x ) const override { 
    const auto& xk = cast(x);
    Real result = 0;
    Kokkos::parallel_reduce( data_.size(), KOKKOS_LAMBDA( const int i, Real& lresult ) {
      lresult += data_[i] * xk.data_[i];
    }, result);
    return result;
  } // end dot

  Real norm() const override { 
    Real result = 0;
    Kokkos::parallel_reduce( data_.size(), KOKKOS_LAMBDA( const int i, Real& lresult ) {
      lresult += data_[i] * data_[i];
    }, result);
    return std::sqrt(result);
  } // end norm 

  Ptr<Vector<Real>> clone() const override {
    Ptr<Vector<Real>> xp = makePtr<KokkosVector>(data_.size());
    return xp;
  } // end clone

  Ptr<Vector<Real>> basis( int i )  const override {
    auto bp = makePtr<KokkosVector>(data_.size(),0);
    (*bp)[i] = 1;
    return bp;
  } // end basis

  void randomize( const Real l = 0.0, const Real u = 1.0 ) override {
    for( size_type i=0; i<data_.size(); ++i ) {
      auto x = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
      data_[i] = u*x + l*(1-x);
    }
  } // end randomize

  void print( std::ostream& os ) const override {
    for( size_type i=0; i<data_.size(); ++i ) os << data_[i] << " ";
    os << std::endl;
  } // end print

  static const KokkosVector& cast( const Vector<Real>& x ) {
    return static_cast<const KokkosVector&>(x);
  }

  static KokkosVector& cast( Vector<Real>& x ) {
    return static_cast<KokkosVector&>(x);
  } 

private:

   container_type data_;

}; // KokkosVector

} // namespace ROL



#endif  // ROL_KOKKOSVECTOR_HPP
