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
#ifndef ROL_KOKKOS_REDUCTIONOP_DECL_HPP
#define ROL_KOKKOS_REDUCTIONOP_DECL_HPP

namespace ROL {
namespace Elementwise {

/** \class KokkosReductionOp
 *  \brief Kokkos::Cuda compatible implementation of reduce operations 
 *
 *  ROL::Elementwise::ReductionOp uses virtual inheritance to perform 
 *  specific reduce type operations on ROL::Vector and derived types.
 */

template<typename Real, typename Device>
class KokkosReductionOp {
public:
 
  virtual ~KokkosReductionOp() = default;
  virtual Real reduce( const ::ROL::KokkosVector<Real,Device>& ) const = 0;

  static std::unique_ptr<KokkosReductionOp>
  create( const ::ROL::Elementwise::ReductionOp<Real>& );

private:

  struct Factory : public ReductionOp<Real>::Visitor {
    void visit( const ::ROL::Elementwise::EuclideanNormSquared<Real>& ) override;
    void visit( const ::ROL::Elementwise::ReductionAnd<Real>&         ) override;
    void visit( const ::ROL::Elementwise::ReductionMax<Real>&         ) override;
    void visit( const ::ROL::Elementwise::ReductionMin<Real>&         ) override;
    void visit( const ::ROL::Elementwise::ReductionSum<Real>&         ) override;
    std::unique_ptr<KokkosReductionOp> rop_;
  }; // Factory

}; // KokkosReductionOp


template<typename Real, typename Device>
class KokkosEuclideanNormSquared : public KokkosReductionOp<Real,Device> {
public:
  virtual ~KokkosEuclideanNormSquared() = default;
  Real reduce( const ::ROL::KokkosVector<Real,Device>& ) const override;
};


//template<typename Real, typename Device>
//class KokkosReductionAnd : public KokkosReductionOp<Real,Device> {
//public:
//  virtual ~KokkosReductionAnd() = default;
//  Real reduce( const ::ROL::KokkosVector<Real,Device>& ) const override;
//
//  struct LAndFunctor {
//    Kokkos::View<const Real*,Device> values;
//    KOKKOS_INLINE_FUNCTION
//    void operator()( const int& i, Real& value ) const { 
//      value = value && values(i);
//    }
//  };
//};


template<typename Real, typename Device>
class KokkosReductionMax : public KokkosReductionOp<Real,Device> {
public:
  virtual ~KokkosReductionMax() = default;
  Real reduce( const ::ROL::KokkosVector<Real,Device>& ) const override;

  struct MaxFunctor {
    Kokkos::View<const Real*,Device> values;
    KOKKOS_INLINE_FUNCTION
    void operator()( const int& i, Real& value ) const { 
      if( values(i) > value ) value = values(i);
    }
  };
};

template<typename Real, typename Device>
class KokkosReductionMin : public KokkosReductionOp<Real,Device> {
public:
  virtual ~KokkosReductionMin() = default;
  Real reduce( const ::ROL::KokkosVector<Real,Device>& ) const override;

  struct MinFunctor {
    Kokkos::View<const Real*,Device> values;
    KOKKOS_INLINE_FUNCTION
    void operator()( const int& i, Real& value ) const { 
      if( values(i) < value ) value = values(i);
    }
  };
};

template<typename Real, typename Device>
class KokkosReductionSum : public KokkosReductionOp<Real,Device> {
public:
  virtual ~KokkosReductionSum() = default;
  Real reduce( const ::ROL::KokkosVector<Real,Device>& ) const override;
};

} // namespace Elementwise
} // namespace ROL 

#endif // ROL_KOKKOS_REDUCTIONOP_DECL_HPP
