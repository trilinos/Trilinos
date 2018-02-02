
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

#include "XROL_Core.hpp"

namespace XROL {

namespace details {

template<class U, class Z> Vector_SimOpt;

template<class U, class Z>
struct IndexType<Vector_SimOpt<U,Z>> {
  using type = common_type_t<index_t<U>,index_t<Z>>;
};

template<class U, class Z>
struct ElementType<Vector_SimOpt<U,Z>> {
  using type = common_type_t<element_t<U>,element_t<Z>>;
};

template<class U, class Z>
struct DualType<Vector_SimOpt<U,Z>> {
  using type = Vector_SimOpt<dual_t<U>,dual_t<Z>>;
};

tempate<class U, class Z>
struct SelfDual<Vector_SimOpt<U,Z>> {
  static constexpr auto value = 
    SelfDual<U>::value and SelfDual<Z>::value;
};



template<class U, class Z>
class Vector_SimOpt : public Vector<U,Z> {
  using IndexT   = index_t<Vector_SimOpt>;
  using ElementT = element_t<Vector_SimOpt>;
  using NormT    = norm_t<Vector_SimOpt>;
  using DualT    = dual_t<Vector_SimOpt>;

private:

  template<class X>
  using ptr_t = typename PointerType<X,is_same<X,dual_t<X>>>::type;

  ptr_t<U> sim_;
  ptr_t<Z> opt_;
  
  mutable unique_ptr<DualT> dual_vec_;

public:

  Vector_SimOpt( ptr_t<U> sim, ptr_t<Z> opt ) {
    if(!SelfDual<Vector_SimOpt>::value) {
      dual_vec_ = make_unique<DualT>(move(sim->dual().clone()),
                                     move(opt->dual().clone());
    }
  }

  void plus( const Vector_SimOpt& x ) {
    sim_->plus(x.get_1());
    opt_->plus(x.get_2());
  }

  void set( const Vector_SimOpt& x ) {
    sim_->set(x.get_1());
    opt_->set(x.get_2());
  }

  NormT dot( const Vector_SimOpt& x ) {
    return sim_->dot(x.get_1()) + opt_->dot(x.get_2());
  }

  NormT norm() const {
    return sqrt(sim_->dot(*sim_) + opt_->dot(*opt_));
  }

  unique_ptr<Vector<Vector_SimOpt>> clone() const {
    auto x = make_unique<Vector_SimOpt>(sim_->clone(),opt_->clone());
    return move(x);
  }

  void axpy( const ElementT alpha, const Vector_SimOpt& x ) {
    sim_->axpy(static_cast<element_t<U>>(alpha),x.get_1());
    opt_->axpy(static_cast<element_t<Z>>(alpha),x.get_2());
  }

  void fill( const ElementT alpha ) {
    sim_->fill(static_cast<element_t<U>>(alpha));
    opt_->fill(static_cast<element_t<Z>>(alpha));
  }

  void scale( const ElementT alpha ) {
    sim_->scale(static_cast<element_t<U>>(alpha));
    opt_->scale(static_cast<element_t<Z>>(alpha));
  }

  unique_ptr<Vector<Vector_SimOpt>> basis( IndexT i ) const {
    IndexT n1 = sim_->dimension();
    IndexT n2 = opt_->dimension();
  }

  IndexT dimension() const { 
    return sim_->dimension() + opt_->dimension();
  }

  enable_if_t<SelfDual<Vector_SimOpt>,const DualT&>
  dual() const {
    return *this;
  }

  enable_if_t<!SelfDual<Vector_SimOpt>,const DualT&>
  dual() const {
    return *dual_vec_;
  }

  void print( ostream &os, const string& delimiter = " " ) const {
    os << sim_->print(os) << delimiter;
    os << opt_->print(os) << delimiter;
  }

  template<class R>
  NormT reduce( R&& r ) const {
    norm_t<U> result1 = sim_->reduce(forward<R>(r));
    norm_t<Z> result2 = opt_->reduce(forward<R>(r));
    return r(static_cast<NormT>(result1), static_cast<NormT>(result2));
  }

  template<class F, class Vs...>
  void applyFunction( F&& f, const Vs&... vs ) {
    sim_->applyFunction(forward<F>(f), make_tuple(vs.get_1()...));
    opt_->applyFunction(forward<F>(f), make_tuple(vs.get_2()...));
  }

  template<class F, class R, class Vs...>
  void applyFunctionAndReduce( F&& f, R&& r, const Vs&... vs ) {
    norm_t<U> result1 = sum_->applyFunctionAndReduce(forward<F>(f),
                                                     forward<R>(r),
                                                     make_tuple(vs.get_1()...));
    norm_t<Z> result2 = opt_->applyFunctionAndReduce(forward<F>(f),
                                                     forward<R>(r),
                                                     make_tuple(vs.get_2()...));
    return r(static_cast<NormT>(result1), static_cast<NormT>(result2));
  }


  U& get_1() { return *sim_; }
  Z& get_2() { return *opt_;}
  const U& get_1() const { return *sim_; }
  const Z& get_2() const { return *opt_; }
  void set_1( const U& sim ) { sim_->set(sim); }
  void set_2( const Z& opt ) { opt_->set(opt); }

};



} // namespace XROL

