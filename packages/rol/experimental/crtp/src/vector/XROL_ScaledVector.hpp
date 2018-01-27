
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


namespace XROL {

namespace details {

template<class V> class PrimalScaledVector;
template<class V> class DualScaledVector;

template<class V> 
struct IndexType<PrimalScaledVector> {
  using type = index_t<V>;
};

template<class V>
struct IndexType<DualScaledVector> {
  using type = index_t<V>;
};

template<class V>
struct ElementType<PrimalScaledVector> {
  using type = element_t<V>;
};

template<class V>
struct ElementType<DualScaledVector> {
  using type = element_t<V>;
};

template<class V>
struct DualType<PrimalScaledVector> {
  using type = DualScaledVector<V>;
};

template<class V>
struct DualType<DualScaledVector> {
  using type = PrimalScaledVector<V>;
};


template<class V>
class PrimalScaledVector : public Vector<PrimalScaledVector<V>> {

  using IndexT   = index_t<PrimalScaledVector>;
  using ElementT = element_t<PrimalScaledVector>;
  using NormT    = norm_t<PrimalScaledVector>;
  using DualT    = dual_t<PrimalScaledVector>;

private:

  unique_ptr<Vector<V>>                    vec_;
  shared_ptr<Vector<V>>                    scaling_vec_;
  unique_ptr<Vector<DualScaledVector<V>>>  dual_vec_; 
  bool isInitialized_;

public:

  template<class X> friend class DualScaledVector;

  PrimalScaledVector( unique_ptr<Vector<V>> vec,
                      shared_ptr<Vector<V>> scaling_vec ) :
    vec_(move(vec)), scaling_vec_(scaling_vec),
    isInitialized_(false) {}   

  void plus( const PrimalScaledVector& x ) {
    vec_->plus(*(x.vec_));
  }  

  void set( const PrimalScaledVector& x ) {
    vec_->set(*(x.vec_));
  }  

  NormT dot( const PrimalScaledVector& x ) const {
    return vec_->applyFunctionAndReduce(product,sum,*vec_,*scaling_vec_,*(x.vec_));
  }  
  
  NormT norm() const {
    return sqrt(this->dot(*this));
  } 

  unique_ptr<Vector<PrimalScaledVector>> 
  clone() const {
    auto psv = make_unique<PrimalScaledVector>(move(vec_->clone()),scaling_vec_);
    return move(psv);
  }

  void axpy( const ElementT alpha, const PrimalScaledVector& x ) {
    vec_->axpy(alpha,*(x.vec_));
  }

  void fill( const ElementT alpha ) {
    vec_->fill(alpha);
  }

  void scale( const ElementT alpha ) {
    vec_->scale(alpha);
  }

  unique_ptr<Vector<PrimalScaledVector>>
  basis( IndexT i ) const {
    auto psv = make_unique<PrimalScaledVector>(move(vec_->basis(i)),scaling_vec_);
    return move(psv);
  }

  IndexT dimension() const { return vec_->dimension(); }

  const DualT& dual() const {
    if( !isInitialized_ ) {
      dual_vec_ = make_unique<DualScaledVector>( move(vec_->clone()), scaling_vec_ );
      isInitialized_ = true;
    }
    dual_vec_->vec_->applyFunction(product,*vec_,*scaling_vec_);   
    return *dual_vec_;
  }

  void print( ostream &os, const string& delimiter=" " ) const {
    vec_->print(os,delimiter);
  }

  template<class F, class... Vs>
  void applyFunction( F&& f, const Vs&... vs ) {
    vec_->applyFunction(forward<F>(f), make_tuple(*(vs.vec_)...));
  }

  template<class R, class... Vs>
  NormT reduce( R&& r ) const {
    return vec_->reduce(forward<R>(r));
  }

  template<class F, class R, class... Vs>
  NormT applyFunctionAndReduce( F&& f, R&& r, const Vs&... vs ) {
    return vec_->applyFunctionAndReduce(forward<F>(f),forward<R>(r),make_tuple((*vs.vec_)...));  
  }
}; // class PrimalScaledVector





template<class V>
class DualScaledVector : public Vector<DualScaledVector<V>> {

  using IndexT   = index_t<DualScaledVector>;
  using ElementT = element_t<DualScaledVector>;
  using NormT    = norm_t<DualScaledVector>;
  using DualT    = dual_t<DualScaledVector>;

private:

  unique_ptr<Vector<V>>                      vec_;
  shared_ptr<Vector<V>>                      scaling_vec_;
  unique_ptr<Vector<PrimalScaledVector<V>>>  primal_vec_; 
  bool isInitialized_;

public:

  template<class X> friend class PrimalScaledVector;

  DualScaledVector( unique_ptr<Vector<V>> vec,
                    shared_ptr<Vector<V>> scaling_vec ) :
    vec_(move(vec)), scaling_vec_(scaling_vec),
    isInitialized_(false) {}   

  void plus( const DualScaledVector& x ) {
    vec_->plus(*(x.vec_));
  }  

  void set( const DualScaledVector& x ) {
    vec_->set(*(x.vec_));
  }  

  NormT dot( const DualScaledVector& x ) const {
    auto f = [](auto num1, auto num2, auto den) { return num1*num2/den; }; 
    return vec_->applyFunctionAndReduce(f,sum,*vec_,*scaling_vec_,*(x.vec_));
  }  
  
  NormT norm() const {
    return sqrt(this->dot(*this));
  } 

  unique_ptr<Vector<DualScaledVector>> 
  clone() const {
    auto psv = make_unique<DualScaledVector>(move(vec_->clone()),scaling_vec_);
    return move(psv);
  }

  void axpy( const ElementT alpha, const DualScaledVector& x ) {
    vec_->axpy(alpha,*(x.vec_));
  }

  void fill( const ElementT alpha ) {
    vec_->fill(alpha);
  }

  void scale( const ElementT alpha ) {
    vec_->scale(alpha);
  }

  unique_ptr<Vector<DualScaledVector>>
  basis( IndexT i ) const {
    auto psv = make_unique<DualScaledVector>(move(vec_->basis(i)),scaling_vec_);
    return move(psv);
  }

  IndexT dimension() const { return vec_->dimension(); }

  const DualT& dual() const {
    if( !isInitialized_ ) {
      primal_vec_ = make_unique<PrimalScaledVector>( move(vec_->clone()), scaling_vec_ );
      isInitialized_ = true;
    }
    primal_vec_->vec_->applyFunction(divide,*vec_,*scaling_vec_);   
    return *dual_vec_;
  }

  void print( ostream &os, const string& delimiter=" " ) const {
    vec_->print(os,delimiter);
  }

  template<class F, class... Vs>
  void applyFunction( F&& f, const Vs&... vs ) {
    vec_->applyFunction(forward<F>(f), make_tuple(*(vs.vec_)...));
  }

  template<class R, class... Vs>
  NormT reduce( R&& r ) const {
    return vec_->reduce(forward<R>(r));
  }

  template<class F, class R, class... Vs>
  NormT applyFunctionAndReduce( F&& f, R&& r, const Vs&... vs ) {
    return vec_->applyFunctionAndReduce(forward<F>(f),forward<R>(r),make_tuple((*vs.vec_)...));  
  }
}; // class DualScaledVector


} // namespace details

template<class V> using PrimalScaledVector = details::PrimalScaledVector<V>;
template<class V> using DualScaledVector = details::DualScaledVector<V>;


} // namespace XROL














