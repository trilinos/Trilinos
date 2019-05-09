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

#include "ROL_StdVector.hpp"
#include <complex>

#ifndef ROL_COMPLEXSTDVECTOR_H
#define ROL_COMPLEXSTDVECTOR_H


// Specialize the StdVector class for complex-valued elements
namespace ROL {

template<class Real>
class StdVector<Real, std::complex<Real> > : public Vector<Real> {

  typedef std::complex<Real>         Complex;
  typedef std::vector<Complex>       vector;
  typedef Vector<Real>               V;
  typedef StdVector                  SV;
  typedef typename vector::size_type uint;

private:

  ROL::Ptr<vector> std_vec_;

public:
  StdVector( const ROL::Ptr<vector> &std_vec) : std_vec_(std_vec) {}

  void set( const V &x ) {
    const SV &ex = dynamic_cast<const SV&>(x);
    const vector& xval = *ex.getVector();
    std::copy(xval.begin(),xval.end(),std_vec_->begin());  
  }

  void plus( const V& x ) {
    const SV &ex = dynamic_cast<const SV&>(x);
    const vector& xval = *ex.getVector();
    uint dim = std_vec_->size();
    for( uint i=0; i<dim; ++i ) {
      (*std_vec_)[i] = (*std_vec_)[i] + xval[i];
    }
  }  

  void axpy( const Real alpha, const V& x ) {
    const SV &ex = dynamic_cast<const SV&>(x);
    const vector& xval = *ex.getVector();
    uint dim = std_vec_->size();
    for( uint i=0; i<dim; ++i ) {
      (*std_vec_)[i] = (*std_vec_)[i] + alpha*xval[i];
    }
  }  

  void scale( const Real alpha ) {
    uint dim = std_vec_->size();
    for( uint i=0; i<dim; ++i ) {
      (*std_vec_)[i] = alpha*(*std_vec_)[i]; 
    }
  }

  Real dot( const V &x ) const {
  
    const SV &ex = dynamic_cast<const SV&>(x);
    const vector& xval = *ex.getVector();
    uint dim = std_vec_->size();
    Real val = 0;
    for(uint i=0; i<dim; ++i) {
      val += rol_cast<Complex,Real>((*std_vec_)[i]*std::conj(xval[i]));
    }
    return val;
  }
 
  Real norm() const {
    Real val = 0;
    uint dim = std_vec_->size();
    for(uint i=0; i<dim; ++i) {
      val += std::norm((*std_vec_)[i]);
    }
    return std::sqrt(val);
  }

  virtual ROL::Ptr<Vector<Real> > clone() const {
    return ROL::makePtr<StdVector>( ROL::makePtr<vector>(std_vec_->size()));
  }

  ROL::Ptr<const vector> getVector() const {
    return std_vec_;
  }

  ROL::Ptr<vector> getVector() {
    return std_vec_;
  }

  ROL::Ptr<Vector<Real> > basis( const int i ) const {

    ROL::Ptr<StdVector> e = ROL::makePtr<StdVector>( ROL::makePtr<vector>(std_vec_->size(), 0.0));
    (*e->getVector())[i] = 1.0;
    return e;
  }

  int dimension() const {
    return static_cast<int>(std_vec_->size());
  }

  void applyUnary( const Elementwise::UnaryFunction<Complex> &f ) {
    uint dim  = std_vec_->size();
    for(uint i=0; i<dim; ++i) {
      (*std_vec_)[i] = f.apply((*std_vec_)[i]);
    }
  }

  void applyBinary( const Elementwise::BinaryFunction<Complex> &f, const V &x ) {

    const SV & ex = dynamic_cast<const SV&>(x);
    const vector& xval = *ex.getVector();
    uint dim  = std_vec_->size();
    for (uint i=0; i<dim; i++) {
      (*std_vec_)[i] = f.apply((*std_vec_)[i],xval[i]);
    }

  }

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const {
    Real result = r.initialValue();
    throw Exception::NotImplemented(">>> StdVector<Real,complex<Real>>::reduce(r) is not implemented!");
    return result;
  }  
}; 

}

#endif // ROL_COMPLEXSTDVECTOR_H


