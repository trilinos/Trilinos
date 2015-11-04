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

#ifndef ROL_CPLXSTDVECTOR_H
#define ROL_CPLXSTDVECTOR_H

#include <algorithm>
#include <complex>
#include <cstdlib>

#include "ROL_Vector.hpp"

/** \class ROL::CplxStdVector 
    \brief Variant on ROL::StdVector with complex entries. Perhaps this can
           be accomplished with partial template specialization of 
           ROL::StdVector instead?
*/
           

namespace ROL {

template<class Real, template<class> class complex>
class CplxStdVector : public Vector<Real> { 

  typedef std::vector<complex<Real> > vector;
  typedef Vector<Real>                V; 
  typedef typename vector::size_type  uint;

private:
  
  Teuchos::RCP<vector> std_vec_;

public:

  CplxStdVector(const Teuchos::RCP<vector> &std_vec) : std_vec_(std_vec) {}  
 
  void set( const V &x ) {
    const CplxStdVector &ex = Teuchos::dyn_cast<const CplxStdVector>(x);
    const vector &xval = *ex.getVector();
    std::copy(xval.begin(),xval.end(),std_vec_->begin()); 
  }

  void plus( const V &x ) {
    const CplxStdVector &ex = Teuchos::dyn_cast<const CplxStdVector>(x);
    const vector& xval = *ex.getVector();
    uint dimension  = std_vec_->size();
    for (uint i=0; i<dimension; i++) {
      (*std_vec_)[i] += xval[i];
    }
  }

  void axpy( const Real alpha, const V &x ) {
    const CplxStdVector &ex = Teuchos::dyn_cast<const CplxStdVector>(x);
    const vector& xval = *ex.getVector();
    uint dimension  = std_vec_->size();
    for (uint i=0; i<dimension; i++) {
      (*std_vec_)[i] += alpha*xval[i];
    }
  }

  void scale( const Real alpha ) {
    uint dimension = std_vec_->size();
    for (uint i=0; i<dimension; i++) {
      (*std_vec_)[i] *= alpha;
    }
  }
 
  complex<Real> dot( const V &x ) const {
    const CplxStdVector & ex = Teuchos::dyn_cast<const CplxStdVector>(x);
    const vector& xval = *ex.getVector();
    uint dimension  = std_vec_->size();
    complex<Real> val = 0;
    for (uint i=0; i<dimension; i++) {
      val += (*std_vec_)[i]*std::conj(xval[i]);
    }
    return val;
  }  

  Real norm() const {
    Real val = 0;
    val = std::sqrt(std::abs(dot(*this))); 
  } 

  Teuchos::RCP<V> clone() const {
    return Teuchos::rcp( new StdVector( Teuchos::rcp(new std::vector<Element>(std_vec_->size())) ));
  }

  Teuchos::RCP<const vector> getVector() const {
    return std_vec_;
  }

  Teuchos::RCP<vector> getVector() {
    return std_vec_;
  }

  Teuchos::RCP<V> basis( const int i ) const {
    Teuchos::RCP<CplxStdVector> e = Teuchos::rcp( new CplxStdVector( Teuchos::rcp(new vector(std_vec_->size(), 0.0)) ));
    (*e->getVector())[i] = 1.0;
    return e;
  }

  int dimension() const {
    return static_cast<int>(std_vec_->size());
  }

  Teuchos::RCP<V> dual() const {
    uint dimension = std_vec_->size(); 
    Teuchos::RCP<vector> x_rcp = Teuchos::rcp( new vector(dimension,0.0) );
    for( uint i=0; i<dimension; ++i) {
      (*x_rcp)[i] = std::conj( (*std_vec_)[i] );
    }

  }


  // Elementwise operations probably do not make sense in general here

};

} // namespace ROL

#endif // ROL_CPLXSTDVECTOR_H
