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

#ifndef ROL2_STDVECTOR_DEF_H
#define ROL2_STDVECTOR_DEF_H

/** \class ROL2::StdVector
    \brief Provides the std::vector implementation of the ROL2::Vector interface.
*/

namespace ROL2 {

template<class Real>
StdVector<Real>::StdVector( const Ptr<std::vector<Real>>& std_vec ) : 
  std_vec_(std_vec) {
}

template<class Real>
StdVector<Real>::StdVector( int dim, Real val ) {
  std_vec_ = makePtr<std::vector<Real>>(dim,val);
}

template<class Real>
StdVector<Real>::
StdVector( std::initializer_list<Real> ilist ) : 
  std_vec_( makePtr<std::vector<Real>>(ilist) ) {}

template<class Real>
Real& StdVector<Real>::operator[] ( int i ) { return (*std_vec_)[i]; }

template<class Real>
const Real& StdVector<Real>::operator[] ( int i ) const { return (*std_vec_)[i]; }

template<class Real>
void StdVector<Real>::set( const Vector<Real> &x ) {
  ROL_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                          std::invalid_argument,
                          "Error: Vectors must have the same dimension." );
  const auto& xval = getData(x);
  std::copy(xval.begin(),xval.end(),std_vec_->begin());
}

template<class Real>
void StdVector<Real>::plus( const Vector<Real> &x ) {
  ROL_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                          std::invalid_argument,
                          "Error: Vectors must have the same dimension." );

  const auto& xval = getData(x);
  size_type dim  = std_vec_->size();
  for (size_type i=0; i<dim; i++) (*std_vec_)[i] += xval[i];
}

template<class Real>
void StdVector<Real>::axpy( Real alpha, const Vector<Real> &x ) {
  ROL_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                          std::invalid_argument,
                          "Error: Vectors must have the same dimension." );
  const auto& xval = getData(x);
  size_type dim  = std_vec_->size();
  for (size_type i=0; i<dim; i++) (*std_vec_)[i] += alpha*xval[i];
}

template<class Real>
void StdVector<Real>::scale( Real alpha ) {
  for( auto& e : *std_vec_ ) e *= alpha;
}

template<class Real>
Real StdVector<Real>::dot( const Vector<Real> &x ) const {
  ROL_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                          std::invalid_argument,
                          "Error: Vectors must have the same dimension." );
  const auto& xval = getData(x);
  return std::inner_product(std_vec_->begin(), std_vec_->end(), xval.begin(), Real(0));
}

template<class Real>
Real StdVector<Real>::norm() const {
  Real val = 0;
  val = std::sqrt( dot(*this) );
  return val;
}

template<class Real>
Ptr<Vector<Real>> StdVector<Real>::clone() const {
  return makePtr<StdVector>( makePtr<std::vector<Real>>(std_vec_->size(), 0) );
}

template<class Real>
Ptr<const std::vector<Real>> StdVector<Real>::getVector() const {
  return std_vec_;
}

template<class Real>
Ptr<std::vector<Real>> StdVector<Real>::getVector() {
  return std_vec_;
}

template<class Real>
Ptr<Vector<Real>> StdVector<Real>::basis( int i ) const {
  auto bptr = clone();
  auto& b = getData(*bptr);
  for( int j=0; j<dimension(); ++j ) b[j] = static_cast<Real>(j==i);
  return bptr;
}

template<class Real>
int StdVector<Real>::dimension() const {
 return static_cast<int>(std_vec_->size());
}

template<class Real>
void StdVector<Real>::applyUnary( const Elementwise::UnaryFunction<Real> &f ) {
 for( auto& e : *std_vec_ ) e = f.apply(e);
}

template<class Real>
void StdVector<Real>::applyBinary( const Elementwise::BinaryFunction<Real> &f, 
                                   const Vector<Real> &x ) {
  ROL_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                              std::invalid_argument,
                              "Error: Vectors must have the same dimension." );

  const auto& xval = getData(x);
  size_type dim  = std_vec_->size();
  for (size_type i=0; i<dim; i++) {
    (*std_vec_)[i] = f.apply((*std_vec_)[i],xval[i]);
  }
}

template<class Real>
Real StdVector<Real>::reduce( const Elementwise::ReductionOp<Real> &r ) const {
  Real result = r.initialValue();
  size_type dim  = std_vec_->size();
  for(size_type i=0; i<dim; ++i) {
    r.reduce((*std_vec_)[i],result);
  }
  return result;
}

template<class Real>
void StdVector<Real>::setScalar( Real C ) {
  size_type dim = std_vec_->size();
  std_vec_->assign(dim,C);
}

template<class Real>
void StdVector<Real>::randomize( Real l, Real u ) {
  Real a = (u-l);
  Real b = l;
  auto get_rand = [a,b]( Real& e ) { 
    auto x = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX); 
    e = a*x+b;
  };
  std::for_each( std_vec_->begin(), std_vec_->end(), get_rand );
}

template<class Real>
void StdVector<Real>::print( std::ostream &outStream ) const {
 for( auto e : *std_vec_ ) outStream << e << " ";
  outStream << std::endl;
}


} // namespace ROL2

#endif // ROL2_STDVECTOR_DEF_HPP
