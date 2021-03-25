// @HEADER
// ************************************************************************
//
//               Rapid Optimi*thisation Library (ROL) Package
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
//              Denis Rid*thisal (drid*thisal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL2_VECTOR_DEF_H
#define ROL2_VECTOR_DEF_H

namespace ROL2 {

template<class Real>
void Vector<Real>::axpy( Real alpha, const Vector<Real>& x ) {
  auto ax = x.clone();
  ax->set(x);
  ax->scale(alpha);
  plus(*ax);
}

template<class Real>
void Vector<Real>::zero() { 
  scale( static_cast<Real>(0) );
}

template<class Real>
Ptr<Vector<Real>> Vector<Real>::basis( int i ) const {
  return nullPtr;
}

template<class Real>
void Vector<Real>::set( const Vector<Real>& x ) {
  zero();
  plus(x);
}

template<class Real>
Real Vector<Real>::apply(const Vector<Real>& x) const {
  return dot(x.dual());
}

template<class Real>
void Vector<Real>::applyUnary( const Elementwise::UnaryFunction<Real>& f ) {
  std::ignore = f;
  ROL_TEST_FOR_EXCEPTION( true, std::logic_error,
    "The method applyUnary wass called, but not implemented" << std::endl); 
}


template<class Real>
void Vector<Real>::applyBinary( const Elementwise::BinaryFunction<Real>& f, const Vector& x ) {
  std::ignore = f;
  std::ignore = x;
  ROL_TEST_FOR_EXCEPTION( true, std::logic_error,
    "The method applyBinary wass called, but not implemented" << std::endl); 
}

template<class Real>
Real Vector<Real>::reduce( const Elementwise::ReductionOp<Real>& r ) const {
  std::ignore = r;
  ROL_TEST_FOR_EXCEPTION( true, std::logic_error,
    "The method reduce was called, but not implemented" << std::endl); 
}

template<class Real>
void Vector<Real>::print( std::ostream& os ) const {
  os << "The method print was called, but not implemented" << std::endl;
}

template<class Real>
void Vector<Real>::setScalar( Real C ) {
  applyUnary(Elementwise::Fill<Real>(C));
}

template<class Real>
void Vector<Real>::randomize( Real l, Real u ) {
  Elementwise::UniformlyRandom<Real> ur(l,u);
  applyUnary(ur);
}

template<class Real>
int Vector<Real>::checkVector( const Vector<Real>& x,
                               const Vector<Real>& y,
                                     bool          printToStream,
                                     std::ostream& os ) const {

  auto table = Table( "Verification of Linear Algebra",
                      {"Check","Consistency Error"}, true, 8 );

  Real tol = std::sqrt(ROL_EPSILON<Real>);

  int error = 0;
  Real one  =  1.0;
  Real a    =  1.234;
  Real b    = -0.4321;

  std::vector<Real> vCheck;

  NullStream no_output;

  Ptr<std::ostream> pStream;
  if( printToStream ) pStream = makePtrFromRef(os);
  else                pStream = makePtrFromRef(no_output);

  // Save the original format of the pStream
  NullStream oldFormatState, headerFormatState;

  auto v_ptr    = clone();
  auto vtmp_ptr = clone();
  auto xtmp_ptr = x.clone();
  auto ytmp_ptr = y.clone();

  auto& v    = *v_ptr;
  auto& vtmp = *vtmp_ptr;
  auto& xtmp = *xtmp_ptr;
  auto& ytmp = *ytmp_ptr;

  auto update_check = [v_ptr,&vCheck,&table,tol]( std::string text ) {
    auto vnorm = v_ptr->norm();
    table.add_row(text,vnorm);
    vCheck.push_back(vnorm);
    return vnorm>tol;
  };

  auto reset_vxy = [&]() {
    v.set(*this);
    xtmp.set(x);
    ytmp.set(y);
  };

  headerFormatState.copyfmt(*pStream);

  // Commutativity of addition
  reset_vxy();
  v.plus(x);
  xtmp.plus(*this);
  v.axpy(-one,xtmp);
  error += update_check("Commutivity of addition");

  // Associativity of addition
  reset_vxy();
  ytmp.plus(x);
  v.plus(ytmp);
  xtmp.plus(*this);
  xtmp.plus(y);
  v.axpy(-one,xtmp);
  error += update_check("Associativity of addition");

  // Identity element of addition
  reset_vxy();
  v.zero();
  v.plus(x);
  v.axpy(-one,x);
  error += update_check("Identity element of addition");

  // Inverse elements of addition
  v.set(*this);
  v.scale(-one);
  v.plus(*this);
  error += update_check("Inverse element of addition");

  // Identity element of scalar multiplication
  v.set(*this);
  v.scale(one);
  v.axpy(-one,*this);
  error += update_check("Identity element of scalar multiplication");

  // Consistency of scalar multiplication with field multiplication
  v.set(*this);
  vtmp.set(*this);
  v.scale(a);
  v.scale(b);
  vtmp.scale(a*b);
  v.axpy(-one,vtmp);
  error += update_check("Consistency of scalar multiplication with field multiplication");

  // Distributivity of scalar multiplication with respect to field addition
  v.set(*this);
  vtmp.set(*this);
  v.scale(a+b);
  vtmp.scale(a);
  vtmp.axpy(b,*this);
  v.axpy(-one,vtmp);
  error += update_check("Distributivity of scalar multiplication with respect to field addition");

  // Distributivity of scalar multiplication with respect to vector addition
  reset_vxy();
  v.plus(x);
  v.scale(a);
  xtmp.scale(a);
  xtmp.axpy(a,*this);
  v.axpy(-one,xtmp);
  error += update_check("Distributivity of scalar multiplication with respect to vector addition");

  // Commutativity of dot (inner) product over the field of reals
  v.set(*this);
  v.scale(x.dot(y));
  vtmp.set(*this);
  vtmp.scale(y.dot(x));
  v.axpy(-one,vtmp);
  error += update_check("Commutativity of dot (inner) product over the field of reals");

  // Additivity of inner product
  v.set(x);
  v.plus(y);
  auto z_dot_v = this->dot(v);
  auto z_dot_x = this->dot(x);
  auto z_dot_y = this->dot(y);
  auto diff    = z_dot_v - z_dot_x - z_dot_y;
  v.scale(diff);
  error += update_check("Additivity of inner product");

  *pStream << table;

  // Restore format state of pStream used for the header info.
  pStream->copyfmt(oldFormatState);

  return error;
} // Vector::checkVector

} // namespace ROL2

#endif // ROL2_VECTOR_DEF_H
