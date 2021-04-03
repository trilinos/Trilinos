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

#ifndef ROL_SLACKLESSCONSTRAINT_DEF_HPP
#define ROL_SLACKLESSCONSTRAINT_DEF_HPP

namespace ROL {

template<typename Real> 
SlacklessConstraint<Real>::SlacklessConstraint( const Ptr<Constraint<Real>> &con ) : con_(con) {}
 
template<typename Real> 
void SlacklessConstraint<Real>::update( const Vector<Real> &x, UpdateType type, int iter ) {
  con_->update( *getOpt(x), type, iter );
}
 
template<typename Real> 
void SlacklessConstraint<Real>::update( const Vector<Real> &x, bool flag, int iter ) {
  con_->update( *getOpt(x), flag, iter );
}

template<typename Real> 
void SlacklessConstraint<Real>::value(Vector<Real> &c, const Vector<Real> &x, Real &tol ) {
  con_->value( c, *getOpt(x), tol );
}

template<typename Real> 
void SlacklessConstraint<Real>::applyJacobian( Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  con_->applyJacobian( jv, *getOpt(v), *getOpt(x), tol );
}

template<typename Real> 
void SlacklessConstraint<Real>::applyAdjointJacobian( Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &dualv, Real &tol ) {
  zeroSlack(ajv);
  con_->applyAdjointJacobian( *getOpt(ajv), v, *getOpt(x), dualv, tol );
}

template<typename Real> 
void SlacklessConstraint<Real>::applyAdjointHessian( Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  zeroSlack(ahuv);
  con_->applyAdjointHessian( *getOpt(ahuv), u, *getOpt(v), *getOpt(x), tol );     
}

template<typename Real> 
void SlacklessConstraint<Real>::setParameter(const std::vector<Real> &param) {
  Constraint<Real>::setParameter(param);
  con_->setParameter(param);
}

template<typename Real> 
Ptr<Vector<Real>> SlacklessConstraint<Real>::getOpt( Vector<Real> &xs ) const {
  return dynamic_cast<PartitionedVector<Real>&>(xs).get(0);
}

template<typename Real> 
Ptr<const Vector<Real>> SlacklessConstraint<Real>::getOpt( const Vector<Real> &xs ) const {
  return dynamic_cast<const PartitionedVector<Real>&>(xs).get(0);
}

template<typename Real> 
void SlacklessConstraint<Real>::zeroSlack( Vector<Real> &x ) const {
  PartitionedVector<Real> &xpv
    = dynamic_cast<PartitionedVector<Real>&>(x);
  const int nvec = static_cast<int>(xpv.numVectors());
  for (int i = 1; i < nvec; ++i) {
    xpv.get(i)->zero();
  }
} 

} // namespace ROL

#endif // ROL__SLACKLESSCONSTRAINT_HPP

