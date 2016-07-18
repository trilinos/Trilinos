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

#ifndef ROL_COMPOSITECONSTRAINT_H
#define ROL_COMPOSITECONSTRAINT_H

#include "ROL_PartitionedVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_InequalityConstraint.hpp"

namespace ROL {

/** @ingroup func_group
 *  \class ROL::CompositeConstraint
 *  \brief Has both inequality and equality constraints.
 *        Treat inequality constraint as equality with slack variable
 */

template<class Real>
class CompositeConstraint : public EqualityConstraint<Real> {
private:

  typedef Vector<Real>            V;
  typedef PartitionedVector<Real> PV;
  typedef typename PV::size_type  size_type;

  const static size_type OPT   = 0;
  const static size_type SLACK = 1;

  const static size_type INEQ  = 0;
  const static size_type EQUAL = 1;

  Teuchos::RCP<InequalityConstraint<Real> > incon_;
  Teuchos::RCP<EqualityConstraint<Real> >   eqcon_;

  bool hasEquality_;         // True if an equality constraint is present
  int  ncval_;               // Number of constraint evaluations


public:

  // Constructor with inequality and equality constraints
  CompositeConstraint( const Teuchos::RCP<InequalityConstraint<Real> > &incon,
                       const Teuchos::RCP<EqualityConstraint<Real> > &eqcon ) :
                       incon_(incon), eqcon_(eqcon),
                       hasEquality_(true), ncval_(0) { }

  // Constructor with inequality constraint only
  CompositeConstraint( const Teuchos::RCP<InequalityConstraint<Real> > &incon ) :
                       incon_(incon), eqcon_(Teuchos::null),
                       hasEquality_(false), ncval_(0) { }


  int getNumberConstraintEvaluations(void) {
    return ncval_;
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {

    const PV &xpv = Teuchos::dyn_cast<const PV>(x);

    Teuchos::RCP<const V> xo = xpv.get(OPT);
    Teuchos::RCP<const V> xs = xpv.get(SLACK);

    incon_->update(*xo,flag,iter);

    if( hasEquality_ ) {
      eqcon_->update(*xo,flag,iter);
    }

  }

  void value( Vector<Real> &c, const Vector<Real> &x, Real &tol ) {

    PV &cpv = Teuchos::dyn_cast<PV>(c);
    const PV &xpv = Teuchos::dyn_cast<const PV>(x);

    Teuchos::RCP<const V> xo = xpv.get(OPT);
    Teuchos::RCP<const V> xs = xpv.get(SLACK);

    Teuchos::RCP<V> ci = cpv.get(INEQ);
    Teuchos::RCP<V> ce;

    incon_->value(*ci, *xo, tol);
    ci->axpy(-1.0,*xs);

    if(hasEquality_) {
      ce = cpv.get(EQUAL);
      eqcon_->value(*ce, *xo, tol);
    }

    ++ncval_;

  }

  void applyJacobian( Vector<Real> &jv,
                      const Vector<Real> &v,
                      const Vector<Real> &x,
                      Real &tol ) {

    using Teuchos::RCP;  using Teuchos::dyn_cast;

    // Partition vectors and extract subvectors
    const PV &xpv = dyn_cast<const PV>(x);
    const PV &vpv = dyn_cast<const PV>(v);

    RCP<const V> xo = xpv.get(OPT);
    RCP<const V> xs = xpv.get(SLACK);

    RCP<const V> vo = vpv.get(OPT);
    RCP<const V> vs = vpv.get(SLACK);

    PV &jvpv = dyn_cast<PV>(jv);

    RCP<V> jvi = jvpv.get(INEQ);
    incon_->applyJacobian(*jvi, *vo, *xo, tol);
    jvi->axpy(-1.0,*vs);

    if(hasEquality_) {
      RCP<V> jve = jvpv.get(EQUAL);
      eqcon_->applyJacobian(*jve, *vo, *xo, tol);
    }

  }

  void applyAdjointJacobian( Vector<Real> &ajv,
                             const Vector<Real> &v,
                             const Vector<Real> &x,
                             Real &tol ) {

    using Teuchos::RCP;  using Teuchos::dyn_cast;

    // Partition vectors and extract subvectors
    const PV &xpv = dyn_cast<const PV>(x);
    PV &ajvpv = dyn_cast<PV>(ajv);

    RCP<const V> xo = xpv.get(OPT);
    RCP<const V> xs = xpv.get(SLACK);

    RCP<V> ajvo = ajvpv.get(OPT);
    RCP<V> ajvs = ajvpv.get(SLACK);

    const PV &vpv = dyn_cast<const PV>(v);

    RCP<const V> vi = vpv.get(INEQ);

    incon_->applyAdjointJacobian(*ajvo,*vi,*xo,tol);

    ajvs->set(*vi);
    ajvs->scale(-1.0);

    if(hasEquality_) {

      RCP<const V> ve = vpv.get(EQUAL);
      RCP<V> temp = ajvo->clone();
      eqcon_->applyAdjointJacobian(*temp,*ve,*xo,tol);
      ajvo->plus(*temp);

    }

  }

  void applyAdjointHessian( Vector<Real> &ahuv,
                            const Vector<Real> &u,
                            const Vector<Real> &v,
                            const Vector<Real> &x,
                            Real &tol ) {

    using Teuchos::RCP;  using Teuchos::dyn_cast;

    const PV &xpv = dyn_cast<const PV>(x);
    const PV &vpv = dyn_cast<const PV>(v);
    PV &ahuvpv = dyn_cast<PV>(ahuv);

    RCP<const V> xo = xpv.get(OPT);
    RCP<const V> xs = xpv.get(SLACK);

    RCP<const V> vo = vpv.get(OPT);

    RCP<V> ahuvo = ahuvpv.get(OPT);
    RCP<V> ahuvs = ahuvpv.get(SLACK);

    RCP<V> temp = ahuvo->clone();

    const PV &upv = dyn_cast<const PV>(u);

    RCP<const V> ui = upv.get(INEQ);

    incon_->applyAdjointHessian(*ahuvo,*ui,*vo,*xo,tol);
    ahuvs->zero();

    if(hasEquality_) {
      RCP<const V> ue   = upv.get(EQUAL);
      eqcon_->applyAdjointHessian(*temp,*ue,*vo,*xo,tol);
      ahuvo->plus(*temp);
    }

  }

}; // class CompositeConstraint

} // namespace ROL

#endif
