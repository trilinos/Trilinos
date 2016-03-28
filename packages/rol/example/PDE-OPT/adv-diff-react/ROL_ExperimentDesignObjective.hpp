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

#ifndef ROL_EXPERIMENTDESIGN_OBJECTIVE_H
#define ROL_EXPERIMENTDESIGN_OBJECTIVE_H

#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_Vector.hpp"

namespace ROL {

template <class Real>
class ExperimentDesignObjective : public Objective<Real> {

  class ReducedHessian : public LinearOperator<Real> {
  public:
    ReducedHessian(ExperimentDesignObjective<Real> & edObj) : edObj_(edObj) {}

    void apply( Vector<Real> &hv, const Vector<Real> &v, Real &tol ) const {
      Real mytol(1e-8);
      const Teuchos::RCP<Vector<Real> > Bv = edObj_.constraint_->clone();
      const Teuchos::RCP<Vector<Real> > AiBv = edObj_.state_->clone();
      const Teuchos::RCP<Vector<Real> > QAiBv = edObj_.observation_->clone();
      const Teuchos::RCP<Vector<Real> > WQAiBv = (edObj_.observation_->dual()).clone();
      const Teuchos::RCP<Vector<Real> > QtWQAiBv = (edObj_.state_->dual()).clone();
      const Teuchos::RCP<Vector<Real> > AitQtWQAiBv = (edObj_.constraint_->dual()).clone();

      const Teuchos::RCP<Vector<Real> > BtAitQtWQAiBv = (edObj_.control_->dual()).clone();

      edObj_.con_->applyJacobian_2(*Bv, v, *(edObj_.state_), *(edObj_.control_), mytol);
      edObj_.con_->applyInverseJacobian_1(*AiBv, *Bv, *(edObj_.state_), *(edObj_.control_), mytol);
      edObj_.applyObserveOp(*QAiBv, *AiBv);
      edObj_.applyWeightOp(*WQAiBv, *QAiBv, *w_);
      edObj_.applyAdjointObserveOp(*QtWQAiBv, *WQAiBv);
      edObj_.con_->applyInverseAdjointJacobian_1(*AitQtWQAiBv, *QtWQAiBv, *(edObj_.state_), *(edObj_.control_), mytol);
      edObj_.con_->applyAdjointJacobian_2(*BtAitQtWQAiBv, *AitQtWQAiBv, *(edObj_.state_), *(edObj_.control_), mytol);

      edObj_.obj_->hessVec_22(hv, v, *(edObj_.state_), *(edObj_.control_), mytol);
      hv.plus(*BtAitQtWQAiBv);
   
    }

    void setWeights(const Teuchos::RCP<const Vector<Real> > &w) {
      w_ = w;
    }

  private:
    ExperimentDesignObjective<Real> & edObj_;
    Teuchos::RCP<const Vector<Real> > w_;

  };

public:
  ExperimentDesignObjective(const Teuchos::RCP<Objective_SimOpt<Real> > &obj,
                            const Teuchos::RCP<EqualityConstraint_SimOpt<Real> > &con,
                            const Teuchos::RCP<const Vector<Real> > &state,
                            const Teuchos::RCP<const Vector<Real> > &adjoint,
                            const Teuchos::RCP<const Vector<Real> > &control,
                            const Teuchos::RCP<const Vector<Real> > &constraint,
                            const Teuchos::RCP<const Vector<Real> > &observation,
                            const std::vector<Teuchos::RCP<const Vector<Real> > > &training,
                            const Teuchos::RCP<const Vector<Real> > &ones) :
    obj_(obj), con_(con), state_(state), adjoint_(adjoint),
    control_(control), constraint_(constraint), observation_(observation),
    training_(training), ones_(ones), reducedHessian_(*this) {
    // Start ExperimentDesignObjective.
    //redobj_ = Teuchos::rcp(new ROL::Reduced_Objective_SimOpt<Real> (obj_, con_, state_, adjoint_));
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    const Teuchos::RCP<Vector<Real> > Mtrain = (control_->dual()).clone();
    const Teuchos::RCP<Vector<Real> > CinvMtrain = control_->clone();
    int numTraining = training_.size();
    Real sumB(0);
    for (int i=0; i<numTraining; ++i) {
      obj_->hessVec_22(*Mtrain, *(training_[i]), *state_, *control_, tol);
      applyInverseReducedHessian(*CinvMtrain, *Mtrain, x);
      sumB += CinvMtrain->dot(*CinvMtrain);
    }
    return sumB;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    g.scale(0.0);
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    hv.zero();
  }

private:
  const Teuchos::RCP<Objective_SimOpt<Real> > obj_;
  const Teuchos::RCP<EqualityConstraint_SimOpt<Real> > con_;
  const Teuchos::RCP<const Vector<Real> > state_;
  const Teuchos::RCP<const Vector<Real> > adjoint_;
  const Teuchos::RCP<const Vector<Real> > control_;
  const Teuchos::RCP<const Vector<Real> > constraint_;
  const Teuchos::RCP<const Vector<Real> > observation_;
  const std::vector<Teuchos::RCP<const Vector<Real> > > training_;
  const Teuchos::RCP<const Vector<Real> > ones_;
  const Teuchos::RCP<Reduced_Objective_SimOpt<Real> > redobj_;
  ReducedHessian reducedHessian_;

  /*void applyReducedHessian(Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &w) const {
    Real tol(1e-8);
    const Teuchos::RCP<Vector<Real> > Bv = constraint_->clone();
    const Teuchos::RCP<Vector<Real> > AiBv = state_->clone();
    const Teuchos::RCP<Vector<Real> > QAiBv = observation_->clone();
    const Teuchos::RCP<Vector<Real> > WQAiBv = (observation_->dual()).clone();
    const Teuchos::RCP<Vector<Real> > QtWQAiBv = (state_->dual()).clone();
    const Teuchos::RCP<Vector<Real> > AitQtWQAiBv = (constraint_->dual()).clone();

    const Teuchos::RCP<Vector<Real> > BtAitQtWQAiBv = (control_->dual()).clone();

    con_->applyJacobian_2(Bv, v, state_, control_, tol);
    con_->applyInverseJacobian_1(AiBv, Bv, state_, control_, tol);
    applyObserveOp(QAiBv, AiBv);
    applyWeightOp(WQAiBv, QAiBv, w);
    applyAdjointObserveOp(QtWQAiBv, WQAiBv);
    con_->applyInverseAdjointJacobian_1(AitQtWQAiBv, QtWQAiBv, state_, control_, tol);
    con_->applyAdjointJacobian_2(BtAitQtWQAiBv, AitQtWQAiBv, state_, control_, tol);

    obj_->hessVec_22(hv, v, state_, control_, tol);
    hv.plus(BtAitQtWQAiBv);
  }*/

  void applyInverseReducedHessian(Vector<Real> &ihv, const Vector<Real> &v, const Vector<Real> &w) const {
    Real abstol(1e-4);
    Real reltol(1e-8);
    int  maxiter(100);
    int iter(0);
    int flag(0);
    ConjugateGradients<Real> cg(abstol, reltol, maxiter, false);
    reducedHessian_.setWeights(Teuchos::rcpFromRef(w));
    cg.run(ihv, reducedHessian_, v, reducedHessian_, iter, flag);
  }

  void applyObserveOp(Vector<Real> &obsv, const Vector<Real> &v) const {
    obsv.set(v);
  }


  void applyAdjointObserveOp(Vector<Real> &aobsv, const Vector<Real> &v) const {
    aobsv.set(v);
  }


  void applyWeightOp(Vector<Real> &weightv, const Vector<Real> &v, const Vector<Real> &w) const {
    weightv.set(v.dual());
  }

  

}; // class ExperimentDesignObjective

} // namespace ROL

#endif
