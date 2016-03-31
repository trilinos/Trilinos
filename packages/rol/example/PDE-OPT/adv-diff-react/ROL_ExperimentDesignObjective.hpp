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
class ExperimentDesignObjective;

template <class Real>
class ReducedHessian : public LinearOperator<Real> {
private:
  const ExperimentDesignObjective<Real> edo_;
  Teuchos::RCP<const Vector<Real> > w_;

public:

  ReducedHessian(const ExperimentDesignObjective<Real> & edo, Teuchos::RCP<const Vector<Real> > w) : edo_(edo), w_(w) {}

  void apply( Vector<Real> &hv, const Vector<Real> &v, Real &tol ) const {
    Real mytol(1e-8);
    Teuchos::RCP<Vector<Real> > Bv = edo_.getConstraintVec()->clone();
    Teuchos::RCP<Vector<Real> > AiBv = edo_.getStateVec()->clone();
    Teuchos::RCP<Vector<Real> > QAiBv = edo_.getObservationVec()->clone();
    Teuchos::RCP<Vector<Real> > WQAiBv = (edo_.getObservationVec()->dual()).clone();
    Teuchos::RCP<Vector<Real> > QtWQAiBv = (edo_.getStateVec()->dual()).clone();
    Teuchos::RCP<Vector<Real> > AitQtWQAiBv = (edo_.getConstraintVec()->dual()).clone();

    Teuchos::RCP<Vector<Real> > BtAitQtWQAiBv = (edo_.getControlVec()->dual()).clone();

    edo_.getConstraint()->applyJacobian_2(*Bv, v, *(edo_.getStateVec()), *(edo_.getControlVec()), mytol);
    edo_.getConstraint()->applyInverseJacobian_1(*AiBv, *Bv, *(edo_.getStateVec()), *(edo_.getControlVec()), mytol);
    edo_.applyObserveOp(*QAiBv, *AiBv);
    edo_.applyWeightOp(*WQAiBv, *QAiBv, *w_);
    edo_.applyAdjointObserveOp(*QtWQAiBv, *WQAiBv);
    edo_.getConstraint()->applyInverseAdjointJacobian_1(*AitQtWQAiBv, *QtWQAiBv, *(edo_.getStateVec()), *(edo_.getControlVec()), mytol);
    edo_.getConstraint()->applyAdjointJacobian_2(*BtAitQtWQAiBv, *AitQtWQAiBv, *(edo_.getStateVec()), *(edo_.getControlVec()), mytol);

    edo_.getObjective()->hessVec_22(hv, v, *(edo_.getStateVec()), *(edo_.getControlVec()), mytol);
    hv.plus(*BtAitQtWQAiBv);
   
  }

};

template <class Real>
class ExperimentDesignObjective : public Objective<Real> {
private:
  Teuchos::RCP<Objective_SimOpt<Real> > obj_;
  Teuchos::RCP<EqualityConstraint_SimOpt<Real> > con_;
  Teuchos::RCP<Vector<Real> > state_;
  Teuchos::RCP<Vector<Real> > adjoint_;
  Teuchos::RCP<Vector<Real> > control_;
  Teuchos::RCP<Vector<Real> > constraint_;
  Teuchos::RCP<Vector<Real> > observation_;
  std::vector<Teuchos::RCP<Vector<Real> > > training_;
  Teuchos::RCP<Vector<Real> > ones_;
  Teuchos::RCP<Reduced_Objective_SimOpt<Real> > redobj_;

  Real cgabstol_;
  Real cgreltol_;
  int  cgmaxiter_;

public:
  ExperimentDesignObjective(const Teuchos::RCP<Objective_SimOpt<Real> > &obj,
                            const Teuchos::RCP<EqualityConstraint_SimOpt<Real> > &con,
                            const Teuchos::RCP<Vector<Real> > &state,
                            const Teuchos::RCP<Vector<Real> > &adjoint,
                            const Teuchos::RCP<Vector<Real> > &control,
                            const Teuchos::RCP<Vector<Real> > &constraint,
                            const Teuchos::RCP<Vector<Real> > &observation,
                            const std::vector<Teuchos::RCP<Vector<Real> > > &training,
                            const Teuchos::RCP<Vector<Real> > &ones,
                            const Teuchos::RCP<Teuchos::ParameterList> &parlist) :
    obj_(obj), con_(con), state_(state), adjoint_(adjoint),
    control_(control), constraint_(constraint), observation_(observation),
    training_(training), ones_(ones) {
    //
    cgabstol_  = parlist->sublist("Problem").get("OED CG Absolute Tolerance", 1e10);
    cgreltol_  = parlist->sublist("Problem").get("OED CG Relative Tolerance", 1e-4);
    cgmaxiter_ = parlist->sublist("Problem").get("OED CG Iteration Limit", 1000);
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    Teuchos::RCP<Vector<Real> > Mtrain = (control_->dual()).clone();
    Teuchos::RCP<Vector<Real> > CinvMtrain = control_->clone();
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

  virtual void applyObserveOp(Vector<Real> &obsv, const Vector<Real> &v) const {
    obsv.set(v);
  }


  virtual void applyAdjointObserveOp(Vector<Real> &aobsv, const Vector<Real> &v) const {
    aobsv.set(v);
  }


  virtual void applyWeightOp(Vector<Real> &weightv, const Vector<Real> &v, const Vector<Real> &w) const {
    weightv.set(v.dual());
  }

  Teuchos::RCP<Vector<Real> > getStateVec() const { return state_; }
  Teuchos::RCP<Vector<Real> > getAdjointVec() const { return adjoint_; }
  Teuchos::RCP<Vector<Real> > getControlVec() const { return control_; }
  Teuchos::RCP<Vector<Real> > getObservationVec() const { return observation_; }
  Teuchos::RCP<Vector<Real> > getConstraintVec() const { return constraint_; }
  Teuchos::RCP<Objective_SimOpt<Real> > getObjective() const { return obj_; }
  Teuchos::RCP<EqualityConstraint_SimOpt<Real> > getConstraint() const { return con_; }

private:

  void applyInverseReducedHessian(Vector<Real> &ihv, const Vector<Real> &v, const Vector<Real> &w) const {
    int iter(0);
    int flag(0);
    ConjugateGradients<Real> cg(cgabstol_, cgreltol_, cgmaxiter_, false);
    ReducedHessian<Real> reducedHessian(*this, Teuchos::rcpFromRef(w));
    cg.run(ihv, reducedHessian, v, reducedHessian, iter, flag);
std::cout << "iter = " << iter << std::endl;
std::cout << "flag = " << flag << std::endl;
  }

}; // class ExperimentDesignObjective

} // namespace ROL

#endif
