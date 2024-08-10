// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_EXPERIMENTDESIGN_OBJECTIVE_H
#define ROL_EXPERIMENTDESIGN_OBJECTIVE_H

#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_Vector.hpp"

namespace ROL {

template <class Real>
class ExperimentDesignObjective;

/*
  Reduced Hessian for the OED problem.
  This is a utility class, called by ExperimentDesignObjective.
*/
template <class Real>
class ReducedHessian : public LinearOperator<Real> {
private:
  const ExperimentDesignObjective<Real> edo_;
  ROL::Ptr<const Vector<Real> > w_;

public:

  ReducedHessian(const ExperimentDesignObjective<Real> & edo, ROL::Ptr<const Vector<Real> > w) : edo_(edo), w_(w) {}

  void apply( Vector<Real> &hv, const Vector<Real> &v, Real &tol ) const {
    Real mytol(1e-8);
    ROL::Ptr<Vector<Real> > Bv = edo_.getConstraintVec()->clone();
    ROL::Ptr<Vector<Real> > AiBv = edo_.getStateVec()->clone();
    ROL::Ptr<Vector<Real> > QAiBv = edo_.getObservationVec()->clone();
    ROL::Ptr<Vector<Real> > WQAiBv = (edo_.getObservationVec()->dual()).clone();
    ROL::Ptr<Vector<Real> > QtWQAiBv = (edo_.getStateVec()->dual()).clone();
    ROL::Ptr<Vector<Real> > AitQtWQAiBv = (edo_.getConstraintVec()->dual()).clone();

    ROL::Ptr<Vector<Real> > BtAitQtWQAiBv = (edo_.getControlVec()->dual()).clone();

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


/** @ingroup func_group
    \class ROL::ExperimentDesignObjective
    \brief Provides the interface to evaluate objective functions
           used for optimal experimental design (OED).

    ROL's objective function used for optimal experimental design
    relies on the SimOpt functional interface and the Elementwise
    vector interface.  It adds four virtual functions:

    \li #applyObserveOp -- apply observation operator;
    \li #applyAdjointObserveOp -- apply adjoint of the observation operator; and
    \li #applyWeightOp -- apply weight operator.

    ---
*/
template <class Real>
class ExperimentDesignObjective : public Objective<Real> {
private:
  ROL::Ptr<Objective_SimOpt<Real> > obj_;            // objective function used for the conventional inverse problem 
  ROL::Ptr<Constraint_SimOpt<Real> > con_;   // constraint function used for the conventional inverse problems
  ROL::Ptr<Vector<Real> > state_;                    // state vector, used for cloning
  ROL::Ptr<Vector<Real> > adjoint_;                  // adjoint vector, used for cloning
  ROL::Ptr<Vector<Real> > control_;                  // control vector, used for cloning
  ROL::Ptr<Vector<Real> > constraint_;               // constraint vector, used for cloning
  ROL::Ptr<Vector<Real> > observation_;              // observation vector, used for cloning
  std::vector<ROL::Ptr<Vector<Real> > > training_;   // training-set vectors used in OED
  ROL::Ptr<Vector<Real> > rand01_;                   // a vector of 0 and 1 entries occurring with probability 1/2

  Real cgabstol_;   // CG absolute tolerance to solve reduced-Hessian subproblems
  Real cgreltol_;   // CG relative tolerance to solve reduced-Hessian subproblems
  int  cgmaxiter_;  // max number of CG iterations to solve reduced-Hessian subproblems

  Real sigma_;      // standard deviation of the noise in data
  Real beta_;       // sparsity regularization factor

public:
  ExperimentDesignObjective(const ROL::Ptr<Objective_SimOpt<Real> > &obj,
                            const ROL::Ptr<Constraint_SimOpt<Real> > &con,
                            const ROL::Ptr<Vector<Real> > &state,
                            const ROL::Ptr<Vector<Real> > &adjoint,
                            const ROL::Ptr<Vector<Real> > &control,
                            const ROL::Ptr<Vector<Real> > &constraint,
                            const ROL::Ptr<Vector<Real> > &observation,
                            const std::vector<ROL::Ptr<Vector<Real> > > &training,
                            const ROL::Ptr<Vector<Real> > &rand01,
                            const Teuchos::RCP<Teuchos::ParameterList> &parlist) :
    obj_(obj), con_(con), state_(state), adjoint_(adjoint),
    control_(control), constraint_(constraint), observation_(observation),
    training_(training), rand01_(rand01) {
    // get problem parameters
    cgabstol_  = parlist->sublist("Problem").get("OED CG Absolute Tolerance", 1e10);
    cgreltol_  = parlist->sublist("Problem").get("OED CG Relative Tolerance", 1e-4);
    cgmaxiter_ = parlist->sublist("Problem").get("OED CG Iteration Limit", 1000);
    sigma_     = parlist->sublist("Problem").get("OED Noise Standard Deviation", 1e-2);
    beta_      = parlist->sublist("Problem").get("OED Sparsity Regularization", 1e-2);
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    ROL::Ptr<Vector<Real> > Mtrain = (control_->dual()).clone();
    ROL::Ptr<Vector<Real> > CinvMtrain = control_->clone();
    ROL::Ptr<Vector<Real> > Vx = observation_->dual().clone();
    ROL::Ptr<Vector<Real> > QtVx = state_->dual().clone();
    ROL::Ptr<Vector<Real> > AitQtVx = constraint_->dual().clone();
    ROL::Ptr<Vector<Real> > BtAitQtVx = control_->dual().clone();
    ROL::Ptr<Vector<Real> > CinvBtAitQtVx = control_->clone();

    Real mytol(1e-8);
    // Initialize sum of bias, variance and sparse regularization.
    Real sumBVR(0);
    // Norm computation for the bias term.
    // Comment out, for now.
    /*
    int numTraining = static_cast<int>(training_.size());
    for (int i=0; i<numTraining; ++i) {
      obj_->hessVec_22(*Mtrain, *(training_[i]), *state_, *control_, tol);
      applyInverseReducedHessian(*CinvMtrain, *Mtrain, x);
      sumBVR += CinvMtrain->dot(*CinvMtrain);
    }
    */
    // Trace estimation for the variance term.
    Vx->set(*rand01_);
    Vx->applyBinary(Elementwise::Multiply<Real>(), x);
    applyAdjointObserveOp(*QtVx, *Vx);
    con_->applyInverseAdjointJacobian_1(*AitQtVx, *QtVx, *state_, *control_, mytol);
    con_->applyAdjointJacobian_2(*BtAitQtVx, *AitQtVx, *state_, *control_, mytol);
    applyInverseReducedHessian(*CinvBtAitQtVx, *BtAitQtVx, x);
    sumBVR += (sigma_*sigma_)*CinvBtAitQtVx->dot(*CinvBtAitQtVx);
    // Sparse regularization term.
    sumBVR += beta_*x.reduce(ROL::Elementwise::ReductionSum<Real>());
    return sumBVR;
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
    weightv.applyBinary(Elementwise::Multiply<Real>(), w);
  }

  // Access functions.
  ROL::Ptr<Vector<Real> > getStateVec() const { return state_; }
  ROL::Ptr<Vector<Real> > getAdjointVec() const { return adjoint_; }
  ROL::Ptr<Vector<Real> > getControlVec() const { return control_; }
  ROL::Ptr<Vector<Real> > getObservationVec() const { return observation_; }
  ROL::Ptr<Vector<Real> > getConstraintVec() const { return constraint_; }
  ROL::Ptr<Objective_SimOpt<Real> > getObjective() const { return obj_; }
  ROL::Ptr<Constraint_SimOpt<Real> > getConstraint() const { return con_; }

private:

  void applyInverseReducedHessian(Vector<Real> &ihv, const Vector<Real> &v, const Vector<Real> &w) const {
    int iter(0);
    int flag(0);
    ConjugateGradients<Real> cg(cgabstol_, cgreltol_, cgmaxiter_, false);
    ReducedHessian<Real> reducedHessian(*this, ROL::makePtrFromRef(w));
    cg.run(ihv, reducedHessian, v, reducedHessian, iter, flag);
std::cout << "iter = " << iter << std::endl;
std::cout << "flag = " << flag << std::endl;
  }

}; // class ExperimentDesignObjective

} // namespace ROL

#endif
