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

#include "ROL_ConjugateGradients.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_Vector.hpp"

namespace ROL {

/** @ingroup func_group
    \class ROL::ExperimentDesignInterface
    \brief Provides the interface for optimal experimental design (OED),
           which plugs into ROL::ExperimentDesignObjective.

    ROL's interface class for optimal experimental design; it relies on
    the SimOpt functional interface and the Elementwise
    vector interface.  It consists of three virtual functions:

    \li #applyObserveOp -- apply observation operator;
    \li #applyAdjointObserveOp -- apply adjoint of the observation operator; and
    \li #applyWeightOp -- apply weight operator.

    ---
*/
template <class Real>
class ExperimentDesignInterface {
private:
  ROL::Ptr<Objective_SimOpt<Real> > obj_;            // objective function used for the conventional inverse problem 
  ROL::Ptr<Constraint_SimOpt<Real> > con_;   // constraint function used for the conventional inverse problems
  ROL::Ptr<Vector<Real> > state_;                    // state vector, used for cloning
  ROL::Ptr<Vector<Real> > stateDual_;                // state dual vector, used for cloning
  ROL::Ptr<Vector<Real> > control_;                  // control vector, used for cloning
  ROL::Ptr<Vector<Real> > controlDual_;              // control dual vector, used for cloning
  ROL::Ptr<Vector<Real> > constraint_;               // constraint vector, used for cloning
  ROL::Ptr<Vector<Real> > constraintDual_;           // constraint dual vector, used for cloning
  ROL::Ptr<Vector<Real> > observation_;              // observation vector, used for cloning
  ROL::Ptr<Vector<Real> > observationDual_;          // observation dual vector, used for cloning
  std::vector<ROL::Ptr<Vector<Real> > > randvecs_;   // a set of vectors of -1 and 1 entries occurring with probability 1/2
  std::vector<ROL::Ptr<Vector<Real> > > training_;   // training-set vectors used in OED

public:
  ExperimentDesignInterface(const ROL::Ptr<Objective_SimOpt<Real> > &obj,
                            const ROL::Ptr<Constraint_SimOpt<Real> > &con,
                            const ROL::Ptr<Vector<Real> > &state,
                            const ROL::Ptr<Vector<Real> > &stateDual,
                            const ROL::Ptr<Vector<Real> > &control,
                            const ROL::Ptr<Vector<Real> > &controlDual,
                            const ROL::Ptr<Vector<Real> > &constraint,
                            const ROL::Ptr<Vector<Real> > &constraintDual,
                            const ROL::Ptr<Vector<Real> > &observation,
                            const ROL::Ptr<Vector<Real> > &observationDual,
                            const std::vector<ROL::Ptr<Vector<Real> > > &randvecs,
                            const std::vector<ROL::Ptr<Vector<Real> > > &training) :
    obj_(obj), con_(con),
    state_(state), stateDual_(stateDual), 
    control_(control), controlDual_(controlDual),
    constraint_(constraint), constraintDual_(constraintDual),
    observation_(observation), observationDual_(observationDual), training_(training) {

    for (unsigned i=0; i < randvecs.size(); ++i) {
      randvecs_.push_back(randvecs[i]->clone());
      randvecs_[i]->set(*randvecs[i]);  // deep copy the random -1/1 vector
    }
  }

  virtual ~ExperimentDesignInterface() {}

  /** \brief Apply the observation operator \f$Q \in L(\mathcal{U}, \mathcal{O})\f$,
             to vector \f$v\f$.

             @param[out]      obsv is the result of applying the observation operator to @b v ; an observation-space vector
             @param[in]       v is a simulation-space vector

             On return, \f$\mathsf{obsv} = Qv\f$, where
             \f$v \in \mathcal{U}\f$, \f$\mathsf{obsv} \in \mathcal{O}\f$. \n\n
             The default implementation is to set @b obsv to @b v; note that this is only possible
             if the spaces \f$\mathcal{U}\f$ and \f$\mathcal{O}\f$ are the same.

             ---
  */
  virtual void applyObserveOp(Vector<Real> &obsv, const Vector<Real> &v) const {
    obsv.set(v);
  }

  /** \brief Apply the adjoint observation operator \f$Q^* \in L(\mathcal{O}^*, \mathcal{U}^*)\f$,
             to vector \f$v\f$.

             @param[out]      aobsv is the result of applying the observation operator to @b v ; a dual simulation-space vector
             @param[in]       v is a dual observation-space vector

             On return, \f$\mathsf{aobsv} = Q^*v\f$, where
             \f$v \in \mathcal{O}^*\f$, \f$\mathsf{aobsv} \in \mathcal{U}^*\f$. \n\n
             The default implementation is to set @b aobsv to @b v; note that this is only possible
             if the spaces \f$\mathcal{U}^*\f$ and \f$\mathcal{O}^*\f$ are the same.

             ---
  */
  virtual void applyAdjointObserveOp(Vector<Real> &aobsv, const Vector<Real> &v) const {
    aobsv.set(v);
  }

  /** \brief Apply the weight operator \f$W(w) \in L(\mathcal{O}, \mathcal{O}^*)\f$, evaluated at \f$w \in \mathcal{O}^*\f$,
             to vector \f$v\f$.

             @param[out]      weightv is the result of applying the weight operator at @b w to @b v ; a dual observation-space vector
             @param[in]       v is an observation-space vector

             On return, \f$\mathsf{weightv} = W(w)*v\f$, where
             \f$v \in \mathcal{O}\f$, \f$\mathsf{weightv} \in \mathcal{O}^*\f$. \n\n
             The default implementation assumes that the spaces \f$\mathcal{O}\f$ and \f$\mathcal{O}^*\f$ are the same.

             ---
  */
  virtual void applyWeightOp(Vector<Real> &weightv, const Vector<Real> &v, const Vector<Real> &w) const {
    weightv.set(v.dual());
    ROL::Ptr<ROL::Vector<Real> > wDual = weightv.clone();
    wDual->set(w.dual());
    weightv.applyBinary(Elementwise::Multiply<Real>(), *wDual);
  }

  // Access functions.
  ROL::Ptr<Vector<Real> > getStateVec() const { return state_; }
  ROL::Ptr<Vector<Real> > getStateDualVec() const { return stateDual_; }
  ROL::Ptr<Vector<Real> > getControlVec() const { return control_; }
  ROL::Ptr<Vector<Real> > getControlDualVec() const { return controlDual_; }
  ROL::Ptr<Vector<Real> > getObservationVec() const { return observation_; }
  ROL::Ptr<Vector<Real> > getObservationDualVec() const { return observationDual_; }
  ROL::Ptr<Vector<Real> > getConstraintVec() const { return constraint_; }
  ROL::Ptr<Vector<Real> > getConstraintDualVec() const { return constraintDual_; }
  std::vector<ROL::Ptr<Vector<Real> > > getRandVecs() const { return randvecs_; }
  std::vector<ROL::Ptr<Vector<Real> > > getTrainingVecs() const {return training_; }
  ROL::Ptr<Objective_SimOpt<Real> > getObjective() const { return obj_; }
  ROL::Ptr<Constraint_SimOpt<Real> > getConstraint() const { return con_; }

}; // class ExperimentDesignInterface



/** @ingroup func_group
    \class ROL::ExperimentDesignObjective
    \brief Implements objective functions used for optimal experimental design (OED),
           based on ROL::ExperimentDesignInterface.

    ROL's objective function used for optimal experimental design;
    relies on ROL::ExperimentDesignInterface and ROL::Vector::applyUnary, ROL::Vector::applyBinary and ROL::Vector::reduce functions.

    ---
*/
template <class Real>
class ExperimentDesignObjective : public Objective<Real> {
private:
  const ROL::Ptr<ExperimentDesignInterface<Real> > edi_;   // experiment design interface

  Real cgabstol_;   // CG absolute tolerance to solve reduced-Hessian subproblems
  Real cgreltol_;   // CG relative tolerance to solve reduced-Hessian subproblems
  int  cgmaxiter_;  // max number of CG iterations to solve reduced-Hessian subproblems

  Real sigma_;      // standard deviation of the noise in data
  Real beta_;       // sparsity regularization factor

  int  nrandvecs_;  // number of random vectors

public:
  /*
    Constructor.
  */
  ExperimentDesignObjective(const ROL::Ptr<ExperimentDesignInterface<Real> > &edi,
                            const Teuchos::RCP<Teuchos::ParameterList> &parlist) : edi_(edi) {
    // get problem parameters
    cgabstol_  = parlist->sublist("Problem").get("OED CG Absolute Tolerance", 1e10);
    cgreltol_  = parlist->sublist("Problem").get("OED CG Relative Tolerance", 1e-10);
    cgmaxiter_ = parlist->sublist("Problem").get("OED CG Iteration Limit", 1000);
    sigma_     = parlist->sublist("Problem").get("OED Noise Standard Deviation", 1e-2);
    beta_      = parlist->sublist("Problem").get("OED Sparsity Regularization", 1e-2);
    nrandvecs_ = (edi_->getRandVecs()).size();
  }

  /*
    Reduced Hessian for the OED problem.
    This is a utility class, called by ExperimentDesignObjective.
  */
  class ReducedHessian : public LinearOperator<Real> {
  private:
    const ROL::Ptr<ExperimentDesignInterface<Real> > edi_;  // experiment design interface
    const ROL::Ptr<const Vector<Real> > w_;                 // weight vector

  public:

    ReducedHessian(const ROL::Ptr<ExperimentDesignInterface<Real> > & edi, ROL::Ptr<const Vector<Real> > w) : edi_(edi), w_(w) {}

    void apply( Vector<Real> &hv, const Vector<Real> &v, Real &tol ) const {
      Real mytol(1e-8);
      ROL::Ptr<Vector<Real> > Bv = (edi_->getConstraintVec())->clone();
      ROL::Ptr<Vector<Real> > AiBv = (edi_->getStateVec())->clone();
      ROL::Ptr<Vector<Real> > QAiBv = (edi_->getObservationVec())->clone();
      ROL::Ptr<Vector<Real> > WQAiBv = (edi_->getObservationDualVec())->clone();
      ROL::Ptr<Vector<Real> > QtWQAiBv = (edi_->getStateDualVec())->clone();
      ROL::Ptr<Vector<Real> > AitQtWQAiBv = (edi_->getConstraintDualVec())->clone();
      ROL::Ptr<Vector<Real> > BtAitQtWQAiBv = (edi_->getControlDualVec())->clone();

      (edi_->getConstraint())->applyJacobian_2(*Bv, v, *(edi_->getStateVec()), *(edi_->getControlVec()), mytol);
      (edi_->getConstraint())->applyInverseJacobian_1(*AiBv, *Bv, *(edi_->getStateVec()), *(edi_->getControlVec()), mytol);
      edi_->applyObserveOp(*QAiBv, *AiBv);
      edi_->applyWeightOp(*WQAiBv, *QAiBv, *w_);
      edi_->applyAdjointObserveOp(*QtWQAiBv, *WQAiBv);
      (edi_->getConstraint())->applyInverseAdjointJacobian_1(*AitQtWQAiBv, *QtWQAiBv, *(edi_->getStateVec()), *(edi_->getControlVec()), mytol);
      (edi_->getConstraint())->applyAdjointJacobian_2(*BtAitQtWQAiBv, *AitQtWQAiBv, *(edi_->getStateVec()), *(edi_->getControlVec()), mytol);

      (edi_->getObjective())->hessVec_22(hv, v, *(edi_->getStateVec()), *(edi_->getControlVec()), mytol);
      hv.plus(*BtAitQtWQAiBv);
   
    }
  }; // class ReducedHessian


  Real value( const Vector<Real> &x, Real &tol ) {

    ROL::Ptr<Vector<Real> > Vx = (edi_->getControlDualVec())->clone();
    ROL::Ptr<Vector<Real> > CinvVx = (edi_->getControlVec())->clone();

    // Initialize objective value.
    Real val(0);

    // Trace estimation term.
    for (int i=0; i<nrandvecs_; ++i) {
      Vx->set( *((edi_->getRandVecs())[i]) );
      applyInverseReducedHessian(*CinvVx, *Vx, x);
      val += (sigma_)*Vx->dot(*CinvVx);
    }

    // Sparse regularization term.
    val += beta_*x.reduce(ROL::Elementwise::ReductionSum<Real>());
    return val;
  }


  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    ROL::Ptr<Vector<Real> > v = (edi_->getControlDualVec())->clone();
    ROL::Ptr<Vector<Real> > Civ = (edi_->getControlVec())->clone();
    ROL::Ptr<Vector<Real> > BCiv = (edi_->getConstraintVec())->clone();
    ROL::Ptr<Vector<Real> > AiBCiv = (edi_->getStateVec())->clone();
    ROL::Ptr<Vector<Real> > QAiBCiv = (edi_->getObservationVec())->clone();
    ROL::Ptr<Vector<Real> > gtmp = g.clone();
    ROL::Ptr<Vector<Real> > vecones = g.clone();

    Real mytol(1e-8);
    g.zero();

    // Trace estimation term.
    for (int i=0; i<nrandvecs_; ++i) {
      v->set( *((edi_->getRandVecs())[i]) );
      applyInverseReducedHessian(*Civ, *v, x);
      (edi_->getConstraint())->applyJacobian_2(*BCiv, *Civ, *(edi_->getStateVec()), *(edi_->getControlVec()), mytol);
      (edi_->getConstraint())->applyInverseJacobian_1(*AiBCiv, *BCiv, *(edi_->getStateVec()), *(edi_->getControlVec()), mytol);
      edi_->applyObserveOp(*QAiBCiv, *AiBCiv);
      gtmp->set(*QAiBCiv);
      gtmp->applyBinary(Elementwise::Multiply<Real>(), *QAiBCiv);
      g.axpy(-sigma_, *gtmp);
    }

    // Sparse regularization term.
    vecones->applyUnary(Elementwise::Fill<Real>(1.0));
    g.axpy(beta_, *vecones);
  }


  /*void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    hv.zero();
  }*/

private:

  void applyInverseReducedHessian(Vector<Real> &ihv, const Vector<Real> &v, const Vector<Real> &w) const {
    int iter(0);
    int flag(0);
    ConjugateGradients<Real> cg(cgabstol_, cgreltol_, cgmaxiter_, false);
    ReducedHessian reducedHessian(edi_, ROL::makePtrFromRef(w));
    cg.run(ihv, reducedHessian, v, reducedHessian, iter, flag);
    //std::cout << "iter = " << iter << std::endl;
    //std::cout << "flag = " << flag << std::endl;
    if (flag != 0) {
      std::cout << std::endl << "The inner CG loop in ExperimentDesignObjective is not converged!" << std::endl;
    }
  }

}; // class ExperimentDesignObjective

} // namespace ROL

#endif
