//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperOperatorSplit_decl_hpp
#define Tempus_StepperOperatorSplit_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_Stepper.hpp"
#include "Tempus_StepperOperatorSplitAppAction.hpp"

namespace Tempus {

/** \brief OperatorSplit stepper loops through the Stepper list.
 *
 *  OperatorSplit stepper loops through the provided list of SubSteppers,
 *  and passes the SolutionHistory sequentially between them.  This is
 *  simply a first-order splitting.  It should be noted that specially
 *  constructed sequence of SubSteppers could obtain higher orders.
 *
 *  The OperatorSplit Stepper does not have any model, but the SubSteppers
 *  do.  The OperatorSplit Stepper does not have a solver either, but the
 *  SubSteppers may or may not have a solver depending if they are implicit
 *  or explicit.
 *
 *  Operator Split is only defined for one-step methods, so multi-step
 *  methods (e.g., BDF) should not be used with StepperOperatorSplit.
 *
 *  <b> Algorithm </b>
 *  The algorithm for operator-split stepper is
 *
 *  \f{center}{
 *    \parbox{5in}{
 *    \rule{5in}{0.4pt} \\
 *    {\bf Algorithm} Operator Split \\
 *    \rule{5in}{0.4pt} \vspace{-15pt}
 *    \begin{enumerate}
 *      \setlength{\itemsep}{0pt} \setlength{\parskip}{0pt}
 * \setlength{\parsep}{0pt} \item {\it appAction.execute(solutionHistory,
 * stepper, BEGIN\_STEP)} \item $x^\ast = x^n$ \hfill {\it * Initialize
 * operator-split solution.} \item {\bf for (each subStepper)} \item \quad  {\it
 * appAction.execute(solutionHistory, stepper, BEFORE\_STEPPER)} \item \quad
 * {\bf subStepper take time step.} \hfill {\it * Evolve solution, $x^\ast
 * \rightarrow x^{\ast\ast}$.} \item \quad  {\it
 * appAction.execute(solutionHistory, stepper, AFTER\_STEPPER)} \item \quad {\bf
 * if (subStep failed) then break.} \item \quad  {\bf Promote solution, $ x^\ast
 * \leftarrow x^{\ast\ast}$} \item {\bf end for} \item {\it
 * appAction.execute(solutionHistory, stepper, END\_STEP)} \end{enumerate}
 *    \vspace{-10pt} \rule{5in}{0.4pt}
 *    }
 *  \f}
 *
 *  Note that steppers in general can not use FSAL (useFSAL=true) with
 *  operator splitting as \f$\dot{x}_{n-1}\f$ will usually be modified
 *  by other operators.
 */
template <class Scalar>
class StepperOperatorSplit : virtual public Tempus::Stepper<Scalar> {
 public:
  /** \brief Default constructor.
   *
   *  Requires subsequent setModel(), setSolver() and initialize()
   *  calls before calling takeStep().
   */
  StepperOperatorSplit();

  /// Constructor
  StepperOperatorSplit(
      std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>>> appModels,
      std::vector<Teuchos::RCP<Stepper<Scalar>>> subStepperList, bool useFSAL,
      std::string ICConsistency, bool ICConsistencyCheck, int order,
      int orderMin, int orderMax,
      const Teuchos::RCP<StepperOperatorSplitAppAction<Scalar>>&
          stepperOSAppAction);

  /// \name Basic stepper methods
  //@{
  virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>>& appModel);

  virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>> getModel() const;

  virtual void setSolver(
      Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar>> solver);

  virtual Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar>> getSolver() const
  {
    return Teuchos::null;
  }

  virtual void setAppAction(
      Teuchos::RCP<StepperOperatorSplitAppAction<Scalar>> appAction);

  virtual Teuchos::RCP<StepperOperatorSplitAppAction<Scalar>> getAppAction()
      const
  {
    return stepperOSAppAction_;
  }

  virtual void setTempState(Teuchos::RCP<Tempus::SolutionState<Scalar>> state)
  {
    tempState_ = state;
  }

  /// Initialize during construction and after changing input parameters.
  virtual void initialize();

  /// Set the initial conditions and make them consistent.
  virtual void setInitialConditions(
      const Teuchos::RCP<SolutionHistory<Scalar>>& solutionHistory);

  /// Take the specified timestep, dt, and return true if successful.
  virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar>>& solutionHistory);

  /// Pass initial guess to Newton solver
  virtual void setInitialGuess(
      Teuchos::RCP<const Thyra::VectorBase<Scalar>> /* initial_guess */)
  {
  }

  /// Get a default (initial) StepperState
  virtual Teuchos::RCP<Tempus::StepperState<Scalar>> getDefaultStepperState();
  virtual Scalar getOrder() const { return order_; }
  virtual Scalar getOrderMin() const { return orderMin_; }
  virtual Scalar getOrderMax() const { return orderMax_; }

  virtual Scalar getInitTimeStep(
      const Teuchos::RCP<SolutionHistory<Scalar>>& /* solutionHistory */) const
  {
    return Scalar(1.0e+99);
  }
  virtual void setOrder(Scalar o) { order_ = o; }
  virtual void setOrderMin(Scalar o) { orderMin_ = o; }
  virtual void setOrderMax(Scalar o) { orderMax_ = o; }

  virtual bool isExplicit() const
  {
    bool isExplicit = false;
    typename std::vector<Teuchos::RCP<Stepper<Scalar>>>::const_iterator
        subStepperIter = subStepperList_.begin();
    for (; subStepperIter < subStepperList_.end(); subStepperIter++) {
      if ((*subStepperIter)->isExplicit()) isExplicit = true;
    }
    return isExplicit;
  }
  virtual bool isImplicit() const
  {
    bool isImplicit = false;
    typename std::vector<Teuchos::RCP<Stepper<Scalar>>>::const_iterator
        subStepperIter = subStepperList_.begin();
    for (; subStepperIter < subStepperList_.end(); subStepperIter++) {
      if ((*subStepperIter)->isImplicit()) isImplicit = true;
    }
    return isImplicit;
  }
  virtual bool isExplicitImplicit() const
  {
    return isExplicit() && isImplicit();
  }
  virtual bool isOneStepMethod() const
  {
    bool isOneStepMethod = true;
    typename std::vector<Teuchos::RCP<Stepper<Scalar>>>::const_iterator
        subStepperIter = subStepperList_.begin();
    for (; subStepperIter < subStepperList_.end(); subStepperIter++) {
      if (!(*subStepperIter)->isOneStepMethod()) isOneStepMethod = false;
    }
    return isOneStepMethod;
  }
  virtual bool isMultiStepMethod() const { return !isOneStepMethod(); }
  virtual void setUseFSAL(bool a)
  {
    this->useFSAL_       = a;
    this->isInitialized_ = false;
  }
  virtual OrderODE getOrderODE() const { return FIRST_ORDER_ODE; }
  //@}

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  void createSubSteppers(
      std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>>> appModels,
      Teuchos::RCP<Teuchos::ParameterList> pl);

  /// \name Overridden from Teuchos::Describable
  //@{
  virtual void describe(Teuchos::FancyOStream& out,
                        const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  virtual bool isValidSetup(Teuchos::FancyOStream& out) const;

  virtual std::vector<Teuchos::RCP<Stepper<Scalar>>> getStepperList() const
  {
    return subStepperList_;
  }
  virtual void setStepperList(std::vector<Teuchos::RCP<Stepper<Scalar>>> sl)
  {
    subStepperList_ = sl;
  }

  /** \brief Add Stepper to subStepper list.
   *  In most cases, subSteppers cannot use xDotOld (thus the default),
   *  but in some cases, the xDotOld can be used and save compute cycles.
   *  The user can set this when adding to the subStepper list.
   */
  virtual void addStepper(Teuchos::RCP<Stepper<Scalar>> stepper,
                          bool useFSAL = false);

  virtual void setSubStepperList(
      std::vector<Teuchos::RCP<Stepper<Scalar>>> subStepperList);

  virtual void clearSubStepperList() { subStepperList_.clear(); }

  virtual void setModels(
      std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>>> appModels);

 protected:
  Scalar order_;
  Scalar orderMin_;
  Scalar orderMax_;

  std::vector<Teuchos::RCP<Stepper<Scalar>>> subStepperList_;
  Teuchos::RCP<SolutionHistory<Scalar>> OpSpSolnHistory_;
  Teuchos::RCP<SolutionState<Scalar>> tempState_;
  Teuchos::RCP<StepperOperatorSplitAppAction<Scalar>> stepperOSAppAction_;
};

/// Nonmember constructor - ModelEvaluator and ParameterList
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<StepperOperatorSplit<Scalar>> createStepperOperatorSplit(
    std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>>> appModels,
    Teuchos::RCP<Teuchos::ParameterList> pl);

}  // namespace Tempus

#endif  // Tempus_StepperOperatorSplit_decl_hpp
