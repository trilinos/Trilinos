// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperOperatorSplit_decl_hpp
#define Tempus_StepperOperatorSplit_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_Stepper.hpp"
#include "Tempus_StepperOperatorSplitObserver.hpp"


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
 *  Note that steppers in general can not use FSAL (useFSAL=true) with
 *  operator splitting as \f$\dot{x}_{n-1}\f$ will usually be modified
 *  by other operators.
 */
template<class Scalar>
class StepperOperatorSplit : virtual public Tempus::Stepper<Scalar>
{
public:

  /** \brief Default constructor.
   *
   *  - Constructs with a default ParameterList.
   *  - Can reset ParameterList with setParameterList().
   *  - Requires subsequent addStepper() and initialize() calls before calling
   *    takeStep().
  */
  StepperOperatorSplit();

  /// Constructor
  StepperOperatorSplit(
    std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > appModels,
    Teuchos::RCP<Teuchos::ParameterList> pList);

  /// \name Basic stepper methods
  //@{
    virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel);
    virtual void setNonConstModel(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& appModel);
    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
      getModel();

    virtual void setSolver(std::string solverName);
    virtual void setSolver(
      Teuchos::RCP<Teuchos::ParameterList> solverPL=Teuchos::null);
    virtual void setSolver(
        Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver);
    virtual Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > getSolver() const
      { return Teuchos::null; }
    virtual void setObserver(
      Teuchos::RCP<StepperObserver<Scalar> > obs = Teuchos::null);
    virtual void setTempState(Teuchos::RCP<Tempus::SolutionState<Scalar>> state)
      { tempState_ = state; }

    /// Initialize during construction and after changing input parameters.
    virtual void initialize();

    /// Set the initial conditions and make them consistent.
    virtual void setInitialConditions (
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    /// Pass initial guess to Newton solver
    virtual void setInitialGuess(
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > /* initial_guess */){}

   virtual std::string getStepperType() const
     { return stepperPL_->get<std::string>("Stepper Type"); }

    /// Get a default (initial) StepperState
    virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState();
    virtual Scalar getOrder()    const
      {return stepperPL_->get<int>("Order");}
    virtual Scalar getOrderMin() const
      {return stepperPL_->get<int>("Minimum Order");}
    virtual Scalar getOrderMax() const
      {return stepperPL_->get<int>("Maximum Order");}
    virtual Scalar getInitTimeStep(
        const Teuchos::RCP<SolutionHistory<Scalar> >& /* solutionHistory */) const
      {return Scalar(1.0e+99);}
    virtual void setOrder   (Scalar ord)
      {stepperPL_->set<int>("Order", ord);}
    virtual void setOrderMin(Scalar ord)
      {stepperPL_->set<int>("Minimum Order", ord);}
    virtual void setOrderMax(Scalar ord)
      {stepperPL_->set<int>("Maximum Order", ord);}

    virtual bool isExplicit() const
    {
      bool isExplicit = false;
      typename std::vector<Teuchos::RCP<Stepper<Scalar> > >::const_iterator
        subStepperIter = subStepperList_.begin();
      for (; subStepperIter < subStepperList_.end(); subStepperIter++) {
        if ( (*subStepperIter)->isExplicit() ) isExplicit = true;
      }
      return isExplicit;
    }
    virtual bool isImplicit() const
    {
      bool isImplicit = false;
      typename std::vector<Teuchos::RCP<Stepper<Scalar> > >::const_iterator
        subStepperIter = subStepperList_.begin();
      for (; subStepperIter < subStepperList_.end(); subStepperIter++) {
        if ( (*subStepperIter)->isImplicit() ) isImplicit = true;
      }
      return isImplicit;
    }
    virtual bool isExplicitImplicit() const
      {return isExplicit() and isImplicit();}
    virtual bool isOneStepMethod()   const
    {
      bool isOneStepMethod = true;
      typename std::vector<Teuchos::RCP<Stepper<Scalar> > >::const_iterator
        subStepperIter = subStepperList_.begin();
      for (; subStepperIter < subStepperList_.end(); subStepperIter++) {
        if ( !(*subStepperIter)->isOneStepMethod() ) isOneStepMethod = false;
      }
      return isOneStepMethod;
    }
    virtual bool isMultiStepMethod() const {return !isOneStepMethod();}

    virtual OrderODE getOrderODE()   const {return FIRST_ORDER_ODE;}

    virtual void setUseFSAL(bool a) {stepperPL_->set<bool>("Use FSAL", a);}
    virtual bool getUseFSAL() const
      {return stepperPL_->get<bool>("Use FSAL", false);}

    virtual void setICConsistency(std::string s)
      {stepperPL_->set<std::string>("Initial Condition Consistency", s);}
    virtual std::string getICConsistency() const
      {return stepperPL_->get<std::string>("Initial Condition Consistency",
                                           "None");}

    virtual void setICConsistencyCheck(bool c)
      {stepperPL_->set<bool>("Initial Condition Consistency Check", c);}
    virtual bool getICConsistencyCheck() const
      {return stepperPL_->get<bool>("Initial Condition Consistency Check",
                                    false);}
   //@}

  /// \name ParameterList methods
  //@{
    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & pl);
    Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
    Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
    Teuchos::RCP<Teuchos::ParameterList> getDefaultParameters() const;
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual std::string description() const;
    virtual void describe(Teuchos::FancyOStream        & out,
                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  virtual std::vector<Teuchos::RCP<Stepper<Scalar> > > getStepperList() const
    { return subStepperList_; }
  virtual void setStepperList(std::vector<Teuchos::RCP<Stepper<Scalar> > > sl)
    { subStepperList_ = sl; }
  /** \brief Add Stepper to subStepper list.
   *  In most cases, subSteppers cannot use xDotOld (thus the default),
   *  but in some cases, the xDotOld can be used and save compute cycles.
   *  The user can set this when adding to the subStepper list.
   */
  virtual void addStepper(Teuchos::RCP<Stepper<Scalar> > stepper,
                          bool useFSAL = false)
  {
    stepper->setUseFSAL(useFSAL);
    subStepperList_.push_back(stepper);
  }
  virtual void clearStepperList() { subStepperList_.clear(); }
  /// Take models and ParameterList and create subSteppers
  virtual void createSubSteppers(
    std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > appModels);

protected:

  Teuchos::RCP<Teuchos::ParameterList>                stepperPL_;
  std::vector<Teuchos::RCP<Stepper<Scalar> > >        subStepperList_;
  Teuchos::RCP<SolutionHistory<Scalar> >              OpSpSolnHistory_;
  Teuchos::RCP<SolutionState<Scalar> >                tempState_;
  Teuchos::RCP<StepperOperatorSplitObserver<Scalar> > stepperOSObserver_;

};

} // namespace Tempus

#endif // Tempus_StepperOperatorSplit_decl_hpp
