// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperSubcycling_decl_hpp
#define Tempus_StepperSubcycling_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperExplicit.hpp"
#include "Tempus_StepperSubcyclingObserver.hpp"
#include "Tempus_IntegratorBasic.hpp"


namespace Tempus {

/** \brief Subcycling time stepper.
 *
 *  This stepper wraps an IntegratorBasic object to perform the
 *  subcycling, thus it has all the capbilities of an IntegratorBasic
 *  with the following specializations and defaults:
 *   - Main Integrator operates independently from subcycling Stepper,
 *     and can be driven by another physics.
 *   - Need separate TimeStepControl for the subcycling so it can be
 *     driven by the physics (e.g., variable time stepping).
 *   - Finish the subcycling exactly on the full timestep.
 *   - No solution I/O within the subcycling.
 *   - No restart capability within subcycling, but still have restart
 *     capability from the full timestep.
 *   - Do not need to keep a solution history of the subcycling.
 */
template<class Scalar>
class StepperSubcycling : virtual public Tempus::Stepper<Scalar>
{
public:

  /** \brief Default constructor.
   *
   *  - Requires subsequent setModel() and initialize() call before
   *    calling takeStep().
  */
  StepperSubcycling();

  /// Constructor
  StepperSubcycling(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperSubcyclingObserver<Scalar> >& obs,
    const Teuchos::RCP<IntegratorBasic<Scalar> >& integrator,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck);

  /// \name Basic stepper methods
  //@{
    virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel);

    virtual void setNonConstModel(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& appModel);

    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
      getModel(){return scIntegrator_->getStepper()->getModel();}

    virtual void setObserver(
      Teuchos::RCP<StepperObserver<Scalar> > obs = Teuchos::null);

    virtual Teuchos::RCP<StepperObserver<Scalar> > getObserver() const;

    /// Initialize during construction and after changing input parameters.
    virtual void initialize();

    virtual Scalar getInitTimeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) const;

    /// Set the initial conditions, make them consistent, and set needed memory.
    virtual void setInitialConditions (
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    /// Pass initial guess to Newton solver (only relevant for implicit solvers)
    //  thus a no-op for explicit steppers.
    virtual void setInitialGuess(
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > initial_guess);

    virtual void setSolver(
      Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver = Teuchos::null);

    virtual Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > getSolver() const;

    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    /// Get a default (initial) StepperState
    virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState();

    virtual bool isExplicit()         const;
    virtual bool isImplicit()         const;
    virtual bool isExplicitImplicit() const;
    virtual bool isOneStepMethod()    const;
    virtual bool isMultiStepMethod()  const;

    virtual Scalar getOrder() const;
    virtual Scalar getOrderMin() const;
    virtual Scalar getOrderMax() const;

    virtual OrderODE getOrderODE()   const;
  //@}

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual void describe(Teuchos::FancyOStream        & out,
                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  /// \name Functions to set the subcycling stepper values.
  //@{
    virtual void setSubcyclingStepper(Teuchos::RCP<Stepper<Scalar> > stepper);
    virtual void setSubcyclingMinTimeStep(Scalar MinTimeStep);
    virtual void setSubcyclingInitTimeStep(Scalar InitTimeStep);
    virtual void setSubcyclingMaxTimeStep(Scalar MaxTimeStep);
    virtual void setSubcyclingStepType(std::string StepType);
    virtual void setSubcyclingMaxFailures(int MaxFailures);
    virtual void setSubcyclingMaxConsecFailures(int MaxConsecFailures);
    virtual void setSubcyclingScreenOutputIndexInterval(int i);
    virtual void setSubcyclingScreenOutputIndexList(std::string s);
    virtual void setSubcyclingTimeStepControlStrategy(
      Teuchos::RCP<TimeStepControlStrategy<Scalar> > tscs);
    virtual void setSubcyclingIntegratorObserver(
      Teuchos::RCP<IntegratorObserver<Scalar> > obs);
    virtual void setSubcyclingPrintDtChanges(bool printDtChanges);
  //@}

  /// \name Functions to get the subcycling stepper values.
  //@{
    virtual Teuchos::RCP<const Stepper<Scalar> > getSubcyclingStepper() const;
    virtual Scalar getSubcyclingMinTimeStep() const;
    virtual Scalar getSubcyclingInitTimeStep() const;
    virtual Scalar getSubcyclingMaxTimeStep() const;
    virtual std::string getSubcyclingStepType() const;
    virtual int getSubcyclingMaxFailures() const;
    virtual int getSubcyclingMaxConsecFailures() const;
    virtual int getSubcyclingScreenOutputIndexInterval() const;
    virtual std::string getSubcyclingScreenOutputIndexList() const;
    virtual Teuchos::RCP<TimeStepControlStrategy<Scalar> >
      getSubcyclingTimeStepControlStrategy() const;
    virtual Teuchos::RCP<IntegratorObserver<Scalar> >
      getSubcyclingIntegratorObserver() const;
    virtual bool getSubcyclingPrintDtChanges() const;
  //@}

protected:

  Teuchos::RCP<StepperSubcyclingObserver<Scalar> >  stepperSCObserver_;
  Teuchos::RCP<IntegratorBasic<Scalar> >            scIntegrator_;

};

} // namespace Tempus

#endif // Tempus_StepperSubcycling_decl_hpp
