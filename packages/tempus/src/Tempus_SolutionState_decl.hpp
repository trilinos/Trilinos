// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_SolutionState_decl_hpp
#define Tempus_SolutionState_decl_hpp

// Teuchos
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
// Thrya
#include "Thyra_VectorBase.hpp"
#include "Thyra_ModelEvaluator.hpp"
// Tempus
#include "Tempus_config.hpp"
#include "Tempus_SolutionStateMetaData.hpp"
#include "Tempus_StepperState.hpp"
#include "Tempus_PhysicsState.hpp"

namespace Tempus {

/** \brief Solution state for integrators and steppers.
 *  SolutionState contains the metadata for solutions and the solutions
 *  themselves.
 *
 *  For simple time integration, the SolutionState is sufficient for
 *  checkpointing, restart and undo operations (i.e., it is the Memento
 *  object).
 *
 *  For more complex time integration where the physics has additional
 *  state information or the time integrator is not a one-step method
 *  (i.e., can not accurately start from a single time step), this class
 *  can be inherited and the physics state or additional time-integration
 *  parameters can be managed.
 *
 *  SolutionStates can be interpolated to generate solutions at various
 *  times (see SolutionHistory).  However not all metadata or state
 *  information can be interpolated.  Thus interpolated solutions may not
 *  be suitable for checkpointing, restart and undo operations, but may
 *  be useful for adjoint sensitivities.
 *
 *  The solution vectors, \f$x\f$, \f$\dot{x}\f$, and \f$\ddot{x}\f$, in
 *  SolutionState can be null pointers.  This indicates that the
 *  application does not need them, so do not storage them.  This can be
 *  a huge savings when saving many states in the solution history.
 *  Some Steppers will need temporary memory to store time derivative(s)
 *  (\f$\dot{x}\f$, or \f$\ddot{x}\f$) for evaluation of the ODE/DAE
 *  (\f$f(x, \dot{x}, \ddot{x},t)\f$), but each individual Stepper will
 *  manage that.
 */
template<class Scalar>
class SolutionState :
  public Teuchos::Describable,
  public Teuchos::VerboseObject<Tempus::SolutionState<Scalar> >
{
public:

  /** \brief Default Constructor -- Not meant for immediate adding to SolutionHistory.
   *  This constructor does not set the solution vectors, x, xdot and xdotdot.
   *  which should be set via setX(), setXDot(), and/or setXDotDot() prior
   *  to being added to SolutionHistory.
   */
  SolutionState();

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  SolutionState(
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xdot    = Teuchos::null,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xddot   = Teuchos::null,
    const Teuchos::RCP<StepperState<Scalar> >& stepperState = Teuchos::null,
    const Teuchos::RCP<PhysicsState<Scalar> >& physicsState = Teuchos::null);

  SolutionState(
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& x,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xdot = Teuchos::null,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xddot= Teuchos::null,
    const Teuchos::RCP<const StepperState<Scalar> >& stepperSt = Teuchos::null,
    const Teuchos::RCP<const PhysicsState<Scalar> >& physicsSt = Teuchos::null);
#endif

  SolutionState(
    const Teuchos::RCP<SolutionStateMetaData<Scalar> > ssmd,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xdot,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xdotdot,
    const Teuchos::RCP<StepperState<Scalar> >& stepperState,
    const Teuchos::RCP<PhysicsState<Scalar> >& physicsState);

  SolutionState(
    const Teuchos::RCP<const SolutionStateMetaData<Scalar> > ssmd,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& x,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xdot,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xdotdot,
    const Teuchos::RCP<const StepperState<Scalar> >& stepperState,
    const Teuchos::RCP<const PhysicsState<Scalar> >& physicsState);

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  SolutionState(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    const Teuchos::RCP<StepperState<Scalar> >& stepperState = Teuchos::null,
    const Teuchos::RCP<PhysicsState<Scalar> >& physicsState = Teuchos::null);
#endif

  /// This is a shallow copy constructor, use clone for a deep copy constructor
  SolutionState(const SolutionState<Scalar>& ss);

  /// This is a deep copy constructor
  virtual Teuchos::RCP<SolutionState<Scalar> > clone() const;

  /// This is a deep copy
  virtual void copy(const Teuchos::RCP<const SolutionState<Scalar> >& ss);

  /// Deep copy solution data, but keep metaData untouched.
  virtual void copySolutionData(
    const Teuchos::RCP<const SolutionState<Scalar> >& s);

  /// Destructor
  virtual ~SolutionState() {}

  /// \name Get MetaData values
  //@{
    virtual Teuchos::RCP<const SolutionStateMetaData<Scalar> >
      getMetaData() const { return metaData_; }
    virtual Teuchos::RCP<SolutionStateMetaData<Scalar> > getMetaData()
      { TEUCHOS_ASSERT(metaData_nc_ != Teuchos::null);
        return metaData_nc_; }

    virtual Scalar getTime()             const {return metaData_->getTime();}
    virtual int    getIndex()            const {return metaData_->getIStep();}
    virtual Scalar getTimeStep()         const {return metaData_->getDt();}
    virtual Scalar getErrorAbs()         const {return metaData_->getErrorAbs();}
    virtual Scalar getErrorRel()         const {return metaData_->getErrorRel();}
    virtual int    getOrder()            const {return metaData_->getOrder();}
    virtual int    getNFailures()        const {return metaData_->getNFailures();}
    virtual int    getNRunningFailures() const {return metaData_->getNRunningFailures();}
    virtual int    getNConsecutiveFailures() const {return metaData_->getNConsecutiveFailures();}
    virtual Scalar getTolAbs()           const {return metaData_->getTolAbs();}
    virtual Scalar getTolRel()           const {return metaData_->getTolRel();}
    virtual Scalar getXNormL2()          const {return metaData_->getXNormL2();}
    virtual Scalar getDxNormL2Abs()      const {return metaData_->getDxNormL2Abs();}
    virtual Scalar getDxNormL2Rel()      const {return metaData_->getDxNormL2Rel();}
    virtual bool   getComputeNorms()     const {return metaData_->getComputeNorms();}
    virtual Status getSolutionStatus()   const {return metaData_->getSolutionStatus();}
    virtual bool   getOutput()           const {return metaData_->getOutput();}
    virtual bool   getOutputScreen()     const {return metaData_->getOutputScreen();}
    virtual bool   getIsSynced()         const {return metaData_->getIsSynced();}
    virtual bool   getIsInterpolated()   const {return metaData_->getIsInterpolated();}
    virtual bool   getAccuracy()         const {return metaData_->getAccuracy();}
  //@}

  /// \name Set MetaData values
  //@{
    virtual void setMetaData(
      Teuchos::RCP<const SolutionStateMetaData<Scalar> > md)
      { metaData_ = md; metaData_nc_ = Teuchos::null; }
    virtual void setMetaData(Teuchos::RCP<SolutionStateMetaData<Scalar> > md)
      { metaData_nc_ = md; metaData_ = metaData_nc_; }

    virtual void setTime(Scalar time)   { TEUCHOS_ASSERT(metaData_nc_ != Teuchos::null);
                                          metaData_nc_->setTime(time);}
    virtual void setIndex(Scalar index) { TEUCHOS_ASSERT(metaData_nc_ != Teuchos::null);
                                          metaData_nc_->setIStep(index);}
    virtual void setTimeStep(Scalar dt) { TEUCHOS_ASSERT(metaData_nc_ != Teuchos::null);
                                          metaData_nc_->setDt(dt);}
    virtual void setErrorAbs(Scalar errorAbs){ TEUCHOS_ASSERT(metaData_nc_ != Teuchos::null);
                                               metaData_nc_->setErrorAbs(errorAbs);}
    virtual void setErrorRel(Scalar errorRel){ TEUCHOS_ASSERT(metaData_nc_ != Teuchos::null);
                                               metaData_nc_->setErrorRel(errorRel);}
    virtual void setOrder(int order) { TEUCHOS_ASSERT(metaData_nc_ != Teuchos::null);
                                          metaData_nc_->setOrder(order); }

    virtual void setNFailures(int nFailures) { TEUCHOS_ASSERT(metaData_nc_ != Teuchos::null);
                                               metaData_nc_->setNFailures(nFailures); }
    virtual void setNRunningFailures(int nFailures) { TEUCHOS_ASSERT(metaData_nc_ != Teuchos::null);
                                          metaData_nc_->setNRunningFailures(nFailures); }
    virtual void setNConsecutiveFailures(int nConsecutiveFailures)
                                        { TEUCHOS_ASSERT(metaData_nc_ != Teuchos::null);
                                          metaData_nc_->setNConsecutiveFailures(nConsecutiveFailures); }
    virtual void setTolRel (Scalar tolRel){ TEUCHOS_ASSERT(metaData_nc_ != Teuchos::null);
                                          metaData_nc_->setTolRel(tolRel); }
    virtual void setTolAbs (Scalar tolAbs){ TEUCHOS_ASSERT(metaData_nc_ != Teuchos::null);
                                          metaData_nc_->setTolAbs(tolAbs); }

    virtual void setXNormL2 (Scalar xNormL2){ TEUCHOS_ASSERT(metaData_nc_ != Teuchos::null);
                                          metaData_nc_->setXNormL2(xNormL2); }
    virtual void setDxNormL2Rel (Scalar dxNormL2Rel){ TEUCHOS_ASSERT(metaData_nc_ != Teuchos::null);
                                          metaData_nc_->setDxNormL2Rel(dxNormL2Rel); }
    virtual void setDxNormL2Abs (Scalar dxNormL2Abs){ TEUCHOS_ASSERT(metaData_nc_ != Teuchos::null);
                                          metaData_nc_->setDxNormL2Abs(dxNormL2Abs); }
    virtual void setComputeNorms(bool computeNorms) { TEUCHOS_ASSERT(metaData_nc_ != Teuchos::null);
                                          metaData_nc_->setComputeNorms(computeNorms); }

    virtual void setSolutionStatus(Status s) { TEUCHOS_ASSERT(metaData_nc_ != Teuchos::null);
                                               metaData_nc_->setSolutionStatus(s); }
    virtual void setSolutionStatus(const Thyra::SolveStatus<Scalar> sStatus)
      {
        if (sStatus.solveStatus == Thyra::SOLVE_STATUS_CONVERGED )
          metaData_nc_->setSolutionStatus(Status::PASSED);
        else
          metaData_nc_->setSolutionStatus(Status::FAILED);
      }

    virtual void setOutput(bool output) { TEUCHOS_ASSERT(metaData_nc_ != Teuchos::null);
                                          metaData_nc_->setOutput(output); }
    virtual void setOutputScreen(bool output) { TEUCHOS_ASSERT(metaData_nc_ != Teuchos::null);
                                                metaData_nc_->setOutputScreen(output); }
    virtual void setIsSynced(bool isSynced) { TEUCHOS_ASSERT(metaData_nc_ != Teuchos::null);
                                              metaData_nc_->setIsSynced(isSynced); }
    virtual void setIsInterpolated(bool isInterpolated) { TEUCHOS_ASSERT(metaData_nc_ != Teuchos::null);
                                              metaData_nc_->setIsInterpolated(isInterpolated); }
    virtual void setAccuracy(bool accuracy) { TEUCHOS_ASSERT(metaData_nc_ != Teuchos::null);
                                              metaData_nc_->setAccuracy(accuracy); }
  //@}

  /// \name Get State Data
  //@{
    virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getX()
      { TEUCHOS_ASSERT(x_nc_ != Teuchos::null);
        return x_nc_; }
    virtual Teuchos::RCP<const Thyra::VectorBase<Scalar> > getX() const
      { return x_; }
    virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getXDot()
      { return xdot_nc_; }
    virtual Teuchos::RCP<const Thyra::VectorBase<Scalar> > getXDot() const
      { return xdot_; }
    virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getXDotDot()
      { return xdotdot_nc_; }
    virtual Teuchos::RCP<const Thyra::VectorBase<Scalar> > getXDotDot() const
      { return xdotdot_; }

    virtual Teuchos::RCP<StepperState<Scalar> > getStepperState()
      { TEUCHOS_ASSERT(stepperState_nc_ != Teuchos::null);
        return stepperState_nc_; }
    virtual Teuchos::RCP<const StepperState<Scalar> > getStepperState() const
      { return stepperState_; }

    virtual Teuchos::RCP<PhysicsState<Scalar> > getPhysicsState()
      { return physicsState_nc_; }
    virtual Teuchos::RCP<const PhysicsState<Scalar> > getPhysicsState() const
      { return physicsState_; }
  //@}

  /// \name Set State Data
  //@{
    virtual void setX(Teuchos::RCP<Thyra::VectorBase<Scalar> > x)
      { x_nc_ = x; x_ = x; }
    virtual void setX(Teuchos::RCP<const Thyra::VectorBase<Scalar> > x)
      { x_nc_ = Teuchos::null; x_ = x; }
    virtual void setXDot(Teuchos::RCP<Thyra::VectorBase<Scalar> > xdot)
      { xdot_nc_ = xdot; xdot_ = xdot; }
    virtual void setXDot(Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdot)
      { xdot_nc_ = Teuchos::null; xdot_ = xdot; }
    virtual void setXDotDot(Teuchos::RCP<Thyra::VectorBase<Scalar> > xdotdot)
      { xdotdot_nc_ = xdotdot; xdotdot_ = xdotdot; }
    virtual void setXDotDot(Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdotdot)
      { xdotdot_nc_ = Teuchos::null; xdotdot_ = xdotdot; }

    virtual void setStepperState(Teuchos::RCP<StepperState<Scalar> >& ss)
      { stepperState_nc_ = ss; stepperState_ = ss; }
    virtual void setStepperState(const Teuchos::RCP<StepperState<Scalar> >& ss)
      { stepperState_nc_ = Teuchos::null; stepperState_ = ss; }

    virtual void setPhysicsState(Teuchos::RCP<PhysicsState<Scalar> >& ps)
      { physicsState_nc_ = ps; physicsState_ = ps; }
    virtual void setPhysicsState(const Teuchos::RCP<PhysicsState<Scalar> >& ps)
      { physicsState_nc_ = Teuchos::null; physicsState_ = ps; }
  //@}


  /// \name Comparison methods
  //@{
    /// Less than comparison for sorting based on time
    bool operator< (const SolutionState<Scalar>& ss) const;

    /// Less than comparison for sorting based on time
    bool operator<= (const SolutionState<Scalar>& ss) const;

    /// Less than comparison for sorting based on time
    bool operator< (const Scalar& t) const;

    /// Less than comparison for sorting based on time
    bool operator<= (const Scalar& t) const;

    /// Less than comparison for sorting based on time
    bool operator> (const SolutionState<Scalar>& ss) const;

    /// Less than comparison for sorting based on time
    bool operator>= (const SolutionState<Scalar>& ss) const;

    /// Less than comparison for sorting based on time
    bool operator> (const Scalar& t) const;

    /// Less than comparison for sorting based on time
    bool operator>= (const Scalar& t) const;

    /// Equality comparison for matching
    bool operator== (const SolutionState<Scalar>& ss) const;

    /// Equality comparison for matching
    bool operator== (const Scalar& t) const;
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual std::string description() const;
    virtual void describe(Teuchos::FancyOStream          &out,
                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  /// Compute the solution norms, and solution change from ssIn, if provided.
  virtual void computeNorms(
    const Teuchos::RCP<const SolutionState<Scalar> >& ssIn = Teuchos::null);

private:
  // Member Data

  /// Meta Data for the solution state
  Teuchos::RCP<const SolutionStateMetaData<Scalar> > metaData_;
  Teuchos::RCP<SolutionStateMetaData<Scalar> > metaData_nc_;

  /// Solution
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > x_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > x_nc_;

  /// Time derivative of the solution
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdot_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > xdot_nc_;

  /// Second time derivative of the solution
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdotdot_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > xdotdot_nc_;

  /// StepperState for this SolutionState
  Teuchos::RCP<const StepperState<Scalar> > stepperState_;
  Teuchos::RCP<StepperState<Scalar> > stepperState_nc_;

  /// PhysicsState for this SolutionState
  Teuchos::RCP<const PhysicsState<Scalar> > physicsState_;
  Teuchos::RCP<PhysicsState<Scalar> > physicsState_nc_;

};


/// Nonmember constructor from non-const solution vectors, x.
template<class Scalar>
Teuchos::RCP<SolutionState<Scalar> >
createSolutionStateX(
  const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xdot  = Teuchos::null,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xdotdot = Teuchos::null);

/// Nonmember constructor from const solution vectors, x.
template<class Scalar>
Teuchos::RCP<SolutionState<Scalar> >
createSolutionStateX(
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& x,
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xdot    =Teuchos::null,
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xdotdot =Teuchos::null);

/// Nonmember constructor from Thyra ModelEvaluator.
template<class Scalar>
Teuchos::RCP<SolutionState<Scalar> >
createSolutionStateME(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
  const Teuchos::RCP<StepperState<Scalar> >& stepperState = Teuchos::null,
  const Teuchos::RCP<PhysicsState<Scalar> >& physicsState = Teuchos::null);



} // namespace Tempus

#endif // Tempus_SolutionState_decl_hpp
