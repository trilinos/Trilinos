#ifndef TEMPUS_INTEGRATOR_HPP
#define TEMPUS_INTEGRATOR_HPP

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include <string>

namespace tempus {

  /** \brief Simple time integrator
   */
  template<class Scalar>
  class IntegratorSimple : public tempus::Integrator {
  public:

    /** \brief Constructor with ParameterList, models, initial conditions
     *  and optional solvers. */
    IntegratorSimple(
      RCP<ParameterList>              parameterList,
      RCP<Thyra::VectorBase<Scalar> > x,
      RCP<Thyra::VectorBase<Scalar> > xdot=Teuchos::null,
      RCP<Thyra::VectorBase<Scalar> > xdotdot=Teuchos::null );

    /// Destructor
    virtual ~IntegratorSimple() {}

    //! Unique name for this integrator.
    virtual std::string name() const = 0;

    /// \name Basic integrator methods
    //@{
    /// Advance the solution to time, and return true if successful.
    virtual bool advanceTime(const Scalar time) = 0;

    /// Get the current solution, x.
    virtual std::vector<RCP<const Thyra::VectorBase<Scalar> > > getXs() = 0;

    /// Get the current time derivative of the solution, xdot.
    virtual std::vector<RCP<const Thyra::VectorBase<Scalar> > > getXDots();

    /// Get the current second time derivative of the solution, xdotdot.
    virtual std::vector<RCP<const Thyra::VectorBase<Scalar> > > getXDotDots();
    //@}

    /// \name ParameterList methods
    //@{
    virtual void setParameterList(RCP<ParameterList> const& pl);
    virtual RCP<ParameterList> getNonconstParameterList();
    virtual RCP<ParameterList> unsetParameterList();
    virtual RCP<const ParameterList> getValidParameters() const;
    //@}

    /// \name Solver methods
    //@{
    virtual void setSolvers(
      const std:vector<RCP<Thyra::NonlinearSolverBase<Scalar> > &> solvers);
    virtual std:vector<RCP<Thyra::NonlinearSolverBase<Scalar> > >
      getNonconstSolvers();
    virtual std:vector<RCP<const Thyra::NonlinearSolverBase<Scalar> >
      getSolvers() const;
    //@}

    /// \name Accessor methods
    //@{
    virtual std::string description() const;
    virtual void describe( Teuchos::FancyOStream        & out,
                           const Teuchos::EVerbosityLevel verbLevel) const;
    /// Get time
    virtual Scalar getTime() const{return time;}
    /// Set time
    virtual setTime(Scalar time_){time = time_;}
    /// Get time step index
    virtual Scalar getIStep() const{return iStep;}
    /// Set time step index
    virtual void setIStep(Scalar iStep_){iStep = iStep_;}
    /// Get time step size
    virtual Scalar getDt() const{return dt;}
    /// Set time step size
    virtual void setDt(Scalar dt_){dt = dt_;}
    /// Get order of integrator
    virtual int getOrder() const{return order;}
    /// Get error of integrator
    virtual Scalar getError() const{return error;}
    virtual bool isInitialized() const;
    //@}

    /// \name Time step control methods
    //@{

    /// Error control methods
    /// Ramping control methods
    /// Stability control methods

    //@}

    /// \name Observer methods
    //@{

    //@}

    /// \name Adjoint methods
    //@{
    virtual Scalar advanceAdjointTime();
    //@}

    /// \name Undo type capabilities
    //@{
    /// Only accept step after meeting time step criteria.
    virtual bool acceptStep() {return false;}
    //@}

    /// \name Solution history methods
    //@{

    /// Functionality like InterpolationBuffer for checkpointing and restart

    //@}

  protected:

    Scalar time_;
    Scalar dt_;
    int    iStep_;
    int    order_;
    Scalar error_;
    bool   haveParameterList_;
    bool   haveModel_;
    bool   haveSolution_;
    bool   haveSolver_;
    bool   isInitialized_;
    RCP<const ParameterList>                   parameterList_;
    RCP<const Thyra::ModelEvaluator<Scalar> >& model_;
    RCP<const Thyra::VectorBase<Scalar> >&     solution_;
    RCP<Thyra::NonlinearSolverBase<Scalar> >   solver_;

};
} // namespace tempus
#endif // TEMPUS_INTEGRATOR_HPP
