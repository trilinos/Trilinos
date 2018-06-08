// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperExplicitRK_decl_hpp
#define Tempus_StepperExplicitRK_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_Stepper.hpp"
#include "Tempus_RKButcherTableau.hpp"
#include "Tempus_StepperExplicitRKObserver.hpp"


namespace Tempus {

/** \brief Explicit Runge-Kutta time stepper.
 *
 *  For the explicit ODE system,
 *  \f[
 *    \dot{x} = \bar{f}(x,t),
 *  \f]
 *  the general explicit Runge-Kutta method for \f$s\f$-stages can be
 *  written as
 *  \f[
 *    X_{i} = x_{n-1}
 *    + \Delta t\,\sum_{j=1}^{i-1} a_{ij}\,\bar{f}(X_{j},t_{n-1}+c_{j}\Delta t)
 *  \f]
 *  \f[
 *    x_{n} = x_{n-1}
 *    + \Delta t\,\sum_{i=1}^{s}b_{i}\,\bar{f}(X_{i},t_{n-1}+c_{i}\Delta t)
 *  \f]
 *  where \f$X_{i}\f$ are intermediate approximations to the solution
 *  at times, \f$t_{n-1}+c_{i}\Delta t\f$, (stage solutions) which may
 *  be correct to a lower order of accuracy than the solution, \f$x_{n}\f$.
 *  We should note that these lower-order approximations are combined
 *  through \f$b_{i}\f$ so that error terms cancel out and produce a
 *  more accurate solution. Note for explicit RK that \f$a_{ij}=0\f$ for
 *  \f$j \leq i\f$ and does not require any solves.
 *  Note that the stage time derivatives are
 *  \f[
 *    \dot{X}_{i} = \bar{f}(X_{i},t_{n-1}+c_{i}\Delta t),
 *  \f]
 *  and the time derivative by definition is
 *  \f[
 *    \dot{x}_{n} = \bar{f}(x_{n},t_{n}),
 *  \f]
 *
 *  <b> Algorithm </b>
 *  The single-timestep algorithm for Explicit RK is simply,
 *   - for \f$i = 1 \ldots s\f$ do
 *     - \f$X_i \leftarrow x_{n-1}
 *              + \Delta t\,\sum_{j=1}^{i-1} a_{ij}\,\dot{X}_j\f$
 *     - Evaluate \f$\bar{f}(X_{i},t_{n-1}+c_{i}\Delta t)\f$
 *     - \f$\dot{X}_i \leftarrow \bar{f}(X_i,t_{n-1}+c_i\Delta t)\f$
 *   - end for
 *   - \f$x_n \leftarrow x_{n-1} + \Delta t\,\sum_{i=1}^{s}b_i\,\dot{X}_i\f$
 *   - \f$\dot{x}_n \leftarrow \bar{f}(x_{n},t_{n})\f$ [Optional]
 */
template<class Scalar>
class StepperExplicitRK : virtual public Tempus::Stepper<Scalar>
{
public:

  /// Constructor to use default Stepper parameters.
  StepperExplicitRK(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    std::string stepperType = "RK Explicit 4 Stage");

  /// Constructor to specialize Stepper parameters.
  StepperExplicitRK(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    Teuchos::RCP<Teuchos::ParameterList> pList);

  /// Constructor for StepperFactory.
  StepperExplicitRK(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    std::string stepperType, Teuchos::RCP<Teuchos::ParameterList> pList);

  /// \name Basic stepper methods
  //@{
    virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel);
    virtual void setNonConstModel(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& appModel);
    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
      getModel(){return appModel_;}

    virtual void setSolver(std::string solverName);
    virtual void setSolver(
      Teuchos::RCP<Teuchos::ParameterList> solverPL=Teuchos::null);
    virtual void setSolver(
        Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver);
    virtual Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > getSolver() const
    { return Teuchos::null; }
    virtual void setObserver(
      Teuchos::RCP<StepperObserver<Scalar> > obs = Teuchos::null);

    void setTableau(
      Teuchos::RCP<Teuchos::ParameterList> pList,
      std::string stepperType = "");

    /// Initialize during construction and after changing input parameters.
    virtual void initialize();

    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    /// Get a default (initial) StepperState
    virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState();
    virtual Scalar getOrder() const {return ERK_ButcherTableau_->order();}
    virtual Scalar getOrderMin() const {return ERK_ButcherTableau_->orderMin();}
    virtual Scalar getOrderMax() const {return ERK_ButcherTableau_->orderMax();}
    virtual Scalar getInitTimeStep(
        const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) const;

    virtual bool isExplicit()         const {return true;}
    virtual bool isImplicit()         const {return false;}
    virtual bool isExplicitImplicit() const
      {return isExplicit() and isImplicit();}
    virtual bool isOneStepMethod()   const {return true;}
    virtual bool isMultiStepMethod() const {return !isOneStepMethod();}
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

private:

  /// Default Constructor -- not allowed
  StepperExplicitRK();

protected:

  std::string                                            description_;
  Teuchos::RCP<Teuchos::ParameterList>                   stepperPL_;
  /// Explicit ODE ModelEvaluator
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >     appModel_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar>              inArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar>             outArgs_;

  Teuchos::RCP<const RKButcherTableau<Scalar> >          ERK_ButcherTableau_;

  std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar> > > stageXDot_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >               stageX_;

  Teuchos::RCP<StepperObserver<Scalar> >            stepperObserver_;
  Teuchos::RCP<StepperExplicitRKObserver<Scalar> >  stepperExplicitRKObserver_;

  // For Embedded RK
  Teuchos::RCP<Thyra::VectorBase<Scalar> >               ee_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >               abs_u0;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >               abs_u;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >               sc;
};

} // namespace Tempus

#endif // Tempus_StepperExplicitRK_decl_hpp
