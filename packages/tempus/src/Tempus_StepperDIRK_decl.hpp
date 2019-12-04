// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperDIRK_decl_hpp
#define Tempus_StepperDIRK_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_RKButcherTableau.hpp"
#include "Tempus_StepperImplicit.hpp"
#include "Tempus_WrapperModelEvaluator.hpp"
#include "Tempus_StepperRKObserverComposite.hpp"


namespace Tempus {

/** \brief Diagonally Implicit Runge-Kutta (DIRK) time stepper.
 *
 *  For the implicit ODE system, \f$\mathcal{F}(\dot{x},x,t) = 0\f$,
 *  the general DIRK method for \f$s\f$-stages, can be written as
 *  \f[
 *    X_{i} = x_{n-1}
 *    + \Delta t\,\sum_{j=1}^{i-1} a_{ij}\,\bar{f}(X_{j},t_{n-1}+c_{j}\Delta t)
 *    + \Delta t\, a_{ii}\,\bar{f}(X_{i},t_{n-1}+c_{i}\Delta t)
 *  \f]
 *  \f[
 *    x_{n} = x_{n-1}
 *    + \Delta t\,\sum_{i=1}^{s}b_{i}\,\bar{f}(X_{i},t_{n-1}+c_{i}\Delta t)
 *  \f]
 *  where \f$\dot{x}=\bar{f}(x,t)\f$ is the explicit form of the
 *  ODE, \f$X_{i}\f$ are intermediate approximations to the solution
 *  at times, \f$t_{n-1}+c_{i}\Delta t\f$, (stage solutions) which may
 *  be correct to a lower order of accuracy than the solution, \f$x_{n}\f$.
 *  We should note that these lower-order approximations are combined
 *  through \f$b_{i}\f$ so that error terms cancel out and produce a
 *  more accurate solution. Note for DIRK methods that \f$a_{ii}=a\f$
 *  for all the diagonal components is referred to as
 *  Singly Diagonal Implicit Runge-Kutta (SDIRK) methods.
 *
 *  Note that the stage time derivatives,
 *  \f[
 *    \dot{X}_{i} = \bar{f}(X_{i},t_{n-1}+c_{i}\Delta t),
 *  \f]
 *  can be found via
 *  \f{eqnarray*}{
 *    \dot{X}_{i} & = & \frac{1}{a_{ii} \Delta t} [ X_{i} - x_{n-1}
 *                     - \Delta t\,\sum_{j=1}^{i-1} a_{ij}\,\dot{X}_{j} ] \\
 *    \dot{X}_{i} & = & \frac{X_{i} - \tilde{X}}{a_{ii} \Delta t}
 *  \f}
 *  where
 *  \f[
 *    \tilde{X} = x_{n-1} + \Delta t \sum_{j=1}^{i-1} a_{ij}\, \dot{X}_{j}
 *  \f]
 *  Recalling that the definition for a DIRK is that for \f$j>i\f$,
 *  \f$a_{ij} = 0\f$ and \f$a_{ii} \neq 0\f$ for at least one \f$i\f$.
 *  Thus for stages where \f$a_{ii} = 0\f$, we can use the explicit RK
 *  methods (see StepperExplicitRK for additional details).
 *
 *  <b> Algorithm </b>
 *  The single-timestep algorithm for DIRK is,
 *   - for \f$i = 1 \ldots s\f$ do
 *     - if \f$a_{ii} = 0\f$
 *       - \f$X_i \leftarrow x_{n-1}
 *                + \Delta t\,\sum_{j=1}^{i-1} a_{ij}\,\dot{X}_j\f$
 *       - Evaluate \f$\bar{f}(X_{i},t_{n-1}+c_{i}\Delta t)\f$
 *       - \f$\dot{X}_i \leftarrow \bar{f}(X_i,t_{n-1}+c_i\Delta t)\f$
 *     - else
 *       - \f$\tilde{X} =
 *           x_{n-1} +\Delta t \sum_{j=1}^{i-1} a_{ij}\,\dot{X}_{j}\f$
 *       - Define \f$\dot{X}_i \leftarrow
 *                              \frac{X_{i} - \tilde{X}}{a_{ii} \Delta t}\f$
 *       - Solve \f$f(\dot{x} =
 *           \dot{X}_i,X_i,t_{n-1}+c_{i}\Delta t)=0\f$ for \f$X_i\f$
 *       - \f$\dot{X}_i \leftarrow \frac{X_{i} - \tilde{X}}{a_{ii} \Delta t}\f$
 *   - end for
 *   - \f$x_n \leftarrow x_{n-1} + \Delta t\,\sum_{i=1}^{s}b_i\,\dot{X}_i\f$
 *
 *  The First-Step-As-Last (FSAL) principle is not needed with DIRK, but
 *  maybe useful if the first stage is explicit (i.e, EDIRK).
 *  The default is to set useFSAL=false.
 */
template<class Scalar>
class StepperDIRK : virtual public Tempus::StepperImplicit<Scalar>
{
public:

  /// \name Basic stepper methods
  //@{
    virtual void setObserver(
      Teuchos::RCP<StepperObserver<Scalar> > obs = Teuchos::null);

    virtual Teuchos::RCP<const RKButcherTableau<Scalar> > getTableau()
    { return tableau_; }

    virtual Teuchos::RCP<StepperObserver<Scalar> > getObserver() const
    { return this->stepperObserver_; }

    /// Initialize during construction and after changing input parameters.
    virtual void initialize();

    /// Set the initial conditions and make them consistent.
    virtual void setInitialConditions (
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    /// Set parameter so that the initial guess is reset at the beginning of each timestep.
    virtual void setResetInitialGuess(bool reset_guess)
      { resetGuess_ = reset_guess; }
    virtual bool getResetInitialGuess() const
      { return resetGuess_; }

    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    /// Get a default (initial) StepperState
    virtual Teuchos::RCP<Tempus::StepperState<Scalar> >getDefaultStepperState();
    virtual Scalar getOrder()    const{return tableau_->order();}
    virtual Scalar getOrderMin() const{return tableau_->orderMin();}
    virtual Scalar getOrderMax() const{return tableau_->orderMax();}

    virtual bool isExplicit() const
    {
      const int numStages = tableau_->numStages();
      Teuchos::SerialDenseMatrix<int,Scalar> A = tableau_->A();
      bool isExplicit = false;
      for (int i=0; i<numStages; ++i) if (A(i,i) == 0.0) isExplicit = true;
      return isExplicit;
    }
    virtual bool isImplicit()         const {return true;}
    virtual bool isExplicitImplicit() const
      {return isExplicit() and isImplicit();}
    virtual bool isOneStepMethod()   const {return true;}
    virtual bool isMultiStepMethod() const {return !isOneStepMethod();}

    virtual OrderODE getOrderODE()   const {return FIRST_ORDER_ODE;}

    void getValidParametersBasicDIRK(
      Teuchos::RCP<Teuchos::ParameterList> pl) const;

    virtual std::string getDescription() const = 0;
  //@}

  Teuchos::RCP<Thyra::VectorBase<Scalar> >& getStageX() {return stageX_;};
  std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar> > >& getStageXDot() {return stageXDot_;};
  Teuchos::RCP<Thyra::VectorBase<Scalar> >& getXTilde() {return xTilde_;};

  /// Return alpha = d(xDot)/dx.
  virtual Scalar getAlpha(const Scalar dt) const
  {
    const Teuchos::SerialDenseMatrix<int,Scalar> & A=tableau_->A();
    return Scalar(1.0)/(dt*A(0,0));  // Getting the first diagonal coeff!
  }
  /// Return beta  = d(x)/dx.
  virtual Scalar getBeta (const Scalar   ) const { return Scalar(1.0); }

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual void describe(Teuchos::FancyOStream        & out,
                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  /// \name Accessors methods
  //@{
    /** \brief Use embedded if avialable. */
    virtual void setUseEmbedded(bool a) { useEmbedded_ = a; }
    virtual bool getUseEmbedded() const { return useEmbedded_; }
    virtual bool getUseEmbeddedDefault() const { return false; }
  //@}


protected:

  /// Default setup for constructor.
  virtual void setupDefault();

  /// Setup for constructor.
  virtual void setup(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& wrapperModel,
    const Teuchos::RCP<StepperRKObserver<Scalar> >& obs,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess);

  virtual void setupTableau() = 0;

  Teuchos::RCP<RKButcherTableau<Scalar> >                tableau_;

  std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar> > > stageXDot_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >               stageX_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >               xTilde_;

  Teuchos::RCP<StepperRKObserverComposite<Scalar> >        stepperObserver_;

  // For Embedded RK
  bool useEmbedded_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >               ee_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >               abs_u0;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >               abs_u;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >               sc;

  bool resetGuess_ = true;
};


/** \brief Time-derivative interface for DIRK.
 *
 *  Given the stage state \f$X_i\f$ and
 *  \f[
 *    \tilde{X} = x_{n-1} +\Delta t \sum_{j=1}^{i-1} a_{ij}\,\dot{X}_{j},
 *  \f]
 *  compute the DIRK stage time-derivative,
 *  \f[
 *    \dot{X}_i = \frac{X_{i} - \tilde{X}}{a_{ii} \Delta t}
 *  \f]
 *  \f$\ddot{x}\f$ is not used and set to null.
 */
template <typename Scalar>
class StepperDIRKTimeDerivative
  : virtual public Tempus::TimeDerivative<Scalar>
{
public:

  /// Constructor
  StepperDIRKTimeDerivative(
    Scalar s, Teuchos::RCP<const Thyra::VectorBase<Scalar> > xTilde)
  { initialize(s, xTilde); }

  /// Destructor
  virtual ~StepperDIRKTimeDerivative() {}

  /// Compute the time derivative.
  virtual void compute(
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > x,
    Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDot,
    Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDotDot = Teuchos::null)
  {
    xDotDot = Teuchos::null;
    Thyra::V_StVpStV(xDot.ptr(),s_,*x,-s_,*xTilde_);
  }

  virtual void initialize(Scalar s,
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xTilde)
  { s_ = s; xTilde_ = xTilde; }

private:

  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xTilde_;
  Scalar                                         s_;      // = 1/(dt*a_ii)
};


} // namespace Tempus

#endif // Tempus_StepperDIRK_decl_hpp
