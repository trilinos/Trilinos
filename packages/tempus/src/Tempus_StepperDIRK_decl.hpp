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
#include "Tempus_StepperRKBase.hpp"
#include "Tempus_StepperImplicit.hpp"
#include "Tempus_WrapperModelEvaluator.hpp"
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  #include "Tempus_StepperRKObserverComposite.hpp"
#endif


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
 *  The single-timestep algorithm for DIRK is
 *
 *  \f{algorithm}{
 *  \renewcommand{\thealgorithm}{}
 *  \caption{DIRK with the application-action locations indicated.}
 *  \begin{algorithmic}[1]
 *    \State {\it appAction.execute(solutionHistory, stepper, BEGIN\_STEP)}
 *    \If {``Reset initial guess.''}
 *      \State $X \leftarrow x_{n-1}$
 *        \Comment{Reset initial guess to last timestep.}
 *    \EndIf
 *    \For {$i = 0 \ldots s-1$}
 *      \If { $a_{k,i} = 0 \;\forall k = (i+1,\ldots, s-1)$, $b(i) = 0$, $b^\ast(i) = 0$}
 *        \State $\dot{X}_i \leftarrow 0$
 *          \Comment{Not needed for later calculations.}
 *        \State {\bf continue}
 *      \EndIf
 *      \State $\tilde{X} \leftarrow
 *                      x_{n-1} +\Delta t \sum_{j=1}^{i-1} a_{ij}\,\dot{X}_{j}$
 *      \State {\it appAction.execute(solutionHistory, stepper, BEGIN\_STAGE)}
 *      \If {$a_{ii} = 0$}             \Comment{Explicit stage.}
 *        \If {$i=0$ and ``Use FSAL''} \Comment{Save an evaluation?}
 *          \State $\dot{X}_0 \leftarrow \dot{X}_{s-1}$
 *            \Comment{Use $\dot{X}_{s-1}$ from $n-1$ time step.}
 *        \Else
 *          \State $\dot{X}_i \leftarrow \bar{f}(\tilde{X},t_{n-1}+c_i\Delta t)$
 *        \EndIf
 *      \Else                          \Comment{Implicit stage.}
 *        \State {\it appAction.execute(solutionHistory, stepper, BEFORE\_SOLVE)}
 *        \If {``Zero initial guess.''}
 *          \State $X \leftarrow 0$
 *            \Comment{Else use previous stage value as initial guess.}
 *        \EndIf
 *        \State Solve $\mathcal{F}_i(
 *                      \dot{X}_i = \frac{X - \tilde{X}}{a_{ii} \Delta t},
 *                      X, t_{n-1}+c_{i}\Delta t) = 0$ for $X$
 *        \State {\it appAction.execute(solutionHistory, stepper, AFTER\_SOLVE)}
 *        \State $\dot{X}_i \leftarrow \frac{X - \tilde{X}}{a_{ii} \Delta t}$
 *      \EndIf
 *      \State {\it appAction.execute(solutionHistory, stepper, BEFORE\_EXPLICIT\_EVAL)}
 *      \State {\it appAction.execute(solutionHistory, stepper, END\_STAGE)}
 *    \EndFor
 *    \State $x_n \leftarrow x_{n-1} + \Delta t\,\sum_{i=0}^{s-1}b_i\,\dot{X}_i$
 *    \If {``Embedded''}  \Comment{Compute the local truncation error estimate.}
 *      \State $\mathbf{e} \leftarrow
 *                         \sum_{i=0}^{s-1} (b_i-b^\ast_i)\Delta t\,\dot{X}_i$
 *      \State $\tilde{\mathbf{e}} \leftarrow
 *             \mathbf{e}/(a_{tol} + \max(\|x_n\|, \|x_{n-1}\|)r_{tol})$
 *      \State $e_n \leftarrow \|\tilde{\mathbf{e}}\|_\infty$
 *    \EndIf
 *    \State {\it appAction.execute(solutionHistory, stepper, END\_STEP)}
 *  \end{algorithmic}
 *  \f}
 *
 *  The First-Step-As-Last (FSAL) principle is not needed with DIRK, but
 *  maybe useful if the first stage is explicit (EDIRK) (e.g., Trapezoidal
 *  Method).  The default is to set useFSAL=false.
 *
 *  <b> Iteration Matrix, \f$W\f$.</b>
 *  Recalling that the definition of the iteration matrix, \f$W\f$, is
 *  \f[
 *    W = \alpha \frac{\partial \mathcal{F}_n}{\partial \dot{x}_n}
 *      + \beta  \frac{\partial \mathcal{F}_n}{\partial x_n},
 *  \f]
 *  where \f$ \alpha \equiv \frac{\partial \dot{x}_n(x_n) }{\partial x_n}, \f$
 *  and \f$ \beta \equiv \frac{\partial x_n}{\partial x_n} = 1\f$. For the stage
 *  solutions, we have
 *  \f[
 *    \mathcal{F}_i = \dot{X}_{i} - \bar{f}(X_{i},t_{n-1}+c_{i}\Delta t) =0.
 *  \f]
 *  where \f$\mathcal{F}_n \rightarrow \mathcal{F}_i\f$,
 *  \f$x_n \rightarrow X_{i}\f$, and
 *  \f$\dot{x}_n(x_n) \rightarrow \dot{X}_{i}(X_{i})\f$.
 *  The time derivative for the DIRK stages is
 *  \f[
 *    \dot{X}_{i}(X_{i}) = \frac{X_{i} - \tilde{X}}{a_{ii} \Delta t},
 *  \f]
 *  and we can determine that
 *  \f$ \alpha = \frac{1}{a_{ii} \Delta t} \f$
 *  and \f$ \beta = 1 \f$, and therefore write
 *  \f[
 *    W = \frac{1}{a_{ii} \Delta t}
 *        \frac{\partial \mathcal{F}_i}{\partial \dot{X}_i}
 *      + \frac{\partial \mathcal{F}_i}{\partial X_i}.
 *  \f]
 */
template<class Scalar>
class StepperDIRK : virtual public Tempus::StepperImplicit<Scalar>,
                    virtual public Tempus::StepperRKBase<Scalar>
{
public:

  /// \name Basic stepper methods
  //@{
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    virtual void setObserver(
      Teuchos::RCP<StepperObserver<Scalar> > obs = Teuchos::null);

    virtual Teuchos::RCP<StepperObserver<Scalar> > getObserver() const
    { return this->stepperObserver_; }
#endif
    /// Initialize after construction and changing input parameters.
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

    virtual bool isExplicit() const
    {
      const int numStages = this->tableau_->numStages();
      Teuchos::SerialDenseMatrix<int,Scalar> A = this->tableau_->A();
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

  std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar> > >& getStageXDot() {return stageXDot_;};
  Teuchos::RCP<Thyra::VectorBase<Scalar> >& getXTilde() {return xTilde_;};

  /// Return alpha = d(xDot)/dx.
  virtual Scalar getAlpha(const Scalar dt) const
  {
    const Teuchos::SerialDenseMatrix<int,Scalar> & A=this->tableau_->A();
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

  virtual bool isValidSetup(Teuchos::FancyOStream & out) const;

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
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  virtual void setup(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& wrapperModel,
    const Teuchos::RCP<StepperRKObserver<Scalar> >& obs,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess);
#endif
  virtual void setup(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& wrapperModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction);

  virtual void setupTableau() = 0;

  std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar> > > stageXDot_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >               xTilde_;

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  Teuchos::RCP<StepperRKObserverComposite<Scalar> >      stepperObserver_;
#endif

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
