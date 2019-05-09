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
#include "Tempus_StepperExplicit.hpp"
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
 *
 *  When using the First-Step-As-Last (FSAL) priniciple, where one can
 *  reuse the last function evaulation as the first evaluation of the next
 *  time step, the algorithm is only slightly more complicated.
 *   - for \f$i = 1 \ldots s\f$ do
 *     - if ( i==1 && useFSAL && (previous step not failed) )
 *       - tmp = \f$\dot{X}_1\f$
 *       - \f$\dot{X}_1 = \dot{X}_s\f$
 *       - \f$\dot{X}_s\f$ = tmp
 *     - else
 *       - \f$X_i \leftarrow x_{n-1}
 *                + \Delta t\,\sum_{j=1}^{i-1} a_{ij}\,\dot{X}_j\f$
 *       - Evaluate \f$\bar{f}(X_{i},t_{n-1}+c_{i}\Delta t)\f$
 *       - \f$\dot{X}_i \leftarrow \bar{f}(X_i,t_{n-1}+c_i\Delta t)\f$
 *   - end for
 *   - \f$x_n \leftarrow x_{n-1} + \Delta t\,\sum_{i=1}^{s}b_i\,\dot{X}_i\f$
 *
 *   For Explicit RK, FSAL requires \f$c_1 = 0\f$, \f$c_s = 1\f$, and
 *   be stiffly accurate (\f$a_{sj} = b_j\f$).  An example of this is
 *   the Bogacki-Shampine 3(2) method.
 *   \f[
 *   \begin{array}{c|cccc}  0  & 0    &     &     &   \\
 *                         1/3 & 1/2  & 0   &     &   \\
 *                         2/3 & 0    & 3/4 & 0   &   \\
 *                          1  & 2/9  & 1/3 & 4/9 & 0 \\ \hline
 *                             & 2/9  & 1/3 & 4/9 & 0
 *                             & 7/24 & 1/4 & 1/3 & 1/8 \end{array}
 *   \f]
 */
template<class Scalar>
class StepperExplicitRK : virtual public Tempus::StepperExplicit<Scalar>
{
public:

  /** \brief Default constructor.
   *
   *  - Constructs with a default ParameterList.
   *  - Can reset ParameterList with setParameterList().
   *  - Requires subsequent setModel() and initialize() calls before calling
   *    takeStep().
  */
  StepperExplicitRK();

  /// Constructor to specialize Stepper parameters.
  StepperExplicitRK(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    Teuchos::RCP<Teuchos::ParameterList> pList);

  /// Constructor to use default Stepper parameters.
  StepperExplicitRK(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    std::string stepperType = "RK Explicit 4 Stage");

  /// Constructor for StepperFactory.
  StepperExplicitRK(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    std::string stepperType, Teuchos::RCP<Teuchos::ParameterList> pList);

  /// \name Basic stepper methods
  //@{
    virtual void setObserver(
      Teuchos::RCP<StepperObserver<Scalar> > obs = Teuchos::null);

    virtual void setTableau(std::string stepperType);

    virtual void setTableau(
      Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null);

    virtual void setTableau(
      Teuchos::RCP<const RKButcherTableau<Scalar> > ERK_ButcherTableau);

    /// Initialize during construction and after changing input parameters.
    virtual void initialize();

    /// Set the initial conditions and make them consistent.
    virtual void setInitialConditions (
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    virtual std::string getStepperType() const
     { return this->stepperPL_->template get<std::string>("Stepper Type"); }

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

    virtual OrderODE getOrderODE()   const {return FIRST_ORDER_ODE;}
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

protected:

  Teuchos::RCP<const RKButcherTableau<Scalar> >          ERK_ButcherTableau_;

  std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar> > > stageXDot_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >               stageX_;

  Teuchos::RCP<StepperExplicitRKObserver<Scalar> >  stepperExplicitRKObserver_;

  // For Embedded RK
  Teuchos::RCP<Thyra::VectorBase<Scalar> >               ee_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >               abs_u0;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >               abs_u;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >               sc;

};

} // namespace Tempus

#endif // Tempus_StepperExplicitRK_decl_hpp
