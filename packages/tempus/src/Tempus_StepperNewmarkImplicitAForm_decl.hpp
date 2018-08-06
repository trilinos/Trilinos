// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperNewmarkImplicitAForm_decl_hpp
#define Tempus_StepperNewmarkImplicitAForm_decl_hpp

#include "Tempus_StepperImplicit.hpp"
#include "Tempus_WrapperModelEvaluatorSecondOrder.hpp"

namespace Tempus {


/** \brief Newmark time stepper in acceleration form (a-form).
 *
 *  Here, we implement the Newmark scheme in predictor/corrector form;
 *  see equations (34)-(35) in: A. Mota, W. Klug, M. Ortiz,  "Finite element
 *  simulation of firearm injury to the human cranium",  Computational
 *  Mechanics 31(1) 115-121 (2003).
 *
 *  Newmark has two parameters: \f$\beta\f$
 *  and \f$\gamma\f$, both of which need to be in the range \f$[0,1]\f$.
 *  Newmark can be an explicit or implicit method, depending on
 *  the value of the \f$\beta\f$ parameter. If \f$\beta = 0\f$, the method
 *  is explicit.  Regardless of whether the method is implicit
 *  or explicit, a linear solve is required.  This linear solve can be
 *  optimized, however, for the explicit case by lumping the mass matrix.
 *  This optimization can be invoked by running "Newmark Explicit d-Form"
 *  Stepper through the Piro::TempusSolver class.
 *
 *  Newmark is second order accurate if \f$\gamma =  0.5\f$; otherwise it
 *  is first order accurate.  Some additional properties about the Newmark
 *  scheme can be found
 *<a href="http://opensees.berkeley.edu/wiki/index.php/Newmark_Method">here</a>.
 */
template<class Scalar>
class StepperNewmarkImplicitAForm
 : virtual public Tempus::StepperImplicit<Scalar>
{
public:

  /// Constructor
  StepperNewmarkImplicitAForm(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null);

  /// \name Basic stepper methods
  //@{
    virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel);

    virtual void setObserver(
      Teuchos::RCP<StepperObserver<Scalar> > obs = Teuchos::null){}

    /// Initialize during construction and after changing input parameters.
    virtual void initialize();

    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    /// Get a default (initial) StepperState
    virtual Teuchos::RCP<Tempus::StepperState<Scalar> >
      getDefaultStepperState();
    virtual Scalar getOrder() const {
      if (gamma_ == 0.5) return 2.0;
      else return 1.0;
    }
    virtual Scalar getOrderMin() const {return 1.0;}
    virtual Scalar getOrderMax() const {return 2.0;}

    virtual bool isExplicit()         const {return false;}
    virtual bool isImplicit()         const {return true;}
    virtual bool isExplicitImplicit() const
      {return isExplicit() and isImplicit();}
    virtual bool isOneStepMethod()   const {return true;}
    virtual bool isMultiStepMethod() const {return !isOneStepMethod();}
  //@}
  
  /// Pass initial guess to Newton solver  
  virtual void setInitialGuess(Teuchos::RCP<const Thyra::VectorBase<Scalar> > initial_guess)
       {initial_guess_ = initial_guess;}

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

    void predictVelocity(Thyra::VectorBase<Scalar>& vPred,
                             const Thyra::VectorBase<Scalar>& v,
                             const Thyra::VectorBase<Scalar>& a,
                             const Scalar dt) const;

    void predictDisplacement(Thyra::VectorBase<Scalar>& dPred,
                               const Thyra::VectorBase<Scalar>& d,
                               const Thyra::VectorBase<Scalar>& v,
                               const Thyra::VectorBase<Scalar>& a,
                               const Scalar dt) const;

    void correctVelocity(Thyra::VectorBase<Scalar>& v,
                             const Thyra::VectorBase<Scalar>& vPred,
                             const Thyra::VectorBase<Scalar>& a,
                             const Scalar dt) const;

    void correctDisplacement(Thyra::VectorBase<Scalar>& d,
                               const Thyra::VectorBase<Scalar>& dPred,
                               const Thyra::VectorBase<Scalar>& a,
                               const Scalar dt) const;

private:

  /// Default Constructor -- not allowed
  StepperNewmarkImplicitAForm();

private:

  Thyra::ModelEvaluatorBase::InArgs<Scalar>                inArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar>               outArgs_;

  Scalar beta_;
  Scalar gamma_;


  Teuchos::RCP<Teuchos::FancyOStream> out_;

  Teuchos::RCP<const Thyra::VectorBase<Scalar> >      initial_guess_;


};
} // namespace Tempus

#endif // Tempus_StepperNewmarkImplicitAForm_decl_hpp
