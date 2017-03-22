// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperNewmarkImplicit_decl_hpp
#define Tempus_StepperNewmarkImplicit_decl_hpp

#include "Tempus_StepperImplicit.hpp"
#include "Tempus_SecondOrderResidualModelEvaluator.hpp"

namespace Tempus {


/** \brief Backward Euler time stepper.
 *  Backward Euler is an implicit time stepper (i.e., with solver).
 */
template<class Scalar>
class StepperNewmarkImplicit : virtual public Tempus::StepperImplicit<Scalar>
{
public:

  /// Constructor
  StepperNewmarkImplicit(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& transientModel,
    Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null);

  /// \name Basic stepper methods
  //@{
    virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& transientModel);
    virtual void setNonConstModel(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& transientModel);
    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
      getModel(){return residualModel_->getTransientModel();}

    /// Set the solver
    virtual void setSolver(std::string solverName);
    virtual void setSolver(
      Teuchos::RCP<Teuchos::ParameterList> solverPL=Teuchos::null);
    virtual void setSolver(
      Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver);


    /// Initialize during construction and after changing input parameters.
    virtual void initialize();

    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    /// Get a default (initial) StepperState
    virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState();
    virtual Scalar getOrder() const {
      if (gamma_ == 0.5) return 2.0;
      else return 1.0;
    }
    virtual Scalar getOrderMin() const {return 1.0;}
    virtual Scalar getOrderMax() const {return 2.0;}
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
  StepperNewmarkImplicit();

private:

  Teuchos::RCP<Teuchos::ParameterList>                     stepperPL_;
  Teuchos::RCP<SecondOrderResidualModelEvaluator<Scalar> > residualModel_;
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >        solver_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar>                inArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar>               outArgs_;

  Scalar beta_;
  Scalar gamma_;


  Teuchos::RCP<Teuchos::FancyOStream> out_;

};
} // namespace Tempus

#endif // Tempus_StepperNewmarkImplicit_decl_hpp
