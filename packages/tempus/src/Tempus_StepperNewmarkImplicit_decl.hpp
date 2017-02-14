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
    Teuchos::RCP<Teuchos::ParameterList>                      pList,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& transientModel);

  /// \name Basic stepper methods
  //@{
    virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& transientModel);
    virtual void setNonConstModel(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& transientModel);
    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
      getModel(){return residualModel_->getTransientModel();}

    /// Set the solver
    void setSolver(std::string solverName);
    void setSolver(Teuchos::RCP<Teuchos::ParameterList> solverPL=Teuchos::null);
    /// Set the predictor
    void setPredictor(std::string predictorName);
    void setPredictor(Teuchos::RCP<Teuchos::ParameterList> predPL=Teuchos::null);

    /// Initialize during construction and after changing input parameters.
    virtual void initialize();

    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    /// Get a default (initial) StepperState
    virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState();
    virtual Scalar getOrder() const {return 1.0;}
    virtual Scalar getOrderMin() const {return 1.0;}
    virtual Scalar getOrderMax() const {return 1.0;}
  //@}

  /// Compute predictor given the supplied stepper
  virtual void computePredictor(
    const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

  /// \name ParameterList methods
  //@{
    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & pl);
    Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
    Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
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

  Teuchos::RCP<Teuchos::ParameterList>                         pList_;
  Teuchos::RCP<SecondOrderResidualModelEvaluator<Scalar> >     residualModel_;
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >            solver_;
  Teuchos::RCP<Stepper<Scalar> >                               predictorStepper_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar>                    inArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar>                   outArgs_;

  // Compute the balancing time derivative as a function of x
  std::function<void (const Thyra::VectorBase<Scalar> &,
                            Thyra::VectorBase<Scalar> &)>
  xDotFunction(Scalar dt,Teuchos::RCP<const Thyra::VectorBase<Scalar> > x_old);

  Scalar beta_; 
  Scalar gamma_; 


};
} // namespace Tempus

#endif // Tempus_StepperNewmarkImplicit_decl_hpp
