// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperImplicit_decl_hpp
#define Tempus_StepperImplicit_decl_hpp

// Tempus
#include "Tempus_Stepper.hpp"
#include "Tempus_WrapperModelEvaluatorBasic.hpp"


namespace Tempus {


/** \brief Thyra Base interface for implicit time steppers.
 *
 */
template<class Scalar>
class StepperImplicit : virtual public Tempus::Stepper<Scalar>
{
public:

  /// \name Basic implicit stepper methods
  //@{
    virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel);
    virtual void setNonConstModel(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& appModel);
    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
      getModel(){return wrapperModel_->getAppModel();}

    /// Set solver via ParameterList solver name.
    virtual void setSolver(std::string solverName);
    /// Set solver via solver ParameterList.
    virtual void setSolver(
      Teuchos::RCP<Teuchos::ParameterList> solverPL=Teuchos::null);
    /// Set solver.
    virtual void setSolver(
      Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver);
    virtual Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > getSolver() const
      { return solver_; }

    /// Solve problem using x in-place.
    const Thyra::SolveStatus<Scalar> solveImplicitODE(
      const Teuchos::RCP<Thyra::VectorBase<Scalar> > & x);

    /// Set parameter so that the initial guess is set to zero (=True) or use last timestep (=False).
    virtual void setZeroInitialGuess(bool zIG)
      { stepperPL_->set<bool>("Zero Initial Guess", zIG); }
    virtual bool getZeroInitialGuess() const
      { return stepperPL_->get<bool>("Zero Initial Guess", false); }
    virtual Scalar getInitTimeStep(
        const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) const
      {return std::numeric_limits<Scalar>::max();}
  //@}

protected:

  Teuchos::RCP<Teuchos::ParameterList>                stepperPL_;
  Teuchos::RCP<WrapperModelEvaluator<Scalar> >        wrapperModel_;
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >   solver_;
};

} // namespace Tempus
#endif // Tempus_StepperImplicit_decl_hpp
