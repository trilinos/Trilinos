// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperImplicit_hpp
#define Tempus_StepperImplicit_hpp

// Tempus
#include "Tempus_Stepper.hpp"
#include "Tempus_TimeDerivative.hpp"
// Thrya
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"

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
    /// Solve non-linear problem using x in-place.
    const Thyra::SolveStatus<Scalar> solveNonLinear(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > & model,
      Thyra::NonlinearSolverBase<Scalar> & solver,
      const Teuchos::RCP<Thyra::VectorBase<Scalar> > & x);

    /// Compute non-linear solve preserving x and returning solution.
    const Thyra::SolveStatus<Scalar> solveNonLinear(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > & model,
      Thyra::NonlinearSolverBase<Scalar> & solver,
      const Teuchos::RCP<Thyra::VectorBase<Scalar> > & x0,
      const Teuchos::RCP<Thyra::VectorBase<Scalar> > & solution_vec);
  //@}

};

template <typename Scalar>
const Thyra::SolveStatus<Scalar> StepperImplicit<Scalar>::solveNonLinear(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > & model,
  Thyra::NonlinearSolverBase<Scalar> & solver,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & x)
{
  // Set the model to use
  solver.setModel(model);

  // Solve
  const Thyra::SolveStatus<Scalar> solve_status = solver.solve(&*x, NULL, NULL);

  return solve_status;
}


template <typename Scalar>
const Thyra::SolveStatus<Scalar> StepperImplicit<Scalar>::solveNonLinear(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > & model,
  Thyra::NonlinearSolverBase<Scalar> & solver,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & x0,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & solution_vec)
{
  // Set the model to use
  solver.setModel(model);

  // Allocate solution vector
  Thyra::SolveCriteria<Scalar> solve_criteria; // this object is ignored
  if (solution_vec == Teuchos::null)
    solution_vec = Thyra::createMember(model->get_x_space());

  // Set initial guess
  Thyra::assign(solution_vec.ptr(),x0);

  // Solve
  const Thyra::SolveStatus<Scalar> solve_status =
    solver.solve(&*solution_vec, &solve_criteria, NULL);

  return solve_status;
}


} // namespace Tempus
#endif // Tempus_StepperImplicit_hpp
