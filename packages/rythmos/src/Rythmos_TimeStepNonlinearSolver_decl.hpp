//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER


#ifndef RYTHMOS_TIME_STEP_NONLINEAR_SOLVER_DECL_HPP
#define RYTHMOS_TIME_STEP_NONLINEAR_SOLVER_DECL_HPP

#include "Rythmos_Types.hpp"
#include "Thyra_NonlinearSolverBase.hpp"

namespace Rythmos {


/** \brief Simple undampended Newton solver designed to solve time step
 * equations in accurate times-tepping methods.
 * 
 * ToDo: Finish documentation.
 *
 * 2007/05/18: rabartl: ToDo: Derive NonlinearSolverBase from
 * ParameterListAcceptor and accept options through a validated
 * parameter list!  Then remove these STANDARD_MEMBER_COMPOSITION_MEMBERS()
 * macros.
 */
template <class Scalar>
class TimeStepNonlinearSolver : public Thyra::NonlinearSolverBase<Scalar> {
public:

  /** \brief. */
  typedef Teuchos::ScalarTraits<Scalar> ST;
  /** \brief. */
  typedef typename ST::magnitudeType ScalarMag;
  /** \brief. */
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;

  /** @name Constructors/Intializers/Misc */
  //@{

  /** \brief Sets parameter defaults . */
  TimeStepNonlinearSolver();

  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);
  /** \brief . */
  RCP<ParameterList> getNonconstParameterList();
  /** \brief . */
  RCP<ParameterList> unsetParameterList();
  /** \brief . */
  RCP<const ParameterList> getParameterList() const;
  /** \brief . */
  RCP<const ParameterList> getValidParameters() const;

  //@}

  /** @name Overridden from NonlinearSolverBase */
  //@{

  /** \brief . */
  void setModel(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &model
    );
  /** \brief . */
  RCP<const Thyra::ModelEvaluator<Scalar> > getModel() const;
  /** \brief . */
  Thyra::SolveStatus<Scalar> solve(
    Thyra::VectorBase<Scalar> *x,
    const Thyra::SolveCriteria<Scalar> *solveCriteria,
    Thyra::VectorBase<Scalar> *delta = NULL
    );
  /** \brief . */
  bool supportsCloning() const;
  /** \brief . */
  RCP<Thyra::NonlinearSolverBase<Scalar> >
  cloneNonlinearSolver() const;  
  /** \brief . */
  RCP<const Thyra::VectorBase<Scalar> > get_current_x() const;
  /** \brief . */
  bool is_W_current() const;
  /** \brief . */
  RCP<Thyra::LinearOpWithSolveBase<Scalar> >
  get_nonconst_W(const bool forceUpToDate);
  /** \brief . */
  RCP<const Thyra::LinearOpWithSolveBase<Scalar> > get_W() const;
  /** \brief . */
  void set_W_is_current(bool W_is_current);

  //@}

private:

  // private object data members

  RCP<ParameterList> paramList_;
  RCP<const Thyra::ModelEvaluator<Scalar> > model_;
  RCP<Thyra::LinearOpWithSolveBase<Scalar> > J_;
  RCP<Thyra::VectorBase<Scalar> > current_x_;
  bool J_is_current_;

  double defaultTol_;
  int defaultMaxIters_;
  double nonlinearSafetyFactor_;
  double linearSafetyFactor_;
  double RMinFraction_;
  bool throwOnLinearSolveFailure_;

  // static class data members

  static const std::string DefaultTol_name_;
  static const double DefaultTol_default_;

  static const std::string DefaultMaxIters_name_;
  static const int DefaultMaxIters_default_;

  static const std::string NonlinearSafetyFactor_name_;
  static const double NonlinearSafetyFactor_default_;

  static const std::string LinearSafetyFactor_name_;
  static const double LinearSafetyFactor_default_;

  static const std::string RMinFraction_name_;
  static const double RMinFraction_default_;

  static const std::string ThrownOnLinearSolveFailure_name_;
  static const bool ThrownOnLinearSolveFailure_default_;

};


/** \brief Nonmember constructor.
 *
 * \relates TimeStepNonlinearSolver
 */
template <class Scalar>
RCP<TimeStepNonlinearSolver<Scalar> > timeStepNonlinearSolver();


/** \brief Nonmember constructor.
 *
 * \relates TimeStepNonlinearSolver
 */
template <class Scalar>
RCP<TimeStepNonlinearSolver<Scalar> >
timeStepNonlinearSolver(const RCP<ParameterList> &pl);

} // namespace Rythmos


#endif // RYTHMOS_TIME_STEP_NONLINEAR_SOLVER_DECL_HPP
