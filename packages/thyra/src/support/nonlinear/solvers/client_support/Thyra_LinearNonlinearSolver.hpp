// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_LINEAR_NONLINEAR_SOLVER_BASE_HPP
#define THYRA_LINEAR_NONLINEAR_SOLVER_BASE_HPP


#include "Thyra_NonlinearSolverBase.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_as.hpp"


namespace Thyra {


/** \brief Concrete nonlinear solver for linear equations.
 *
 * This class basically implements a Newton method with one iteration and
 * never checks the final tolerence.  Otherwise, it is identical to a Newton
 * method with one iteration.
 */
template <class Scalar>
class LinearNonlinearSolver : public NonlinearSolverBase<Scalar> {
public:

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> getParameterList();
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** @name Overridden from NonlinearSolverBase */
  //@{

  /** \brief . */
  void setModel(
    const Teuchos::RCP<const ModelEvaluator<Scalar> > &model
    );
  /** \brief . */
  Teuchos::RCP<const ModelEvaluator<Scalar> > getModel() const;
  /** \brief . */
  SolveStatus<Scalar> solve(
    VectorBase<Scalar> *x,
    const SolveCriteria<Scalar> *solveCriteria,
    VectorBase<Scalar> *delta
    );
  /** \brief . */
  Teuchos::RCP<LinearOpWithSolveBase<Scalar> > get_nonconst_W();
  /** \brief . */
  Teuchos::RCP<const LinearOpWithSolveBase<Scalar> > get_W() const;

  //@}

private:

  Teuchos::RCP<Teuchos::ParameterList> paramList_;
  Teuchos::RCP<const ModelEvaluator<Scalar> > model_;
  Teuchos::RCP<LinearOpWithSolveBase<Scalar> > J_;

};


// ////////////////////////
// Defintions


// Overridden from Teuchos::ParameterListAcceptor


template<class Scalar>
void LinearNonlinearSolver<Scalar>::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& paramList
  )
{
  using Teuchos::get;
  TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParametersAndSetDefaults(*getValidParameters(),0);
  paramList_ = paramList;
  // ToDo: Accept some parameters if this makes sense!
  Teuchos::readVerboseObjectSublist(&*paramList_,this);
#ifdef TEUCHOS_DEBUG
  paramList_->validateParameters(*getValidParameters(),0);
#endif // TEUCHOS_DEBUG
}


template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
LinearNonlinearSolver<Scalar>::getParameterList()
{
  return paramList_;
}


template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
LinearNonlinearSolver<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
LinearNonlinearSolver<Scalar>::getParameterList() const
{
  return paramList_;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
LinearNonlinearSolver<Scalar>::getValidParameters() const
{
  using Teuchos::setDoubleParameter; using Teuchos::setIntParameter;
  static Teuchos::RCP<const Teuchos::ParameterList> validPL;
  if (is_null(validPL)) {
    Teuchos::RCP<Teuchos::ParameterList>
      pl = Teuchos::parameterList();
    // ToDo: Set up some parameters when needed!
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return validPL;
}


// Overridden from NonlinearSolverBase


template <class Scalar>
void LinearNonlinearSolver<Scalar>::setModel(
  const Teuchos::RCP<const ModelEvaluator<Scalar> > &model
  )
{
  TEST_FOR_EXCEPT(model.get()==NULL);
  model_ = model;
  J_ = Teuchos::null;
}


template <class Scalar>
Teuchos::RCP<const ModelEvaluator<Scalar> >
LinearNonlinearSolver<Scalar>::getModel() const
{
  return model_;
}


template <class Scalar>
SolveStatus<Scalar> LinearNonlinearSolver<Scalar>::solve(
  VectorBase<Scalar> *x,
  const SolveCriteria<Scalar> *solveCriteria,
  VectorBase<Scalar> *delta
  )
{

  using std::endl;
  using Teuchos::incrVerbLevel;
  using Teuchos::describe;
  using Teuchos::as;
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::OSTab;
  using Teuchos::getFancyOStream;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::VerboseObjectTempState<MEB> VOTSME;
  typedef Thyra::LinearOpWithSolveBase<Scalar> LOWSB;
  typedef Teuchos::VerboseObjectTempState<LOWSB> VOTSLOWSB;

#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(0==x);
  THYRA_ASSERT_VEC_SPACES(
    "TimeStepNonlinearSolver<Scalar>::solve(...)",
    *x->space(),*model_->get_x_space() );
  TEST_FOR_EXCEPT(
    0!=solveCriteria && "ToDo: Support passed in solve criteria!" );
#endif
  
  const Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  const bool showTrace = (as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW));
  const bool dumpAll = (as<int>(verbLevel) == as<int>(Teuchos::VERB_EXTREME)); 
  TEUCHOS_OSTAB;
  VOTSME stateModel_outputTempState(model_,out,incrVerbLevel(verbLevel,-1));
  if(out.get() && showTrace)
    *out
      << "\nEntering LinearNonlinearSolver::solve(...) ...\n"
      << "\nmodel = " << describe(*model_,verbLevel);

  if(out.get() && dumpAll) {
    *out << "\nInitial guess:\n";
    *out << "\nx = " << *x;
  }

  // Compute the Jacobian and the residual at the input point!
  if(!J_.get()) J_ = model_->create_W();
  Teuchos::RCP<VectorBase<Scalar> >
    f = createMember(model_->get_f_space());
  if(out.get() && showTrace)
    *out << "\nEvaluating the model f and W ...\n";
  eval_f_W( *model_, *x,  &*f, &*J_ );

  // Solve the system: J*dx = -f
  Teuchos::RCP<VectorBase<Scalar> >
    dx = createMember(model_->get_x_space());
  if(out.get() && showTrace)
    *out << "\nSolving the system J*dx = -f ...\n";
  VOTSLOWSB J_outputTempState(J_,out,incrVerbLevel(verbLevel,-1));
  assign( &*dx, ST::zero() );
  Thyra::SolveStatus<Scalar>
    linearSolveStatus = Thyra::solve( *J_, NOTRANS, *f, &*dx );
  if(out.get() && showTrace)
    *out << "\nLinear solve status:\n" << linearSolveStatus;
  Vt_S( &*dx, Scalar(-ST::one()) );
  if(out.get() && dumpAll)
    *out << "\ndx = " << describe(*dx,verbLevel);
  if (delta != NULL) {
    Thyra::assign( delta, *dx );
    if(out.get() && dumpAll)
      *out << "\ndelta = " << describe(*delta,verbLevel);
  }

  // Update the solution: x += dx
  Vp_V( x, *dx );
  if(out.get() && dumpAll)
    *out << "\nUpdated solution x = " << describe(*x,verbLevel);
  
  if(out.get() && showTrace)
    *out << "\nLeaving LinearNonlinearSolver::solve(...) ...\n";
  
  // Return default status
  return SolveStatus<Scalar>();

}


template <class Scalar>
Teuchos::RCP<LinearOpWithSolveBase<Scalar> >
LinearNonlinearSolver<Scalar>::get_nonconst_W()
{
  return J_;
}


template <class Scalar>
Teuchos::RCP<const LinearOpWithSolveBase<Scalar> >
LinearNonlinearSolver<Scalar>::get_W() const
{
  return J_;
}


} // namespace Thyra


#endif // THYRA_LINEAR_NONLINEAR_SOLVER_BASE_HPP
