// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
 *
 * \ingroup Thyra_Nonlin_ME_solvers_grp
 */
template <class Scalar>
class LinearNonlinearSolver : public NonlinearSolverBase<Scalar> {
public:

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(RCP<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  RCP<Teuchos::ParameterList> getNonconstParameterList();
  /** \brief . */
  RCP<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  RCP<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  RCP<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** @name Overridden from NonlinearSolverBase */
  //@{

  /** \brief . */
  void setModel(
    const RCP<const ModelEvaluator<Scalar> > &model
    );
  /** \brief . */
  RCP<const ModelEvaluator<Scalar> > getModel() const;
  /** \brief . */
  SolveStatus<Scalar> solve(
    VectorBase<Scalar> *x,
    const SolveCriteria<Scalar> *solveCriteria,
    VectorBase<Scalar> *delta
    );
  /** \brief . */
  RCP<LinearOpWithSolveBase<Scalar> > get_nonconst_W(const bool forceUpToDate);
  /** \brief . */
  RCP<const LinearOpWithSolveBase<Scalar> > get_W() const;

  //@}

private:

  RCP<Teuchos::ParameterList> paramList_;
  RCP<const ModelEvaluator<Scalar> > model_;
  RCP<LinearOpWithSolveBase<Scalar> > J_;

};


/** \biref Nonmember constructor.
 *
 * \relates LinearNonlinearSolver
 */
template <class Scalar>
RCP<LinearNonlinearSolver<Scalar> > linearNonlinearSolver()
{
  return Teuchos::rcp(new LinearNonlinearSolver<Scalar>());
}


// ////////////////////////
// Defintions


// Overridden from Teuchos::ParameterListAcceptor


template<class Scalar>
void LinearNonlinearSolver<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList
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
RCP<Teuchos::ParameterList>
LinearNonlinearSolver<Scalar>::getNonconstParameterList()
{
  return paramList_;
}


template<class Scalar>
RCP<Teuchos::ParameterList>
LinearNonlinearSolver<Scalar>::unsetParameterList()
{
  RCP<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}


template<class Scalar>
RCP<const Teuchos::ParameterList>
LinearNonlinearSolver<Scalar>::getParameterList() const
{
  return paramList_;
}


template<class Scalar>
RCP<const Teuchos::ParameterList>
LinearNonlinearSolver<Scalar>::getValidParameters() const
{
  using Teuchos::setDoubleParameter; using Teuchos::setIntParameter;
  static RCP<const Teuchos::ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<Teuchos::ParameterList>
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
  const RCP<const ModelEvaluator<Scalar> > &model
  )
{
  TEST_FOR_EXCEPT(model.get()==NULL);
  model_ = model;
  J_ = Teuchos::null;
}


template <class Scalar>
RCP<const ModelEvaluator<Scalar> >
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
  
  const RCP<Teuchos::FancyOStream> out = this->getOStream();
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
  RCP<VectorBase<Scalar> >
    f = createMember(model_->get_f_space());
  if(out.get() && showTrace)
    *out << "\nEvaluating the model f and W ...\n";
  eval_f_W( *model_, *x,  &*f, &*J_ );

  // Solve the system: J*dx = -f
  RCP<VectorBase<Scalar> >
    dx = createMember(model_->get_x_space());
  if(out.get() && showTrace)
    *out << "\nSolving the system J*dx = -f ...\n";
  VOTSLOWSB J_outputTempState(J_,out,incrVerbLevel(verbLevel,-1));
  assign( dx.ptr(), ST::zero() );
  Thyra::SolveStatus<Scalar>
    linearSolveStatus = J_->solve(NOTRANS, *f, dx.ptr() );
  if(out.get() && showTrace)
    *out << "\nLinear solve status:\n" << linearSolveStatus;
  Vt_S( dx.ptr(), Scalar(-ST::one()) );
  if(out.get() && dumpAll)
    *out << "\ndx = " << Teuchos::describe(*dx,verbLevel);
  if (delta != NULL) {
    Thyra::assign( ptr(delta), *dx );
    if(out.get() && dumpAll)
      *out << "\ndelta = " << Teuchos::describe(*delta,verbLevel);
  }

  // Update the solution: x += dx
  Vp_V( ptr(x), *dx );
  if(out.get() && dumpAll)
    *out << "\nUpdated solution x = " << Teuchos::describe(*x,verbLevel);
  
  if(out.get() && showTrace)
    *out << "\nLeaving LinearNonlinearSolver::solve(...) ...\n";
  
  // Return default status
  return SolveStatus<Scalar>();

}


template <class Scalar>
RCP<LinearOpWithSolveBase<Scalar> >
LinearNonlinearSolver<Scalar>::get_nonconst_W(const bool forceUpToDate)
{
  if (forceUpToDate) {
    TEST_FOR_EXCEPT(forceUpToDate);
  }
  return J_;
}


template <class Scalar>
RCP<const LinearOpWithSolveBase<Scalar> >
LinearNonlinearSolver<Scalar>::get_W() const
{
  return J_;
}


} // namespace Thyra


#endif // THYRA_LINEAR_NONLINEAR_SOLVER_BASE_HPP
