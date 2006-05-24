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

namespace Thyra {

/** \brief Concrete nonlinear solver for linear equations :-)
 * 
 * ToDo: Finish documentation.
 */
template <class Scalar>
class LinearNonlinearSolver : public NonlinearSolverBase<Scalar> {
public:

  /** @name Overridden from NonlinearSolverBase */
  //@{

  /** \brief . */
  void setModel(
    const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> > &model
    );
  /** \brief . */
  Teuchos::RefCountPtr<const ModelEvaluator<Scalar> > getModel() const;
  /** \brief . */
  SolveStatus<Scalar> solve(
    VectorBase<Scalar>              *x
    ,const SolveCriteria<Scalar>    *solveCriteria
    ,VectorBase<Scalar>             *delta
    );
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > get_nonconst_W();
  /** \brief . */
  Teuchos::RefCountPtr<const LinearOpWithSolveBase<Scalar> > get_W() const;

  //@}

private:

  Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >    model_;
  Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> >   J_;

};

// ////////////////////////
// Defintions

// Overridden from NonlinearSolverBase

template <class Scalar>
void LinearNonlinearSolver<Scalar>::setModel(
  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> > &model
  )
{
  TEST_FOR_EXCEPT(model.get()==NULL);
  model_ = model;
  J_ = Teuchos::null;
}

template <class Scalar>
Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >
LinearNonlinearSolver<Scalar>::getModel() const
{
  return model_;
}

template <class Scalar>
SolveStatus<Scalar> LinearNonlinearSolver<Scalar>::solve(
  VectorBase<Scalar>             *x
  ,const SolveCriteria<Scalar>   *solveCriteria
  ,VectorBase<Scalar>            *delta = NULL
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  TEST_FOR_EXCEPT(solveCriteria!=NULL); // ToDo: Pass to linear solver?
  // Compute the Jacobian and the residual at the input point!
  if(!J_.get())                               J_ = model_->create_W();
  Teuchos::RefCountPtr<VectorBase<Scalar> >   f = createMember(model_->get_f_space());
  eval_f_W( *model_, *x, ST::one(), &*f, &*J_ );
  // Solve the system: J*m_dx = f
  Teuchos::RefCountPtr<VectorBase<Scalar> > m_dx = createMember(model_->get_x_space());
  Thyra::solve( *J_, NOTRANS, *f, &*m_dx );
  // Set the solution: x = x - m_dx
  Vt_S( &*m_dx, Scalar(-ST::one()) );
  Vp_V( x, *m_dx );
  if (delta != NULL) assign( delta, *m_dx );
  // Return default status
  return SolveStatus<Scalar>();
}

template <class Scalar>
Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> >
LinearNonlinearSolver<Scalar>::get_nonconst_W()
{
  return J_;
}

template <class Scalar>
Teuchos::RefCountPtr<const LinearOpWithSolveBase<Scalar> >
LinearNonlinearSolver<Scalar>::get_W() const
{
  return J_;
}

} // namespace Thyra

#endif // THYRA_LINEAR_NONLINEAR_SOLVER_BASE_HPP
