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

#ifndef THYRA_NONLINEAR_SOLVER_BASE_HPP
#define THYRA_NONLINEAR_SOLVER_BASE_HPP

#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"


namespace Thyra {


/** \brief Base class for all nonlinear equation solvers.
 *
 * <b>Warning!</b> This interface is highly experimental and general
 * developers should not even consider using it in any way if there is any
 * expectation of code stability!
 * 
 * ToDo: Finish documentation.
 *
 * ToDo:<ul>
 * <li> Add a supportsSolveCritiera(...) function that informs client if
 *      a solve criteria is supported or not.
 * <li> Add a getDefaultSolveCriteria() function to return the default
 *      solve criteria.  This will make it very easy to set additional parameters
 *      like the maximum number of iterations etc. in the extraParameters
 *      parameter list field.
 * </ul>
 *
 * \ingroup Thyra_Nonlinear_solver_interfaces_code_grp
 */
template <class Scalar>
class NonlinearSolverBase
  : virtual public Teuchos::Describable
  , virtual public Teuchos::VerboseObject<NonlinearSolverBase<Scalar> >
  , virtual public Teuchos::ParameterListAcceptor
{
public:
  
  /** @name Pure virtual functions that must be overridden in subclasses */
  //@{

  /** \brief Set the model that defines the nonlinear equations.
   *
   * After the model is set, only the residual <tt>f</tt> can change between
   * solves and not the structure of the Jacobian <tt>W</tt>.  If a more
   * significant change to <tt>*model</tt> occurs, then this function must be
   * called again to reset the model and reinitialize.
   */
  virtual void setModel(
    const RCP<const ModelEvaluator<Scalar> > &model
    ) = 0;
  
  /** \brief Get the model that defines the nonlinear equations. */
  virtual RCP<const ModelEvaluator<Scalar> > getModel() const = 0;
  
  /** \brief Solve a set of nonlinear equations from a given starting
   * point.
   *
   * \param x [in/out] On input, <tt>*x</tt> contains the guess for the
   * solution of <tt>f(x)=0</tt> as defined by <tt>*this->getModel()</tt> On
   * output, <tt>*x</tt> will contain the best estimate of the solution.
   *
   * \param solveCriteria [in] If <tt>solveCriteria!=NULL</tt> then
   * <tt>*solveCriteria</tt> will contain the criteria for the nonlinear
   * solve.
   *
   * \return The returned <tt>SolveStatus</tt> object gives the status of the
   * returned solution <tt>*x</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li>this->getModel().get()!=NULL</tt>
   * <li>x->space()->isCompatible(*this->getModel()->get_x_space())==true</tt>
   * </ul>
   */
  virtual SolveStatus<Scalar> solve(
    VectorBase<Scalar> *x,
    const SolveCriteria<Scalar> *solveCriteria = NULL,
    VectorBase<Scalar> *delta = NULL
    ) = 0;
  
  //@}
  
  /** @name Virtual functions with default implementation */
  //@{

  /** \brief Return if this solver object supports cloning or not.
   *
   * The default implementation returns false.
   */
  virtual bool supportsCloning() const;

  /** \brief Clone the solver algorithm if supported.
   *
   * <b>Postconditions:</b><ul>
   * <li>[<tt>supportsCloning()==true</tt>] <tt>returnVal != Teuchos::null</tt>
   * <li>[<tt>supportsCloning()==false</tt>] <tt>returnVal == Teuchos::null</tt>
   * </ul>
   *
   * Note that cloning a nonlinear solver in this case does not imply that the
   * Jacobian state will be copied as well, shallow or deep.  Instead, here
   * cloning means to just clone the solver algorithm and it will do a
   * showllow of the model as well if a model is set.  Since the model is
   * stateless, this is okay.  Therefore, do not assume that the state of
   * <tt>*returnValue</tt> is exactly the same as the state of <tt>*this</tt>.
   * You have been warned!
   * 
   * The default implementation returns <tt>Teuchos::null</tt> which is
   * consistent with the default implementation of <tt>supportsCloning()</tt>.
   * If this function is overridden in a base class to support cloning, then
   * <tt>supportsCloning()</tt> must be overridden to return <tt>true</tt>.
   */
  virtual RCP<NonlinearSolverBase<Scalar> > cloneNonlinearSolver() const;
  
  /** \brief Return the current value of the solution <tt>x</tt> as computed
   * in the last <tt>solve()</tt> operation if supported.
   *
   * The default implementation returns <tt>return.get()==NULL</tt>.
   */
  virtual RCP<const VectorBase<Scalar> > get_current_x() const;

  /** \brief Returns <tt>true</tt> if <tt>*get_W()</tt> is current with
   * respect to <tt>*get_current_x()</tt>.
   *
   * The default implementation returns <tt>false</tt>.
   */
  virtual bool is_W_current() const;

  /** \brief Get a nonconst RCP to the Jacobian if available.
   *
   * \param forceUpToDate [in] If <tt>true</tt>, then the underlying W will be
   * forced to be up to date w.r.t the current x if a Jacobian exists.
   *
   * <b>Postconditions:</b><ul>
   * <li>[<tt>forceUpToDate==true</tt>] <tt>this->is_W_current() == true</tt>
   * </ul>
   *
   * Through this the RCP returned from this function, a client can change the
   * <tt>W</tt> object held internally.  If the object gets changed the client
   * should call <tt>set_W_is_current(false)</tt>.
   *
   * The default implementation returns <tt>return.get()==NULL</tt>.
   */
  virtual RCP<LinearOpWithSolveBase<Scalar> >
  get_nonconst_W( const bool forceUpToDate = false );

  /** \brief Get a const RCP to the Jacobian if available.
   *
   * Through this interface the client should not change the object
   * <tt>W</tt>.
   *
   * The default implementation returns <tt>return.get()==NULL</tt>.
   */
  virtual RCP<const LinearOpWithSolveBase<Scalar> > get_W() const;

  /** \brief Set if <tt>*get_W()</tt> is current with respect to
   * <tt>*get_current_x()</tt>.
   *
   * <b>Preconditions:</b></ul>
   * <li> <tt>this->get_W().get()!=NULL</tt>
   * </ul>
   *
   * The default implementation throwns an exception.
   */
  virtual void set_W_is_current(bool W_is_current);

  //@}

private:
  
  // Not defined and not to be called
  NonlinearSolverBase<Scalar>&
  operator=(const NonlinearSolverBase<Scalar>&);

};


/** \brief . 
 *
 * \relates NonlinearSolverBase
 */
template <class Scalar>
const SolveStatus<Scalar> solve(
  NonlinearSolverBase<Scalar> &nonlinearSolver,
    VectorBase<Scalar> *x,
    const SolveCriteria<Scalar> *solveCriteria = NULL,
    VectorBase<Scalar> *delta = NULL
    )
{
  return nonlinearSolver.solve(x,solveCriteria,delta);
}


// ///////////////////////////////
// Implementations

template <class Scalar>
bool NonlinearSolverBase<Scalar>::supportsCloning() const
{
  return false;
}

template <class Scalar>
RCP<NonlinearSolverBase<Scalar> >
NonlinearSolverBase<Scalar>::cloneNonlinearSolver() const
{
  return Teuchos::null;
}

template <class Scalar>
RCP<const VectorBase<Scalar> >
NonlinearSolverBase<Scalar>::get_current_x() const
{
  return Teuchos::null;
}

template <class Scalar>
bool NonlinearSolverBase<Scalar>::is_W_current() const
{
  return false;
}

template <class Scalar>
RCP<LinearOpWithSolveBase<Scalar> >
NonlinearSolverBase<Scalar>::get_nonconst_W(const bool forceUpToDate)
{
  return Teuchos::null;
}

template <class Scalar>
RCP<const LinearOpWithSolveBase<Scalar> >
NonlinearSolverBase<Scalar>::get_W() const
{
  return Teuchos::null;
}

template <class Scalar>
void NonlinearSolverBase<Scalar>::set_W_is_current(bool W_is_current)
{
  TEST_FOR_EXCEPTION(
    true, std::logic_error,
    "Error, the subclass object described as " << this->description() << " did not"
    " override this function!"
    );
}


} // namespace Thyra


#endif // THYRA_NONLINEAR_SOLVER_BASE_HPP
