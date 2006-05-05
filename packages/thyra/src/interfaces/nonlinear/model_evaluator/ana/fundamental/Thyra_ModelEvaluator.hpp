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

#ifndef THYRA_MODEL_EVALUATOR_HPP
#define THYRA_MODEL_EVALUATOR_HPP

#include "Thyra_ModelEvaluatorBase.hpp"

namespace Thyra {

/** \brief Pure abstract base interface for evaluating a stateless "model"
 * that can be mapped into a number of different types of problems.
 *
 * \section Thyra_ME_outline_sec Outline
 *
 * <ul>
 * <li>\ref Thyra_ME_intro_sec
 * <li>\ref Thyra_ME_problem_types_sec
 *     <ul>
 *     <li>\ref Thyra_ME_nonlinear_equations_sec
 *     <li>\ref Thyra_ME_explicit_ode_sec
 *     <li>\ref Thyra_ME_implicit_dae_sec
 *     <li>\ref Thyra_ME_unconstrained_optimization_sec
 *     <li>\ref Thyra_ME_equality_constrained_optimization_sec
 *     </ul>
 * <li>\ref Thyra_ME_derivatives_sec
 * <li>\ref Thyra_ME_dev_notes_sec
 * </ul>
 *
 * \section Thyra_ME_intro_sec Introduction
 *
 * The model represented by this interface is composed of the following
 * general functions:
 *
 * <ul>
 *
 * <li><b>State vector function:</b>
 *
 * <tt>(x_dot,x,{p(l)},t}) -> f</tt>
 *
 * <li><b>Auxiliary response vector functions:</b>
 *
 * <tt>(x_dot,x,{p(l)},t}) -> g(j)</tt>,
 *
 * for <tt>j=0...Ng-1</tt>
 *
 * </ul>
 *
 * given the general input variables/parameters:
 *
 * <ul>
 *
 * <li><b>State variables vector:</b>
 *
 * <tt>x</tt>
 *
 * <li><b>State variables derivative w.r.t. <tt>t</tt> vector:</b>
 *
 * <tt>x_dot</tt>
 *
 * <li><b>Auxiliary parameter vectors:</b>
 *
 * <tt>p(l)</tt>,
 *
 * for <tt>l=0...Np-1</tt>
 *
 * <li>Time point (or some other independent variable):
 *
 * <tt>t</tt>
 *
 * </ul>
 *
 * Above, the notation <tt>{p(l)}</tt> is shorthand for the set of parameter
 * vectors <tt>{ p(0), p(1), ..., p(Np-1) }</tt>.
 *
 * All of the above variables/parameters and functions are represented as
 * abstract <tt>Thyra::VectorBase</tt> objects.  The vector spaces associated
 * with these vector quantities are returned by <tt>get_space_x()</tt>,
 * <tt>get_p_space(int)</tt>, <tt>get_space_f()</tt>, and
 * <tt>get_g_space(int)</tt>.
 *
 * These functions all get computed at the same time in the function
 * <tt>evalModel()</tt>
 *
 * \section Thyra_ME_problem_types_sec Problem Types
 *
 * There are a number of different types of mathematical problems that can be
 * formulated using this interface.
 *
 * \subsection Thyra_ME_nonlinear_equations_sec Nonlinear Equations
 *
 
 \verbatim

  f(x) = 0.

 \endverbatim 
 
 * \subsection Thyra_ME_explicit_ode_sec Explicit ODEs
 *
 
 \verbatim

  x_dot = f(x,t).

 \endverbatim 

 * Here the argument <tt>t</tt> may or may not be accepted by <tt>*this</tt>.
 *
 * \subsection Thyra_ME_implicit_dae_sec Implicit ODEs or DAEs
 
 \verbatim

  f(x_dot,x,t) = 0.

 \endverbatim
 
 * Whether the problem is an implicit ODE or DAE is determined by the nature
 * of the derivative matrix <tt>d(f)/d(x_dot)</tt>:
 *
 * <ul>
 * <li> ODE: <tt>d(f)/d(x_dot)</tt> is full rank
 * <li> DAE: <tt>d(f)/d(x_dot)</tt> is not full rank
 * </ul>
 *
 * This information is given by the function ??? (ToDo: Add this function!)
 *
 * Here the argument <tt>t</tt> may or may not be accepted by <tt>*this</tt>.
 *
 * \subsection Thyra_ME_unconstrained_optimization_sec Unconstrained optimization
 
 \verbatim

  min g(x,{p(l)})

 \endverbatim

 * where the objective function <tt>g(x,{p(l)})</tt> is some aggregated
 * function built from some subset of the the auxiliary response functions
 * <tt>g(j)(x,{p(l)})</tt>, for <tt>j=1...Ng</tt>.
 *
 * \subsection Thyra_ME_equality_constrained_optimization_sec Equality constrained optimization
 
 \verbatim

  min g(x,{p(l)})

  s.t. f(x,{p(l)}) = 0

 \endverbatim 

 * where the objective function <tt>g(x,{p(l)})</tt> is some aggregated
 * function built from some subset of the the auxiliary response functions
 * <tt>g(j)(x,{p(l)})</tt>, for <tt>j=0...Ng-1</tt>.
 *
 * \subsection Thyra_ME_general_constrained_optimization_sec Equality constrained optimization
 
 \verbatim

  min g(x,{p(l)})

  s.t. f(x,{p(l)}) = 0
       r(x,{p(l)}) = 0
       hL <= h(x,{p(l)}) <= hU
       xL <= x <= xU
       pL(l) <= p(l) <= pU(l)

 \endverbatim 

 * where the objective function <tt>g(x,{p(l)})</tt> and the auxiliary
 * equality <tt>r(x,{p(l)})</tt> and inequality <tt>h(x,{p(l)})</tt>
 * constraint functions are aggregated functions built from some subset of the
 * the auxiliary response functions <tt>g(j)(x,{p(l)})</tt>, for
 * <tt>j=0...Ng-1</tt>.  The auxiliary response functions for a particualar
 * model can be interpreted in a wide variety of ways and can be mapped into a
 * number of different optimization problems.
 *
 * \section Thyra_ME_derivatives_sec Function derivatives and sensitivities
 *
 * A model can also optionally support various derivatives of the underlying
 * model functions.  The primary use for these derivatives is in the
 * computation of varyious types of sensitivities.  Specifically, direct and
 * adjoint sensitivities will be considered.
 *
 * To illustrate the issues involved, consider a single auxiliary parameter
 * <tt>p</tt> and a single auxiliary response function <tt>g</tt> of the form
 * <tt>(x,p) => g</tt>.  Assuming that <tt>(x,p) ==> f</tt> defines the state
 * equation <tt>f(x,p)=0</tt> and that <tt>D(f)/D(x)</tt> is full rank, then
 * <tt>f(x,p)=0</tt> defines the implicit function <tt>p ==> x(p)</tt>.  Given
 * this implicit function, the reduced auxiliary function is
 *
 * <tt>g_hat(p) = g(x(p),p)</tt>
 *
 * The reduced derivative <tt>D(g_hat)/D(p)</tt> is given as:
 *
 * <tt>D(g_hat)/D(p) = D(g)/D(x) * D(x)/D(p) + D(g)/D(p)</tt>
 *
 * where <tt>D(x)/D(p) = - [D(f)/D(x)]^{-1} * [D(f)/D(p)]</tt>
 *
 * Restated, the reduced derivative <tt>D(g_hat)/D(p)</tt> is given as:
 *
 * <tt>D(g_hat)/D(p) = - [D(g)/D(x)] * [D(f)/D(x)]^{-1} * [D(f)/D(p)] + D(g)/D(p)</tt>
 *
 * The reduced derivative <tt>D(g_hat)/D(p)</tt> can be computed using the
 * direct or the adjoint approahces.
 *
 * <ul>
 *
 * <li><b>State function derivatives</b>
 *
 *     <ul>
 *     
 *     <li><b>State variable derivatives</b>
 *
 *     <tt>W = alpha*D(f)/D(x_dot) + beta*D(f)/D(x)</tt>
 *
 *     This is derivative operator is a special object that is derived from
 *     the <tt>LinearOpWithSolveBase</tt> interface and therefore supports
 *     linear solves.  Objects of this type are created with the function
 *     <tt>create_W()</tt> before they are computed in <tt>evalModel()</tt>.
 *
 *     <li><b>State variable Taylor coefficients</b>
 *
 *     <tt>x_poly =</tt> \f$\sum_{i=0}^d x_i t^i\f$,
 *
 *     <tt>x_dot_poly =</tt> \f$\sum_{i=0}^d-1 (i+1)x_{i+1} t^i\f$,
 *
 *     <tt>f_poly =</tt> \f$\sum_{i=0}^{d-1} f_i t^i\f$.
 *
 *     If <tt>x_poly</tt> is a given polynomial of degree \f$d\f$ and
 *     <tt>x_dot_poly =</tt> \f$d(\f$<tt>x_poly</tt>\f$)/dt\f$, then <tt>f_poly
 *     = f(x_poly, x_dot_poly, t) +</tt> \f$O(t^{d})\f$ where \f$f_i =
 *     \frac{1}{k!} \frac{d^k}{dt^k} f(dx/dt ,x(t),t), i=0,\ldots,d-1\f$ are
 *     the Taylor series coefficient of \f$f\f$.  The polynomials
 *     <tt>x_poly</tt>, <tt>x_dot_poly</tt>, and <tt>f_poly</tt> are
 *     represented by <tt>Teuchos::Polynomial</tt> objects where each
 *     coefficient is a <tt>Thyra::VectorBase</tt> object.  The Taylor series
 *     coefficients of \f$f\f$ can easily be computed using automatic
 *     differentiation, but this is generally the only means of doing so.
 *     
 *     <li><b>Auxiliary parameter derivatives</b>
 *
 *     <tt>DfDp(l) = D(f)/D(p(l))</tt>, for <tt>l=0...Np-1</tt>
 *
 *     These are derivative objects that represent the derivative of the state
 *     function <tt>f</tt> with respect to the auxiliary parameters
 *     <tt>p(l)</tt>.  This object can either be represented and manipulated
 *     in a very abstract fashion through the <tt>LinearOpBase</tt> interface
 *     interface or as a multi-vector object.
 *
 *     </ul>
 *
 * <li><b>Auxiliary response function derivatives</b>
 *
 * </ul>
 *
 * \section Thyra_ME_dev_notes_sec Notes to subclass devleopers
 *
 * Subclass developers should consider deriving from either
 * <tt>StateFuncModelEvaluatorBase</tt> or
 * <tt>ModelEvaluatorDelegatorBase</tt>.  The base class
 * <tt>StateFuncModelEvaluatorBase</tt> makes it easy to create models
 * that start with the state function evaluation <tt>x -> f(x)</tt>.  The
 * <tt>ModelEvaluatorDelegatorBase</tt> base class makes it easy to
 * develop and maintain different types of decorator subcalsses.
 *
 *
 * ToDo: Finish Documentation!
 */
template<class Scalar>
class ModelEvaluator : public ModelEvaluatorBase {
public:

  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \name Basic inforamtion */
  //@{

  /** \brief Return the number of sets of auxiliary parameters.
   *
   * If this function returns 0, then there are no auxiliary parameters.
   */
  virtual int Np() const = 0;

  /** \brief Return the number of sets of auxiliary response functions.
   *
   * If this function returns 0, then there are no auxiliary response
   * functions.
   */
  virtual int Ng() const = 0;

  //@}

  /** \name Vector spaces */
  //@{

  /** \brief Return the vector space for the state variables <tt>x</tt> . */
  virtual Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > get_x_space() const = 0;

  /** \brief Return the vector space for the state function <tt>f(...)</tt>. */
  virtual Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > get_f_space() const = 0;

  /** \brief Return the vector space for the auxiliary parameters
   * <tt>p(l)</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->Np() > 0</tt>
   * <li><tt>0 <= l < this->Np()</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li> <tt>return.get()!=NULL</tt>
   * </ul>
   */
  virtual Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > get_p_space(int l) const = 0;

  /** \brief Return the vector space for the auxiliary response functions
   * <tt>g(j)</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->Ng() > 0</tt>
   * <li><tt>0 <= j < this->Ng()</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li> <tt>return.get()!=NULL</tt>
   * </ul>
   */
  virtual Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > get_g_space(int j) const = 0;

  //@}

  /** \name Initial guess and upper and lower bounds */
  //@{
  
  /** \brief Return the set of nominal values or initial guesses for the input
   * arguments.
   */
  virtual ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const = 0;
  
  /** \brief Return the set of lower bounds for the input arguments.
   */
  virtual ModelEvaluatorBase::InArgs<Scalar> getLowerBounds() const = 0;
  
  /** \brief Return the set of upper bounds for the input arguments.
   */
  virtual ModelEvaluatorBase::InArgs<Scalar> getUpperBounds() const = 0;

  //@}

  /** \name Factory functions for creating derivative objects */
  //@{

  /** \brief If supported, create a <tt>LinearOpWithSolveBase</tt> object for
   * <tt>W</tt> to be evaluated.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->createOutArgs().supports(OUT_ARG_W)==true</tt>
   * </ul>
   */
  virtual Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > create_W() const = 0;

  /** \brief If supported, create a <tt>LinearOpBase</tt> object for
   * <tt>W</tt> to be evaluated.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->createOutArgs().supports(OUT_ARG_W_op)==true</tt>
   * </ul>
   */
  virtual Teuchos::RefCountPtr<LinearOpBase<Scalar> > create_W_op() const = 0;

  /** \brief If supported, create a linear operator derivative object for
   * <tt>D(f)/D(p(l))</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->Np() > 0</tt>
   * <li><tt>0 <= l < this->Np()</tt>
   * <li><tt>outArgs.supports_DfDp(l).supports(DERIV_LINEAR_OP)==true</tt>,
   *     where <tt>outArgs = this->createOutArgs()</tt>
   * </ul>
   */
  virtual Teuchos::RefCountPtr<LinearOpBase<Scalar> > create_DfDp_op(int l) const = 0;

  // ToDo: Add functions for creating D(g(j))/D(x_dot) if needed!

  /** \brief If supported, create a linear operator derivative object for
   * <tt>D(g(j))/D(x)</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->Ng() > 0</tt>
   * <li><tt>0 <= j < this->Ng()</tt>
   * <li><tt>outArgs.supports_DgDx(j).supports(DERIV_LINEAR_OP)==true</tt>,
   *     where <tt>outArgs = this->createOutArgs()</tt>
   * </ul>
   */
  virtual Teuchos::RefCountPtr<LinearOpBase<Scalar> > create_DgDx_op(int j) const = 0;

  /** \brief If supported, create a linear operator derivative object for
   * <tt>D(g(j))/D(p(l))</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->Ng() > 0</tt>
   * <li><tt>this->Np() > 0</tt>
   * <li><tt>0 <= j < this->Ng()</tt>
   * <li><tt>0 <= l < this->Np()</tt>
   * <li><tt>outArgs.supports_DgDp(j,l).supports(DERIV_LINEAR_OP)==true</tt>,
   *     where <tt>outArgs = this->createOutArgs()</tt>
   * </ul>
   */
  virtual Teuchos::RefCountPtr<LinearOpBase<Scalar> > create_DgDp_op( int j, int l ) const = 0;
  
  //@}

  /** \name Computational functions */
  //@{

  /** \brief . */
  virtual ModelEvaluatorBase::InArgs<Scalar> createInArgs() const = 0;

  /** \brief . */
  virtual ModelEvaluatorBase::OutArgs<Scalar> createOutArgs() const = 0;

  /** \brief . */
  virtual void evalModel(
    const ModelEvaluatorBase::InArgs<Scalar>       &inArgs
    ,const ModelEvaluatorBase::OutArgs<Scalar>     &outArgs
    ) const = 0;

  //@}

  /** \name Reporting functions */
  //@{

  /** \brief Report the final point and whether the problem was considered
   * solved or not.
   *
   * ToDo: Add a PL to InArgs to allow extra data like Lagrange multipliers to
   * be passed back as well.
   */
  virtual void reportFinalPoint(
    const ModelEvaluatorBase::InArgs<Scalar>      &finalPoint
    ,const bool                                   wasSolved
    ) = 0;

  //@}

};

} // namespace Thyra


#endif // THYRA_MODEL_EVALUATOR_HPP
