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
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"


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
 * <li>\ref Thyra_ME_nominal_values_sec
 * <li>\ref Thyra_ME_bounds_sec
 * <li>\ref Thyra_ME_checking_sec
 * <li>\ref Thyra_ME_dev_notes_sec
 * </ul>
 *
 * \section Thyra_ME_intro_sec Introduction
 *
 * This interface defines a very loosely mathematically typed interface to a
 * very wide variety of simulation-based models that can support a very wide
 * range of simulation-based numerical algorithms.
 *
 * For the most part, a model represented by this interface is composed:
 *
 * <ul>
 *
 * <li><b>State vector function:</b>
 *
 * <tt>(x_dot,x,{p(l)},t,...}) -> f <: f_space</tt>
 *
 * <li><b>Auxiliary response vector functions:</b>
 *
 * <tt>(x_dot,x,{p(l)},t,...}) -> g(j) <: g_space(j)</tt>, for <tt>j=0...Ng-1</tt>
 *
 * <li><b>Other outputs:</b>
 *
 * A model can compute other objects as well (see derivatives below).
 *
 * </ul>
 *
 * given the general input variables/parameters:
 *
 * <ul>
 *
 * <li><b>State variables vector:</b>
 *
 * <tt>x <: x_space</tt>
 *
 * <li><b>State variables derivative w.r.t. <tt>t</tt> vector:</b>
 *
 * <tt>x_dot <: x_space</tt>
 *
 * <li><b>Auxiliary parameter vectors:</b>
 *
 * <tt>p(l) <: p_space(l)</tt>, for <tt>l=0...Np-1</tt>
 *
 * <li><b>Time point (or some other independent variable):</b>
 *
 * <tt>t <: Scalar</tt>
 *
 * <li><b>Other inputs:</b>
 *
 * A model can accept additional input objects as well (see below).
 *
 * </ul>
 *
 * where <tt>x_space <: RE^n_x</tt>, <tt>f_space <: RE^n_x</tt>,
 * <tt>p_space(l) <: RE^n_p_l</tt> (for <tt>l=0...Np-1</tt>), and
 * <tt>g_space(j) <: RE^n_g_j</tt> (for <tt>j=0...Ng-1</tt>) are
 * <tt>%Thyra</tt> vector spaces of the given dimensions.
 *
 * Above, the notation <tt>{p(l)}</tt> is shorthand for the set of parameter
 * vectors <tt>{ p(0), p(1), ..., p(Np-1) }</tt>.
 *
 * All of the above variables/parameters and functions are represented as
 * abstract <tt>Thyra::VectorBase</tt> objects.  The vector spaces associated
 * with these vector quantities are returned by <tt>get_x_space()</tt>,
 * <tt>get_p_space()</tt>, <tt>get_f_space()</tt>, and
 * <tt>get_g_space()</tt>.
 *
 * All of the input variables/parameters are specified as a
 * <tt>ModelEvaluatorBase::InArgs</tt> object, all functions to be computed
 * are specified as a <tt>ModelEvaluatorBase::OutArgs</tt> object, and
 * evaluations of all function at a single set of variable values is performed in a single
 * call to <tt>evalModel()</tt>.
 *
 * A particular <tt>%ModelEvaluator</tt> subclass object can support any
 * subset of these inputs and outputs and it is up to the client to map these
 * variables/parameters and functions into abstract mathematical problems.
 * Some of the different types of abstract mathematical problems that can be
 * represented through this interface are given in the next section.
 *
 * This interface can also support the computation of various derivatives of
 * these functions w.r.t. the input arguments (see the section \ref
 * Thyra_ME_derivatives_sec below).
 *
 * \section Thyra_ME_problem_types_sec Examples of Abstract Problem Types
 *
 * There are a number of different types of mathematical problems that can be
 * formulated using this interface.  In the following subsections, a few
 * different examples of specific abstract problems types are given.
 *
 * \subsection Thyra_ME_nonlinear_equations_sec Nonlinear Equations
 *
 
 \verbatim
    f(x) = 0
 \endverbatim 

 * Here it is assumed that <tt>D(f)/D(x)</tt> is nonsingular in general but
 * this is not strictly required.  If <tt>W=D(f)/D(x)</tt> is supported, the
 * nature of <tt>D(f)/D(x)</tt> may be given by
 * <tt>this->createOutArgs().get_W_properties()</tt>.
 
 * \subsection Thyra_ME_explicit_ode_sec Explicit ODEs
 *
 
 \verbatim
    x_dot = f(x,t)
 \endverbatim 

 * Here it is assumed that <tt>D(f)/D(x)</tt> is nonsingular in general but
 * this is not strictly required.  If <tt>W=D(f)/D(x)</tt> is supported, the
 * nature of <tt>D(f)/D(x)</tt> may be given by
 * <tt>this->createOutArgs().get_W_properties()</tt>.

 * Above, the argument <tt>t</tt> may or may not be accepted by the model
 * (i.e. <tt>createInArgs().supports(IN_ARG_t)</tt> may return
 * <tt>false</tt>).
 *
 * \subsection Thyra_ME_implicit_dae_sec Implicit ODEs or DAEs
 
 \verbatim
    f(x_dot,x,t) = 0
 \endverbatim

 * Here it is assumed that <tt>D(f)/D(x)</tt> is nonsingular in general but
 * this is not strictly required.
 *
 * The problem is either an implicit ODE or DAE depending on the nature of the
 * derivative matrix <tt>D(f)/D(x_dot)</tt>:
 *
 * <ul>
 * <li> ODE: <tt>D(f)/D(x_dot)</tt> is full rank
 * <li> DAE: <tt>D(f)/D(x_dot)</tt> is <em>not</em> full rank
 * </ul>
 *
 * If supported, the nature of <tt>W=alpha*D(f)/D(x_dot)+beta*D(f)/D(x)</tt>
 * may be given by <tt>this->createOutArgs().get_W_properties()</tt>.
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
 * In general, it would be assumed that the Hessian <tt>D^2(g)/D(x^2)</tt> is
 * symmetric semidefinite but this is not strictly required.
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
 * Here it is assumed that <tt>D(f)/D(x)</tt> is nonsingular in general but
 * this is not strictly required.  If <tt>W=D(f)/D(x)</tt> is supported, the
 * nature of <tt>D(f)/D(x)</tt> may be given by
 * <tt>this->createOutArgs().get_W_properties()</tt>.
 *
 * \subsection Thyra_ME_general_constrained_optimization_sec General constrained optimization
 
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
 * <tt>j=0...Ng-1</tt>.  The auxiliary response functions for a particular
 * model can be interpreted in a wide variety of ways and can be mapped into a
 * number of different optimization problems.
 *
 * Here it is assumed that <tt>D(f)/D(x)</tt> is nonsingular in general but
 * this is not strictly required.  If <tt>W=D(f)/D(x)</tt> is supported, the
 * nature of <tt>D(f)/D(x)</tt> may be given by
 * <tt>this->createOutArgs().get_W_properties()</tt>.
 *
 * \section Thyra_ME_derivatives_sec Function derivatives and sensitivities
 *
 * A model can also optionally support the computation of various derivatives
 * of the underlying model functions.  The primary use for these derivatives
 * is in the computation of various types of sensitivities.  Specifically,
 * direct and adjoint sensitivities will be considered.
 *
 * To illustrate the issues involved, consider a single auxiliary parameter
 * <tt>p</tt> and a single auxiliary response function <tt>g</tt> of the form
 * <tt>(x,p) => g</tt>.  Assuming that <tt>(x,p) ==> f</tt> defines the state
 * equation <tt>f(x,p)=0</tt> and that <tt>D(f)/D(x)</tt> is full rank, then
 * <tt>f(x,p)=0</tt> defines the implicit function <tt>p ==> x(p)</tt>.  Given
 * this implicit function, the reduced auxiliary function is
 
 \verbatim
    g_hat(p) = g(x(p),p)
 \endverbatim
 
 * The reduced derivative <tt>D(g_hat)/D(p)</tt> is given as:
 
 \verbatim
    D(g_hat)/D(p) = D(g)/D(x) * D(x)/D(p) + D(g)/D(p)
 \endverbatim
 
 * where

 \verbatim
    D(x)/D(p) = - [D(f)/D(x)]^{-1} * [D(f)/D(p)]
 \endverbatim
 
 * Restated, the reduced derivative <tt>D(g_hat)/D(p)</tt> is given as:
 
 \verbatim
    D(g_hat)/D(p) = - [D(g)/D(x)] * [D(f)/D(x)]^{-1} * [D(f)/D(p)] + D(g)/D(p)
 \endverbatim
 
 * The reduced derivative <tt>D(g_hat)/D(p)</tt> can be computed using the
 * direct or the adjoint approaches.
 *
 * The <b>direct sensitivity approach</b> first solves for

 \verbatim
    D(x)/D(p) = - [D(f)/D(x)]^{-1} * [D(f)/D(p)]
 \endverbatim

 * explicitly, then computes

 \verbatim
    D(g_hat)/D(p) = D(g)/D(x) * D(x)/D(p) + D(g)/D(p).
 \endverbatim

 * In this case, <tt>D(f)/D(p)</tt> is needed as a multivector since it forms
 * the RHS for a set of linear equations.  However, only the action of
 * <tt>D(g)/D(x)</tt> on the multivector <tt>D(x)/D(p)</tt> is needed and
 * therefore <tt>D(g)/D(x)</tt> can be returned as only a linear operator
 * (i.e. a <tt>LinearOpBase</tt> object).  Note that in Thyra a multivector
 * <em>is a</em> linear operator and therefore every derivative object
 * returned as a multivector automatically implements the forward and adjoint
 * linear operators for the derivative operator.
 *
 * The final derivative <tt>D(g)/D(p)</tt> should be returned as a multivector
 * that can be added to the multivector <tt>D(g)/D(x)*D(x)/D(p)</tt>.
 *
 * The <b>adjoint sensitivity approach</b> computes

 \verbatim
    D(g_hat)/D(p)^T =  [D(f)/D(p)]^T * ( - [D(f)/D(x)]^{-T} * [D(g)/D(x)]^T ) + [D(g)/D(p)]^T
 \endverbatim

 * by first solving the adjoint system 

 \verbatim
    Lambda = - [D(f)/D(x)]^{-T} * [D(g)/D(x)]^T
 \endverbatim

 * for the multivector <tt>Lambda</tt> and then computes

 \verbatim
    D(g_hat)/D(p)^T =  [D(f)/D(p)]^T * Lambda + [D(g)/D(p)]^T.
 \endverbatim

 * In this case, <tt>[D(g)/D(x)]^T</tt> is needed as an explicit multivector
 * since it forms the RHS for the linear adjoint equations.  Also, only the
 * adjoint operator application[D(f)/D(p)]^T is needed.  And in this case, the
 * multivector form of the adjoint <tt>[D(g)/D(p)]^T</tt> is required.
 *
 * As demonstrated above, general derivative objects (e.g. <tt>D(f)/D(p)</tt>,
 * <tt>D(g)/D(x)</tt>, and <tt>D(g)/D(p)</tt>) may be needed as either only a
 * linear operator (where it's forward or adjoint application is required) or
 * as a multivector for its forward or adjoint forms.  A derivative
 * <tt>D(h)/D(z)</tt> for some function <tt>h(z)</tt> can be supported in any
 * of the following forms:
 *
 * <ul>
 *
 * <li> <tt><b>D(h)/D(z)</b></tt> as a <tt>LinearOpBase</tt> object where the
 * forward and/or adjoint operator applications are supported
 *
 * <li> <tt><b>D(h)/D(z)</b></tt> as a <tt>MultiVectorBase</tt> object where
 * each column in the multivector <tt>i</tt> represents <tt>D(h)/D(z(i))</tt>,
 * the derivatives for all of the functions <tt>h</tt> for the single variable
 * <tt>z(i)</tt>
 *
 * <li> <tt><b>[D(h)/D(z)]^T</b></tt> as a <tt>MultiVectorBase</tt> object
 * where each column in the multivector <tt>k</tt> represents
 * <tt>[D(h(k))/D(z)]^T</tt>, the derivatives for the function <tt>h(k)</tt>
 * for all of the variables <tt>z</tt>
 *
 * </ul>
 *
 * A model can sign up to compute any, or all, or none of these forms of a
 * derivative and this information is returned from
 * <tt>this->createOutArgs().supports(OUT_ARG_blah,...)</tt> as a
 * <tt>ModelEvaluatorBase::DerivativeSupport</tt> object, where <tt>blah</tt>
 * is either <tt>DfDp</tt>, <tt>DgDx_dot</tt>,  <tt>DgDx</tt>, or <tt>DgDp</tt>.  The
 * <tt>LinearOpBase</tt> form of a derivative is supported if
 * <tt>this->createOutArgs().supports(OUT_ARG_blah,...).supports(DERIV_LINEAR_OP)==true</tt>.
 * The forward <tt>MultiVectorBase</tt> form of the derivative is supported if
 * <tt>this->createOutArgs().supports(OUT_ARG_blah,...).supports(DERIV_MV_BY_COL)==true</tt>
 * while the adjoint form of the derivative is supported if
 * <tt>this->createOutArgs().supports(OUT_ARG_blah,...).supports(DERIV_TRANS_MV_BY_ROW)==true</tt>.
 * 
 * In order to accommodate these different forms of a derivative, the simple
 * class <tt>ModelEvaluatorBase::Derivative</tt> is defined that can store
 * either a <tt>LinearOpBase</tt> or one of two <tt>MultiVectorBase</tt> forms
 * (i.e. forward or adjoint) of the derivative.  A
 * <tt>%ModelEvaluatorBase::%Derivative</tt> object can only store one (or
 * zero) forms of a derivative object at one time.
 *
 * We now describe each of these derivative objects in more detail:
 *
 * <ul>
 *
 * <li><b>State function <tt>f(x_dot,x,{p(l)},t,...)</tt> derivatives</b>
 *
 *     <ul>
 *     
 *     <li><b>State variable derivatives</b>

       \verbatim
       W = alpha*D(f)/D(x_dot) + beta*D(f)/D(x)
       \endverbatim

 *     This derivative operator is a special object that is derived from the
 *     <tt>LinearOpWithSolveBase</tt> interface and therefore supports linear
 *     solves.  Objects of this type are created with the function
 *     <tt>create_W()</tt> and set on the <tt>InArgs</tt> object before they
 *     are computed in <tt>evalModel()</tt>.  Note that if the model does not
 *     define or support an <tt>x_dot</tt> vector then the scalars
 *     <tt>alpha</tt> and <tt>beta</tt> need not be supported.
 *
 *     The <tt>LinearOpWithSolveBase</tt> form of this derivative is supported
 *     if <tt>this->createOutArgs().supports(OUT_ARG_W)==true</tt>.  The
 *     <tt>LinearOpBase</tt>-only form (i.e. no solve operation is given) of
 *     this derivative is supported if
 *     <tt>this->createOutArgs().supports(OUT_ARG_W_op)==true</tt>.  The
 *     <tt>W_op</tt> form of <tt>W</tt> is to be preferred when no linear solve
 *     with <tt>W</tt> will ever be needed.  A valid implementation may
 *     support none, both, or either of these forms <tt>LOWSB</tt> and/or
 *     <tt>LOB</tt> of <tt>W</tt>.
 *
 *     Also note that an underlying model may only support a single copy of
 *     <tt>W</tt> or <tt>W_op</tt> at one time.  This is required to
 *     accommodate some types of underlying applications that are less flexible
 *     in how they maintain their memory and how they deal with their objects.
 *     Therefore, to accommodate these types of less ideal application
 *     implementations, if a client tries to create and maintain more than one
 *     <tt>W</tt> and/or <tt>W_op</tt> object at one time, then
 *     <tt>create_W()</tt> and <tt>create_W_op()</tt> may return <tt>null</tt>
 *     (or may throw an undetermined exception).  However, every "good"
 *     implementation of this interface should support the creation and
 *     maintenance of as many <tt>W</tt> and/or <tt>W_op</tt> objects at one
 *     time as will fit into memory.
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

       \verbatim
       DfDp(l) = D(f)/D(p(l))
       \endverbatim

 *     for <tt>l=0...Np-1</tt>.

 *     These are derivative objects that represent the derivative of the state
 *     function <tt>f</tt> with respect to the auxiliary parameters
 *     <tt>p(l)</tt>.  This derivative is manipulated as a
 *     <tt>ModelEvaluatorBase::Derivative</tt> object.
 *
 *     </ul>
 *
 * <li><b>Auxiliary response function <tt>g(j)(x,{p(l)},t,...)</tt> derivatives</b>
 *
 *     <ul>
 *     
 *     <li><b>State variable derivatives</b>
 
       \verbatim
       DgDx_dot(j) = D(g(j))/D(x_dot)
       \endverbatim
 
 *     for <tt>j=0...Ng-1</tt>.
 *
 *     These are derivative objects that represent the derivative of the
 *     axillary function <tt>g(j)</tt> with respect to the state variables
 *     derivative <tt>x_dot</tt>.  This derivative is manipulated as a
 *     <tt>ModelEvaluatorBase::Derivative</tt> object.
 
       \verbatim
       DgDx(j) = D(g(j))/D(x)
       \endverbatim
 
 *     for <tt>j=0...Ng-1</tt>.
 *
 *     These are derivative objects that represent the derivative of the
 *     axillary function <tt>g(j)</tt> with respect to the state variables
 *     <tt>x</tt>.  This derivative is manipulated as a
 *     <tt>ModelEvaluatorBase::Derivative</tt> object.
 *     
 *     <li><b>Auxiliary parameter derivatives</b>
 *
       \verbatim
       DgDp(j,l) = D(g(j))/D(p(l))
       \endverbatim

 *     for <tt>j=0...Ng-1</tt>, <tt>l=0...Np-1</tt>.
 *
 *     These are derivative objects that represent the derivative of the
 *     axillary function <tt>g(j)</tt> with respect to the auxiliary
 *     parameters <tt>p(l)</tt>.  This derivative is manipulated as a
 *     <tt>ModelEvaluatorBase::Derivative</tt> object.
 *
 *     </ul>
 *
 * </ul>
 *
 * \section Thyra_ME_nominal_values_sec Nominal values
 *
 * A model can optionally define a nominal value for any of the input
 * arguments and these are returned from the <tt>getNominalValues()</tt>
 * function as a <tt>const ModelEvaluatorBase::InArgs</tt> object.  These
 * nominal values can be used as initial guesses, as typical values (e.g. for
 * scaling), or for other purposes.  See <tt>evalModel()</tt> a discussion of
 * how nonminal values are treated in an evaluation where the client does not
 * pass in the values explicitly.
 *
 * \section Thyra_ME_bounds_sec Variable and Function Bounds
 *
 * A model can optionally define a set of upper and/or lower bounds for each
 * of the input variables/parameters.  These bounds are returned as
 * <tt>const</tt> <tt>ModelEvaluatorBase::InArgs</tt> objects from the
 * functions <tt>getLowerBounds()</tt> and <tt>getUpperBounds()</tt>.  These
 * bounds are typically used to define regions in space where the model
 * functions are well defined.  A client algorithm is free to ignore these
 * bounds if they can not handle these types of constraints.
 *
 * That fact that a model can defined reasonable bounds but most numerical
 * algorithms can not handle bounds is no reason to leave bounds out of this
 * interface.  Again, if the client algorithm can not handle bounds then then
 * can be simply ignored and there is not harm done (except the client
 * algorithm might run into lots of trouble computing functions with undefined
 * values)..
 *
 * \section Thyra_ME_param_subvectors_sec Significance of Parameter Subvectors
 *
 * The parameters for any particular model are partitioned into different
 * subvectors <tt>p(l)</tt> for several different reasons:
 *
 * <ul>
 *
 * <li>Parameters are grouped together into a subvectors <tt>p(l)</tt> to
 * allow an ANA to manipulate an entire set at a time for different purposes.
 * It is up to someone to select which parameters from a model will be exposed
 * in a parameter subvector and how they are partitioned into subvectors.  For
 * example, one parameter subvector may be used a design parameters, while
 * another subvector may be used for uncertain parameters, while still another
 * may be used as continuation parameters.
 *
 * <li>Parameters are grouped together into subvectors <tt>p(l)</tt> implies
 * that certain derivatives will be supplied for all the parameters in the
 * subvector or none.  If an ANA wants to flexibility to get the derivative
 * for any individual scalar parameter by itself, then the paramters must
 * be segregated into different parameter subvectors <tt>p(l)</tt> with
 * one component each (e.g. <tt>get_p_space(l)->dim()==1</tt>).
 *
 * </ul>
 *
 * \section Thyra_ME_failed_evals_sec Failed evaluations
 *
 * The way for a ModelEvalutor object to return a failed evaluation is to set
 * NaN in one or more of the output objects.  If an algebraic model executes
 * and happens to pass a negative number to a squareroot or something, then a
 * NaN is what will get created anyway (even if the ModelEvaluator object does
 * not detect this).  Therefore, clients of the ME interface really need to be
 * searching for a NaN to see if an evaluation has faield.  Also, a ME object
 * can set the isFailed flag on the outArgs object on the return from the
 * <tt>evalModel()</tt> function if it knows that the evaluation has failed
 * for some reason.
 *
 * \section Thyra_ME_checking_sec Compile-Time and Run-Time Safety and Checking
 *
 * The <tt>%ModelEvaluator</tt> interface is designed to allow for great
 * flexibility in how models are defined.  The idea is that, at runtime, a
 * model can decide what input and output arguments it will support and client
 * algorithms must decide how to interpret what the model provides in order to
 * form an abstract problem to solve.  As a result, the
 * <tt>%ModelEvaluator</tt> interface is weakly typed mathematically.
 * However, the interface is strongly typed in terms of the types of objects
 * involved.  For example, while a single <tt>%ModelEvaluator</tt> object can
 * represent anything and everything from a set of nonlinear equations to a
 * full-blown constrained transient optimization problem, if the state vector
 * <tt>x</tt> is supported, it must be manipulated as a <tt>VectorBase</tt>
 * object, and this is checked at comple time.
 *
 * In order for highly dynamically configurable software to be safe and
 * usable, a great deal of runtime checking and good error reporting is
 * required.  Practically all of the runtime checking and error reporting work
 * associated with a <tt>%ModelEvaluator</tt> object is handled by the
 * concrete <tt>ModelEvaluatorBase::InArgs</tt> and
 * <tt>ModelEvaluatorBase::OutArgs</tt> classes.  Once a concrete
 * <tt>%ModelEvaluator</tt> object has setup what input and output arguments
 * the model supports, as returned by the <tt>createInArgs()</tt> and
 * <tt>createOutArgs()</tt> functions, all runtime checking for the proper
 * setting of input and output arguments and error reporting is automatically
 * handled.
 *
 * \section Thyra_ME_dev_notes_sec Notes to subclass developers
 *
 * Nearly every subclass should directly or indirectly derive from the node
 * subclass <tt>ModelEvaluatorDefaultBase</tt> since provides checking for
 * correct specification of a model evaluator subclass and provides default
 * implementations for various features.
 *
 * Subclass developers should consider deriving from one (but not more than
 * one) of the following subclasses rather than directly deriving from
 * <tt>ModelEvaluatorDefaultBase</tt>:
 *
 * <ul>
 *
 * <li><tt>StateFuncModelEvaluatorBase</tt> makes it easy to create models that
 * start with the state function evaluation <tt>x -> f(x)</tt>.
 *
 * <li><tt>ResponseOnlyModelEvaluatorBase</tt> makes it easy to create models
 * that start with the non-state response evaluation <tt>p -> g(p)</tt>.
 *
 * <li><tt>ModelEvaluatorDelegatorBase</tt> makes it easy to develop and
 * maintain different types of decorator subclasses.
 *
 * </ul>
 *
 * When deriving from any of the above intermediate base classes, the subclass
 * can override any of virtual functions in any way that it would like.
 * Therefore these subclasses just make the job of creating concrete
 * subclasses easier without removing any of the flexibility of creating a
 * subclass.
 *
 * ToDo: Finish Documentation!
 *
 * \ingroup Thyra_Nonlinear_model_evaluator_interfaces_code_grp
 */
template<class Scalar>
class ModelEvaluator : public ModelEvaluatorBase {
public:

  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \name Basic information */
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

  /** \brief Return the vector space for the state variables <tt>x <: RE^n_x</tt>. */
  virtual RCP<const VectorSpaceBase<Scalar> > get_x_space() const = 0;

  /** \brief Return the vector space for the state function <tt>f(...) <: RE^n_x</tt>. */
  virtual RCP<const VectorSpaceBase<Scalar> > get_f_space() const = 0;

  /** \brief Return the vector space for the auxiliary parameters
   * <tt>p(l) <: RE^n_p_l</tt>.
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
  virtual RCP<const VectorSpaceBase<Scalar> > get_p_space(int l) const = 0;

  /** \brief Get the names of the parameters associated with parameter
   * subvector l if available.
   *
   * \return Returns an RCP to a Teuchos::Array<std::string> object that
   * contains the names of the parameters.   If returnVal == Teuchos::null,
   * then there are no names available for the parameter subvector p(l).
   * If returnVal->size() == 1, then a single name is given to the entire
   * parameter subvector.  If returnVal->size() == get_p_space(l)->dim(),
   * then a name is given to every parameter scalar entry.
   *
   * The default implementation return returnVal==Teuchos::null which means
   * by default, parameters have no names associated with them.
   */
  virtual RCP<const Teuchos::Array<std::string> > get_p_names(int l) const = 0;

  /** \brief Return the vector space for the auxiliary response functions
   * <tt>g(j) <: RE^n_g_j</tt>.
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
  virtual RCP<const VectorSpaceBase<Scalar> > get_g_space(int j) const = 0;

  //@}

  /** \name Initial guess and upper and lower bounds */
  //@{
  
  /** \brief Return the set of nominal values or initial guesses for the
   * supported the input arguments.
   *
   * In most cases, when a supported input argument is not specified in an
   * <tt>InArgs</tt> object passed to <tt>evalModel()</tt>, the nominal value
   * is assumed (see <tt>evalModel()</tt> for more details).
   *
   * A model does not have to return nominal values for every supported input
   * argument.  Therefore, just because
   * <tt>returnVal.supports(IN_ARG_blah)==true</tt>, the client should not
   * assume that <tt>returnVal.get_blah()!=null</tt>.
   *
   * <b>Warning!</b> Clients should not try to modify the state of the
   * contained objects.  Doing so might invalidate the state of <tt>*this</tt>
   * model or of other clients.  This is a hole in the const support and the
   * use of InArgs.  If a client wants to modify a value, then a clone of that
   * value should be created and reset on the returned InArgs object.  Do
   * *not* modify the values in place!
   */
  virtual ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const = 0;
  
  /** \brief Return the set of lower bounds for the input arguments.
   *
   * It is not required for the client to supply lower bounds for every input
   * argument.  If a lower bound object for an input argument is null, then
   * the bounds of negative infinity should be assumed by the client.
   *
   * <b>Warning!</b> Clients should not try to modify the state of the
   * contained objects.  Doing so might invalidate the state of <tt>*this</tt>
   * model or of other clients.  This is a hole in the const support and the
   * use of InArgs.  If a client wants to modify a value, then a clone of that
   * value should be created and reset on the returned InArgs object.  Do
   * *not* modify the values in place!
   */
  virtual ModelEvaluatorBase::InArgs<Scalar> getLowerBounds() const = 0;
  
  /** \brief Return the set of upper bounds for the input arguments.
   *
   * It is not required for the client to supply upper bounds for every input
   * argument.  If an upper bound object for an input argument is null, then
   * the bounds of positive infinity should be assumed by the client.
   *
   * <b>Warning!</b> Clients should not try to modify the state of the
   * contained objects.  Doing so might invalidate the state of <tt>*this</tt>
   * model or of other clients.  This is a hole in the const support and the
   * use of InArgs.  If a client wants to modify a value, then a clone of that
   * value should be created and reset on the returned InArgs object.  Do
   * *not* modify the values in place!
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
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>!is_null(returnVal)</tt>
   * </ul>
   *
   * Note that a model is only required to support a single <tt>W</tt>
   * object if the precondition below is satisfied and if the client asks for
   * more than one <tt>W</tt> object, the response should be to return
   * <tt>null</tt> from this function.
   */
  virtual RCP<LinearOpWithSolveBase<Scalar> > create_W() const = 0;

  /** \brief If supported, create a <tt>LinearOpBase</tt> object for
   * <tt>W</tt> to be evaluated.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->createOutArgs().supports(OUT_ARG_W_op)==true</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>!is_null(returnVal)</tt>
   * <li><tt>isPartiallyInitialized(*returnVal) || isFullyInitialized(*returnVal)</tt>
   * </ul>
   *
   * Note that a model is only required to support a single <tt>W_op</tt>
   * object if the precondition below is satisfied and if the client asks for
   * more than one <tt>W_op</tt> object, the response should be to return
   * <tt>null</tt> from this function.
   *
   * Also note the above post-condition that requires that a created
   * <tt>W_op</tt> object also needs to be at least partially initialized.
   * This means that its range and domain spaces must be fully formed and
   * accessible.  This greatly simplifies the creation of composite structures
   * and simplifies the implementation of certain types of implicit
   * ModelEvaluator subclasses.
   */
  virtual RCP<LinearOpBase<Scalar> > create_W_op() const = 0;

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
  virtual RCP<LinearOpBase<Scalar> > create_DfDp_op(int l) const = 0;

  /** \brief If supported, create a linear operator derivative object for
   * <tt>D(g(j))/D(x_dot)</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->Ng() > 0</tt>
   * <li><tt>0 <= j < this->Ng()</tt>
   * <li><tt>outArgs.supports_DgDx_dot(j).supports(DERIV_LINEAR_OP)==true</tt>,
   *     where <tt>outArgs = this->createOutArgs()</tt>
   * </ul>
   */
  virtual RCP<LinearOpBase<Scalar> > create_DgDx_dot_op(int j) const = 0;

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
  virtual RCP<LinearOpBase<Scalar> > create_DgDx_op(int j) const = 0;

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
  virtual RCP<LinearOpBase<Scalar> > create_DgDp_op( int j, int l ) const = 0;
  
  //@}

  /** \name Linear solver factory for W */
  //@{

  /** \brief If supported, return a <tt>LinearOpWithSolveFactoryBase</tt>
   * object that can be used to initialize a <tt>LOWSB</tt> object for
   * <tt>W</tt> given a <tt>LOB</tt> object for <tt>W_op</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->createOutArgs().supports(OUT_ARG_W) ||
   *     this->createOutArgs().supports(OUT_ARG_W_op)</tt>
   * </ul>
   *
   * By returning this factory object, a model evaluator allows a client to
   * compute the <tt>LOB</tt>-only version <tt>W_op</tt> and then allow client
   * to create a linear solver associated through the <tt>LOWSB</tt> interface
   * at any time.  Note that this allow allows access to the underlying
   * precondioner factory object
   * (i.e. <tt>this->get_W_factory()->getPreconditionerFactory()</tt.) if such
   * a factory object is supported.
   *
   * This function will return <tt>retalVal==Teuchos::null</tt> if no such
   * factory object is supported.
   */
  virtual RCP<const LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const = 0;

  //@}

  /** \name Computational functions */
  //@{

  /** \brief Create an empty input arguments object that can be set up and
   * passed to <tt>evalModel()</tt>.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>returnVal.supports(IN_ARG_blah)</tt> gives the variables that are
   * supported or not supported by the underlying model.
   * <li><tt>returnVal.get_blah(...)==null</tt> for all variables/parameters
   * that are supported.
   * </ul>
   * 
   */
  virtual ModelEvaluatorBase::InArgs<Scalar> createInArgs() const = 0;

  /** \brief Create an empty output functions/derivatives object that can be
   * set up and passed to <tt>evalModel()</tt>.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>returnVal.supports(OUT_ARG_blah...)</tt> gives the functions/derivatives that
   * are supported or not supported by the underlying model.
   * <li><tt>returnVal.get_blah(...) == null</tt> for all functions/derivatives
   * that are supported.
   * </ul>
   *
   */
  virtual ModelEvaluatorBase::OutArgs<Scalar> createOutArgs() const = 0;

  /** \brief Compute all of the requested functions/derivatives at the given
   * point.
   *
   * \param inArgs [in] Gives the values of the input variables and parameters
   * that are supported by the model.  This object must have been initially
   * created by <tt>this->createInArgs()</tt> before being set up by the
   * client.  All supported variables/parameters that are not set by the
   * client (i.e. <tt>inArgs.get_blah(...)==null</tt>) are assumed to be at
   * the nominal value.  If a particular input variable/parameter is not set
   * and does not have a nominal value
   * (i.e. <tt>this->getNominalValues().get_blah(...)==null</tt>, then the
   * evaluation is undefined and an exception should be thrown by a good
   * implementation.  The one exception this rule is support for
   * <tt>x_dot</tt>.  If <tt>x_dot</tt> is supported but not specified in
   * <tt>inArgs</tt>, then it is implicitly assumed to be zero.
   *
   * \param outArgs [out] Gives the objects for the supported functions and
   * derivatives that are to be computed at the given point.  This object must
   * have been created by <tt>this->createOutArgs()</tt> before being set up
   * by the client.  Only functions and derivatives that are set will be
   * computed.  The client should check for NaN in any of objects computed to
   * see if the evaluation failed.  Also, if the underlying object knows that
   * the evaluation failed, then <tt>outArgs.isFailed()</tt> will be
   * <tt>true</tt> on output.  The client needs to check for both NaN and
   * <tt>outArgs.isFailed()</tt> to be sure.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->createInArgs().isCompatible(inArgs)==true</tt>
   * <li><tt>this->creatOutnArgs().isCompatible(outArgs)==true</tt>
   * <li>[<tt>inArgs.supports(IN_ARG_blah)==true</tt>] <tt>inArgs.get_blah(...)!=null
   * || this->getNominalValues().get_blah(...)!=null</tt>.
   * <li><tt>inArgs.get_p(l)!=null ||
   * this->getNominalValues().get_pl(l)!=null</tt>, for
   * <tt>l=0...inArgs.Np()-1</tt>.
   * </ul>
   *
   * This function is the real meat of the <tt>%ModelEvaluator</tt> interface.
   * A whole set of supported functions and/or their various supported
   * derivatives are computed all in one shot at a given point (i.e. set of
   * values for input variables and parameters).  This allows a model to be
   * stateless in the sense that, in general, a model's behavior is assumed to
   * be unaffected by evaluations at previous points.  This greatly simplifies
   * software maintenance and makes data dependences explicit.
   *
   * <b>WARNING:</b> The implementation of this function must not create any
   * persisting associations involving any of the input or output objects
   * (even though these are returned through RCP objects).  This includes not
   * embedding any of the input objects in the created output objects.  One
   * reason for this requirement is that the RCP objects may not be strong
   * owning RCPs in which case the objects may not stay around.  This will
   * likely result in a dangling reference exception in a debug-mode bulid (or
   * a segfault in a non-debug build).  Also, even if the objects do stay
   * around, the implementation of <tt>evalModel(...)</tt> can not assume that
   * the objects will not be further modified by the client after
   * <tt>evalModel(...)</tt> is called.  For example, the value of the vector
   * returned from <tt>*inArgs.get_x()</tt> may changed by the client just
   * after this function finishes which would invalidate any objects that
   * might have expected it to not change.
   *
   * TODO: Define behavior for state-full models!
   */
  virtual void evalModel(
    const ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const ModelEvaluatorBase::OutArgs<Scalar> &outArgs
    ) const = 0;

  //@}

  /** \name Reporting functions */
  //@{

  /** \brief Report the final point and whether the problem was considered
   * solved or not.
   *
   * ToDo: Add a PL to InArgs to allow extra data like Lagrange multipliers to
   * be passed back as well?
   */
  virtual void reportFinalPoint(
    const ModelEvaluatorBase::InArgs<Scalar> &finalPoint,
    const bool wasSolved
    ) = 0;

  //@}

private:
  
  // Not defined and not to be called
  ModelEvaluator<Scalar>&
  operator=(const ModelEvaluator<Scalar>&);

};

} // namespace Thyra


#endif // THYRA_MODEL_EVALUATOR_HPP
