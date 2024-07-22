// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_LINEAR_OP_WITH_SOLVE_BASE_DECL_HPP
#define THYRA_LINEAR_OP_WITH_SOLVE_BASE_DECL_HPP

#include "Thyra_OperatorSolveTypes.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_SolveSupportTypes.hpp"
#include "Teuchos_VerboseObject.hpp"


namespace Thyra {


/** \brief Base class for all linear operators that can support a high-level
 * solve operation.
 *
 * \ingroup Thyra_Op_Solve_fundamental_interfaces_code_grp
 *
 * \section LOWSB_outline_sec Outline
 *
 * <ul>
 * <li>\ref LOWSB_intro_sec
 * <li>\ref LOWSB_solve_criteria_sec
 * <li>\ref LOWSB_solve_status_sec
 * <li>\ref LOWSB_use_cases_sec 
 * <li>\ref LOWSB_developer_notes_sec
 * </ul>
 * 
 * \section LOWSB_intro_sec Introduction
 * 
 * This interface supports linear operators (with potentially different range
 * and domain scalar types) that can also support a forward solve operation
 * (using <tt>solve()</tt>) of the form:

  \f[
    A X = B 
  \f]
 
 * and/or a transpose solve operation (using <tt>A_trans==TRANS</tt>) of the
 * form:

  \f[
    A^T X = B 
  \f]

 * and/or an adjoint solve operation (using <tt>A_trans==CONJTRANS()</tt>) of
 * the form:

  \f[
    A^H X = B 
  \f]

 * where \f$A\f$ is <tt>*this</tt> linear operator, \f$B\f$ is an appropriate
 * RHS multi-vector with \f$m\f$ columns, and \f$X\f$ is a LHS multi-vector
 * with \f$m\f$ columns that is computed by this interface.  Note that if the
 * underlying operator has real-valued entries then the transpose \f$A^T\f$
 * and the adjoint \f$A^H\f$ are the same.

 * Note that this interface does not assume that the linear operator itself is
 * nonsingular or invertible in the classic sense (i.e. an inverse operator
 * may not exist).

 * Let \f$op(A)\f$ signify either the forward operator \f$A\f$, the transpose
 * operator \f$A^T\f$ or the adjoint operator \f$A^H\f$.  What this interface
 * assumes is that for any appropriately selected consistent multi-vector RHS
 * <tt>B</tt> that a solve of \f$A X = B\f$ will yield an approximate solution
 * LHS multi-vector <tt>X</tt> such that <tt>A X == B</tt>.  Note that this
 * interface does not assume that a solution \f$X\f$ can be computed for any
 * random RHS multi-vector \f$B\f$.  Solutions for any random RHS can on be
 * expected for relatively well conditioned non-singular operators.

 * <b>Note:</b> It is recommended that clients use the non-member helper
 * functions defined \ref Thyra_LinearOpWithSolveBase_helper_grp "here" rather
 * than call these member functions directly as they support a number of other
 * simpler use cases.
 * 
 * \section LOWSB_solve_criteria_sec Solve Criteria
 * 
 * This interface potentially allows clients to specify a relative tolerance
 * on either the relative residual norm or the relative norm of the solution
 * error and can target different solution criteria to different blocks of
 * linear systems.  This interface tries to allow for mathematically rigorous
 * solution tolerances that are not based only any implementation-dependent
 * features like the number of iterations of some solver algorithm.  This
 * interface, however, allows <tt>*this</tt> operator to exclude support
 * certain types of solve measures (see the functions
 * <tt>solveSupportsSolveMeasureType()</tt> and
 * <tt>solveTransposeSupportsSolveMeasureType()</tt>).  Also, this interface
 * assumes that all implementations can support a "default" solve criteria
 * that is determined internally to <tt>*this</tt>.
 * 
 * This interface is meant to support direct and iterative linear solvers as
 * well as combinations of the two in a variety of configurations.  Because of
 * the almost infinite number of types of linear solver configurations
 * possible, this interface tries not to specify any particular
 * solver-specific solution control options.  The one exception is a maximum
 * number of iterations which is totally implementation defined.  These types
 * of control options are better specified in lower lever implementations and
 * should be kept out of an interface such as this.
 *
 * Let <tt>solveCriteria</tt> be a <tt>SolveCriteria</tt> object setup by a
 * client to be passed into a solve operation.  This object can be set up in a
 * variety of ways to support several different use cases which are described
 * below:
 * 
 * <ul>
 * 
 * <li> <b>Unspecified (Default) solution criteria type</b> [
 * <tt>solveCriteria.solveMeasureType.useDefault()==true</tt> ]: In this mode,
 * criteria for the solution tolerance is determined internally by
 * <tt>*this</tt> object.  Usually, it would be assumed that <tt>*this</tt>
 * would solve the linear systems to a sufficient tolerance for the given ANA
 * client.  This is the mode that many ANAs not designed to control
 * inexactness would work with.  In this mode, the solution criteria will be
 * determined in the end by the application or the user in an "appropriate"
 * manner.  In this case, the value of <tt>solveCriteria.requestedTol</tt> is
 * ignored and no meaningful solve status can be returned to the client.
 * 
 * <li> <b>Relative residual tolerance solution criteria</b> [
 * <tt>solveCriteria.solveMeasureType.numerator==SOLVE_MEASURE_NORM_RESIDUAL
 * && solveCriteria.solveMeasureType.denominator==SOLVE_MEASURE_NORM_RHS &&
 * solveCriteria.requestedTol!=SolveCriteria::unspecifiedTolerance()</tt> ]:
 * In this mode, the solution algorithm would be requested to solve the block
 * system to the relative residual tolerance of

  \f[
     \frac{|| op(A) X_{(:,i)} - B_{(:,i)} ||}{ ||B_{(:,i)}||} \le \mu_r,
  \f]

 * for \f$i=0 {}\ldots m-1\f$, where \f$\mu_r =\f$
 * <tt>solveCriteria.requestedTol</tt> and where the norm \f$||.||\f$ is
 * given by the natural norm defined by the range space of \f$op(A)\f$ and
 * computed from <tt>Thyra::norm()</tt>.  Many linear solvers should be able
 * to monitor this tolerance and be able to achieve it, or on failure be able
 * to report the actual tolerance achieved.
 *
 * <li> ToDo: Add examples of other types of solve measures when they are
 * needed!
 * 
 * </ul>
 * 
 * \section LOWSB_solve_status_sec Solve Status
 * 
 * After the <tt>solve()</tt> and <tt>solveTranspose()</tt> functions return,
 * the client can optionally get back a solution status for each block of
 * linear systems for of block solve criteria.  Specifically, for each block
 * of linear systems
  
  \f[
    A  X_{(:,i_{j}+1:i_j)} = B_{(:,i_{j}+1:i_j)} 
  \f]

 * whose solution criteria is specified by a <tt>SolveCriteria</tt> object, a
 * <tt>SolveStatus</tt> object can optionally be returned that lets the client
 * know the status of the linear solve.
 *
 * A note about direct solvers is in order.  The "inexact" solve features of
 * this interface are primarily designed to support "loose" solve tolerances
 * that exploit the properties of iterative linear solvers.  With that said,
 * any decent direct solver can assume that it has met the convergence
 * criteria as requested by the client but does not have to return an estimate
 * of the actual tolerance achieved.
 *
 * If <tt>solveStatus</tt> is a <tt>SolveStatus</tt> object returned for the
 * above block linear system the the following return status are significant:
 * 
 * <ul>
 * 
 * <li><b>Converged</b> [
 * <tt>solveStatus.solveStatus==SOLVE_STATUS_CONVERGED</tt> ]: This status is
 * returned by the linear solver if the solution criteria was likely achieved
 * (within some acceptable cushion cased by round-off etc.).  This should
 * almost always be the return value for a direct linear solver.  The maximum
 * actual tolerance achieved may or may not be returned in the field
 * <tt>solveStatus.achievedTol</tt>.  The two sub-cases are:
 * 
 *   <ul>
 * 
 *   <li><b>Known tolerance</b> [
 *   <tt>solveStatus.achievedTol!=SolveStatus::unknownTolerance()</tt> ] : The
 *   linear solver knows the approximate order-of-magnitude estimate of the
 *   maximum tolerance achieved.  An order-of-magnitude (or so) estimate of
 *   the achieved tolerance would likely be known by any iterative linear
 *   solver when
 *   <tt>solveCriteria.solveMeasureType.numerator==SOLVE_MEASURE_NORM_RESIDUAL
 *   &&
 *   solveCriteria.solveMeasureType.denominator==SOLVE_MEASURE_NORM_RHS</tt>.
 *   Most direct linear solvers will not return a known tolerance except for
 *   the case where
 *   <tt>solveCriteria.solveMeasureType.numerator==SOLVE_MEASURE_NORM_RESIDUAL
 *   &&
 *   solveCriteria.solveMeasureType.denominator==SOLVE_MEASURE_NORM_RHS</tt>
 *   and where iterative refinement is used.

 *   <li><b>Unknown tolerance</b> [
 *   <tt>solveStatus.achievedTol==SolveStatus::unknownTolerance()</tt> ] : The
 *   linear solver does not know the tolerance that was achieved but the
 *   achieved tolerance should be very close to the requested tolerance.  This
 *   would be the most likely return status for a direct linear solver.
 * 
 *   </ul>
 * 
 * <li><b>Unconverged</b> [
 * <tt>solveStatus.solveStatus==SOLVE_STATUS_UNCONVERGED</tt> ]: The linear
 * solver was most likely not able to achieve the requested tolerance.  A
 * direct linear solver should almost never return this status except for in
 * extreme cases (e.g. highly ill conditioned matrix and a tight requested
 * tolerance).  The linear solver may not be able to return the actual
 * tolerance achieved and the same two cases as for the <it>unconverged</it>
 * case are possible and the two subdcases are:
 * 
 *   <ul>
 * 
 *   <li><b>Known tolerance</b> [
 *   <tt>solveStatus.achievedTol!===SolveStatus::unknownTolerance()0</tt> ] :
 *   The linear solver knows the approximate order-of-magnitude estimate of
 *   the maximum tolerance achieved.
 * 
 *   <li><b>Unknown tolerance</b> [
 *   <tt>solveStatus.achievedTol==SolveStatus::unknownTolerance()</tt> ] : The
 *   linear solver does not know the tolerance that was achieved but the
 *   achieved tolerance is most likely significantly greater than the
 *   requested tolerance.
 * 
 *   </ul>
 * 
 * <li><b>Unknown status</b> [
 * <tt>solveStatus.solveStatus==SOLVE_STATUS_UNKNOWN</tt> ]: The linear solver
 * does not know if the solution status was achieved or not.  In this case,
 * the value of
 * <tt>solveStatus.achievedTol==SolveStatus::unknownTolerance()</tt> may or
 * many not be returned.  This may be the return value when there is no
 * reasonable way that the linear solver algorithm can know now to compute or
 * estimate the requested tolerance.  This will also always be the return
 * status when <tt>solveCriteria.solveMeasureType.useDefault()==true</tt> since
 * the client would have no way to interpret this tolerance.  The value of
 * <tt>solveStatus.achievedTol!=SolveStatus::unknownTolerance()</tt> in this
 * case should only be returned when
 * <tt>solveCriteria.solveMeasureType.useDefault()==true</tt> and therefore the
 * client would have no way to interpret this tolerance as a residual or an
 * solution norm.
 * 
 * </ul>
 *
 * The implementation of the function <tt>accumulateSolveStatus()</tt> defines
 * how to accumulate the individual solve status for each RHS in a block into
 * the overall solve status for an entire block returned in the
 * <tt>SolveStatus</tt> object from teh <tt>solve()</tt>.
 * 
 * \section LOWSB_use_cases_sec Use cases
 *
 * This interface supports a variety of use cases where where described, more
 * or less, in the above sections.  Here, we give specific examples for a
 * number of important use cases and show how to use the non-member helper
 * functions defined \ref Thyra_LinearOpWithSolveBase_helper_grp "here".
 *
 * ToDo: Finish documentation!
 * 
 * \section LOWSB_developer_notes_sec Notes to subclass developers
 * 
 * This interface assumes, by default, that subclasses will only support the
 * forward solve operation.  See <tt>LinearOpBase</tt> for what other virtual
 * functions must be overridden to completely define a concrete subclass.
 */
template<class Scalar>
class LinearOpWithSolveBase
  : virtual public LinearOpBase<Scalar>
  , virtual public Teuchos::VerboseObject<LinearOpWithSolveBase<Scalar> >
{
public:

  /** \name Public interface funtions. */
  //@{

  // 2010/08/22: rabartl: To properly handle the new SolveCriteria struct with
  // reduction functionals (bug 4915) the function solveSupports() must be
  // refactored.  Here is how this refactoring can be done incrementally and
  // safely:
  //
  // (*) Create new override solveSupports(transp, solveCriteria) that calls
  // virtual solveSupportsNewImpl(transp, solveCriteria).
  //
  // (*) One by one, refactor existing LOWSB subclasses to implement
  // solveSupportsNewImpl(transp, solveCriteria).  This can be done by
  // basically copying the existing solveSupportsSolveMeasureTypeImpl()
  // override.  Then have each of the existing
  // solveSupportsSolveMeasureTypeImpl() overrides call
  // solveSupportsNewImpl(transp, solveCriteria) to make sure that
  // solveSupportsNewImpl() is getting tested right away.  Also, have the
  // existing solveSupportsImpl(...) overrides call
  // solveSupportsNewImpl(transp, null).  This will make sure that all
  // functionality is now going through solveSupportsNewImpl(...) and is
  // getting tested.
  //
  // (*) Refactor Teko software.
  //
  // (*) Once all LOWSB subclasses implement solveSupportsNewImpl(transp,
  // solveCriteria), finish off the refactoring in one shot:
  //
  //   (-) Remove the function solveSupports(transp), give solveCriteria a
  //   default null in solveSupports(transp, solveCriteria).
  //
  //   (-) Run all tests.
  //
  //   (-) Remove all of the solveSupportsImpl(transp) overrides, rename solve
  //   solveSupportsNewImpl() to solveSupportsImpl(), and make
  //   solveSupportsImpl(...) pure virtual.
  //
  //   (-) Run all tests.
  //
  //   (-) Change solveSupportsSolveMeasureType(transp, solveMeasureType) to
  //   call solveSupportsImpl(transp, solveCriteria) by setting
  //   solveMeasureType on a temp SolveCriteria object.  Also, deprecate the
  //   function solveSupportsSolveMeasureType(...).
  //
  //   (-) Run all tests.
  //
  //   (-) Remove all of the existing solveSupportsSolveMeasureTypeImpl()
  //   overrides.
  //
  //   (-) Run all tests.
  //
  //   (-) Clean up all deprecated warnings about calling
  //   solveSupportsSolveMeasureType() and instead have them call
  //   solveSupports(...) with a SolveCritera object.
  //
  // (*) Enter a release notes item about this breaking backward compatiblilty
  // for existing subclasses of LOWSB.
  //
  // This refactoring will be done by and by has bug 4915 is implemented.
  // 

  /** \brief Return if <tt>solve()</tt> supports the argument <tt>transp</tt>.
   *
   * The default implementation returns <tt>true</tt> for non-transposed,
   * non-conjugate solves..
   */
  bool solveSupports(EOpTransp transp) const
    { return solveSupportsImpl(transp); }

  /** \brief Return if <tt>solve()</tt> supports a given transpose and solve
   * criteria specification.
   *
   */
  bool solveSupports(EOpTransp transp,
    const Ptr<const SolveCriteria<Scalar> > solveCriteria) const
    { return solveSupportsNewImpl(transp, solveCriteria); }

  /** \brief Return if <tt>solve()</tt> supports the given the solve measure
   * type.
   *
   * The default implementation returns <tt>true</tt> for
   * <tt>solveMeasureType.useDefault()==true</tt>.
   */
  bool solveSupportsSolveMeasureType(EOpTransp transp,
    const SolveMeasureType& solveMeasureType
    ) const
    { return solveSupportsSolveMeasureTypeImpl(transp, solveMeasureType); }

  /** \brief Request the solution of a block linear system.
   *
   * \param A_trans [in] Determines if the elements are non-conjugate
   * non-transposed (<tt>NONTRANS</tt>) or conjugate transposed
   * (<tt>CONJTRANS</tt>).
   *
   * \param B [in] The RHS multi-vector with <tt>m = B.domain()->dim()</tt>
   * columns.
   *
   * \param X [in/out] The LHS multi-vector with with <tt>m =
   * X->domain()->dim()</tt> columns.  On input, contains the initial guess
   * for the solution (only significant for iterative solvers) and on output
   * contains an estimate of the solution.
   *
   * \param solveCriteria [in] Gives the desired solution criteria for linear
   * systems.  A value of <tt>solveCriteria==null</tt> means to use the
   * default solve criteria.
   *
   * \returns Return the solve status if any program has been made.
   *
   * <b>Preconditions:</b><ul>
   *
   * <li><tt>this->solveSupports(transp, solveCriteria)==true</tt>
   *
   * <li><tt>nonnull(X)==true</tt>
   *
   * <li><tt>op(this)->range()->isCompatible(*B.range())==true</tt>
   *
   * <li><tt>op(this)->domain()->isCompatible(*X->range())==true</tt>
   *
   * <li><tt>B->domain()->isCompatible(*X->domain())==true</tt>
   *
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   *
   * <li>If any progress on solving the linear systems could be achieved, then
   * the function will return normally with the solve status.  However, if no
   * progress could be made in solving the linear systems, then an exception
   * of type <tt>CatastrophicSolveFailure</tt> will be thrown.  If the
   * function returns normally, then the following postconditions apply.
   *
   * </ul>
   *
   * See the above introduction for a more complete description of how this
   * function behaves and the meaning of the argument <tt>solveCriteria</tt>
   * and the return value.
   */
  SolveStatus<Scalar> solve(
    const EOpTransp A_trans,
    const MultiVectorBase<Scalar> &B,
    const Ptr<MultiVectorBase<Scalar> > &X,
    const Ptr<const SolveCriteria<Scalar> > solveCriteria = Teuchos::null
    ) const 
    { return solveImpl(A_trans, B, X, solveCriteria); }

  //@}

protected:

  /** \name Protected virtual functions to be overridden by subclasses. */
  //@{

  /** \brief Virtual implementation for solveSupports(). */
  virtual bool solveSupportsImpl(EOpTransp transp) const;

  /** \brief Virtual implementation of <tt>solveSupports()</tt>. */
  virtual bool solveSupportsNewImpl(EOpTransp /* transp */,
    const Ptr<const SolveCriteria<Scalar> > /* solveCriteria */
    ) const
    {
      TEUCHOS_TEST_FOR_EXCEPT(true);
      TEUCHOS_UNREACHABLE_RETURN(false);
    }

  /** \brief Virtual implementation for solveSupportsSolveMeasureType(). */
  virtual bool solveSupportsSolveMeasureTypeImpl(EOpTransp transp,
    const SolveMeasureType& solveMeasureType) const;

  /** \brief Virtual implementation for solve(). */
  virtual SolveStatus<Scalar> solveImpl(
    const EOpTransp transp,
    const MultiVectorBase<Scalar> &B,
    const Ptr<MultiVectorBase<Scalar> > &X,
    const Ptr<const SolveCriteria<Scalar> > solveCriteria
    ) const = 0;

  //@}

private:

private:
  
  // Not defined and not to be called
  LinearOpWithSolveBase<Scalar>&
  operator=(const LinearOpWithSolveBase<Scalar>&);

};


/** \brief Call <tt>solveSupports()</tt> as a non-member function.
 *
 * \relates LinearOpWithSolveBase
 */
template<class Scalar>
inline
bool solveSupports(const LinearOpWithSolveBase<Scalar> &A, const EOpTransp transp)
{
  return A.solveSupports(transp);
}


/** \brief Call <tt>solveSupports()</tt> as a non-member function.
 *
 * \relates LinearOpWithSolveBase
 */
template<class Scalar>
inline
bool solveSupports(
  const LinearOpWithSolveBase<Scalar> &A,
  const EOpTransp transp,
  const Ptr<const SolveCriteria<Scalar> > solveCriteria
  )
{
  return A.solveSupports(transp, solveCriteria);
}


/** \brief Call <tt>solveSupportsSolveMeasureType()</tt> as a non-member
 * function.
 *
 * \relates LinearOpWithSolveBase
 */
template<class Scalar>
inline
bool solveSupportsSolveMeasureType(
  const LinearOpWithSolveBase<Scalar> &A,
  const EOpTransp transp,
  const SolveMeasureType &solveMeasureType
  )
{
  return A.solveSupportsSolveMeasureType(transp, solveMeasureType);
}


/** \brief Call <tt>solve()</tt> as a non-member function
 *
 * \relates LinearOpWithSolveBase
 */
template<class Scalar>
inline
SolveStatus<Scalar> solve(
  const LinearOpWithSolveBase<Scalar> &A,
  const EOpTransp A_trans,
  const MultiVectorBase<Scalar> &B,
  const Ptr<MultiVectorBase<Scalar> > &X,
  const Ptr<const SolveCriteria<Scalar> > solveCriteria = Teuchos::null
  )
{
  return A.solve(A_trans, B, X, solveCriteria);
}



} // namespace Thyra


#endif // THYRA_LINEAR_OP_WITH_SOLVE_BASE_DECL_HPP
