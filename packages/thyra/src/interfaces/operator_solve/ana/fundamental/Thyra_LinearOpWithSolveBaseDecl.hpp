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

#ifndef THYRA_LINEAR_OP_WITH_SOLVE_BASE_DECL_HPP
#define THYRA_LINEAR_OP_WITH_SOLVE_BASE_DECL_HPP

#include "Thyra_LinearOpBaseDecl.hpp"
#include "Thyra_SolveSupportTypes.hpp"
#include "Teuchos_VerboseObject.hpp"

namespace Thyra {

/** \brief Base class for all linear operators that can support a high-level
 * solve operation.
 * 
 * \ingroup Thyra_Op_Vec_Interoperability_Extended_Interfaces_grp
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
 
 * and/or a transpose solve operation (using <tt>solveTranspose()</tt>) of the
 * form:

  \f[
    A^T X = B 
  \f]

 * and/or an adjoint solve operation (using <tt>solveTranspose()</tt>) of the
 * form:

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
 * The functions <tt>solve()</tt> and <tt>solveTranspose()</tt> both take the
 * arguments:

 \code
    ,const int                                    numBlocks             = 0
    ,const BlockSolveCriteria<PromotedScalar>     blockSolveCriteria[]  = NULL
    ,SolveStatus<PromotedScalar>                  blockSolveStatus[]    = NULL
 \endcode
 *
 * The array arguments <tt>blockSolveCriteria[]</tt> and
 * <tt>blockSolveStatus[]</tt> specify different blocks of solution criteria
 * and the corresponding solve return statuses for a partitioned set of linear
 * systems.  Assuming that the client passes in arrays of dimensions
 * \f$N=\f$<tt>numBlocks</tt>, these tolerances define the solution criteria
 * for the block systems:
   
 \f[
   op(A) \left[ \begin{array}{ccccc} X_{(:,0:i_1)} & X_{(:,i_1+1:i_2)} & \ldots & X_{(:,i_{N-1}+1:i_N)} \end{array} \right]
   = \left[ \begin{array}{ccccc} B_{(:,0:i_1)} & B_{(:,i_1+1:i_2)} & \ldots & B_{(:,i_{N-1}+1:i_N)} \end{array} \right]
 \f]

 * where the column indexes are given by \f$i_j = \left( \sum_{k=0}^{j}
 * \mbox{blockSolveCriteria[k].numRhs} \right)\f$, for \f$j = 0 \ldots N-1\f$.
 *
 * The solve criteria for the \f$j^{\mbox{th}}\f$ block system
   
 \f[
   op(A)  X_{(:,i_{j}+1:i_j)} = B_{(:,i_{j}+1:i_j)} 
 \f]

 * is given by <tt>blockSolveCriteria[j]</tt> (if
 * <tt>blockSolveCriteria!=NULL</tt>) and the solution status after return is
 * given by <tt>blockSolveStatus[j]</tt> (if
 * <tt>blockSolveStatus!=NULL</tt>).
 *
 * By specifying solution criteria in blocks and then only requesting basic
 * tolerances, we allow linear solver implementations every opportunity to
 * perform as many optimizations as possible in solving the linear systems.
 * For example, SVD could be performed on each block of RHSs and then a
 * reduced set of linear systems could be solved.
 *
 * For the remainder of this discussion we will focus on how the solution
 * criteria for single block of linear systems specified by a single
 * <tt>BlockSolveCriteria</tt> object is specified and how the status of this
 * block linear solve is reported in a <tt>SolveStatus</tt> object.
 *
 * The struct <tt>BlockSolveCriteria</tt> contains a <tt>SolveCriteria</tt>
 * member and a number of RHSs that it applies to.  It is the
 * <tt>SolveCriteria</tt> object that determines the type and tolerance of a
 * block solve request.
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
 *   <tt>solveCriteria.solveMeasureType.numerator==SOLVE_MEASURE_NORM_RESIDUAL&&solveCriteria.solveMeasureType.denominator==SOLVE_MEASURE_NORM_RHS</tt>.
 *   Most direct linear solvers will not return a known tolerance except for
 *   the case where
 *   <tt>solveCriteria.solveMeasureType.numerator==SOLVE_MEASURE_NORM_RESIDUAL&&solveCriteria.solveMeasureType.denominator==SOLVE_MEASURE_NORM_RHS</tt>
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
 * the overall solve status for an entire block returned by
 * <tt>blockSolveStatus[]</tt>.
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
 * forward solve operation in which case only a single virtual function
 * <tt>solve()</tt> must be overridden.  See <tt>LinearOpBase</tt> for what
 * other virtual functions must be overridden to define a concrete subclass.
 */
template <class RangeScalar, class DomainScalar = RangeScalar>
class LinearOpWithSolveBase
  : virtual public LinearOpBase<RangeScalar,DomainScalar>
  , virtual public Teuchos::VerboseObject<LinearOpWithSolveBase<RangeScalar,DomainScalar> >
{
public:

  /** \brief Local typedef for promoted scalar type.*/
  typedef typename Teuchos::PromotionTraits<RangeScalar,DomainScalar>::promote PromotedScalar;
  
  /** @name Pure virtual functions that must be overridden in subclasses */
  //@{

  /** \brief Request the forward solution of a block system with different targeted
   * solution criteria.
   *
   * \param  conj  [in] Determines if the elements are non-conjugate (<tt>NONCONJ_ELE</tt>) or
   *               conjugate (<tt>CONJ_ELE</tt>).  For real valued operator, this argument is meaningless.
   *               Most ANAs will request <tt>NONCONJ_ELE</tt>.
   * \param  B     [in] The RHS multi-vector with <tt>m = B.domain()->dim()</tt> columns.
   * \param  X     [in/out] The LHS multi-vector with with <tt>m = X->domain()->dim()</tt> columns.
   *               On input, contains the initial guess for the solution (only significant for
   *               iterative solvers) and on output contains an estimate of the solution.
   * \param  numBlocks
   *               [in] The number of blocks for which solve tolerances will be specified for.
   *               If <tt>numBlocks==0</tt> then this is a flag that a default set of tolerances
   *               should be used for all the linear systems. Default <tt>numBlocks=0</tt>.
   * \param  blockSolveCriteria
   *               [in] Array (length <tt>numBlocks</tt>) which gives the desired solution criteria
   *               for each of the <tt>numBlocks</tt> blocks of RHS.  If <tt>numBlocks>0</tt> then 
   *               this argument must be non-<tt>NULL</tt> and point to a valid array of <tt>numBlocks</tt>
   *               entries.
   * \param  blockSolveStatus
   *               [out] Array (length <tt>numBlocks</tt>) which gives the status of each set of 
   *               block systems.  A value of <tt>blockSolveStatus==NULL</tt> is allowed and
   *               means that the client is not interested in solution status of the linear systems.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->solveSupportsConj(conj)==true</tt>
   * <li><tt>X!=NULL</tt>
   * <li><tt>this->range()->isCompatible(*B.range())==true</tt>
   * <li><tt>this->domain()->isCompatible(*X->range())==true</tt>
   * <li><tt>B->domain()->isCompatible(*X->domain())==true</tt>
   * <li><tt>numBlocks >= 0</tt>
   * <li>[<tt>numBlocks == 0</tt>] <tt>blockSolveCriteria==NULL && blockSolveStatus==NULL</tt>.
   * <li>[<tt>numBlocks > 0</tt>] <tt>blockSolveCriteria!=NULL</tt> and points to an array of length
   *     at least <tt>numBlocks</tt>.
   * <li>[<tt>numBlocks > 0</tt>] <tt>this->solveSupportsSolveMeasureType(conj,blockSolveCriteria[k].solveMeasureType)==true</tt>,
   *     for <tt>k = 0...numBlocks-1</tt>.
   * <li>[<tt>blockSolveStatus!=NULL</tt>] <tt>blockSolveStatus</tt> and points to an array of length
   *     at least <tt>numBlocks</tt>.
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li>If any progress on solving the linear systems could be achieved, then 
   *     the function will return normally.  However, if no progress could be
   *     made in solving the linear systems, then an exception of type
   *     <tt>CatastrophicSolveFailure</tt> will be thrown.  If the function
   *     returns normally, then the following postconditions apply.
   * <li>If <tt>blockSolveStatus!=NULL</tt> then <tt>blockSolveStatus[k]</tt> gives the solution
   *     status of the block of linear systems specified by <tt>blockSolveCriteria[k]</tt>
   *     where, for <tt>k=0...numBlocks-1</tt>.
   * </ul>
   *
   * See the above introduction for a more complete description of how this
   * function behaves and the meaning of the arguments
   * <tt>blockSolveCriteria[]</tt> and <tt>blockSolveStatus[]</tt>.
   */
  virtual void solve(
    const EConj                                   conj
    ,const MultiVectorBase<RangeScalar>           &B
    ,MultiVectorBase<DomainScalar>                *X
    ,const int                                    numBlocks             = 0
    ,const BlockSolveCriteria<PromotedScalar>     blockSolveCriteria[]  = NULL
    ,SolveStatus<PromotedScalar>                  blockSolveStatus[]    = NULL
    ) const = 0;

  //@}

  /** @name Virtual functions with default implementations */
  //@{

  /** \brief Return if <tt>solve()</tt> supports the argument <tt>conj</tt>.
   *
   * The default implementation returns <tt>true</tt> for real valued scalar types
   * or when <tt>conj==NONCONJ_ELE</tt> for complex valued types.
   */
  virtual bool solveSupportsConj(EConj conj) const;

  /** \brief Return if <tt>solveTranspose()</tt> supports the argument <tt>conj</tt>.
   *
   * The default implementation returns <tt>false</tt>.
   */
  virtual bool solveTransposeSupportsConj(EConj conj) const;

  /** \brief Return if <tt>solve()</tt> supports the given the solve measure
   * type.
   *
   * The default implementation returns <tt>true</tt> for
   * <tt>solveMeasureType.inNone()</tt>.
   */
  virtual bool solveSupportsSolveMeasureType(EConj conj, const SolveMeasureType& solveMeasureType) const;

  /** \brief Return if <tt>solveTranspose()</tt> supports the given the solve
   * measure type.
   *
   * The default implementation returns <tt>true</tt> for
   * <tt>solveMeasureType.inNone()</tt>.
   */
  virtual bool solveTransposeSupportsSolveMeasureType(EConj conj, const SolveMeasureType& solveMeasureType) const;

  /** \brief Request the transpose (or adjoint) solution of a block system
   * with different targeted solution criteria.
   *
   * \param  conj  [in] Determines if the elements are non-conjugate (<tt>NONCONJ_ELE</tt>) or
   *               conjugate (<tt>CONJ_ELE</tt>).  For real valued operator, this argument is meaningless.
   *               The transpose solve is requested with <tt>conj==NONCONJ_ELE</tt> and the adjoint
   *               solve is requested with <tt>conj==CONJ_ELE</tt>.
   * \param  B     [in] The RHS multi-vector with <tt>m = B.domain()->dim()</tt> columns.
   * \param  X     [in/out] The LHS multi-vector with with <tt>m = X->domain()->dim()</tt> columns.
   *               On input, contains the initial guess for the solution (only significant for
   *               iterative solvers) and on output contains an estimate of the solution.
   * \param  numBlocks
   *               [in] The number of blocks for which solve tolerances will be specified for.
   *               If <tt>numBlocks==0</tt> then this is a flag that a default set of tolerances
   *               should be used for all the linear systems. Default <tt>numBlocks=0</tt>.
   * \param  blockSolveCriteria
   *               [in] Array (length <tt>numBlocks</tt>) which gives the desired solution criteria
   *               for each of the <tt>numBlocks</tt> blocks of RHS.  If <tt>numBlocks>0</tt> then 
   *               this argument must be non-<tt>NULL</tt> and point to a valid array of <tt>numBlocks</tt>
   *               entries.
   * \param  blockSolveStatus
   *               [out] Array (length <tt>numBlocks</tt>) which gives the status of each set of 
   *               block systems.  A value of <tt>blockSolveStatus==NULL</tt> is allowed and
   *               means that the client is not interested in solution status of the linear systems.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->solveTransposeSupportsConj(conj)==true</tt>
   * <li><tt>X!=NULL</tt>
   * <li><tt>this->domain()->isCompatible(*B.range())==true</tt>
   * <li><tt>this->range()->isCompatible(*X->range())==true</tt>
   * <li><tt>B->domain()->isCompatible(*X->domain())==true</tt>
   * <li><tt>numBlocks >= 0</tt>
   * <li>[<tt>numBlocks == 0</tt>] <tt>blockSolveCriteria==NULL && blockSolveStatus==NULL</tt>.
   * <li>[<tt>numBlocks > 0</tt>] <tt>blockSolveCriteria!=NULL</tt> and points to an array of length
   *     at least <tt>numBlocks</tt>.
   * <li>[<tt>numBlocks > 0</tt>] <tt>this->solveTransposeSupportsSolveMeasureType(conj,blockSolveCriteria[k].solveMeasureType)==true</tt>,
   *     for <tt>k = 0...numBlocks-1</tt>.
   * <li>[<tt>blockSolveStatus!=NULL</tt>] <tt>blockSolveStatus</tt> and points to an array of length
   *     at least <tt>numBlocks</tt>.
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li>If any progress on solving the linear systems could be achieved, then 
   *     the function will return normally.  However, if no progress could be
   *     made in solving the linear systems, then an exception of type
   *     <tt>CatastrophicSolveFailure</tt> will be thrown.  If the function
   *     returns normally, then the following postconditions apply.
   * <li>If <tt>blockSolveStatus!=NULL</tt> then <tt>blockSolveStatus[k]</tt> gives the solution
   *     status of the block of linear systems specified by <tt>blockSolveCriteria[k]</tt>
   *     where, for <tt>k=0...numBlocks-1</tt>.
   * </ul>
   *
   * See the above introduction for a more complete description of how this
   * function behaves and the meaning of the arguments
   * <tt>blockSolveCriteria[]</tt> and <tt>blockSolveStatus[]</tt>.
   */
  virtual void solveTranspose(
    const EConj                                   conj
    ,const MultiVectorBase<DomainScalar>          &B
    ,MultiVectorBase<RangeScalar>                 *X
    ,const int                                    numBlocks             = 0
    ,const BlockSolveCriteria<PromotedScalar>     blockSolveCriteria[]  = NULL
    ,SolveStatus<PromotedScalar>                  blockSolveStatus[]    = NULL
    ) const;

  //@}

};

/** \defgroup Thyra_LinearOpWithSolveBase_helper_grp Non-member LinearOpWithSolveBase helper functions.
 *
 * These functions allow for simpler calling sequences for solving linear
 * systems given a <tt>LinearOpWithSolveBase</tt> object.  In fact, these
 * functions help to document the various use cases associated with a
 * <tt>LinearOpWithSolveBase</tt> object.
 *
 * \ingroup Thyra_Op_Vec_Interoperability_Extended_Interfaces_grp
 */
//@{

/** \brief Determine if a <tt>LinearOpWithSolveBase<Scalar></tt> object
 * supports a particular type of solve or not..
 *
 * \ingroup Thyra_LinearOpWithSolveBase_helper_grp
 */
template<class Scalar>
bool solveSupports(
  const LinearOpWithSolveBase<Scalar>   &A
  ,const ETransp                        A_trans
  )
{
  if( real_trans(A_trans) == NOTRANS ) {
    return A.solveSupportsConj(
      A_trans == NOTRANS ? NONCONJ_ELE : CONJ_ELE
      );
  }
  return A.solveTransposeSupportsConj(
    A_trans == TRANS ? NONCONJ_ELE : CONJ_ELE
    );
}

/** \brief Solve a set of forward linear systems with a single set of
 * tolerances and a single scalar type.
 *
 * This function allows a single interface for accessing
 * <tt>LinearOpWithSolveBase::solve()</tt> or
 * <tt>LinearOpWithSolveBase::solveTranspose()</tt> for operators only
 * templated on a single scalar type for the domain and the range.
 *
 * See the implementation of this function for details.
 *
 * \ingroup Thyra_LinearOpWithSolveBase_helper_grp
 */
template<class Scalar>
SolveStatus<Scalar>
solve(
  const LinearOpWithSolveBase<Scalar>   &A
  ,const ETransp                        A_trans
  ,const MultiVectorBase<Scalar>        &B
  ,MultiVectorBase<Scalar>              *X
  ,const SolveCriteria<Scalar>          *solveCriteria
#ifndef __sun
                                                        = NULL
#endif
  )
{
  if(real_trans(A_trans)==NOTRANS)
    return solve(A,transToConj(A_trans),B,X,solveCriteria);
  return solveTranspose(A,transToConj(A_trans),B,X,solveCriteria);
}

#ifdef __sun

template<class Scalar>
SolveStatus<Scalar>
solve(
  const LinearOpWithSolveBase<Scalar>   &A
  ,const ETransp                        A_trans
  ,const MultiVectorBase<Scalar>        &B
  ,MultiVectorBase<Scalar>              *X
  )
{
  return solve(A,A_trans,B,X,NULL);
}

#endif // __sun

/** \brief Solve a set of forward linear systems with a single set of
 * tolerances.
 *
 * See the implementation of this function for details.
 *
 * \ingroup Thyra_LinearOpWithSolveBase_helper_grp
 */
template<class RangeScalar, class DomainScalar>
SolveStatus<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::PromotedScalar>
solve(
  const LinearOpWithSolveBase<RangeScalar,DomainScalar>   &A
  ,const EConj                                            conj
  ,const MultiVectorBase<RangeScalar>                     &B
  ,MultiVectorBase<DomainScalar>                          *X
  ,const SolveCriteria<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::PromotedScalar>
                                                          *solveCriteria
#ifndef __sun
                                                                          = NULL
#endif
  )
{
  typedef SolveCriteria<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::PromotedScalar>       SC;
  typedef BlockSolveCriteria<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::PromotedScalar>  BSC;
  typedef SolveStatus<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::PromotedScalar>         BSS;
  SC  defaultSolveCriteria;
  BSC blockSolveCriteria[1];
  BSS blockSolveStatus[1];
  blockSolveCriteria[0] = BSC(solveCriteria?*solveCriteria:defaultSolveCriteria,B.domain()->dim());
  A.solve(
    conj,B,X,1
    ,blockSolveCriteria
    ,blockSolveStatus
    );
  return blockSolveStatus[0];
}

#ifdef __sun

template<class RangeScalar, class DomainScalar>
SolveStatus<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::PromotedScalar>
solve(
  const LinearOpWithSolveBase<RangeScalar,DomainScalar>   &A
  ,const EConj                                            conj
  ,const MultiVectorBase<RangeScalar>                     &B
  ,MultiVectorBase<DomainScalar>                          *X
  )
{
  return solve(A,conj,B,X,NULL);
}

#endif // __sun

/** \brief Solve a set of transpose linear systems with a single set of
 * tolerances.
 *
 * See the implementation of this function for details.
 *
 * \ingroup Thyra_LinearOpWithSolveBase_helper_grp
 */
template <class RangeScalar, class DomainScalar>
SolveStatus<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::PromotedScalar>
solveTranspose(
  const LinearOpWithSolveBase<RangeScalar,DomainScalar>   &A
  ,const EConj                                            conj
  ,const MultiVectorBase<DomainScalar>                    &B
  ,MultiVectorBase<RangeScalar>                           *X
  ,const SolveCriteria<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::PromotedScalar>
                                                          *solveCriteria
#ifndef __sun
                                                                          = NULL
#endif
  )
{
  typedef SolveCriteria<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::PromotedScalar>       SC;
  typedef BlockSolveCriteria<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::PromotedScalar>  BSC;
  typedef SolveStatus<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::PromotedScalar>         BSS;
  SC  defaultSolveCriteria;
  BSC blockSolveCriteria[1];
  BSS blockSolveStatus[1];
  blockSolveCriteria[0] = BSC(solveCriteria?*solveCriteria:defaultSolveCriteria,B.domain()->dim());
  A.solveTranspose(
    conj,B,X,1
    ,blockSolveCriteria
    ,blockSolveStatus
    );
  return blockSolveStatus[0];
}

#ifdef __sun

template <class RangeScalar, class DomainScalar>
SolveStatus<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::PromotedScalar>
solveTranspose(
  const LinearOpWithSolveBase<RangeScalar,DomainScalar>   &A
  ,const EConj                                            conj
  ,const MultiVectorBase<DomainScalar>                    &B
  ,MultiVectorBase<RangeScalar>                           *X
  )
{
  return solveTranspose(A,conj,B,X,NULL);
}

#endif // __sun

/** \brief Solve a set of forward linear systems with two or more sets of
 * tolerances.
 *
 * See the implementation of this function for details.
 *
 * \ingroup Thyra_LinearOpWithSolveBase_helper_grp
 */
template<class RangeScalar, class DomainScalar>
void solve(
  const LinearOpWithSolveBase<RangeScalar,DomainScalar>   &A
  ,const EConj                                            conj
  ,const MultiVectorBase<RangeScalar>                     &B
  ,MultiVectorBase<DomainScalar>                          *X
  ,const int                                              numBlocks
  ,const BlockSolveCriteria<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::PromotedScalar>
                                                          blockSolveCriteria[]
#ifndef __sun
                                                                                = NULL
#endif
  ,SolveStatus<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::PromotedScalar>
                                                          blockSolveStatus[]
#ifndef __sun
                                                                                = NULL
#endif
  )
{
  A.solve(conj,B,X,numBlocks,blockSolveCriteria,blockSolveStatus);
}

#ifdef __sun

template<class RangeScalar, class DomainScalar>
void solve(
  const LinearOpWithSolveBase<RangeScalar,DomainScalar>   &A
  ,const EConj                                            conj
  ,const MultiVectorBase<RangeScalar>                     &B
  ,MultiVectorBase<DomainScalar>                          *X
  ,const int                                              numBlocks
  ,const BlockSolveCriteria<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::PromotedScalar>
                                                          blockSolveCriteria[]
  )
{
  solve(A,conj,B,X,numBlcoks,blockSolveCriteria,NULL);
}

template<class RangeScalar, class DomainScalar>
void solve(
  const LinearOpWithSolveBase<RangeScalar,DomainScalar>   &A
  ,const EConj                                            conj
  ,const MultiVectorBase<RangeScalar>                     &B
  ,MultiVectorBase<DomainScalar>                          *X
  ,const int                                              numBlocks
  )
{
  solve(A,conj,B,X,numBlcoks,NULL,NULL);
}

#endif // __sun

/** \brief Solve a set of transpose linear systems with two or more sets of
 * tolerances.
 *
 * See the implementation of this function for details.
 *
 * \ingroup Thyra_LinearOpWithSolveBase_helper_grp
 */
template <class RangeScalar, class DomainScalar>
void solveTranspose(
  const LinearOpWithSolveBase<RangeScalar,DomainScalar>   &A
  ,const EConj                                            conj
  ,const MultiVectorBase<DomainScalar>                    &B
  ,MultiVectorBase<RangeScalar>                           *X
  ,const int                                              numBlocks
  ,const BlockSolveCriteria<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::PromotedScalar>
                                                          blockSolveCriteria[]
#ifndef __sun
                                                                                = NULL
#endif
  ,SolveStatus<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::PromotedScalar>
                                                          blockSolveStatus[]
#ifndef __sun
                                                                                = NULL
#endif
  )
{
  A.solveTranspose(conj,B,X,numBlocks,blockSolveCriteria,blockSolveStatus);
}

#ifdef __sun

template <class RangeScalar, class DomainScalar>
void solveTranspose(
  const LinearOpWithSolveBase<RangeScalar,DomainScalar>   &A
  ,const EConj                                            conj
  ,const MultiVectorBase<DomainScalar>                    &B
  ,MultiVectorBase<RangeScalar>                           *X
  ,const int                                              numBlocks
  ,const BlockSolveCriteria<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::PromotedScalar>
                                                          blockSolveCriteria[]
  )
{
  solveTranspose(A,conj,B,X,numBlocks,blockSolveCriteria,NULL);
}

template <class RangeScalar, class DomainScalar>
void solveTranspose(
  const LinearOpWithSolveBase<RangeScalar,DomainScalar>   &A
  ,const EConj                                            conj
  ,const MultiVectorBase<DomainScalar>                    &B
  ,MultiVectorBase<RangeScalar>                           *X
  ,const int                                              numBlocks
  )
{
  solveTranspose(A,conj,B,X,numBlocks,NULL,NULL);
}

#endif // __sun

//@}

} // namespace Thyra

#endif // THYRA_LINEAR_OP_WITH_SOLVE_BASE_DECL_HPP





