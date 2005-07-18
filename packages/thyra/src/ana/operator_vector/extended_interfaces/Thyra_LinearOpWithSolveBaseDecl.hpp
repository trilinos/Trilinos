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

namespace Thyra {

/** \brief Base class for all linear operators that can support a high-level
 * solve operation.
 * 
 * <b>Introduction</b>
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
 * assumes is that for any appropriately selected random multi-vector
 * <tt>X</tt> which gives the operator application of

  \f[
    B \leftarrow op(A) X
  \f]

 * then a solve operation should be able to be performed that can recover
 * \f$X\f$ to some tolerance.  Note that this interface does not assume that a
 * solution can be achieved for any random RHS multi-vector \f$X\f$ but this
 * should be the case for non-singular operators.

 * It is recommended that clients use the non-member helper functions defined
 * \ref Thyra_LinearOpWithSolveBase_helper_grp "here" rather than call these
 * member functions directly as they support a number of other simpler use
 * cases.
 * 
 * <b>Solve Criteria</b>
 * 
 * This interface potentially allows clients to specify a relative tolerance
 * on either the relative residual norm or the relative norm of the solution
 * error and can target different solution criteria to different blocks of
 * linear system.  This interface tries to allow for mathematically rigorous
 * solution tolerances that are not based only any implementation-dependent
 * features like the number of iterations of some solver algorithm.  This
 * interface, however, allows <tt>*this</tt> operator to exclude support for
 * either or both types of solution tolerances (see the functions
 * <tt>solveSupportsSolveTolType()</tt> and
 * <tt>solveTransposeSupportsSolveTolType()</tt>).
 * 
 * This interface is meant to support direct and iterative linear solvers as
 * well as combinations of the two in a variety of configurations.  Because of
 * the almost infinite possible types of linear solver configurations
 * possible, this interface does not specify any particular solver-specific
 * types of control options, like maximum number of iterations.  These types
 * of control options can be specified in lower lever implementations but have
 * not place in this interface.
 *
 * The functions <tt>solve()</tt> and <tt>solveTranspose()</tt> both take the
 * arguments:

 \code
    ,const int                            numBlocks             = 0
    ,const BlockSolveCriteria<Scalar>     blockSolveCriteria[]  = NULL
    ,SolveStatus<Scalar>                  blockSolveStatus[]    = NULL
 \endcode
 *
 * The array arguments <tt>blockSolveCriteria[]</tt> and
 * <tt>blockSolveStatus[]</tt> specify different blocks of solution criteria
 * and the corresponding solve return statuses for a partitioned set of linear
 * systems.  Assuming that the client passes in arrays of dimensions
 * \f$N=\f$<tt>numBlocks</tt>, these tolerances define the solution criteria
 * for the block systems:
   
 \f[
   op(A) \left[ \begin{array}{ccccc} X_{(:,1:i_1)} & X_{(:,i_1+1:i_2)} & \ldots & X_{(:,i_{N-1}+1:i_N)} \end{array} \right]
   = \left[ \begin{array}{ccccc} B_{(:,1:i_1)} & B_{(:,i_1+1:i_2)} & \ldots & B_{(:,i_{N-1}+1:i_N)} \end{array} \right]
 \f]

 * where the column indexes are given by \f$i_j = \left( \sum_{k=1}^{j}
 * \mbox{blockSolveCriteria[k-1].numRhs} \right)\f$, for \f$j = 1 \ldots N\f$.
 *
 * The solve criteria for the \f$j^{\mbox{th}}\f$ system
   
 \f[
   op(A)  X_{(:,i_{j-1}+1:i_j)} = B_{(:,i_{j-1}+1:i_j)} 
 \f]

 * is given by <tt>blockSolveCriteria[j-1]</tt> (if
 * <tt>blockSolveCriteria!=NULL</tt>) and the solution status after return is
 * given by <tt>blockSolveStatus[j-1]</tt> (if
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
 * member and number of RHSs that it applies to.  It is the
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
 * <li> <b>Default tolerance</b> [
 * <tt>solveCriteria.requestedTol==SolveCriteria::defaultTolerance()</tt> ]:
 * In this mode, criteria for the solution tolerance is determined internally
 * by <tt>*this</tt> object.  Usually, it would be assumed that <tt>*this</tt>
 * would solve the linear systems to a sufficient tolerance for the given ANA
 * client.  This is the mode that many ANAs not designed to control
 * inexactness would work with.  In this mode, the solution criteria will be
 * determined in the end by the application or the user in an "appropriate"
 * manner.  In this case no meaningful solve status can be returned to the
 * client.
 * 
 * <li> <b>Residual tolerance</b> [
 * <tt>solveCriteria.requestedTol!=SolveCriteria::defaultTolerance() &&
 * solveCriteria.solveTolType==SOLVE_TOL_REL_RESIDUAL_NORM</tt> ]: In this
 * mode, the solution algorithm would be requested to solve the block system
 * to the relative residual tolerance of

  \f[
     \frac{|| op(A) X_{(:,i)} - B_{(:,i)} ||}{ ||B_{(:,i)}||} \le \mu_r,
  \f]

 * for \f$i=1 {}\ldots m\f$, where \f$\mu_r =\f$
 * <tt>solveCriteria.requestedTol</tt> and where the norm \f$||.||\f$ is
 * given by the natural norm defined by the range space of \f$op(A)\f$ and
 * computed from <tt>Thyra::norm()</tt>.  Many linear solvers should be able
 * to monitor this tolerance and be able to achieve it, or on failure be able
 * to report the actual tolerance achieved.
 * 
 * <li> <b>Solution error tolerance</b> [
 * <tt>solveCriteria.requestedTol!=SolveCriteria::defaultTolerance() &&
 * solveCriteria.solveTolType==SOLVE_TOL_REL_SOLUTION_ERR_NORM</tt> ]: In this
 * mode, the solution algorithm would be requested to solve the system to the
 * relative solution error tolerance of

  \f[
    \frac{|| X^*_{(:,i)} - X_{(:,i)} ||}{ ||X^*_{(:,i)}||} \le \mu_e,
  \f]
 
 * for \f$i=1 {}\ldots m\f$, where \f$||.||\f$, where \f$\mu_e =\f$
 * <tt>solveCriteria.requestedTol</tt> is the natural norm defined by the
 * domain space of \f$op(A)\f$ computed using <tt>Thyra::norm()</tt> and
 * \f$X^*\f$ is the true solution of the set of linear systems.  This is a
 * more difficult tolerance to monitor and achieve since the true solution in
 * almost never known.  Even in the base case, usually only an
 * order-of-magnitude estimate of this error will be known and can be reported
 * by the linear solver.
 * 
 * </ul>
 * 
 * <b>Solve Status</b>
 * 
 * After the <tt>solve()</tt> and <tt>solveTranspose()</tt> functions return,
 * the client can optionally get back a solution status for each block of
 * linear systems for of block solve criteria.  Specifically, for each block
 * of linear systems
  
  \f[
    A  X_{(:,i_{j-1}+1:i_j)} = B_{(:,i_{j-1}+1:i_j)} 
  \f]

 * whose solution criteria is specified by a <tt>SolveCriteria</tt> object, a
 * <tt>SolveStatus</tt> object can optionally be returned that lets the client
 * know the status of the linear solve.  If <tt>solveStatus</tt> is a
 * <tt>SolveStatus</tt> object returned for the above block linear system the
 * the following return status are significant:
 * 
 * <ul>
 * 
 * <li><b>Converged</b> [
 * <tt>solveStatus.solveStatus==SOLVE_STATUS_CONVERGED</tt> ]: This status is
 * returned by the linear solver if the solution criteria was likely achieved.
 * The maximum actual tolerance achieved may or may not be returned in the
 * field <tt>solveStatus.achievedTol</tt>.  The two sub-cases are:
 * 
 *   <ul>
 * 
 *   <li><b>Known tolerance</b> [ <tt>solveStatus.achievedTol >= 0</tt> ] :
 *   The linear solver knows the approximate order-of-magnitude estimate of
 *   the maximum tolerance achieved.  An order-of-magnitude (or so) estimate
 *   of the achieved tolerance would likely be known by any iterative linear
 *   solver where when
 *   <tt>solveCriteria.solveTolType==SOLVE_TOL_REL_RESIDUAL_NORM</tt>.
 * 
 *   <li><b>Unknown tolerance</b> [
 *   <tt>solveStatus.achievedTol==SolveStatus::unknownTolerance()</tt> ] : The
 *   linear solver does not know the tolerance that was achieved but the
 *   achieved tolerance should be very close to the requested tolerance.  This
 *   would be the most likely return status for a direct linear solver or for
 *   any linear solver where
 *   <tt>solveCriteria.solveTolType==SOLVE_TOL_REL_SOLUTION_ERR_NORM</tt>.
 * 
 *   </ul>
 * 
 * <li><b>Unconverged</b> [
 * <tt>solveStatus.solveStatus==SOLVE_STATUS_UNCONVERGED</tt> ]: The linear
 * solver was most likely not able to achieve the requested tolerance.  The
 * linear solver may not be able to return the actual tolerance achieved and
 * the same to cases as for the <it>unconverged</it> case are possible:
 * 
 *   <ul>
 * 
 *   <li><b>Known tolerance</b> [ <tt>solveStatus.achievedTol >= 0</tt> ] :
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
 * <tt>solveStatus.achievedTol==SolveStatus::unknownTolerance()</tt> will
 * always be returned.  This may be the return value when there is no
 * reasonable way that the linear solver algorithm can know now to compute or
 * estimate the requested tolerance.  This will also always be the return
 * status when
 * <tt>solveCriteria.requestedTol==SolveCriteria::defaultTolerance()</tt>
 * since the client would have no way to interpret this tolerance.
 * 
 * </ul>
 *
 * The implementation of the function <tt>accumulateSolveStatus()</tt> defines
 * how to accumulate the individual solve status for each RHS in a block into
 * the overall solve status for a block returned by
 * <tt>blockSolveStatus[]</tt>.
 * 
 * <b>Notes to subclass developers</b>
 * 
 * This interface assumes, by default, that subclasses will only support the
 * forward solve operation in which case only a single virtual function
 * <tt>solve()</tt> must be overridden.  See <tt>LinearOpBase</tt> for what
 * other virtual functions must be overridden to define a concrete subclass.
 * 
 * \ingroup Thyra_Op_Vec_Interoperability_Extended_Interfaces_grp
 */
template <class RangeScalar, class DomainScalar = RangeScalar>
class LinearOpWithSolveBase : virtual public LinearOpBase<RangeScalar,DomainScalar> {
public:

  /** \brief .*/
  typedef typename Teuchos::PromotionTraits<RangeScalar,DomainScalar>::promote  Scalar;
  
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
   * <li>[<tt>numBlocks > 0</tt>] <tt>this->solveSupportsSolveTolType(conj,blockSolveCriteria[k].solveTolType)==true</tt>,
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
    const EConj                           conj
    ,const MultiVectorBase<RangeScalar>   &B
    ,MultiVectorBase<DomainScalar>        *X
    ,const int                            numBlocks             = 0
    ,const BlockSolveCriteria<Scalar>     blockSolveCriteria[]  = NULL
    ,SolveStatus<Scalar>                  blockSolveStatus[]    = NULL
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

  /** \brief Return if <tt>solve()</tt> supports the solve tolerance type
   * <tt>ESolveTolType</tt>.
   *
   * The default implementation returns <tt>true</tt> for
   * <tt>solveTolType==SOLVE_TOL_REL_RESIDUAL_NORM</tt>.  Therefore, it is
   * assumed by default that the solver implementation will only be able to
   * check for and enforce a tolerance on the residual.
   */
  virtual bool solveSupportsSolveTolType(EConj conj, ESolveTolType solveTolType) const;

  /** \brief Return if <tt>solveTranspose()</tt> supports the solve tolerance
   * type <tt>ESolveTolType</tt>.
   *
   * The default implementation returns <tt>true</tt> for
   * <tt>solveTolType==SOLVE_TOL_REL_RESIDUAL_NORM</tt>.  Therefore, it is
   * assumed by default that the solver implementation will only be able to
   * check for and enforce a tolerance on the residual.
   */
  virtual bool solveTransposeSupportsSolveTolType(EConj conj, ESolveTolType solveTolType) const;

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
   * <li>[<tt>numBlocks > 0</tt>] <tt>this->solveTransposeSupportsSolveTolType(conj,blockSolveCriteria[k].solveTolType)==true</tt>,
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
    const EConj                           conj
    ,const MultiVectorBase<DomainScalar>  &B
    ,MultiVectorBase<RangeScalar>         *X
    ,const int                            numBlocks             = 0
    ,const BlockSolveCriteria<Scalar>     blockSolveCriteria[]  = NULL
    ,SolveStatus<Scalar>                  blockSolveStatus[]    = NULL
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
  ,const SolveCriteria<Scalar>          *solveCriteria = NULL
  )
{
  if(real_trans(A_trans)==NOTRANS)
    return solve(A,transToConj(A_trans),B,X,solveCriteria);
  return solveTranspose(A,transToConj(A_trans),B,X,solveCriteria);
}

/** \brief Solve a set of forward linear systems with a single set of
 * tolerances.
 *
 * See the implementation of this function for details.
 *
 * \ingroup Thyra_LinearOpWithSolveBase_helper_grp
 */
template<class RangeScalar, class DomainScalar>
SolveStatus<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::Scalar>
solve(
  const LinearOpWithSolveBase<RangeScalar,DomainScalar>   &A
  ,const EConj                                            conj
  ,const MultiVectorBase<RangeScalar>                     &B
  ,MultiVectorBase<DomainScalar>                          *X
  ,const SolveCriteria<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::Scalar>
                                                          *solveCriteria = NULL
  )
{
  typedef BlockSolveCriteria<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::Scalar>  BSC;
  typedef SolveStatus<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::Scalar>         BSS;
  BSC blockSolveCriteria[1];
  BSS blockSolveStatus[1];
  if(solveCriteria)
    blockSolveCriteria[0] = BSC(*solveCriteria,B.domain()->dim());
    A.solve(
      conj,B,X,solveCriteria?1:0
      ,solveCriteria ? blockSolveCriteria : NULL
      ,solveCriteria ? blockSolveStatus   : NULL
      );
  return blockSolveStatus[0];
}

/** \brief Solve a set of transpose linear systems with a single set of
 * tolerances.
 *
 * See the implementation of this function for details.
 *
 * \ingroup Thyra_LinearOpWithSolveBase_helper_grp
 */
template <class RangeScalar, class DomainScalar>
SolveStatus<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::Scalar>
solveTranspose(
  const LinearOpWithSolveBase<RangeScalar,DomainScalar>   &A
  ,const EConj                                            conj
  ,const MultiVectorBase<DomainScalar>                    &B
  ,MultiVectorBase<RangeScalar>                           *X
  ,const SolveCriteria<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::Scalar>
                                                          *solveCriteria = NULL
  )
{
  typedef BlockSolveCriteria<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::Scalar>  BSC;
  typedef SolveStatus<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::Scalar>         BSS;
  BSC blockSolveCriteria[1];
  BSS blockSolveStatus[1];
  if(solveCriteria)
    blockSolveCriteria[0] = BSC(*solveCriteria,B.domain()->dim());
    A.solveTranspose(
      conj,B,X,solveCriteria?1:0
      ,solveCriteria ? blockSolveCriteria : NULL
      ,solveCriteria ? blockSolveStatus   : NULL
      );
  return blockSolveStatus[0];
}

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
  ,const BlockSolveCriteria<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::Scalar>
                                                          blockSolveCriteria[]  = NULL
  ,SolveStatus<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::Scalar>
                                                          blockSolveStatus[]    = NULL
  )
{
  A.solve(conj,B,X,numBlocks,blockSolveCriteria,blockSolveStatus);
}

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
  ,const BlockSolveCriteria<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::Scalar>
                                                          blockSolveCriteria[]  = NULL
  ,SolveStatus<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::Scalar>
                                                          blockSolveStatus[]    = NULL
  )
{
  A.solveTranspose(conj,B,X,numBlocks,blockSolveCriteria,blockSolveStatus);
}

//@}

} // namespace Thyra

#endif // THYRA_LINEAR_OP_WITH_SOLVE_BASE_DECL_HPP
