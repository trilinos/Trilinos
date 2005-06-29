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

namespace Thyra {

/** \defgroup Equation_solve_foundation_code_grp
 *
 * \ingroup Thyra_Op_Vec_Interoperability_Extended_Interfaces_grp
 */

/** \brief The type of the tolerance for the solve
 *
 * \ingroup Equation_solve_foundation_code_grp
 */
enum ESolveTolType {
  SOLVE_TOL_REL_RESIDUAL_NORM           ///< Enforce the tolerance on the relative natural norm in the residual vector
  ,SOLVE_TOL_REL_SOLUTION_ERR_NORM      ///< Enforce the tolerance on the relative natural norm in the error in the solution vector
};

/** \brief Simple struct for a solve tolerance for a requested solve.
 *
 * \ingroup Equation_solve_foundation_code_grp
 */
template <class Scalar>
struct SolveTolerance {

  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \brief . */
  static const ScalarMag  DEFAULT_SOLVE_TOLERANCE;

  /** \brief The requested solve tolerance (what the client would like to see).
   *
   * A value of DEFAULT_SOLVE_TOLERANCE means that the solver implementation can
   * define convergence any way it sees fit.
   */
  ScalarMag      requestedTol;  

  /** \brief The type of solve tolerance.
   */
  ESolveTolType    solveTolType;

  /** \brief . */
  SolveTolerance()
    : requestedTol(DEFAULT_SOLVE_TOLERANCE), solveTolType(SOLVE_TOL_REL_RESIDUAL_NORM)
    {}

  /** \brief . */
  SolveTolerance(ScalarMag _requestedTol, ESolveTolType _solveTolType)
    : requestedTol(_requestedTol), solveTolType(_solveTolType)
    {}

};

/** \brief Simple struct for a solve tolerance of a requested block solve.
 *
 * \ingroup Equation_solve_foundation_code_grp
 */
template <class Scalar>
struct BlockSolveTolerance {

  /** \brief Solve tolerance struct */
  SolveTolerance<Scalar>   solveTolerance;

  /** \brief Number of RHS that solve tolerance applies to. */
  int                      numRhs;

  /** \brief . */
  BlockSolveTolerance()
    : solveTolerance(), numRhs(1)
    {}

  /** \brief . */
  BlockSolveTolerance( const SolveTolerance<Scalar> &_solveTolerance, int _numRhs )
    : solveTolerance(_solveTolerance), numRhs(_numRhs)
    {}
  
};

/** \brief The type of the solution returned.
 *
 * \ingroup Equation_solve_foundation_code_grp
 */
enum ESolveReturnStatus {
  SOLVE_STATUS_CONVERGED        ///< The requested tolerance has been achieved
  ,SOLVE_STATUS_UNCONVERGED     ///< The requested tolerance has not been achieved
  ,SOLVE_STATUS_UNKNOWN         ///< The final solution tolerance is unknown
};

/** \brief Simple struct for the return status from a solve.
 *
 * In the future, more fields may be added to aid in user diagnostics.
 *
 * \ingroup Equation_solve_foundation_code_grp
 */
template <class Scalar>
struct SolveReturn {

  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \brief . */
  static const ScalarMag  UNKNOWN_SOLVE_TOLERANCE;

  /** \brief The actual tolerance achieved.
   *
   * The value of UNKNOWN_SOLVE_TOLERANCE means that even an estimate of the
   * the final value of the tolerance is unknown.
   */
  ScalarMag             achievedTol;

  /** \brief The status of the solve. */
  ESolveReturnStatus    solveReturnStatus;

  /** \brief The number of iterations taken.
   *
   * This number is total implementation dependent and it only to be used for
   * user diagnostics and not for any algorithmic purpose.
   */
  int                   numIterations;

  SolveReturn()
    :achievedTol(UNKNOWN_SOLVE_TOLERANCE),solveReturnStatus(SOLVE_STATUS_UNKNOWN)
     ,numIterations(0)
    {}

};

/** \brief Exception type thrown on an catastrophic solve failure
 *
 * \ingroup Equation_solve_foundation_code_grp
 */
class CatastrophicSolveFailure;


/** \brief Base class for all linear operators that can support a high-level
 * solve operation.
 *
 * <b>Introduction</b>
 *
 * This interface supports operators with potentially different range and
 * domain scalar types that can support a solve operation of the form:
 
 \verbatim

   A*X = B 

 \endverbatim

 * where <tt>M</tt> is <tt>*this</tt> linear operator, <tt>B</tt> is an
 * appropriate multi-vector and <tt>X</tt> is a multi-vector that is computed
 * by this interface.
 *
 * Note that this interface does not assume that the linear operator itself is
 * nonsingular or invertable in the classic sence (i.e. an inverse operator
 * may not exist).
 *
 * What this interface assumes is that for any appropriately selected random
 * multi-vector <tt>X</tt> which gives the forward application of

 \verbatim

   B = A*X 

 \endverbatim

 * then a solve operation should be able to be performed what can recover
 * <tt>X</tt> to some tolerance.  Note that this interface does not assume
 * that a solution can be achieved for any random RHS multi-vector <tt>B</tt>.
 *
 * It is recommended that clients use the non-member helper functions defined
 * \ref Thyra_LinearOpWithSolveBase_helper_grp "here" rather than call these
 * member functions directly as they support a number of other simplier use
 * cases.
 *
 * Forward solves are performed using the
 *
 * In addition to the forward solve, this interface may also support the the
 * transpose, and/or the adjoint solve through the function
 * <tt>solveTranspose()</tt>.
 *
 * <b>Solve Tolerances</b>
 *
 * This interface allows clients to specify a relative tolerance on either the
 * relative residual norm or the relative norm of the solution error.
 *
 * ToDo: Finish documentation!
 *
 * <b>Use cases:</b>
 *
 * Here present examples of several different use cases involving this
 * interface.
 *
 * <ul>
 *
 * <li><b>Standard forward solves with single scalar type</b>:
 *
 *     These are examples of some forward solves using non-conjugate elements which
 *     is expected to be the most typical mode of use.
 *
 * <ul>
 *
 * <li><b>Standard forward solve using default tolerance</b>
 
 \code

    template <class Scalar>
    void solveUsingDefaultTol(
      const LinearOpWithSolveBase<Scalar>    &A
      ,const MultiVectorBase<Scalar>         &B
      ,MultiVectorBase<Scalar>               *X
      )
    {
      SolveReturn<Scalar>
        solveReturn = Thyra::solve(A,NONCONJ_ELE,B,X);
      // Check return
      ...
    }

 \endcode
 
 * <li><b>Standard forward solve based on residual norm</b>
 
 \code

    template <class Scalar>
    void solveUsingResidualTol(
      const LinearOpWithSolve<Scalar>                      &A
      ,const MultiVectorBase<Scalar>                       &B
      ,const Teuchos::ScalarTraits<Scalar>::magnitudeType  tol
      ,MultiVectorBase<Scalar>                             *X
      )
    {
      SolveTolerance<Scalar> solveTolerance(tol,SOLVE_TOL_REL_RESIDUAL_NORM);
      SolveReturn<Scalar>
        solveReturn = Thyra::solve(A,NONCONJ_ELE,B,X,&solveTolerance);
      // Check return
      ...
    }

 \endcode

 * <li><b>Standard forward solve based on solution error norm</b>

  \code

    template <class Scalar>
    void solveUsingSolutionTol(
      const LinearOpWithSolve<Scalar>                      &A
      ,const MultiVectorBase<Scalar>                       &B
      ,const Teuchos::ScalarTraits<Scalar>::magnitudeType  tol
      ,MultiVectorBase<Scalar>                             *X
      )
    {
      SolveTolerance<Scalar> solveTolerance(tol,SOLVE_TOL_REL_SOLUTION_ERR_NORM);
      SolveReturn<Scalar>
        solveReturn = Thyra::solve(A,NONCONJ_ELE,B,X,&solveTolerance);
      // Check return
      ...
    }

 \endcode

 * <li><b>Standard forward solve given two sets of tolerances</b>

  \code

    template <class Scalar>
    void solveTwoBlockSystems(
      const LinearOpWithSolve<Scalar>                      &A
      ,const MultiVectorBase<Scalar>                       &B1
      ,const SolveTolerance<Scalar>                        solveTolerance1
      ,const MultiVectorBase<Scalar>                       &B2
      ,const SolveTolerance<Scalar>                        solveTolerance2
      ,MultiVectorBase<Scalar>                             *X1
      ,MultiVectorBase<Scalar>                             *X2
      )
    {
      Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
        B = concatColumns( Teuchos::rcp(&B1,false), Teuchos::rcp(&B2,false) );
      Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
        X = concatColumns( Teuchos::rcp(X1,false), Teuchos::rcp(&X2,false) );
      BlockSolveTolerance<Scalar>
        blockTolerances[2]
          = {
              BlockSolveTolerance<Scalar>(solveTolerance1,B1.domain()->dim()),
              BlockSolveTolerance<Scalar>(solveTolerance2,B2.domain()->dim())
            };
      SolveReturn<Scalar>
        solveReturn = Thyra::solve(A,NONCONJ_ELE,*B,&*X,blockSolveTolerances);
      // Check return
      ...
    }

 \endcode

 * </ul>
 *
 * </ul>
 *
 * <b>Notes to subclass develoeprs</b>
 *
 * This interface assumes that by default, that subclasses will only support
 * the forward solve operation in which case only a single addition virtual
 * function <tt>solve()</tt> must be overridden.  See <tt>LinearOpBase</tt>
 * for what other virtual functions must be overridden to define a concrete
 * subclass.
 *
 * \ingroup Thyra_Op_Vec_Interoperability_Extended_Interfaces_grp
 */
template <class RangeScalar, class DomainScalar = RangeScalar>
class LinearOpWithSolveBase : virtual public LinearOpBase<RangeScalar,DomainScalar> {
public:

  /** \brief .*/
  typedef typename LinearOpBase<RangeScalar,DomainScalar>::Scalar Scalar;
  
  /** @name Pure virtual functions that must be overridden in subclasses */
  //@{

  /** \brief Solve (or try to solve) a block system with different targeted tolerances.
   *
   * \param  conj  [in] Determines if the elements are non-conjugate (NONCONJ_ELE) or
   *               conjugate (CONJ_ELE).  For real valued operator, this argument is meaningless
   * \param  B     [in] Contains the RHS multi-vector.
   * \param  X     [in/out] On input, contains the initial guess for the solution (only significant for
   *               iterative solvers) and on output contains an estimate of the solution.
   * \param  numBlocks
   *               [in] The number of blocks that solve tolerances will be specified for.
   *               Default = 1.
   * \param  blockSolveTolerances
   *               [in] Array (length numBlocks) which gives the desired solution criteria
   *               for each of the <tt>numBlocks</tt> blocks of RHS.  A value of <tt>blockSolveTolerances==NULL</tt>
   *               means that a default set of tolerances will be used.
   *
   * Preconditions:<ul>
   * <li><tt>solveSupports(conj)==true</tt> (thorw <tt>OpNotSupported</tt>
   * <li> ...
   * </ul>
   *
   * Postconditions:<ul>
   * <li> ...
   * </ul>
   *
   * \return Returns a <tt>SolveReturn</tt> struct object that gives the
   * status of the solution.
   *
   * With throw <tt>CatastrophicSolveFailure</tt> on complete failure!
   *
   * See the introduction above for a description of how this function
   * behaves.
   */
  virtual SolveReturn<Scalar> solve(
    const EConj                           conj
    ,const MultiVectorBase<RangeScalar>   &B
    ,MultiVectorBase<DomainScalar>        *X
    ,const int                            numBlocks               = 1
    ,const BlockSolveTolerance<Scalar>    blockSolveTolerances[]  = NULL
    ) const = 0;

  //@}

  /** @name Virtual functions with default implementations */
  //@{

  /** \brief Return if <tt>solve()</tt> supports the argument <tt>conj</tt>.
   *
   * The default implementation returns <tt>true</tt> for real valued scalar types
   * or when <tt>conj==NONCONJ_ELE</tt> for complex valued types.
   */
  virtual bool solveSupports(EConj conj) const;

  /** \brief Return if <tt>solveTranspose()</tt> supports the argument <tt>conj</tt>.
   *
   * The default implementation returns <tt>false</tt>.
   */
  virtual bool solveTransposeSupports(EConj conj) const;

  /** \brief Transpose solve (or try to solve) a block system with different targeted tolerances.
   *
   * \param  conj  [in] Determines if the elements are non-conjugate (NONCONJ_ELE) or
   *               conjugate (CONJ_ELE).  For real valued operator, this argument is meaningless
   * \param  B     [in] Contains the RHS multi-vector.
   * \param  X     [in/out] On input, contains the initial guess for the solution (only significant for
   *               iterative solvers) and on output contains an estimate of the solution.
   * \param  numBlocks
   *               [in] The number of blocks that solve tolerances will be specified for.
   *               Default = 1.
   * \param  blockSolveTolerances
   *               [in] Array (length numBlocks) which gives the desired solution criteria
   *               for each of the <tt>numBlocks</tt> blocks of RHS.  A value of <tt>blockSolveTolerances==NULL</tt>
   *               means that a default set of tolerances will be used.
   *
   * Preconditions:<ul>
   * <li><tt>solveTransposeSupports(conj)==true</tt> (thorw <tt>OpNotSupported</tt>
   * <li> ...
   * </ul>
   *
   * Postconditions:<ul>
   * <li> ...
   * </ul>
   *
   * \return Returns a <tt>SolveReturn</tt> struct object that gives the
   * status of the solution.
   *
   * With throw <tt>CatastrophicSolveFailure</tt> on complete failure!
   *
   * See the introduction above for a description of how this function
   * behaves.
   *
   * With throw <tt>CatastrophicSolveFailure</tt> on complete failure!
   *
   * The default implementation throws an exception with a very good error
   * message.  This is consistent with the default implementation of
   * <tt>solveTransposeSupports()</tt> which returns <tt>false</tt>.
   */
  virtual SolveReturn<Scalar> solveTranspose(
    const EConj                           conj
    ,const MultiVectorBase<DomainScalar>  &B
    ,MultiVectorBase<RangeScalar>         *X
    ,const int                            numBlocks               = 1
    ,const BlockSolveTolerance<Scalar>    blockSolveTolerances[]  = NULL
    ) const;

  //@}

};

/** \defgroup Thyra_LinearOpWithSolveBase_helper_grp Non-member LinearOpWithSolveBase helper functions.
 *
 * These functions allow for simpler calling sequences for solving linear systems given
 * a <tt>LinearOpWithSolveBase</tt> object.
 *
 * \ingroup Thyra_Op_Vec_Interoperability_Extended_Interfaces_grp
 */
//@{

/** \brief Solve a set of forward linear systems with a single set of tolerances.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_LinearOpWithSolveBase_helper_grp
 */
template<class RangeScalar, class DomainScalar>
SolveReturn<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::Scalar>
solve(
  const LinearOpWithSolveBase<RangeScalar,DomainScalar>   &A
  ,const EConj                                            conj
  ,const MultiVectorBase<RangeScalar>                     &B
  ,MultiVectorBase<DomainScalar>                          *X
  ,const SolveTolerance<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::Scalar>
                                                          *solveTolerance = NULL
  )
{
  typedef BlockSolveTolerance<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::Scalar>  BST;
  if(solveTolerance) {
    BST blockSolveTolerances[] = { BST(*solveTolerance,B.domain()->dim()) };
    return A.solve(conj,B,X,1,blockSolveTolerances);
  }
  return A.solve(conj,B,X);
}

/** \brief Solve a set of transpose linear systems with a single set of tolerances.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_LinearOpWithSolveBase_helper_grp
 */
template <class RangeScalar, class DomainScalar>
SolveReturn<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::Scalar>
solveTranspose(
  const LinearOpWithSolveBase<RangeScalar,DomainScalar>   &A
  ,const EConj                                            conj
  ,const MultiVectorBase<DomainScalar>                    &B
  ,MultiVectorBase<RangeScalar>                           *X
  ,const SolveTolerance<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::Scalar>
                                                          *solveTolerance = NULL
  )
{
  typedef BlockSolveTolerance<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::Scalar>  BST;
  if(solveTolerance) {
    BST blockSolveTolerances[] = { BST(*solveTolerance,B.domain()->dim()) };
    return A.solveTranspose(conj,B,X,1,blockSolveTolerances);
  }
  return A.solveTranspose(conj,B,X);
}

/** \brief Solve a set of forward linear systems with two or more sets of tolerances.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_LinearOpWithSolveBase_helper_grp
 */
template<class RangeScalar, class DomainScalar>
SolveReturn<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::Scalar>
solve(
  const LinearOpWithSolveBase<RangeScalar,DomainScalar>   &A
  ,const EConj                                            conj
  ,const MultiVectorBase<RangeScalar>                     &B
  ,MultiVectorBase<DomainScalar>                          *X
  ,const int                                              numBlocks
  ,const BlockSolveTolerance<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::Scalar>
                                                          blockSolveTolerances[]  = NULL
  )
{
  A.solve(conj,B,X,numBlocks,blockSolveTolerances);
}

/** \brief Solve a set of transpose linear systems with two or more sets of tolerances.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_LinearOpWithSolveBase_helper_grp
 */
template <class RangeScalar, class DomainScalar>
SolveReturn<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::Scalar>
solveTranspose(
  const LinearOpWithSolveBase<RangeScalar,DomainScalar>   &A
  ,const EConj                                            conj
  ,const MultiVectorBase<DomainScalar>                    &B
  ,MultiVectorBase<RangeScalar>                           *X
  ,const int                                              numBlocks
  ,const BlockSolveTolerance<typename LinearOpWithSolveBase<RangeScalar,DomainScalar>::Scalar>
                                                          blockSolveTolerances[]  = NULL
  )
{
  A.solveTranspose(conj,B,X,numBlocks,blockSolveTolerances);
}

//@}

} // namespace Thyra

#endif // THYRA_LINEAR_OP_WITH_SOLVE_BASE_DECL_HPP
