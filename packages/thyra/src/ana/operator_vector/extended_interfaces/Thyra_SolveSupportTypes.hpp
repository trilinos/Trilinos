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

#ifndef THYRA_SOLVE_SUPPORT_TYPES_HPP
#define THYRA_SOLVE_SUPPORT_TYPES_HPP

#include "Thyra_OperatorVectorTypes.hpp"

namespace Thyra {

/** \defgroup Equation_solve_foundation_code_grp  Equation solve foundational code
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

/** \brief Simple struct that defines the requested solution criteria for a solve.
 *
 * \ingroup Equation_solve_foundation_code_grp
 */
template <class Scalar>
struct SolveCriteria {
  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
  /** \brief . */
  static const ScalarMag defaultTolerance() { return ScalarMag(-1); }
  /** \brief The type of solve tolerance requested. */
  ESolveTolType    solveTolType;
  /** \brief The requested solve tolerance (what the client would like to see).
   *
   * A value of <tt>defaultTolerance()</tt> means that the solver
   * implementation can define convergence any way it sees fit.
   */
  ScalarMag        requestedTol;  
  /** \brief . */
  SolveCriteria()
    : solveTolType(SOLVE_TOL_REL_RESIDUAL_NORM), requestedTol(defaultTolerance())
    {}
  /** \brief . */
  SolveCriteria(ESolveTolType _solveTolType, ScalarMag _requestedTol)
    : solveTolType(_solveTolType), requestedTol(_requestedTol)
    {}
};

/** \brief Simple struct that defines the requested solution criteria for a block solve.
 *
 * \ingroup Equation_solve_foundation_code_grp
 */
template <class Scalar>
struct BlockSolveCriteria {
  /** \brief Solve tolerance struct */
  SolveCriteria<Scalar>   solveCriteria;
  /** \brief Number of RHS that solve tolerance applies to. */
  int                     numRhs;
  /** \brief . */
  BlockSolveCriteria()
    : solveCriteria(), numRhs(1)
    {}
  /** \brief . */
  BlockSolveCriteria( const SolveCriteria<Scalar> &_solveCriteria, int _numRhs )
    : solveCriteria(_solveCriteria), numRhs(_numRhs)
    {}
};

/** \brief Solution status
 *
 * \ingroup Equation_solve_foundation_code_grp
 */
enum ESolveStatus {
  SOLVE_STATUS_CONVERGED        ///< The requested solution criteria has likely been achieved
  ,SOLVE_STATUS_UNCONVERGED     ///< The requested solution criteria has likely not been achieved
  ,SOLVE_STATUS_UNKNOWN         ///< The final solution status is unknown but he solve did not totally fail
};

/** \brief Simple struct for the return status from a solve.
 *
 * In the future, more fields may be added to aid in user diagnostics.
 *
 * \ingroup Equation_solve_foundation_code_grp
 */
template <class Scalar>
struct SolveStatus {
  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
  /** \brief . */
  static const ScalarMag unknownTolerance() { return ScalarMag(-1); }
  /** \brief The return status of the solve. */
  ESolveStatus solveStatus;
  /** \brief The maximum final tolerance actually achieved by the (block) linear solve.
   *
   * A value of <tt>unknownTolerance()</tt> means that even an estimate of the
   * the final value of the tolerance is unknown.
   */
  ScalarMag achievedTol;
  /** \brief The number of iterations taken.
   *
   * This number is totally implementation dependent and should only be used
   * for user diagnostics and not for any algorithmic purpose.
   */
  int numIterations;
  /** \brief Return message from the linear solver */
  std::string message;
  /** \brief . */
  SolveStatus()
    :solveStatus(SOLVE_STATUS_UNKNOWN), achievedTol(unknownTolerance())
     ,numIterations(0)
    {}
};

/** \brief Exception type thrown on an catastrophic solve failure.
 *
 * \ingroup Equation_solve_foundation_code_grp
 */
class CatastrophicSolveFailure : public std::runtime_error
{public: CatastrophicSolveFailure(const std::string& what_arg) : std::runtime_error(what_arg) {}};

} // namespace Thyra

#endif // THYRA_SOLVE_SUPPORT_TYPES_HPP
