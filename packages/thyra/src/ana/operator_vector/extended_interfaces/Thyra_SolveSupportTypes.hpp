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
  SOLVE_TOL_DEFAULT                     ///< Use some default solve criteria
  ,SOLVE_TOL_REL_RESIDUAL_NORM          ///< Enforce the tolerance on the relative natural norm in the residual vector
  ,SOLVE_TOL_REL_SOLUTION_ERR_NORM      ///< Enforce the tolerance on the relative natural norm in the error in the solution vector
};

/** \brief . */
inline
const char* toString(const ESolveTolType solveTolType)
{
  switch(solveTolType) {
    case SOLVE_TOL_DEFAULT: return"SOLVE_TOL_DEFAULT";
    case SOLVE_TOL_REL_RESIDUAL_NORM: return "SOLVE_TOL_REL_RESIDUAL_NORM";
    case SOLVE_TOL_REL_SOLUTION_ERR_NORM: return "SOLVE_TOL_REL_SOLUTION_ERR_NORM";
    default: TEST_FOR_EXCEPT(true);
  }
  return ""; // Never be called!
}

/** \brief Simple struct that defines the requested solution criteria for a solve.
 *
 * \ingroup Equation_solve_foundation_code_grp
 */
template <class Scalar>
struct SolveCriteria {
  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
  /** \brief . */
  static const ScalarMag unspecifiedTolerance() { return ScalarMag(-1); }
  /** \brief . */
  static const int unspecifiedMaxIterations() { return -1; }
  /** \brief The type of solve tolerance requested as given in
   * <tt>this->requestedTol</tt>.
   */
  ESolveTolType    solveTolType;
  /** \brief The requested solve tolerance (what the client would like to see).
   *
   * Only significant if <tt>this->solveTolType!=SOLVE_TOL_DEFAULT</TT>
   */
  ScalarMag        requestedTol;
  /** \brief The maximum number of iterations the solver is allowed to take.
   *
   * Note that the interpretation of this integer is 100% implementation
   * defined and should not be used in the most general setting.
   *
   * This argument may be significant regardless of the value of <tt>this->solveTolType</tt>.
   */
  int maxIterations;
  /** \brief Default construction is for a default solve. */
  SolveCriteria()
    :solveTolType(SOLVE_TOL_DEFAULT)
    ,requestedTol(unspecifiedTolerance())
     ,maxIterations(unspecifiedMaxIterations())
    {}
  /** \brief Construct with a specified solve criteria. */
  SolveCriteria(ESolveTolType _solveTolType, ScalarMag _requestedTol, int _maxIterations = unspecifiedMaxIterations())
    : solveTolType(_solveTolType), requestedTol(_requestedTol), maxIterations(_maxIterations)
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

/** \brief Exception type thrown on an catastrophic solve failure.
 *
 * \ingroup Equation_solve_foundation_code_grp
 */
class CatastrophicSolveFailure : public std::runtime_error
{public: CatastrophicSolveFailure(const std::string& what_arg) : std::runtime_error(what_arg) {}};

/** \brief Solution status
 *
 * \ingroup Equation_solve_foundation_code_grp
 */
enum ESolveStatus {
  SOLVE_STATUS_CONVERGED        ///< The requested solution criteria has likely been achieved
  ,SOLVE_STATUS_UNCONVERGED     ///< The requested solution criteria has likely not been achieved
  ,SOLVE_STATUS_UNKNOWN         ///< The final solution status is unknown but he solve did not totally fail
};

/** \brief . */
inline
const char* toString(const ESolveStatus solveStatus)
{
  switch(solveStatus) {
    case SOLVE_STATUS_CONVERGED:    return "SOLVE_STATUS_CONVERGED";
    case SOLVE_STATUS_UNCONVERGED:  return "SOLVE_STATUS_UNCONVERGED";
    case SOLVE_STATUS_UNKNOWN:      return "SOLVE_STATUS_UNKNOWN";
    default: TEST_FOR_EXCEPT(true);
  }
  return ""; // Never be called!
}

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
  int iterations;
  /** \brief A simple one-line message (i.e. no newlines) returned from the solver */
  std::string message;
  /** \brief . */
  SolveStatus()
    :solveStatus(SOLVE_STATUS_UNKNOWN), achievedTol(unknownTolerance())
     ,iterations(1)
    {}

  /** \brief Output the achieveTol field.
   */
  static std::string achievedTolToString( const ScalarMag &achievedTol )
    {
      if(achievedTol==unknownTolerance()) return "unknownTolerance()";
      std::ostringstream oss; oss << achievedTol; return oss.str();
    }
};

/** \brief Enum for defining the status of a preconditioner object. */
enum EPreconditionerInputType {
  PRECONDITIONER_INPUT_TYPE_AS_OPERATOR  ///< The input preconditioner should just be applied as an operator
  ,PRECONDITIONER_INPUT_TYPE_AS_MATRIX   ///< The input preconditioner should viewed as a matrix to be factored then backsolved as a preconditioner
};

/** \brief Accumulate solve status objects for solving a block of RHSs is
 * smaller sub-blocks..
 *
 *
 * \param  overallSolveCriteria  [in] The overall solve criteria for the overall blocks.
 * \param  solveStatus           [in] The solve status for a sub-block (or a single RHS) 
 * \param  overallSolveStatus    [in/out] The accumulated solve status for all the
 *                               sub-blocks of RHS.
 *
 * On the first call, set <tt>overallSolveStatus->solveStatus =
 * SOLVE_STATUS_CONVERGED</tt> and <tt>overallSolveStatus->iterations =
 * 0</tt>!
 *
 * \ingroup Equation_solve_foundation_code_grp
 */
template <class Scalar>
void accumulateSolveStatus(
  const SolveCriteria<Scalar>    &overallSolveCriteria
  ,const SolveStatus<Scalar>     &solveStatus
  ,SolveStatus<Scalar>           *overallSolveStatus
  )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(overallSolveStatus==NULL);
#endif
  switch(solveStatus.solveStatus) {
    case SOLVE_STATUS_UNCONVERGED:
    {
      // First, if we see any unconverged solve status, then the entire block is
      // unconverged!
      overallSolveStatus->solveStatus = SOLVE_STATUS_UNCONVERGED;
      overallSolveStatus->message = solveStatus.message;
      break;
    }
    case SOLVE_STATUS_UNKNOWN:
    {
      // Next, if any solve status is unknown, then if the overall solve status
      // says converged, then we have to mark it as unknown.  Note than unknown
      // could mean that the system is actually converged!
      switch(overallSolveStatus->solveStatus) {
        case SOLVE_STATUS_CONVERGED:
          overallSolveStatus->solveStatus = SOLVE_STATUS_UNKNOWN;
          break;
        case SOLVE_STATUS_UNCONVERGED:
        case SOLVE_STATUS_UNKNOWN:
          // If we get here then the overall solve status is either unknown
          // already or says unconverged and this will not change here!
          overallSolveStatus->message = solveStatus.message;
          break;
        default:
          TEST_FOR_EXCEPT(true); // Corrupted enum?
      }
      break;
    }
    case SOLVE_STATUS_CONVERGED:
    {
      // If we get here then the overall solve status is either unknown,
      // unconverged, or converged and this will not change here!
      if(overallSolveStatus->message == "")
        overallSolveStatus->message = solveStatus.message;
      break;
    }
    default:
      TEST_FOR_EXCEPT(true); // Corrupted enum?
  }
  // Update the achieved tolerence to the maximum returned
  if( solveStatus.achievedTol > overallSolveStatus->achievedTol ) {
    overallSolveStatus->achievedTol = solveStatus.achievedTol;
  }
  // Accumulate the total number of iterations
  overallSolveStatus->iterations += solveStatus.iterations;
  // Set a message if non is set
  if(overallSolveStatus->message == "")
    overallSolveStatus->message = solveStatus.message;
}

} // namespace Thyra

#endif // THYRA_SOLVE_SUPPORT_TYPES_HPP
