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
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_FancyOStream.hpp"


namespace Thyra {


/** \brief Type of solve measure norm.
 *
 * For reference we refer to solving a single linear system <tt>A*x=b</tt>.
 *
 * \ingroup Thyra_Op_Solve_fundamental_interfaces_code_grp
 */
enum ESolveMeasureNormType {
  SOLVE_MEASURE_ONE                   ///< No solve measure (i.e. same as 1.0)
  ,SOLVE_MEASURE_NORM_RESIDUAL        ///< Norm of the current residual vector (i.e. <tt>||A*x-b||</tt>)
  ,SOLVE_MEASURE_NORM_SOLUTION        ///< Norm of the current solution vector (i.e. <tt>||x||</tt>)
  ,SOLVE_MEASURE_NORM_INIT_RESIDUAL   ///< Norm of the initial residual vector given a non-zero guess (i.e. <tt>||A*xo-b||</tt>)
  ,SOLVE_MEASURE_NORM_RHS             ///< Norm of the Right-hand side (i.e. <tt>||b||</tt>)
};


/** \brief .
 *
 * \ingroup Thyra_Op_Solve_fundamental_interfaces_code_grp
 */
inline
const std::string toString(const ESolveMeasureNormType solveMeasureNormType)
{
  switch(solveMeasureNormType) {
    case SOLVE_MEASURE_ONE:
      return "SOLVE_MEASURE_ONE";
    case SOLVE_MEASURE_NORM_RESIDUAL:
      return "SOLVE_MEASURE_NORM_RESIDUAL";
    case SOLVE_MEASURE_NORM_SOLUTION:
      return "SOLVE_MEASURE_NORM_SOLUTION";
    case SOLVE_MEASURE_NORM_INIT_RESIDUAL:
      return "SOLVE_MEASURE_NORM_INIT_RESIDUAL";
    case SOLVE_MEASURE_NORM_RHS:
      return "SOLVE_MEASURE_NORM_RHS";
    default:
      TEST_FOR_EXCEPT(true);
  }
  return NULL; // Never be called!
}


/** \brief Solve tolerance type.
 *
 * This represents the solve tolerance measure of the form:
 \verbatim
 (numerator)/(denominator)
 \endverbatim
 *
 * Note that
 * <tt>numerator==SOLVE_MEASURE_ONE&&denominator==SOLVE_MEASURE_ONE</tt>
 * (i.e. 1/1) means that there is no solve measure type specified.
 *
 * \ingroup Thyra_Op_Solve_fundamental_interfaces_code_grp
 */
struct SolveMeasureType {
  /** \brief . */
  ESolveMeasureNormType  numerator;
  /** \brief . */
  ESolveMeasureNormType  denominator;
  /** \brief . */
  SolveMeasureType()
    :numerator(SOLVE_MEASURE_ONE),denominator(SOLVE_MEASURE_ONE)
    {}
  /** \brief . */
  SolveMeasureType(ESolveMeasureNormType _numerator, ESolveMeasureNormType _denominator)
    :numerator(_numerator),denominator(_denominator)
    {}
  /** \brief . */
  void set(ESolveMeasureNormType _numerator, ESolveMeasureNormType _denominator)
    { numerator = _numerator; denominator = _denominator; }
  /** \brief . */
  bool useDefault() const
    { return ( numerator==SOLVE_MEASURE_ONE && denominator==SOLVE_MEASURE_ONE ); }
  /** \brief . */
  bool operator()(ESolveMeasureNormType _numerator, ESolveMeasureNormType _denominator) const
    { return ( numerator==_numerator && denominator==_denominator ); }
};


/** \brief .
 *
 * \ingroup Thyra_Op_Solve_fundamental_interfaces_code_grp
 */
inline
std::string toString(const SolveMeasureType& solveMeasureType)
{
  std::ostringstream oss;
  oss << "("<<toString(solveMeasureType.numerator)<<")/("<<solveMeasureType.denominator<<")";
  return oss.str();
}


/** \brief Simple struct that defines the requested solution criteria for a solve.
 *
 * \ingroup Thyra_Op_Solve_fundamental_interfaces_code_grp
 */
template <class Scalar>
struct SolveCriteria {
  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
  /** \brief . */
  static ScalarMag unspecifiedTolerance() { return ScalarMag(-1); }
  /** \brief The type of solve tolerance requested as given in
   * <tt>this->requestedTol</tt>. */
  SolveMeasureType solveMeasureType;
  /** \brief The requested solve tolerance (what the client would like to see).
   * Only significant if <tt>!this->solveMeasureType.useDefault()</tt> */
  ScalarMag requestedTol;
  /** \brief Any extra control parameters.
   * Note that the contents of this parameter list is totally undefined
   * and any client that uses this does so at their own peril! */
  Teuchos::RCP<Teuchos::ParameterList> extraParameters;
  /** \brief Default construction to use default solve criteria. */
  SolveCriteria()
    :solveMeasureType()
    ,requestedTol(unspecifiedTolerance())
    {}
  /** \brief Construct with a specified solve criteria. */
  SolveCriteria(
    SolveMeasureType _solveMeasureType, ScalarMag _requestedTol
    ,const Teuchos::RCP<Teuchos::ParameterList> &_extraParameters = Teuchos::null
    )
    :solveMeasureType(_solveMeasureType),requestedTol(_requestedTol),extraParameters(_extraParameters)
    {}
};


/** \brief Deprecated.
 *
 * \ingroup Thyra_Op_Solve_fundamental_interfaces_code_grp
 */
template <class Scalar>
struct THYRA_DEPRECATED BlockSolveCriteria {
  /** \brief Solve tolerance struct */
  SolveCriteria<Scalar> solveCriteria;
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
 * \ingroup Thyra_Op_Solve_fundamental_interfaces_code_grp
 */
class CatastrophicSolveFailure : public std::runtime_error
{public: CatastrophicSolveFailure(const std::string& what_arg) : std::runtime_error(what_arg) {}};


/** \brief Solution status
 *
 * \ingroup Thyra_Op_Solve_fundamental_interfaces_code_grp
 */
enum ESolveStatus {
  SOLVE_STATUS_CONVERGED        ///< The requested solution criteria has likely been achieved
  ,SOLVE_STATUS_UNCONVERGED     ///< The requested solution criteria has likely not been achieved
  ,SOLVE_STATUS_UNKNOWN         ///< The final solution status is unknown but he solve did not totally fail
};


/** \brief .
 *
 * \ingroup Thyra_Op_Solve_fundamental_interfaces_code_grp
 */
inline
const std::string toString(const ESolveStatus solveStatus)
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
 * \ingroup Thyra_Op_Solve_fundamental_interfaces_code_grp
 */
template <class Scalar>
struct SolveStatus {
  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
  /** \brief . */
  static ScalarMag unknownTolerance() { return ScalarMag(-1); }
  /** \brief The return status of the solve. */
  ESolveStatus solveStatus;
  /** \brief The maximum final tolerance actually achieved by the (block) linear solve.
   * A value of <tt>unknownTolerance()</tt> means that even an estimate of the
   * the final value of the tolerance is unknown. */
  ScalarMag achievedTol;
  /** \brief A simple one-line message (i.e. no newlines) returned from the solver */
  std::string message;
  /** \brief Any extra status parameters.
   * Note that the contents of this parameter list is totally undefined. */
  Teuchos::RCP<Teuchos::ParameterList> extraParameters;
  /** \brief . */
  SolveStatus()
    :solveStatus(SOLVE_STATUS_UNKNOWN), achievedTol(unknownTolerance())
    {}
  /** \brief Output the achieveTol field.
   */
  static std::string achievedTolToString( const ScalarMag &achievedTol )
    {
      if(achievedTol==unknownTolerance()) return "unknownTolerance()";
      std::ostringstream oss; oss << achievedTol; return oss.str();
    }
};


/** \brief Print the solve status to a stream.
 *
 * \relates SolveStatus
 */
template <class Scalar>
std::ostream& operator<<( std::ostream& out_arg, const SolveStatus<Scalar> &solveStatus )
{
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::getFancyOStream(Teuchos::rcp(&out_arg,false));
  Teuchos::OSTab tab(out);
  *out
    << "solveStatus = " << toString(solveStatus.solveStatus) << std::endl
    << "achievedTol = " << SolveStatus<Scalar>::achievedTolToString(solveStatus.achievedTol) << std::endl;
  *out << "message:";
  if (solveStatus.message.length()) {
    Teuchos::OSTab tab2(out);
    *out << "\n" << solveStatus.message << "\n";
  }
  *out << "extraParameters:";
  if(solveStatus.extraParameters.get()) {
    *out << "\n";
    Teuchos::OSTab tab3(out);
    solveStatus.extraParameters->print(*out, 1000, true);
  }
  else {
    *out << " NONE\n";
  }
  return out_arg;
}


/** \brief Enum that specifies how a <tt>LinearOpWithSolveBase</tt> object
 * will be used for solves after it is constructed.
 *
 * \ingroup Thyra_Op_Solve_fundamental_interfaces_code_grp
 */
enum ESupportSolveUse {
  SUPPORT_SOLVE_UNSPECIFIED  ///< How the output LOWSB object will be useded for solves in unspecified
  ,SUPPORT_SOLVE_FORWARD_ONLY  ///< The output LOWSB object will only be used for forward solves
  ,SUPPORT_SOLVE_TRANSPOSE_ONLY  ///< The output LOWSB object will only be used for transpose solves
  ,SUPPORT_SOLVE_FORWARD_AND_TRANSPOSE  ///< The output LOWSB object will used for forward and transpose solves
};


/** \brief Enum defining the status of a preconditioner object.
 *
 * \ingroup Thyra_Op_Solve_fundamental_interfaces_code_grp
 */
enum EPreconditionerInputType {
  PRECONDITIONER_INPUT_TYPE_AS_OPERATOR  ///< The input preconditioner should just be applied as an operator
  ,PRECONDITIONER_INPUT_TYPE_AS_MATRIX   ///< The input preconditioner should viewed as a matrix to be factored then backsolved as a preconditioner
};


/** \brief Initial overallSolveStatus before calling accumulateSolveStatus().
 *
 * \relates SolveStatus
 */
template <class Scalar>
void accumulateSolveStatusInit(
  const Ptr<SolveStatus<Scalar> > &overallSolveStatus
  )
{
  overallSolveStatus->solveStatus = SOLVE_STATUS_CONVERGED;
}


/** \brief Accumulate solve status objects for solving a block of RHSs is
 * smaller sub-blocks.
 *
 * \param overallSolveCriteria [in] The overall solve criteria for the overall
 * blocks.
 *
 * \param solveStatus [in] The solve status for a sub-block (or a single RHS)
 *
 * \param overallSolveStatus [in/out] The accumulated solve status for all the
 * sub-blocks of RHS.
 *
 * Before the first initialize with
 * <tt>accumulateSolveStatusInit(overallSolveStatus)</tt>.
 *
 * \relates SolveStatus
 */
template <class Scalar>
void accumulateSolveStatus(
  const SolveCriteria<Scalar>, // ToDo: Never used, need to take this out!
  const SolveStatus<Scalar> &solveStatus,
  const Ptr<SolveStatus<Scalar> > &overallSolveStatus
  )
{
  switch(solveStatus.solveStatus) {
    case SOLVE_STATUS_UNCONVERGED:
    {
      // First, if we see any unconverged solve status, then the entire block is
      // unconverged!
      overallSolveStatus->solveStatus = SOLVE_STATUS_UNCONVERGED;
      overallSolveStatus->message = solveStatus.message;
      overallSolveStatus->extraParameters = solveStatus.extraParameters;
      break;
    }
    case SOLVE_STATUS_UNKNOWN:
    {
      // Next, if any solve status is unknown, then if the overall solve
      // status says converged, then we have to mark it as unknown.  Note that
      // unknown could mean that the system is actually converged!
      switch(overallSolveStatus->solveStatus) {
        case SOLVE_STATUS_CONVERGED:
          overallSolveStatus->solveStatus = SOLVE_STATUS_UNKNOWN;
          break;
        case SOLVE_STATUS_UNCONVERGED:
        case SOLVE_STATUS_UNKNOWN:
          // If we get here then the overall solve status is either unknown
          // already or says unconverged and this will not change here!
          overallSolveStatus->message = solveStatus.message;
          overallSolveStatus->extraParameters = solveStatus.extraParameters;
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
  // Set a message if none is set
  if(overallSolveStatus->message == "")
    overallSolveStatus->message = solveStatus.message;
  // Set the extra parameters if none is set
  if(overallSolveStatus->extraParameters.get()==NULL)
    overallSolveStatus->extraParameters = solveStatus.extraParameters;
}


/** \brief Deprecated.
 *
 * \relates SolveStatus
 */
template <class Scalar>
THYRA_DEPRECATED
void accumulateSolveStatus(
  const SolveCriteria<Scalar>, // ToDo: Never used, need to take this out!
  const SolveStatus<Scalar> &solveStatus,
  SolveStatus<Scalar> *overallSolveStatus
  )
{
  accumulateSolveStatus(
    SolveCriteria<Scalar>(),
    solveStatus, Teuchos::ptr(overallSolveStatus)
    );
}


} // namespace Thyra


#endif // THYRA_SOLVE_SUPPORT_TYPES_HPP
