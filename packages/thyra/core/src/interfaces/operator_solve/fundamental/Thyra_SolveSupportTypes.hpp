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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_SOLVE_SUPPORT_TYPES_HPP
#define THYRA_SOLVE_SUPPORT_TYPES_HPP

#include "Thyra_OperatorVectorTypes.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_Describable.hpp"


namespace Thyra {


/** \brief Type of solve measure norm.
 *
 * For reference we refer to solving a single linear system <tt>A*x=b</tt>.
 *
 * \ingroup Thyra_Op_Solve_fundamental_interfaces_code_grp
 */
enum ESolveMeasureNormType {
  /// No solve measure (i.e. same as 1.0)
  SOLVE_MEASURE_ONE,
  /// Norm of the current residual vector (i.e. <tt>||A*x-b||</tt>)
  SOLVE_MEASURE_NORM_RESIDUAL,
  /// Norm of the current solution vector (i.e. <tt>||x||</tt>)
  SOLVE_MEASURE_NORM_SOLUTION,
  /// Norm of the initial residual vector given a non-zero guess (i.e. <tt>||A*xo-b||</tt>)
  SOLVE_MEASURE_NORM_INIT_RESIDUAL,
  /// Norm of the right-hand side (i.e. <tt>||b||</tt>)
  SOLVE_MEASURE_NORM_RHS
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
      TEUCHOS_TEST_FOR_EXCEPT(true);
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
  /** \brief Return if this is a default solve measure (default
   * constructed).
   */
  bool useDefault() const
    { return ( numerator==SOLVE_MEASURE_ONE && denominator==SOLVE_MEASURE_ONE ); }
  /** \brief Return if (numerator,denominataor) matches this. */
  bool operator()(ESolveMeasureNormType numerator_in,
    ESolveMeasureNormType denominator_in
    ) const
    { return ( numerator==numerator_in && denominator==denominator_in ); }
  /** \breif Return if single measure matches numerator or denominator. */
  bool contains(ESolveMeasureNormType measure) const
    { return ( numerator==measure || denominator==measure ); }
};


/** \brief Output operator.
 *
 * \relates SolveMeasureType
 */
inline
std::ostream& operator<<(std::ostream &out, const SolveMeasureType &solveMeasureType)
{
  out << "("<<toString(solveMeasureType.numerator)
      << "/"<<toString(solveMeasureType.denominator)<<")";
  return out;
}


/** \brief A general reduction functional to be used in specialized solve
 * convergence criteria.
 */
template<class Scalar>
class ReductionFunctional : public Teuchos::Describable {
public:

  /** \name Public non-virtual functions. */
  //@{

  /** \brief Compute the reduction over a vector.
   *
   * \param v [in] The vector being reduced into a Scalar.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>this->isCompatible(v) == true</tt>
   * </ul>
   */
  typename ScalarTraits<Scalar>::magnitudeType
  reduce( const VectorBase<Scalar> &v ) const
    {
#ifdef THYRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(!isCompatible(v), Exceptions::IncompatibleVectorSpaces,
        "Error, the vector v="<<v.description()<<" is not compatiable with"
        " *this="<<this->description()<<"!");
#endif
      return reduceImpl(v);
    }      

  /** \brief Returns <tt>true</tt> if <tt>v</tt> is compatible with
   * <tt>*this</tt>.
   */
  bool isCompatible( const VectorBase<Scalar> &v ) const
    { return isCompatibleImpl(v); }

  //@}

protected:

  /** \name Protected virtual functions. */
  //@{

  /** \brief . */
  virtual typename ScalarTraits<Scalar>::magnitudeType
  reduceImpl( const VectorBase<Scalar> &v ) const = 0;

  /** \brief . */
  virtual bool isCompatibleImpl( const VectorBase<Scalar> &v ) const = 0;

  //@}

};


/** \brief Simple struct that defines the requested solution criteria for a solve.
 *
 * A solve criteria defines the solution to a linear (or nonlinear) system of
 * equations in terms of purely mathematical entities.  The form of the linear
 * system is:

 \verbatim

  A * x = b

  r = b - A * x

 \endverbatim

 * with <tt>x0</tt> defining the initial guess for the solution and:

 \verbatim

  r0 = b - A * x0

 \endverbatim

 * The mathematical representation of the solve criteria takes the form:

 \verbatim

    gN(vN) / gD(vD) <= requestedTol

 \endverbatim

 * where <tt>gN(vN)</tt> and <tt>gD(vD)</tt> are defined as <tt>g(v)=</tt>
 * <table>
 * <tr><td><tt>||r||</tt></td>
 *   <td>: if <tt>solveMeasureValue==SOLVE_MEASURE_NORM_RESIDUAL && reductionFunc==null</tt></td></tr>
 * <tr><td><tt>reductionFunc.reduce(r)</tt></td>
 *   <td>: if <tt>solveMeasureValue==SOLVE_MEASURE_NORM_RESIDUAL && reductionFunc!=null</tt></td></tr>
 * <tr><td><tt>||x||</tt></td>
 *   <td>: if <tt>solveMeasureValue==SOLVE_MEASURE_NORM_SOLUTION && reductionFunc==null</tt></td></tr>
 * <tr><td><tt>reductionFunc.reduce(x)</tt></td>
 *   <td>: if <tt>solveMeasureValue==SOLVE_MEASURE_NORM_SOLUTION && reductionFunc!=null</tt></td></tr>
 * <tr><td><tt>||r0||</tt></td>
 *   <td>: if <tt>solveMeasureValue==SOLVE_MEASURE_NORM_INIT_RESIDUAL && reductionFunc==null</tt></td></tr>
 * <tr><td><tt>reductionFunc.reduce(r0)</tt></td>
 *   <td>: if <tt>solveMeasureValue==SOLVE_MEASURE_NORM_INIT_RESIDUAL && reductionFunc!=null</tt></td></tr>
 * <tr><td><tt>||b||</tt></td>
 *   <td>: if <tt>solveMeasureValue==SOLVE_MEASURE_NORM_RHS && reductionFunc==null</tt></td></tr>
 * <tr><td><tt>reductionFunc.reduce(b)</tt></td>
 *   <td>: if <tt>solveMeasureValue==SOLVE_MEASURE_NORM_RHS && reductionFunc!=null</tt></td></tr>
 * <tr><td><tt>1</tt></td>
 *   <td>: if <tt>solveMeasureValue==SOLVE_MEASURE_ONE</tt></td></tr>
 * </table>
 *
 * where <tt>solveMeasureValue = solveMeasure.numerator</tt> and
 * <tt>reductionFunc = numeratorReductionFunc</tt> for <tt>gN(vN)</tt> while
 * <tt>solveMeasureValue = solveMeasure.denominator</tt> and <tt>reductionFunc
 * = denominatorReductionFunc</tt> for <tt>gD(vD)</tt>.
 *
 * For example, for
 * <tt>solveMeasure.numerator==SOLVE_MEASURE_NORM_RESIDUAL</tt> and
 * <tt>solveMeasure.denominator==SOLVE_MEASURE_ONE</tt> we have the solve
 * convergence criteria:
 *

 \verbatim

    ||r|| / 1 <= requestedTol

 \endverbatim

 * For <tt>solveMeasure.numerator==SOLVE_MEASURE_NORM_RESIDUAL</tt> and
 * <tt>solveMeasure.denominator==SOLVE_MEASURE_NORM_INIT_RESIDUAL</tt> we have
 * the solve convergence criteria:
 *

 \verbatim

    ||r|| / ||r0|| <= requestedTol

 \endverbatim

 * The objects <tt>numeratorReductionFunc</tt> and
 * <tt>denominatorReductionFunc</tt> basically override the use of the natural
 * norm <tt>||.||</tt> for the given vector.  This is needed to implement some
 * unusual convergence criteria needed for certain types of nonlinear ANAs
 * (such as the optimization solvers in the Aristos package).
 *
 * There are several reasons for the structure of the solve convergence
 * criteria shown above.  First, we want to give the solver implementation as
 * much information as we can as to the nature of the solve convergence
 * criteria.  That way, the solver implementation can compute the different
 * quantities more efficiently in many cases..  For example, with GMRES no
 * direct estimate of the residual vector <tt>r</tt> is cheaply available but
 * a cheap estimate of the natural norm <tt>||r||</tt> is readily available.
 * Also, while the vectors <tt>r0</tt> and <tt>b</tt> could be computed by the
 * client before the solve, it is potentially more efficient to let the solver
 * do it since it may compute theses quantities as a natural byproduct of the
 * solve process.
 *
 * \ingroup Thyra_Op_Solve_fundamental_interfaces_code_grp
 */
template <class Scalar>
struct SolveCriteria {
  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
  /** \brief . */
  static ScalarMag unspecifiedTolerance() { return ScalarMag(-1.0); }
  /** \brief The type of solve tolerance requested as given in
   * <tt>this->requestedTol</tt>. */
  SolveMeasureType solveMeasureType;
  /** \brief The requested solve tolerance (what the client would like to see).
   * Only significant if <tt>!this->solveMeasureType.useDefault()</tt> */
  ScalarMag requestedTol;
  /** \brief Any extra control parameters (e.g. max iterations).
   *
   * Note that the contents of this parameter list is totally undefined and
   * any client that uses this does so at their own peril!
   */
  RCP<ParameterList> extraParameters;
  /** \brief Reduction function to be used in place of the natural norm of the
   * numerator. */
  RCP<const ReductionFunctional<Scalar> > numeratorReductionFunc;
  /** \brief Reduction function to be used in place of the natural norm of the
   * numerator. */
  RCP<const ReductionFunctional<Scalar> > denominatorReductionFunc;
  /** \brief Default construction to use default solve criteria. */
  SolveCriteria()
    : requestedTol(unspecifiedTolerance())
    {}
  /** \brief Construct with a specified solve criteria. */
  SolveCriteria(
    SolveMeasureType solveMeasureType_in,
    ScalarMag requestedTol_in,
    const RCP<ParameterList> &extraParameters_in = Teuchos::null,
    const RCP<ReductionFunctional<Scalar> > &numeratorReductionFunc_in = Teuchos::null,
    const RCP<ReductionFunctional<Scalar> > &denominatorReductionFunc_in = Teuchos::null
    )
    : solveMeasureType(solveMeasureType_in),
      requestedTol(requestedTol_in), 
      extraParameters(extraParameters_in),
      numeratorReductionFunc(numeratorReductionFunc_in),
      denominatorReductionFunc(denominatorReductionFunc_in)
    {}
};


/** \brief Output operator.
 *
 * \relates SolveCriteria
 */
template<class Scalar>
std::ostream& operator<<(std::ostream &out, const SolveCriteria<Scalar> &solveCriteria)
{
  out << typeName(solveCriteria) << "{";
  out << "solveMeasureType="<<solveCriteria.solveMeasureType;
  out << ", requestedTol="<<solveCriteria.requestedTol;
  if (nonnull(solveCriteria.extraParameters)) {
    out << ", extraParameters="<<solveCriteria.extraParameters;
  }
  if (nonnull(solveCriteria.numeratorReductionFunc)) {
    out << ", numeratorReductionFunc="<<solveCriteria.numeratorReductionFunc->description();
  }
  if (nonnull(solveCriteria.denominatorReductionFunc)) {
    out << ", denominatorReductionFunc="<<solveCriteria.denominatorReductionFunc->description();
  }
  out << "}";
  return out;
}

#ifndef THYRA_HIDE_DEPRECATED_CODE
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
#endif // THYRA_HIDE_DEPRECATED_CODE

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
    default: TEUCHOS_TEST_FOR_EXCEPT(true);
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
  RCP<ParameterList> extraParameters;
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
  RCP<Teuchos::FancyOStream>
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
          TEUCHOS_TEST_FOR_EXCEPT(true); // Corrupted enum?
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
      TEUCHOS_TEST_FOR_EXCEPT(true); // Corrupted enum?
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

#ifndef THYRA_HIDE_DEPRECATED_CODE
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
#endif // THYRA_HIDE_DEPRECATED_CODE

} // namespace Thyra


#endif // THYRA_SOLVE_SUPPORT_TYPES_HPP
