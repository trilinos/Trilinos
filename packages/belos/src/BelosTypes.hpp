//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef BELOS_TYPES_HPP
#define BELOS_TYPES_HPP

/*!
  \file BelosTypes.hpp
  \brief Collection of types and exceptions used within the Belos solvers.
*/

#include "BelosConfigDefs.hpp"
#include "Teuchos_Assert.hpp"

namespace Belos {

  //! @name Belos Exceptions
  //@{

  /// \class BelosError
  /// \brief Parent class to all Belos exceptions.
  class BelosError : public std::logic_error {
  public:
    BelosError (const std::string& what_arg) : std::logic_error(what_arg) {}
  };

  //@}


  /// \enum ETrans
  /// \brief Whether to apply the (conjugate) transpose of an operator.
  ///
  /// The OperatorTraits and Operator abstract interfaces to
  /// multivector-multivector operators have an Apply method, which
  /// implements applying the operator to a multivector.  This method
  /// optionally allows callers to apply the transpose or conjugate
  /// transpose (if applicable) of the operator, by supplying an
  /// ETrans enum value.
  ///
  /// \note This enum is useful to you only if you are specializing
  ///   OperatorTraits, implementing Operator, or implementing a Belos
  ///   solver.
  enum ETrans     {     NOTRANS = 0,  /*!< The operator should not be transposed during this application. */
                        TRANS = 1,    /*!< Apply the transpose of the operator. */
                        CONJTRANS = 2 /*!< Apply the conjugate transpose of the operator. */
  };

  /// \enum NormType
  /// \brief The type of vector norm to compute.
  ///
  /// The MultiVecTraits and MultiVec abstract interfaces to
  /// multivector operations have a method (MvNorm) for computing the
  /// norm of one or more vectors.  This method uses an optional
  /// NormType enum value to decide which kind of norm to compute.
  ///
  /// \note This enum is useful to you only if you are specializing
  ///   MultiVecTraits, implementing MultiVec, or implementing a Belos
  ///   solver.
  enum NormType {   OneNorm,       /*!< Compute the one-norm \f$\sum_{i=1}^{n}(|x_i w_i|)\f$ for each vector. */
                    TwoNorm,       /*!< Compute the two-norm \f$\sqrt(\sum_{i=1}^{n}((x_i w_i)^2))\f$ for each vector. */
                    InfNorm        /*!< Compute the infinity-norm \f$\max_{i=1}^{n}\{|x_i w_i|\}\f$ for each vector. */
  };

  /// \enum ScaleType
  /// \brief The type of scaling to use on the residual norm value.
  ///
  /// For residual norm convergence tests, the residual norm is scaled
  /// (divided) by some value.  This enum specifies which value to use
  /// -- for example, the norm of the right-hand side, the initial
  /// residual vector norm (preconditioned or not), a user-specified
  /// value, or no scaling at all.
  ///
  /// A user-specified value is useful when invoking an iterative
  /// solver as the inner loop of a nonlinear solver or iterative
  /// optimization method.  These may require special scaling factors
  /// (for example, that include the dimension of the linear system
  /// \f$N\f$).
  ///
  /// \note This enum is part of the public interface of many Belos
  ///   solvers, via the input list of parameters.
  enum ScaleType {NormOfRHS,     /*!< Use the norm of the right-hand-side. */
                  NormOfInitRes, /*!< Use the initial residual vector. */
                  NormOfPrecInitRes, /*!< Use the preconditioned initial residual vector. */
                  None,          /*!< Use unscaled residual. */
                  UserProvided,  /*!< User provides an explicit value by which the residual norm will be divided. */
                  NormOfFullInitRes, /*!< Use the full initial residual vector. Only interesting for partial subnorms. For singular fields this is equivalent with NormOfInitRes. */
                  NormOfFullPrecInitRes, /*!< Use the full preconditioned initial residual vector. Only interesting for partial subnorms. For singular fields this is equivalent with NormOfPrecInitRes. */
                  NormOfFullScaledInitRes, /*!< Use the full initial residual vector scaled by the inverse of the length of the subvector compared to the overall residual vector length. Only interesting for partial subnorms. For singular fields this is equivalent with NormOfInitRes. */
                  NormOfFullScaledPrecInitRes /*!< Use the full initial residual vector scaled by the inverse of the length of the subvector compared to the overall residual vector length. Only interesting for partial subnorms. For singular fields this is equivalent with NormOfInitRes. */
  };

  /// \enum OutputType
  /// \brief Style of output used to display status test information.
  ///
  /// Belos solvers optionally display intermediate status output
  /// during iterations.  "General" output may be quite verbose; each
  /// of the stopping criteria gets a chance to print something.
  /// "Brief" was intended to imitate AztecOO's status output.
  ///
  /// \note This enum is part of the public interface of many Belos
  ///   solvers, via the input list of parameters.
  enum OutputType {General,     /*!< Extensive output of status test information. */
                   Brief,       /*!< Simple output of current residual information. */
                   User         /*!< Simple output of user-specified status tests tagged by user */
  };

  /// \enum ReturnType
  /// \brief Whether the Belos solve converged for all linear systems.
  ///
  /// When the solve() method of any Belos::SolverManager is called,
  /// it returns a value of this type.  This indicates whether the
  /// solver manager sucessfully solved the linear system(s).
  ///
  /// \note This enum is part of the public interface of all Belos
  ///   solvers.
  enum ReturnType {
    Converged,  /*!< Convergence was reached for all linear systems. */
    Unconverged /*!< Convergence was not reached for some or all linear systems. */
  };

  //! Convert the given \c ReturnType enum value to its corresponding string.
  std::string
  convertReturnTypeToString (const ReturnType result);

  /// \enum StatusType
  /// \brief Whether the \c StatusTest wants iteration to stop.
  ///
  /// \note Only Belos developers and advanced users need this.
  ///
  /// Belos solvers use subclasses of \c StatusTest as iteration
  /// stopping criteria.  StatusTest's getStatus() and checkStatus()
  /// methods return a StatusType value, indicating whether iteration
  /// should stop.  "Passed" means that according to the particular
  /// test, iteration should stop.  "Failed" means that it should not
  /// yet stop.
  ///
  /// Note that "Passed" does not necessarily mean "Converged."  This
  /// can be a confusing point for new Belos developers.  For most
  /// Belos solvers, the StatusTest instance that controls iteration
  /// is really a combination of two or more StatusTest instances.
  /// These include a test for reaching the maximum iteration count,
  /// and one or more convergence tests.  When the maximum iteration
  /// count test Passes, it means that the maximum iteration count has
  /// been reached.  In this case, the outer test also Passes,
  /// indicating that iteration should stop.  However, this does not
  /// necessarily mean that the solver has converged.
  ///
  /// This enum is useful to you only if you are implementing a
  /// subclass of StatusTest, Iteration, or SolverManager.
  enum StatusType {     Passed = 0x1,      /*!< Some event occured, the iteration needs to stop. */
                        Failed = 0x2,      /*!< No event has occurred requiring the iteration to stop. */
                        Undefined = 0x4    /*!< Status test has not been checked yet. */
  };

  /// \enum ResetType
  /// \brief How to reset the solver.
  ///
  /// \c SolverManager's reset() method accepts an enum value of this
  /// type.  Some Belos solvers have special capabilities, such as
  /// creating and using a recycling subspace.  The RecycleSubspace
  /// enum value allows you to access these capabilities.  For more
  /// details, see the documentation of the specific \c SolverManager
  /// subclass.
  ///
  /// \note This enum is part of the public interface of all Belos
  ///   solvers.
  enum ResetType  {     Problem = 0x1,           /*!< Reset the linear problem inside the solver. */
                        RecycleSubspace = 0x2    /*!< Destroy any existing subspace inside the solver. */
  };

  /// \brief The string name corresponding to the given StatusType enum value.
  ///
  /// This method is DEPRECATED because the generic-sounding name of
  /// this function makes it easy to pass in the wrong enum type.  Use
  /// \c convertStatusTypeToString() instead.
  BELOS_DEPRECATED const char*
  toString (const StatusType status);

  //! The string name corresponding to the given StatusType enum value.
  std::string
  convertStatusTypeToString (const StatusType status);

  //! The StatusType enum value corresponding to the given string name.
  StatusType
  convertStringToStatusType (const std::string& status);

  //! Convert the given string to its \c ScaleType enum value.
  ScaleType
  convertStringToScaleType (const std::string& scaleType);

  //! Convert the given \c ScaleType enum value to its corresponding string.
  std::string
  convertScaleTypeToString (const ScaleType scaleType);

  /// \enum ConjType
  /// \brief Whether or not to conjugate the transpose for block inner products.
  ///
  /// The MultiVecTraits and MultiVec abstract interfaces to linear
  /// algebra operations can compute block inner products using
  /// MvTransMv and MvDot.  This enum specifies whether to use the
  /// transpose or conjugate transpose when computing block inner
  /// products.
  enum ConjType {
    NO_CONJ,      /*!< Not conjugated */
    CONJ          /*!< Conjugated */
  };

  /// \enum MsgType
  /// \brief Available message types recognized by the linear solvers.
  ///
  /// This enum allows users to control the kinds of output that
  /// Belos' solvers generate.  MsgType is not an enum in the
  /// strictest sense, because the values are not exclusive.  Rather,
  /// it is a C-style bit set.  You can bitwise OR together different
  /// enum values to get a combination.
  ///
  /// \note This enum is part of the public interface of many Belos
  ///   solvers, via the input list of parameters.
  enum MsgType {
    Errors= 0,                  /*!< Errors [ always printed ] */
    Warnings = 0x1,             /*!< Internal warnings */
    IterationDetails = 0x2,     /*!< Approximate/exact residuals */
    OrthoDetails = 0x4,         /*!< Orthogonalization/orthonormalization details */
    FinalSummary = 0x8,         /*!< Final computational summary */
    TimingDetails = 0x10,       /*!< Timing details */
    StatusTestDetails = 0x20,   /*!< Status test details */
    Debug = 0x40                /*!< Debugging information */
  };

  /// \brief Show MsgType as a comma-delimited list of names.
  ///
  /// The \c MsgType enum is really a C-style bit set (where you
  /// bitwise OR together different names to get a combination of
  /// values).  This function returns a string representing the given
  /// MsgType (represented here as an int, mainly because
  /// Teuchos::ParameterList seems to prefer storing these kind of
  /// C-style bit sets as int rather than MsgType) as a
  /// comma-delimited, human-readable list of names.  This is useful
  /// for debugging.
  std::string
  convertMsgTypeToString (const MsgType msgType);

  /// \brief Default parameters common to most Belos solvers
  ///
  /// Both Belos solvers and users may override these defaults.  Real
  /// floating-point values are deliberately double, in order to avoid
  /// issues with constexpr construction of certain MagnitudeTypes.
  struct DefaultSolverParameters {
    /// \brief Default convergence tolerance
    ///
    /// This assumes that implicit conversion from double to
    /// Teuchos::ScalarTraits<Scalar>::magnitudeType always works, but
    /// Belos already assumed that.  See discussion starting here:
    ///
    /// https://github.com/trilinos/Trilinos/pull/2677#issuecomment-395453521
    static constexpr double convTol = 1.0e-8;

    //! Relative residual tolerance for matrix polynomial construction
    static constexpr double polyTol = 1.0e-12;

    //! DGKS orthogonalization constant
    static constexpr double orthoKappa = -1.0;

    //! User-defined residual scaling factor
    static constexpr double resScaleFactor = 1.0;

    //! "Implicit Tolerance Scale Factor"
    static constexpr double impTolScale = 10.0;
  };


} // end Belos namespace

#endif /* BELOS_TYPES_HPP */
