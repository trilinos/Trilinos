// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <BelosTypes.hpp>
#include <vector>

namespace Belos {

  namespace {
    const char*
    convertStatusTypeToRawString (const StatusType status)
    {
      if (status == Passed) {
        return "Passed";
      } else if (status == Failed) {
        return "Failed";
      } else if (status == Undefined) {
        return "Undefined";
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Belos::convertStatusTypeToRawString: Invalid StatusType enum value "
          << status << ".");
      }
    }
  } // namespace (anonymous)

  std::string
  convertStatusTypeToString (const StatusType status)
  {
    return convertStatusTypeToRawString (status);
  }

  StatusType
  convertStringToStatusType (const std::string& status)
  {
    if (status == "Passed") {
      return Passed;
    } else if (status == "Failed") {
      return Failed;
    } else if (status == "Undefined") {
      return Undefined;
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Belos::convertStringToStatusType: Invalid string \"" << status
        << "\".");
    }
  }

  NormType
  convertStringToNormType (const std::string& normType)
  {
    if (normType == "OneNorm") {
      return Belos::OneNorm;
    } else if (normType == "TwoNorm") {
      return Belos::TwoNorm;
    } else if (normType == "InfNorm") {
       return Belos::InfNorm;
    } else if (normType == "PreconditionerNorm") {
      return Belos::PreconditionerNorm;
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
        "Belos::convertStringToNormType(): Invalid norm type \""
        << normType << "\".");
    }
  }

  ScaleType
  convertStringToScaleType (const std::string& scaleType)
  {
    if (scaleType == "Norm of Initial Residual") {
      return Belos::NormOfInitRes;
    } else if (scaleType == "Norm of Preconditioned Initial Residual") {
      return Belos::NormOfPrecInitRes;
    } else if (scaleType == "Norm of RHS") {
       return Belos::NormOfRHS;
    } else if (scaleType == "None") {
      return Belos::None;
    } else if (scaleType == "User Provided") {
      return Belos::UserProvided;
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
        "Belos::convertStringToScaleType(): Invalid residual scaling type \""
        << scaleType << "\".");
    }
  }

  std::string
  convertScaleTypeToString (const ScaleType scaleType)
  {
    if (scaleType == Belos::NormOfInitRes) {
      return "Norm of Initial Residual";
    } else if (scaleType == Belos::NormOfPrecInitRes) {
      return "Norm of Preconditioned Initial Residual";
    } else if (scaleType == Belos::NormOfRHS) {
      return "Norm of RHS";
    } else if (scaleType == Belos::None) {
      return "None";
    } else if (scaleType == Belos::UserProvided) {
      return "User Provided";
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
        "Belos::convertScaleTypeToString(): Invalid residual scaling type "
        "value " << scaleType << ".");
    }
  }

  std::string
  convertMsgTypeToString (const MsgType msgType)
  {
    typedef std::vector<int>::size_type size_type;

    // Wouldn't it be nice if C++ enums had introspection and could
    // be enumerated?
    const size_type numValidTypes = 8;
    const int validTypes[] = {
      Belos::Errors,
      Belos::Warnings,
      Belos::IterationDetails,
      Belos::OrthoDetails,
      Belos::FinalSummary,
      Belos::TimingDetails,
      Belos::StatusTestDetails,
      Belos::Debug
    };
    const char* typeNames[] = {
      "Errors",
      "Warnings",
      "IterationDetails",
      "OrthoDetails",
      "FinalSummary",
      "TimingDetails",
      "StatusTestDetails",
      "Debug"
    };

    // We first generate a list, and only then build a single string.
    // This helps us decide where to put the commas.  The list just
    // uses the indices of the valid names, rather than the valid
    // names themselves, in order to save space and time.  We use
    // size_type for the indices to avoid signed/unsigned comparisons.
    std::vector<size_type> theList;
    for (size_type nameIndex = 0; nameIndex < numValidTypes; ++nameIndex) {
      if (msgType & validTypes[nameIndex]) {
        theList.push_back (nameIndex);
      }
    }
    std::ostringstream os;
    for (size_type k = 0; k < theList.size(); ++k) {
      const size_type nameIndex = theList[k];
      os << typeNames[nameIndex];
      if (nameIndex < theList.size() - 1) {
        os << ",";
      }
    }
    return os.str();
  }

  std::string
  convertReturnTypeToString (const ReturnType result)
  {
    if (result == Belos::Converged) {
      return "Converged";
    } else if (result == Belos::Unconverged) {
      return "Unconverged";
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Belos::convertReturnTypeToString: Invalid ReturnType enum value "
        << result << ".");
    }
  }

  // Initialize DefaultSolverParameters.  Has to be done this way because
  // the "static consexpr double blah = ...;" pattern can create ODR-used
  // linking errors (usually in debug builds).
  const double DefaultSolverParameters::convTol = 1.0e-8;
  const double DefaultSolverParameters::polyTol = 1.0e-12;
  const double DefaultSolverParameters::orthoKappa = -1.0;
  const double DefaultSolverParameters::resScaleFactor = 1.0;
  const double DefaultSolverParameters::impTolScale = 10.0;

} // end Belos namespace

