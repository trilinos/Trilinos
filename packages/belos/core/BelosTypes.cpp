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

#include <BelosTypes.hpp>

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

  const char* TEUCHOS_DEPRECATED
  toString (const StatusType status)
  {
    return convertStatusTypeToRawString (status);
  }

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

} // end Belos namespace

