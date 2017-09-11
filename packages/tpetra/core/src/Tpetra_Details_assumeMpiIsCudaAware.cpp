// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

#include "TpetraCore_config.h"
#include "Tpetra_Details_assumeMpiIsCudaAware.hpp"
#include <string>
#include <cctype> // toupper
#include <cstdlib> // getenv

namespace Tpetra {
namespace Details {

bool
assumeMpiIsCudaAware (std::ostream* out)
{
  using std::endl;
#ifdef TPETRA_ASSUME_CUDA_AWARE_MPI
  const bool configureTimeDefault = true;
#else
  const bool configureTimeDefault = false;
#endif // TPETRA_ASSUME_CUDA_AWARE_MPI
  const char envVarName[] = "TPETRA_ASSUME_CUDA_AWARE_MPI";

  if (out != NULL) {
    *out << "Test whether users want Trilinos to assume that its MPI "
      "implementation is CUDA aware." << endl;
  }

  // Environment variable overrides default configure-time value.
  const char* envVarVal = std::getenv (envVarName);
  if (envVarVal == NULL) { // use configure-time default value
    if (out != NULL) {
      *out << "  - Did not find environment variable \"" << envVarName
           << "\"." << endl
           << "  - Return configure-time default: "
           << (configureTimeDefault ? "true" : "false") << "." << endl;
    }
    return configureTimeDefault;
  }
  else {
    if (out != NULL) {
      *out << "  - Found environment variable \"" << envVarName << "\"."
           << endl << "  - Its value is \"" << envVarVal << "\"." << endl;
    }
    const std::string varVal ([=] { // compare to Lisp' LET (yay const!)
        std::string varVal_nc (envVarVal);
        for (auto& c : varVal_nc) {
          c = toupper (c);
        }
        return varVal_nc;
      } ());
    // We include "" with the "true" values rather than the "false"
    // values, since just setting a Boolean environment variable
    // without setting its value suggests a true value.  Note that
    // this requires actually setting the value to "".  For example,
    // in bash, one must do this:
    //
    // export TPETRA_ASSUME_CUDA_AWARE_MPI=""
    //
    // not this:
    //
    // export TPETRA_ASSUME_CUDA_AWARE_MPI
    //
    // If you do the latter, getenv will return NULL.
    const char* falseVals[] = {"0", "NO", "OFF", "FALSE"};
    for (auto falseVal : falseVals) {
      if (varVal == falseVal) {
        if (out != NULL) {
          *out << "  - Report value of environment variable as \"false\"."
               << endl;
        }
        return false;
      }
    }
    const char* trueVals[] = {"", "1", "YES", "ON", "TRUE"};
    for (auto trueVal : trueVals) {
      if (varVal == trueVal) {
        if (out != NULL) {
          *out << "  - Report value of environment variable as \"true\"."
               << endl;
        }
        return true;
      }
    }

    if (out != NULL) {
      *out << "  - Found environment variable, but didn't know how to "
           << "interpret its value." << endl
           << "  - Return configure-time default: "
           << (configureTimeDefault ? "true" : "false") << endl;
    }
    return configureTimeDefault;
  }
}

} // namespace Details
} // namespace Tpetra
