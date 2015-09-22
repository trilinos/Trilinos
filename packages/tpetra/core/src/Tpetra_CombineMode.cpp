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

#include <Tpetra_CombineMode.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

namespace Tpetra {

  void
  setCombineModeParameter (Teuchos::ParameterList& plist,
                           const std::string& paramName)
  {
    typedef Tpetra::CombineMode enum_type;
    typedef Teuchos::StringToIntegralParameterEntryValidator<enum_type>
      validator_type;

    const std::string docString = "Tpetra::CombineMode: rule for combining "
      "entries that overlap across processes, when redistributing data via a "
      "Tpetra::Import or Tpetra::Export";
    const std::string defaultVal = "ADD";
    const bool caseSensitive = false;

    const Teuchos::Array<std::string>::size_type numParams = 5;
    Teuchos::Array<std::string> strs (numParams);
    Teuchos::Array<std::string> docs (numParams);
    Teuchos::Array<enum_type> vals (numParams);

    strs[0] = "ADD";
    strs[1] = "INSERT";
    strs[2] = "REPLACE";
    strs[3] = "ABSMAX";
    strs[4] = "ZERO";

    docs[0] = "Sum new values into existing values";
    docs[1] = "Insert new values that don't currently exist";
    docs[2] = "Replace existing values with new values";
    docs[3] = "Replace old value with maximum of magnitudes of old and new values";
    docs[4] = "Replace old values with zero";

    vals[0] = ADD;
    vals[1] = INSERT;
    vals[2] = REPLACE;
    vals[3] = ABSMAX;
    vals[4] = ZERO;

    plist.set (paramName, defaultVal, docString,
               Teuchos::rcp (new validator_type (strs (), docs (), vals (),
                                                 defaultVal, caseSensitive)));
  }

} // namespace Tpetra
