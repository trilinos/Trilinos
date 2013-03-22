// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER


#include "Teuchos_VerbosityLevelCommandLineProcessorHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_implicit_cast.hpp"


namespace {


using Teuchos::Array;
using Teuchos::tuple;


const Array<Teuchos::EVerbosityLevel>
verbosityLevelValues = tuple<Teuchos::EVerbosityLevel>(
  Teuchos::VERB_DEFAULT,
  Teuchos::VERB_NONE,
  Teuchos::VERB_LOW,
  Teuchos::VERB_MEDIUM,
  Teuchos::VERB_HIGH,
  Teuchos::VERB_EXTREME
  );


const Array<std::string>
verbosityLevelNamesStorage = tuple<std::string>(
  toString(Teuchos::VERB_DEFAULT),
  toString(Teuchos::VERB_NONE),
  toString(Teuchos::VERB_LOW),
  toString(Teuchos::VERB_MEDIUM),
  toString(Teuchos::VERB_HIGH),
  toString(Teuchos::VERB_EXTREME)
  );


Array<const char*> verbosityLevelNames;
// The above variable is Intialized below in order so that if exceptions are
// thrown then they will be caught in main()!


} // namespace


void Teuchos::setVerbosityLevelOption(
  const std::string &optionName,
  EVerbosityLevel *verbLevel,
  const std::string &docString,
  CommandLineProcessor *clp,
  const bool required
  )
{
  const int numVerbLevels = verbosityLevelValues.size();

  if ( !verbosityLevelNames.size() ) {
    verbosityLevelNames = tuple<const char*>(
      verbosityLevelNamesStorage[0].c_str(),
      verbosityLevelNamesStorage[1].c_str(),
      verbosityLevelNamesStorage[2].c_str(),
      verbosityLevelNamesStorage[3].c_str(),
      verbosityLevelNamesStorage[4].c_str(),
      verbosityLevelNamesStorage[5].c_str()
      );
  }
  
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT( optionName.length() );
  TEUCHOS_ASSERT( verbLevel );
  TEUCHOS_ASSERT( clp );
  TEUCHOS_ASSERT( implicit_cast<int>(verbosityLevelNamesStorage.size()) == numVerbLevels );
  TEUCHOS_ASSERT( implicit_cast<int>(verbosityLevelNames.size()) == numVerbLevels );
#endif
  clp->setOption(
    optionName.c_str(), verbLevel,
    numVerbLevels, &verbosityLevelValues[0], &verbosityLevelNames[0],
    docString.c_str(), required
    );
}
