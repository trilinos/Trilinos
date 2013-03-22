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

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"


namespace {


const std::string VerboseObject_name = "VerboseObject";

const std::string OutputFile_name = "Output File";
const std::string OutputFile_default = "none";

const std::string VerbosityLevel_name = "Verbosity Level";
const std::string  VerbosityLevel_default = "default";
Teuchos::RCP<
  Teuchos::StringToIntegralParameterEntryValidator<Teuchos::EVerbosityLevel>
  >
VerbosityLevel_validator;


} // namespace



Teuchos::RCP<const Teuchos::ParameterList>
Teuchos::getValidVerboseObjectSublist()
{
  using Teuchos::rcp_implicit_cast;
  static RCP<const ParameterList> validParams;
  if (is_null(validParams)) {
    RCP<ParameterList>
      pl = rcp(new ParameterList(VerboseObject_name));
    VerbosityLevel_validator = verbosityLevelParameterEntryValidator(VerbosityLevel_name);
    pl->set(
      VerbosityLevel_name, VerbosityLevel_default,
      "The verbosity level to use to override whatever is set in code.\n"
      "The value of \"default\" will allow the level set in code to be used.",
      rcp_implicit_cast<const ParameterEntryValidator>(VerbosityLevel_validator)
      );
    pl->set(
      OutputFile_name, OutputFile_default,
      "The file to send output to.  If the value \"none\" is used, then\n"
      "whatever is set in code will be used.  However, any other std::string value\n"
      "will be used to create an std::ofstream object to a file with the given name.\n"
      "Therefore, any valid file name is a valid std::string value for this parameter."
      );
    validParams = pl;
  }
  return validParams;
}


void Teuchos::setupVerboseObjectSublist( ParameterList* paramList )
{
  TEUCHOS_TEST_FOR_EXCEPT(0==paramList);
  paramList->sublist(VerboseObject_name).setParameters(
    *getValidVerboseObjectSublist()
    ).disableRecursiveValidation();
}


void Teuchos::readVerboseObjectSublist(
  ParameterList* paramList,
  RCP<FancyOStream> *oStream, EVerbosityLevel *verbLevel
  )
{
  // Validate input
  TEUCHOS_TEST_FOR_EXCEPT(0==paramList);
  TEUCHOS_TEST_FOR_EXCEPT(0==oStream);
  TEUCHOS_TEST_FOR_EXCEPT(0==verbLevel);
  ParameterList
    &voSublist = paramList->sublist(VerboseObject_name);
  voSublist.validateParameters(*getValidVerboseObjectSublist());
  const std::string
    outputFileStr = voSublist.get(OutputFile_name,OutputFile_default);
  *verbLevel = VerbosityLevel_validator->getIntegralValue(
    voSublist,VerbosityLevel_name,VerbosityLevel_default
    );
  if (outputFileStr==OutputFile_default) {
    *oStream = null;
  }
  else {
    RCP<std::ofstream>
      oFileStream = rcp(new std::ofstream(outputFileStr.c_str()));
    TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(
      oFileStream->eof(), Exceptions::InvalidParameterValue,
      "Error, the file \"" << outputFileStr << "\n given by the parameter\n"
      "\'" << OutputFile_name << "\' in the sublist\n"
      "\'" << voSublist.name() << "\' count not be opened for output!"
      );
    *oStream = fancyOStream(rcp_implicit_cast<std::ostream>(oFileStream));
  }
#ifdef TEUCHOS_DEBUG
  voSublist.validateParameters(*getValidVerboseObjectSublist());
#endif
}
