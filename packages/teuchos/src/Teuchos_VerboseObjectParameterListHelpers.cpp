// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
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
  static RCP<const ParameterList> validParams;
  if (is_null(validParams)) {
    RCP<ParameterList>
      pl = rcp(new ParameterList(VerboseObject_name));
    pl->set(
      VerbosityLevel_name, VerbosityLevel_default,
      "The verbosity level to use to override whatever is set in code.\n"
      "The value of \"default\" will allow the level set in code to be used.",
      VerbosityLevel_validator = verbosityLevelParameterEntryValidator(
        VerbosityLevel_name
        )
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
  TEST_FOR_EXCEPT(0==paramList);
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
  TEST_FOR_EXCEPT(0==paramList);
  TEST_FOR_EXCEPT(0==oStream);
  TEST_FOR_EXCEPT(0==verbLevel);
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
    TEST_FOR_EXCEPTION_PURE_MSG(
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
