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

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include <fstream>
#include <string>

std::string filename;

namespace Teuchos {

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::UnitTestRepository::getCLP().setOption(
      "filename", &filename, "XML file to parse" );
  }

  TEUCHOS_UNIT_TEST( ParameterList, ExistingSublistIsOkay )
  {
    std::string xmlstring(
        "<ParameterList>                   \n"
        "  <ParameterList name=\"SubList\">\n"
        "  </ParameterList>                \n"
        "</ParameterList>                  \n");
    RCP<ParameterList> plist = getParametersFromXmlString(xmlstring);
    updateParametersFromXmlString( xmlstring, plist() );
  }

  TEUCHOS_UNIT_TEST( ParameterList, XMLDuplicatedSublists )
  {
    ParameterList pl;
    TEST_THROW( updateParametersFromXmlFile(filename, inOutArg(pl) ), DuplicateParameterSublist );
    TEST_THROW( getParametersFromXmlFile(filename), DuplicateParameterSublist );
    TEST_THROW( getParametersFromXmlFile(filename,null), DuplicateParameterSublist );
    //
    std::ifstream fin(filename.c_str());
    std::string xmlstring( (std::istreambuf_iterator<char>(fin)), 
                            std::istreambuf_iterator<char>()      );
    TEST_THROW( updateParametersFromXmlString(xmlstring,inOutArg(pl) ), DuplicateParameterSublist );
    TEST_THROW( getParametersFromXmlString(xmlstring), DuplicateParameterSublist );
    TEST_THROW( getParametersFromXmlString(xmlstring,null), DuplicateParameterSublist );
  }

  TEUCHOS_UNIT_TEST( XMLParameterListReader, XMLDuplicatedSublistsThrowsError )
  {
    FileInputSource xmlFile(filename);
    XMLObject xmlParams = xmlFile.getObject();
    XMLParameterListReader xmlPLReader;
    TEST_EQUALITY_CONST( xmlPLReader.getAllowsDuplicateSublists(), true );
    out << "Changing policy to disallow duplicate sublists" << std::endl;
    xmlPLReader.setAllowsDuplicateSublists( false );
    TEST_EQUALITY_CONST( xmlPLReader.getAllowsDuplicateSublists(), false );
    TEST_THROW( xmlPLReader.toParameterList(xmlParams), DuplicateParameterSublist );
  }

  TEUCHOS_UNIT_TEST( XMLParameterListReader, XMLDuplicatedSublistsBackwardsCompatible )
  {
    FileInputSource xmlFile(filename);
    XMLObject xmlParams = xmlFile.getObject();
    XMLParameterListReader xmlPLReader;
    TEST_EQUALITY_CONST( xmlPLReader.getAllowsDuplicateSublists(), true );
    TEST_NOTHROW( xmlPLReader.toParameterList(xmlParams) );
  }

} // namespace Teuchos
