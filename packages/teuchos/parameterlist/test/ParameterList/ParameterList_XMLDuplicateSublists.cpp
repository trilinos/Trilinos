// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
