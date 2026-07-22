// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \brief test that verifies zoltan2 parameters can be written to
 *       and read from xml files.
 *
 * validators used by zoltan2:
 *   Teuchos::AnyNumberParameterEntryValidator
 *   Teuchos::EnhancedNumberValidator
 *   Teuchos::StringValidator
 *   Teuchos::FileNameValidator
 *   Teuchos::StringToIntegralParameterEntryValidator
 *   Zoltan2::IntegerRangeListValidator
 */

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLObject.hpp>
#include <Teuchos_XMLParameterListWriter.hpp>
#include <Teuchos_ParameterXMLFileReader.hpp>
#include <Teuchos_ValidatorXMLConverterDB.hpp>

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Zoltan2_Parameters.hpp>
#include <Zoltan2_IntegerRangeList.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::tuple;

using std::string;

int main()
{
  // Create a parameter list with each validator type that we use.

  Teuchos::ParameterList pl("pl");

  {
    string parameterName("speed_versus_quality");
    RCP<const Teuchos::StringValidator> strValidatorP =
      rcp(new Teuchos::StringValidator(
        tuple<string>("speed", "balance", "quality")));
    std::ostringstream docString;
    strValidatorP->printDoc(
      "When algorithm choices exist, opt for speed or solution quality?\n"
      "(Default is a balance of speed and quality)\n",
       docString);
    pl.set<string>(parameterName, "balance", docString.str(), strValidatorP);
  }

  {
    string parameterName("debug_output_file");
    RCP<const Teuchos::FileNameValidator > fnameValidatorP =
      fnameValidatorP = rcp(new Teuchos::FileNameValidator(false));
    std::ostringstream docString;
    fnameValidatorP->printDoc(
      "name of file to which debug/status messages should be written\n"
      "(process rank will be included in file name)\n",
       docString);
    pl.set<string>(parameterName, "/dev/null", docString.str(), 
      fnameValidatorP);
  }

  {
    string parameterName("random_seed");
    RCP<const Teuchos::AnyNumberParameterEntryValidator> anyNumValidatorP =
      rcp(new Teuchos::AnyNumberParameterEntryValidator);
    std::ostringstream docString;
    anyNumValidatorP->printDoc("random seed\n", docString);
    pl.set<string>(parameterName, "0.5", docString.str(), anyNumValidatorP);
  }

  {
    string parameterName("debug_level");
    RCP<const Teuchos::StringToIntegralParameterEntryValidator<int> > 
      str2intValidatorP = 
      rcp(new Teuchos::StringToIntegralParameterEntryValidator<int>(
       tuple<string>("no_status",
                    "basic_status",
                    "detailed_status",
                    "verbose_detailed_status"),

       tuple<string>(
        "library outputs no status information",
        "library outputs basic status information (default)",
        "library outputs detailed information",
        "library outputs very detailed information"),

       tuple<int>(0, 1, 2, 3),

       parameterName));

    string info("the amount of status/warning/debugging info printed\n");
    info.append("(If the compile flag Z2_OMIT_ALL_STATUS_MESSAGES was set,\n");
    info.append("then message output code is not executed at runtime.)\n");

    std::ostringstream docString;
    str2intValidatorP->printDoc(info, docString);
    pl.set<string>(parameterName, "basic_status", docString.str(), 
      str2intValidatorP);
  }

  {
    string parameterName("debug_procs");
    typedef Zoltan2::IntegerRangeListValidator<int> irl_t;
    RCP<const irl_t> intRangeValidatorP = rcp(new irl_t);
    std::ostringstream docString;
    intRangeValidatorP->printDoc(
      "list of ranks that output debugging/status messages (default \"0\")\n",
       docString);
    pl.set<string>(parameterName, "0", docString.str(), intRangeValidatorP);

    // An XML converter for irl_t only needs to be added once to the converter db.

    typedef Zoltan2::IntegerRangeListValidatorXMLConverter<int> irlConverter_t;
    RCP<irlConverter_t > converter = rcp(new irlConverter_t);
    Teuchos::ValidatorXMLConverterDB::addConverter(
          intRangeValidatorP,    // can be a dummy of this type
          converter);
  }

  // Write out to XML
  Teuchos::XMLParameterListWriter plw;
  Teuchos::XMLObject obj = plw.toXML(pl);

  std::cout << "Parameter list: " << std::endl;
  std::cout << obj << std::endl;

  std::ofstream of;
  of.open("params.xml");
  of << obj << std::endl;
  of.close();

  // Read parameter list in from XML file.

  Teuchos::ParameterXMLFileReader rdr("params.xml");
  Teuchos::ParameterList newpl = rdr.getParameters();
  Teuchos::XMLObject objnew = plw.toXML(newpl);

  std::cout << "After reading in from XML file: " << std::endl;
  std::cout << objnew << std::endl;
}

