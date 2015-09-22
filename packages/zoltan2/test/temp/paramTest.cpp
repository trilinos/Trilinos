// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
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
using std::ostringstream;
using std::endl;
using std::cout;

int main()
{
  // Create a parameter list with each validator type that we use.

  Teuchos::ParameterList pl("pl");

  {
    string parameterName("speed_versus_quality");
    RCP<const Teuchos::StringValidator> strValidatorP =
      rcp(new Teuchos::StringValidator(
        tuple<string>("speed", "balance", "quality")));
    ostringstream docString;
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
    ostringstream docString;
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
    ostringstream docString;
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

    ostringstream docString;
    str2intValidatorP->printDoc(info, docString);
    pl.set<string>(parameterName, "basic_status", docString.str(), 
      str2intValidatorP);
  }

  {
    string parameterName("bisection_num_test_cuts");
    RCP<const Teuchos::EnhancedNumberValidator<int> > intValidatorP =
      rcp(new Teuchos::EnhancedNumberValidator<int>(1,250,1));
    ostringstream docString;
    intValidatorP->printDoc(
     "Experimental: number of test cuts to do simultaneously (default is 1)\n",
        docString);
    pl.set<int>(parameterName, 1, docString.str(), intValidatorP);
  }

  {
    string parameterName("debug_procs");
    typedef Zoltan2::IntegerRangeListValidator<int> irl_t;
    RCP<const irl_t> intRangeValidatorP = rcp(new irl_t);
    ostringstream docString;
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

  cout << "Parameter list: " << endl;
  cout << obj << endl;

  std::ofstream of;
  of.open("params.xml");
  of << obj << endl;
  of.close();

  // Read parameter list in from XML file.

  Teuchos::ParameterXMLFileReader rdr("params.xml");
  Teuchos::ParameterList newpl = rdr.getParameters();
  Teuchos::XMLObject objnew = plw.toXML(newpl);

  cout << "After reading in from XML file: " << endl;
  cout << objnew << endl;
}

