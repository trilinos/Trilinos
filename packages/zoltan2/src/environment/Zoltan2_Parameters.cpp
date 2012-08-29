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

/*! \file Zoltan2_Parameters.cpp
    \brief Methods that support the Zoltan2 ParameterList
*/

#include <Zoltan2_Parameters.hpp>
#include <Zoltan2_IntegerRangeList.hpp>

#include <Teuchos_StringInputSource.hpp>
#include <Teuchos_XMLParser.hpp>
#include <Teuchos_XMLObject.hpp>
#include <Teuchos_XMLParameterListReader.hpp>
#include <Teuchos_ValidatorXMLConverterDB.hpp>

// A generated header file in {zoltan2-binary-directory}/src
// which defines ZOLTAN2_XML_PARAMETER_STRING.

#include <Zoltan2_XML_Parameters.hpp>  

using namespace std;

namespace Zoltan2 {

/*! \brief  Create a list of all Zoltan2 parameters and validators.
 *
 *  \param pList on return, pList is the parameter list that was created.
 *
 *  This is the validating parameter list that can be
 *  used to process the user's parameter list.
 */

void createAllParameters(Teuchos::ParameterList &pList)
{
  // An XML converter for IntegerRangeListValidator
  // needs to be added to the converter database.

  typedef Zoltan2::IntegerRangeListValidator<int> irl_t;
  typedef Zoltan2::IntegerRangeListValidatorXMLConverter<int> irlConverter_t;

  RCP<const irl_t> intRangeValidatorP = rcp(new irl_t); // dummy
  RCP<irlConverter_t > converter = rcp(new irlConverter_t);
  Teuchos::ValidatorXMLConverterDB::addConverter(
        intRangeValidatorP,
        converter);

  // Create a Teuchos::ParameterList from an XML string.

  std::string xmlParameterString(ZOLTAN2_XML_PARAMETER_STRING);
  Teuchos::StringInputSource src(xmlParameterString);
  Teuchos::XMLParser parser(src.stream());
  Teuchos::XMLObject xmlObj = parser.parse();
  Teuchos::XMLParameterListReader rdr;
  pList = rdr.toParameterList(xmlObj);
}

/*! \brief  Create a parameter list that can validate a
 *          specific list of parameters.
 *
 *   \param plSome   the user's parameter list
 *   \param plAll    a parameter list with all Zoltan2 parameters
 *                   and their validators
 *   \param plVal     on return is a parameter list with all of the
 *                       user's parameters, but with validators.
 *
 *  Environment::commitParameters() calls 
 *  validateParametersAndSetDefaults() on the user's parameter list
 *  rather than validateParameters() because we want the validators' 
 *  validateAndModify() methods to be called rather than the validate()
 *  methods.  But unfortunately, validateAndModify() in addition to
 *  modifying the parameter if necessary also sets it to a default if
 *  the parameter does not appear in the user's parameter list.
 *
 *  We want the user's parameters to be modified, but we do not want unset 
 *  parameters to appear in the validated user's parameter list. To
 *  achieve this, we create a validating list that contains only the
 *  parameters that appear in the user's parameter list.
 *  
 */
static void setValidatorsInList(
  const Teuchos::ParameterList &plSome,   // in: user's parameters
  const Teuchos::ParameterList &plAll,    // in: validators for all params
  Teuchos::ParameterList &plVal)          // out: validators for user's params
{
  ParameterList::ConstIterator next = plSome.begin();

  while (next != plSome.end()){

    const std::string &name = next->first;
    const ParameterEntry &entrySome = plSome.getEntry(name);
    const ParameterEntry &entryAll = plAll.getEntry(name);

    if (entrySome.isList()){
      const ParameterList &sublistSome = plSome.sublist(name);   // get
      const ParameterList &sublistAll = plAll.sublist(name);     // get
      ParameterList &sublistVal = plVal.sublist(name);     // create & get
      setValidatorsInList(sublistSome, sublistAll, sublistVal);
    }
    else{
      plVal.setEntry(name, entryAll);
    }

    ++next;
  }
}

/*! \brief Create a list by adding validators to the users parameter list.
 *  \param plIn  the user's parameter list
 *  \param plOut  a new parameter list which is the user's list with
 *                     our validators added.
 */

void createValidatorList(
   const Teuchos::ParameterList &plIn,
   Teuchos::ParameterList &plOut)
{
  ParameterList allParameters;

  createAllParameters(allParameters);

  setValidatorsInList(plIn, allParameters, plOut);
}

// Why isn't there a Teuchos method that does this?

void printListDocumentation(
  const Teuchos::ParameterList &pl,
  std::ostream &os,
  std::string listNames)
{
  using std::string;

  if (listNames.size() == 0)
    listNames = string("top");

  Array<string> subLists;
  ParameterList::ConstIterator next = pl.begin();

  while (next != pl.end()){
    const string &name = next->first;
    const ParameterEntry &entry = pl.getEntry(name);

    if (entry.isList()){
      subLists.append(name);
    }
    else{
      string doc = entry.docString();
      os << "List: "<< listNames << ", parameter: " << name << "\n";
      if (doc.size())
        os << doc << "\n";
    }

    ++next;
  }

  for (int i=0; i < subLists.size(); i++){
    string newListName = listNames + string("/") + subLists[i];
    const ParameterList &sublist = pl.sublist(subLists[i]);
    printListDocumentation(sublist, os, newListName);
  }
}


}  //namespace Zoltan2
