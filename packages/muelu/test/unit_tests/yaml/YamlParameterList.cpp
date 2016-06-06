// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER


#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <MueLu_Utilities.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Exceptions.hpp>
#include <MueLu_YamlParser_def.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_Exceptions.hpp>

static bool paramListsEquivalent(const ParameterList& lhs, const ParameterList& rhs)
{
  typedef Teuchos::ParameterList::ConstIterator Iter;
  Iter i = lhs.begin();
  Iter j = rhs.begin();
  if(lhs.name() != rhs.name())
    return false;
  for(; i != lhs.end(); i++)
  {
    const std::string& key = lhs.name(i);
    const Teuchos::ParameterEntry& val1 = lhs.entry(i);
    //check that rhs also contains this key
    if(!rhs.isParameter(key))
    {
      return false;
    }
    const Teuchos::ParameterEntry& val2 = rhs.getEntry(key);
    const Teuchos::any& any1 = val1.getAny(false);
    const Teuchos::any& any2 = val2.getAny(false);
    //check that types match
    if(any1.type() != any2.type())
    {
      return false;
    }
    //check for parameter list special case (don't use operator==)
    if(any1.type() == typeid(Teuchos::ParameterList))
    {
      if(!paramListsEquivalent(Teuchos::any_cast<Teuchos::ParameterList>(any1), Teuchos::any_cast<Teuchos::ParameterList>(any2)))
      {
        return false;
      }
    }
    else
    {
      //otherwise, use == to compare the values
      if(!(val1 == val2))
      {
        return false;
      }
    }
    j++;
  }
  //lists must have same # of entries
  if(j != rhs.end())
  {
    return false;
  }
  return true;
}

namespace MueLuTests
{
  TEUCHOS_UNIT_TEST(YAML, XmlEquivalence)
  {
    std::string matchStems[4] = {"Match1", "Match2", "Match3", "Match4"};
    for(int i = 0; i < 4; i++)
    {
      std::string xmlFile = std::string("yaml/") + matchStems[i] + ".xml";
      std::string yamlFile = std::string("yaml/") + matchStems[i] + ".yaml";
      RCP<ParameterList> xmlList = Teuchos::getParametersFromXmlFile(xmlFile);
      RCP<ParameterList> yamlList = MueLu::getParametersFromYamlFile(yamlFile);
      TEST_EQUALITY(paramListsEquivalent(*xmlList, *yamlList), true);
    }
  }
  TEUCHOS_UNIT_TEST(YAML, IntVsDouble)
  {
    //YAML2 has a double param that has the same name/value as an int param in XML2
    //YAML reader should recognize the double and the param lists should not be equivalent
    RCP<ParameterList> xmlList = Teuchos::getParametersFromXmlFile("yaml/IntVsDouble.xml");
    RCP<ParameterList> yamlList = MueLu::getParametersFromYamlFile("yaml/IntVsDouble.yaml");
    TEST_EQUALITY(paramListsEquivalent(*xmlList, *yamlList), false);
  }
  TEUCHOS_UNIT_TEST(YAML, IllegalKeyString)
  {
    TEST_THROW(MueLu::getParametersFromYamlFile("yaml/IllegalKeyString.yaml");, YAML::ParserException);
  }
  TEUCHOS_UNIT_TEST(YAML, IntAndDoubleArray)
  {
    int correctInts[5] = {2, 3, 5, 7, 11};
    //the last number with 10 dec digits of precision should test correct double conversion
    double correctDoubles[5] = {2.718, 3.14159, 1.618, 1.23456789, 42.1337};
    TEST_NOTHROW({
      RCP<Teuchos::ParameterList> params = MueLu::getParametersFromYamlFile("yaml/Arrays.yaml");
      //Retrieve arrays from a specific sublist (tests the mixed nesting of sequence/map)
      ParameterList& sublist = params->get<ParameterList>("smoother: params");
      Teuchos::Array<int>& intArr = sublist.get<Teuchos::Array<int> >("intArray");
      Teuchos::Array<double>& doubleArr = sublist.get<Teuchos::Array<double> >("doubleArray");
      TEST_EQUALITY(intArr.size(), 5);
      TEST_EQUALITY(doubleArr.size(), 5);
      for(int i = 0; i < 5; i++)
      {
        TEST_EQUALITY(correctInts[i], intArr[i]);
        TEST_EQUALITY(correctDoubles[i], doubleArr[i]);
      }
    });
  }
  TEUCHOS_UNIT_TEST(YAML, InconsistentArrayType)
  {
    std::string correctStrings[5] = {"2", "3", "5", "7", "imastring"};
    double correctDoubles[5] = {2, 3, 1.618, 1.23456789, 42.1337};
    TEST_NOTHROW({
      RCP<ParameterList> plist = MueLu::getParametersFromYamlFile("yaml/InconsistentArray.yaml");
      //verify that stringArray and doubleArray have the correct types and the correct values
      const Teuchos::Array<std::string>& stringArr = plist->get<Teuchos::Array<std::string> >("stringArray");
      const Teuchos::Array<double>& doubleArr = plist->get<Teuchos::Array<double> >("doubleArray");
      for(int i = 0; i < 5; i++)
      {
        if(stringArr[i] != correctStrings[i])
        {
          throw std::runtime_error(std::string("stringArray[") + std::to_string(i) + "] is incorrect.");
        }
        if(doubleArr[i] != correctDoubles[i])
        {
          throw std::runtime_error(std::string("doubleArray value [") + std::to_string(i) + "] is incorrect.");
        }
      }
    });
  }
  TEUCHOS_UNIT_TEST(YAML, MPIBroadcast)
  {
    //load Match1.xml and Match1.yaml on proc 0, broadcast them, and make sure it matches on all procs
    //note: the following code is run on all procs
    RCP<ParameterList> xmlList = rcp(new ParameterList);
    RCP<ParameterList> yamlList = rcp(new ParameterList);
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    Teuchos::updateParametersFromXmlFileAndBroadcast("yaml/Match1.xml", xmlList.ptr(), *comm, true);
    MueLu::updateParametersFromYamlFileAndBroadcast("yaml/Match1.yaml", yamlList.ptr(), *comm, true);
    TEST_EQUALITY(paramListsEquivalent(*xmlList, *yamlList), true);
  }
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(YAML, MueLuConfig, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_IFPACK) && defined(HAVE_MUELU_IFPACK2) && defined(HAVE_MUELU_AMESOS) && defined(HAVE_MUELU_AMESOS2)
    
    RCP<ParameterList> testLists[2] = {Teuchos::null, Teuchos::null};
    testLists[0] = Teuchos::getParametersFromXmlFile("yaml/MueLuConfig.xml");
    testLists[1] = MueLu::getParametersFromYamlFile("yaml/MueLuConfig.yaml");

    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(99);
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    for(int i = 0; i < 2; i++)
    {
      ParameterListInterpreter mueluFactory(*testLists[i], comm);
      RCP<Hierarchy> H = mueluFactory.CreateHierarchy();
      H->GetLevel(0)->Set("A", A);
      mueluFactory.SetupHierarchy(*H);
    }
#   else
    out << "Skipping test because some required packages are not enabled (Tpetra, Epetra, EpetraExt, Ifpack, Ifpack2, Amesos, Amesos2)." << std::endl;
#   endif
  }
/* End YAML unit tests */
#define MUELU_ETI_GROUP(Scalar, LocalOrdinal, GlobalOrdinal, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(YAML, MueLuConfig, Scalar, LocalOrdinal, GlobalOrdinal, Node)

#include <MueLu_ETI_4arg.hpp>
} //namespace MueLuTests
