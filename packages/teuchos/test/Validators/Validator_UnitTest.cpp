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
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"


// 2010/07/30: rabartl: Here I just added all the unit tests to the Teuchos
// namespace to remove some clutter.


namespace Teuchos {


/**
 * Tests Number Validators.
 */
TEUCHOS_UNIT_TEST(Teuchos_Validators, numberValidators)
{
	/*
	 * Testing Int Validator.
	 */
	RCP<ParameterList> intList =
    rcp(new ParameterList("Int List"));
	RCP<EnhancedNumberValidator<int> > intVali =
    rcp(new EnhancedNumberValidator<int>(0,10,4));
	TEST_ASSERT(intVali->getMin() == 0);
	TEST_ASSERT(intVali->getMax() == 10);
	TEST_ASSERT(intVali->getStep() == 4);
	TEST_ASSERT(intVali->hasMin());
	TEST_ASSERT(intVali->hasMax());
	RCP<EnhancedNumberValidator<int> > intVali2 =
    rcp(new EnhancedNumberValidator<int>());
	TEST_ASSERT(!intVali2->hasMin());
	TEST_ASSERT(!intVali2->hasMax());
	TEST_ASSERT(intVali2->getMin() == std::numeric_limits<int>::min());
	TEST_ASSERT(intVali2->getMax() == std::numeric_limits<int>::max());
	TEST_ASSERT(intVali2->getStep() == EnhancedNumberTraits<int>::defaultStep());
	intList->set("Int Parameter", 5, "int parameter", intVali);
	TEST_NOTHROW(intList->validateParameters(*intList));
	TEST_THROW(intList->set("Int Parameter", 11),
    Exceptions::InvalidParameterValue);
	TEST_THROW(intList->set("Double Parameter", 5.0, "double parameter", intVali),
    Exceptions::InvalidParameterType);

	/*
	 * Testing Short Validator.
	 */
	RCP<ParameterList> shortList =
    rcp(new ParameterList("Short List"));
	RCP<EnhancedNumberValidator<short> > shortVali =
    rcp(new EnhancedNumberValidator<short>(0,10,4));
	TEST_ASSERT(shortVali->getMin() == 0);
	TEST_ASSERT(shortVali->getMax() == 10);
	TEST_ASSERT(shortVali->getStep() == 4);
	TEST_ASSERT(shortVali->hasMin());
	TEST_ASSERT(shortVali->hasMax());
	RCP<EnhancedNumberValidator<short> > shortVali2 =
    rcp(new EnhancedNumberValidator<short>());
	TEST_ASSERT(!shortVali2->hasMin());
	TEST_ASSERT(!shortVali2->hasMax());
	TEST_ASSERT(shortVali2->getMin() == std::numeric_limits<short>::min());
	TEST_ASSERT(shortVali2->getMax() == std::numeric_limits<short>::max());
	TEST_ASSERT(shortVali2->getStep() == EnhancedNumberTraits<short>::defaultStep());
	shortList->set("Short Parameter", (short)5, "short parameter", shortVali);
	TEST_NOTHROW(shortList->validateParameters(*shortList));
	TEST_THROW(shortList->set("Short Parameter", (short)11),
    Exceptions::InvalidParameterValue);
	TEST_THROW(shortList->set("Double Parameter", 5.0, "double parameter", shortVali),
    Exceptions::InvalidParameterType);

	/*
	 * Testing Float Validator.
	 */
	RCP<ParameterList> floatList = rcp(new ParameterList("Float List"));
	RCP<EnhancedNumberValidator<float> > floatVali =
    rcp(new EnhancedNumberValidator<float>(0,10.0,4.0,6));
	TEST_ASSERT(floatVali->getMin() == 0.0);
	TEST_ASSERT(floatVali->getMax() == 10.0);
	TEST_ASSERT(floatVali->getStep() == 4.0);
	TEST_ASSERT(floatVali->getPrecision() == 6);
	TEST_ASSERT(floatVali->hasMin());
	TEST_ASSERT(floatVali->hasMax());
	RCP<EnhancedNumberValidator<float> > floatVali2 =
    rcp(new EnhancedNumberValidator<float>());
	TEST_ASSERT(!floatVali2->hasMin());
	TEST_ASSERT(!floatVali2->hasMax());
	TEST_ASSERT(floatVali2->getMin() == EnhancedNumberTraits<float>::min());
	TEST_ASSERT(floatVali2->getMax() == EnhancedNumberTraits<float>::max());
	TEST_ASSERT(floatVali2->getStep() == EnhancedNumberTraits<float>::defaultStep());
	TEST_ASSERT(floatVali2->getPrecision() == EnhancedNumberTraits<float>::defaultPrecision());
	floatList->set("Float Parameter", (float)5.0, "float parameter", floatVali);
	TEST_NOTHROW(floatList->validateParameters(*floatList));
	TEST_THROW(floatList->set("Float Parameter", (float)11.0),
    Exceptions::InvalidParameterValue);
	TEST_THROW(floatList->set("Int Parameter", 5, "int parameter", floatVali),
    Exceptions::InvalidParameterType);

	/*
	 * Testing Double Validator.
	 */
	RCP<ParameterList> doubleList = rcp(new ParameterList("Double List"));
	RCP<EnhancedNumberValidator<double> > doubleVali =
    rcp(new EnhancedNumberValidator<double>(0,10.0,4.0,6));
	TEST_ASSERT(doubleVali->getMin() == 0.0);
	TEST_ASSERT(doubleVali->getMax() == 10.0);
	TEST_ASSERT(doubleVali->getStep() == 4.0);
	TEST_ASSERT(doubleVali->getPrecision() == 6);
	TEST_ASSERT(doubleVali->hasMin());
	TEST_ASSERT(doubleVali->hasMax());
	RCP<EnhancedNumberValidator<double> > doubleVali2 =
    rcp(new EnhancedNumberValidator<double>());
	TEST_ASSERT(!doubleVali2->hasMin());
	TEST_ASSERT(!doubleVali2->hasMax());
	TEST_ASSERT(doubleVali2->getMin() == EnhancedNumberTraits<double>::min());
	TEST_ASSERT(doubleVali2->getMax() == EnhancedNumberTraits<double>::max());
	TEST_ASSERT(doubleVali2->getStep() == EnhancedNumberTraits<double>::defaultStep());
	TEST_ASSERT(doubleVali2->getPrecision() == EnhancedNumberTraits<double>::defaultPrecision());
	doubleList->set("Double Parameter", (double)5.0, "double parameter", doubleVali);
	TEST_NOTHROW(doubleList->validateParameters(*doubleList));
	TEST_THROW(doubleList->set("Double Parameter", (double)11.0),
    Exceptions::InvalidParameterValue);
	TEST_THROW(doubleList->set("Int Parameter", 5, "int parameter", doubleVali),
    Exceptions::InvalidParameterType);
}


/*
 * Testing StringValidator.
 */
TEUCHOS_UNIT_TEST(Teuchos_Validators, stringValidator)
{
	RCP<ParameterList> stringList = rcp(new ParameterList("String List"));
	Array<std::string> stringVals = tuple<std::string>("str1", "str2", "str3");
	RCP<StringValidator> stringVali = rcp(new StringValidator(stringVals));
	RCP<const Array<std::string> > valiVals = stringVali->validStringValues();
  /*bool local_success = true;
  for(int i =0; i<valiVals.size() ++i){
	  TEST_ARRAY_ELE_EQUALITY(*valiVals, i, stringVals[i]);
  }
  if (local_success) out << "passed\n";
  else success = false;*/
  TEST_COMPARE_ARRAYS(*valiVals, stringVals);
	TEST_NOTHROW(stringList->set("String param1", "str1", "a string parameter", stringVali));
	TEST_THROW(stringList->set("String param2", "not in list", "a string parameter", stringVali),
    Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("int param", 5, "a int parameter", stringVali),
    Exceptions::InvalidParameterType);
}


/*
 * Testing FileNameValidator.
 */
TEUCHOS_UNIT_TEST(Teuchos_Validators, fileNameValidator)
{
	RCP<ParameterList> fileNameList = rcp(new ParameterList("Filename List"));
	RCP<FileNameValidator> fileNameVali = rcp(new FileNameValidator(true));
	TEST_ASSERT(fileNameVali->fileMustExist());
	fileNameVali->setFileMustExist(false);
	TEST_ASSERT(!fileNameVali->fileMustExist());
	TEST_NOTHROW(fileNameList->set("File name param", "../path", "file name parameter",
      fileNameVali));
	TEST_THROW(fileNameList->set("int param", 5, "int parameter", fileNameVali),
    Exceptions::InvalidParameterType);
	fileNameVali->setFileMustExist(true);
	TEST_NOTHROW(fileNameList->set("file name param", "testFile.txt", "a file name", fileNameVali));
	TEST_THROW(fileNameList->set("file name param", "doesntexist.txt", "a file name", fileNameVali),
    Exceptions::InvalidParameterValue);
}


/*
 * Testing Array Validators.
 */
TEUCHOS_UNIT_TEST(Teuchos_Validators, arrayValidators)
{

	/*
	 * Testing StringArrayValidator.
	 */
	RCP<ParameterList> stringList = rcp(new ParameterList("String List"));
	Array<std::string> stringVals = tuple<std::string>("str1", "str2", "str3");
	RCP<StringValidator> stringVali = rcp(new StringValidator(stringVals));
	RCP<ArrayStringValidator> stringArrayVali = rcp(new ArrayStringValidator(stringVali));
	TEST_ASSERT(stringVali.get() == stringArrayVali->getPrototype().get());
	Array<std::string> stringArray = tuple<std::string>("str2","str3","str1","str3","str2");
	TEST_NOTHROW(stringList->set("String Array Param", stringArray, "string array parameter", stringArrayVali));
	Array<std::string> badStringArray = tuple<std::string>("not valid","str3","str1","str3","str2");
	TEST_THROW(stringList->set("String Array Param", badStringArray, "string array parameter", stringArrayVali),
    Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("Int param", 5, "int parameter", stringArrayVali),
    Exceptions::InvalidParameterType);
	Array<long> longArray = tuple<long>((long)5,(long)5,(long)3);
	TEST_THROW(stringList->set("Long array param", longArray, "long array parameter", stringArrayVali),
    Exceptions::InvalidParameterType);
	
	/*
	 * Testing Int ArrayValidator.
	 */
	RCP<ParameterList> intList = rcp(new ParameterList("Int List"));
	RCP<EnhancedNumberValidator<int> > intVali = rcp(new EnhancedNumberValidator<int>(0, 10));
	RCP<ArrayNumberValidator<int> > intArrayVali = rcp(new ArrayNumberValidator<int>(intVali));
	TEST_ASSERT(intVali.get() == intArrayVali->getPrototype().get());
	Array<int> intArray = tuple<int>(1,4,2,5);
	TEST_NOTHROW(intList->set("int array param", intArray, "int array parameter", intArrayVali));
	Array<int> intBadArray = tuple<int>(11,4,2,5);
	TEST_THROW(intList->set("int bad array param", intBadArray, "int bad array parameter", intArrayVali),
    Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("Long array param", longArray, "long array parameter", intArrayVali),
    Exceptions::InvalidParameterType);

	/*
	 * Testing Short ArrayValidator.
	 */
	RCP<ParameterList> shortList = rcp(new ParameterList("Short List"));
	RCP<EnhancedNumberValidator<short> > shortVali =
    rcp(new EnhancedNumberValidator<short>(0, 10));
	RCP<ArrayNumberValidator<short> > shortArrayVali =
    rcp(new ArrayNumberValidator<short>(shortVali));
	TEST_ASSERT(shortVali.get() == shortArrayVali->getPrototype().get());
	Array<short> shortArray = tuple<short>(1,4,2,5);
	TEST_NOTHROW(shortList->set("short array param", shortArray, "short array parameter", shortArrayVali));
	Array<short> shortBadArray = tuple<short>(11,4,2,5);
	TEST_THROW(shortList->set("short bad array param", shortBadArray, "short bad array parameter", shortArrayVali),
    Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("Long array param", longArray, "long array parameter", shortArrayVali),
    Exceptions::InvalidParameterType);

	/*
	 * Testing Float ArrayValidator.
	 */
	RCP<ParameterList> floatList = rcp(new ParameterList("Float List"));
	RCP<EnhancedNumberValidator<float> > floatVali =
    rcp(new EnhancedNumberValidator<float>(0.0, 10.0));
	RCP<ArrayNumberValidator<float> > floatArrayVali =
    rcp(new ArrayNumberValidator<float>(floatVali));
	TEST_ASSERT(floatVali.get() == floatArrayVali->getPrototype().get());
	Array<float> floatArray = tuple<float>(1.0,4.0,2.0,5.0);
	TEST_NOTHROW(floatList->set("float array param", floatArray, "float array parameter", floatArrayVali));
	Array<float> floatBadArray = tuple<float>(11.0,4.0,2.0,5.0);
	TEST_THROW(floatList->set("float bad array param", floatBadArray, "float bad array parameter", floatArrayVali),
    Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("Long array param", longArray, "long array parameter", floatArrayVali),
    Exceptions::InvalidParameterType);

	/*
	 * Testing Double ArrayValidator.
	 */
	RCP<ParameterList> doubleList = rcp(new ParameterList("Double List"));
	RCP<EnhancedNumberValidator<double> > doubleVali =
    rcp(new EnhancedNumberValidator<double>(0.0, 10.0));
	RCP<ArrayNumberValidator<double> > doubleArrayVali =
    rcp(new ArrayNumberValidator<double>(doubleVali));
	TEST_ASSERT(doubleVali.get() == doubleArrayVali->getPrototype().get());
	Array<double> doubleArray = tuple<double>(1.0,4.0,2.0,5.0);
	TEST_NOTHROW(doubleList->set("double array param", doubleArray, "double array parameter", doubleArrayVali));
	Array<double> doubleBadArray = tuple<double>(11.0,4.0,2.0,5.0);
	TEST_THROW(doubleList->set("double bad array param", doubleBadArray, "double bad array parameter", doubleArrayVali),
    Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("Long array param", longArray, "long array parameter", doubleArrayVali),
    Exceptions::InvalidParameterType);

	/*
	 * Testing FileName ArrayValidator.
	 */
	RCP<ParameterList> fileNameList = rcp(new ParameterList("Filename List"));
	RCP<FileNameValidator> fileNameVali = rcp(new FileNameValidator(true));
	RCP<ArrayFileNameValidator> arrayFileNameVali = rcp(new ArrayFileNameValidator(fileNameVali));
	TEST_ASSERT(arrayFileNameVali->getPrototype().get() == fileNameVali.get());
	Array<std::string> fileNameArray = tuple<std::string>("testFile.txt", "testFile2.txt", "testFile3.txt");
	Array<std::string> fileNameBadArray = tuple<std::string>("doesnexist.txt", "testFile2.txt", "testFile3.txt");
	TEST_NOTHROW(fileNameList->set("File name array", fileNameArray, "file name array parameter", arrayFileNameVali));
	TEST_THROW(fileNameList->set("Bad File name array", fileNameBadArray, "bad file name array parameter", arrayFileNameVali),
    Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("Long array param", longArray, "long array parameter", arrayFileNameVali),
    Exceptions::InvalidParameterType);
}


} // namespace Teuchos

