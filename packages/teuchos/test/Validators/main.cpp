// @HEADER // ***********************************************************************
// 
//         Optika: A Tool For Developing Parameter Obtaining GUIs
//                Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, with Sandia Corporation, the 
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Kurtis Nusbaum (klnusbaum@gmail.com) 
// 
// ***********************************************************************
// @HEADER
#include "Teuchos_LocalTestingHelpers.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

/**
 * Tests Number Validators.
 */
int testNumberValis(Teuchos::FancyOStream &out){
	bool success = true;
	/*
	 * Testing Int Validator.
	 */
	Teuchos::RCP<Teuchos::ParameterList> intList = Teuchos::rcp(new Teuchos::ParameterList("Int List"));
	Teuchos::RCP<Teuchos::EnhancedNumberValidator<int> > intVali = Teuchos::rcp(new Teuchos::EnhancedNumberValidator<int>(0,10,4));
	TEST_ASSERT(intVali->getMin() == 0);
	TEST_ASSERT(intVali->getMax() == 10);
	TEST_ASSERT(intVali->getStep() == 4);
	TEST_ASSERT(intVali->hasMin());
	TEST_ASSERT(intVali->hasMax());
	Teuchos::RCP<Teuchos::EnhancedNumberValidator<int> > intVali2 = Teuchos::rcp(new Teuchos::EnhancedNumberValidator<int>());
	TEST_ASSERT(!intVali2->hasMin());
	TEST_ASSERT(!intVali2->hasMax());
	TEST_ASSERT(intVali2->getMin() == std::numeric_limits<int>::min());
	TEST_ASSERT(intVali2->getMax() == std::numeric_limits<int>::max());
	TEST_ASSERT(intVali2->getStep() == Teuchos::EnhancedNumberValidator<int>::intDefaultStep);
	intList->set("Int Parameter", 5, "int parameter", intVali);
	TEST_NOTHROW(intList->validateParameters(*intList));
	TEST_THROW(intList->set("Int Parameter", 11), Teuchos::Exceptions::InvalidParameterValue);
	TEST_THROW(intList->set("Double Parameter", 5.0, "double parameter", intVali), Teuchos::Exceptions::InvalidParameterType);

	/*
	 * Testing Short Validator.
	 */
	Teuchos::RCP<Teuchos::ParameterList> shortList = Teuchos::rcp(new Teuchos::ParameterList("Short List"));
	Teuchos::RCP<Teuchos::EnhancedNumberValidator<short> > shortVali = Teuchos::rcp(new Teuchos::EnhancedNumberValidator<short>(0,10,4));
	TEST_ASSERT(shortVali->getMin() == 0);
	TEST_ASSERT(shortVali->getMax() == 10);
	TEST_ASSERT(shortVali->getStep() == 4);
	TEST_ASSERT(shortVali->hasMin());
	TEST_ASSERT(shortVali->hasMax());
	Teuchos::RCP<Teuchos::EnhancedNumberValidator<short> > shortVali2 = Teuchos::rcp(new Teuchos::EnhancedNumberValidator<short>());
	TEST_ASSERT(!shortVali2->hasMin());
	TEST_ASSERT(!shortVali2->hasMax());
	TEST_ASSERT(shortVali2->getMin() == std::numeric_limits<short>::min());
	TEST_ASSERT(shortVali2->getMax() == std::numeric_limits<short>::max());
	TEST_ASSERT(shortVali2->getStep() == Teuchos::EnhancedNumberValidator<short>::shortDefaultStep);
	shortList->set("Short Parameter", (short)5, "short parameter", shortVali);
	TEST_NOTHROW(shortList->validateParameters(*shortList));
	TEST_THROW(shortList->set("Short Parameter", (short)11), Teuchos::Exceptions::InvalidParameterValue);
	TEST_THROW(shortList->set("Double Parameter", 5.0, "double parameter", shortVali), Teuchos::Exceptions::InvalidParameterType);

	/*
	 * Testing Float Validator.
	 */
	Teuchos::RCP<Teuchos::ParameterList> floatList = Teuchos::rcp(new Teuchos::ParameterList("Float List"));
	Teuchos::RCP<Teuchos::EnhancedNumberValidator<float> > floatVali = Teuchos::rcp(new Teuchos::EnhancedNumberValidator<float>(0,10.0,4.0,6));
	TEST_ASSERT(floatVali->getMin() == 0.0);
	TEST_ASSERT(floatVali->getMax() == 10.0);
	TEST_ASSERT(floatVali->getStep() == 4.0);
	TEST_ASSERT(floatVali->getPrecision() == 6);
	TEST_ASSERT(floatVali->hasMin());
	TEST_ASSERT(floatVali->hasMax());
	Teuchos::RCP<Teuchos::EnhancedNumberValidator<float> > floatVali2 = Teuchos::rcp(new Teuchos::EnhancedNumberValidator<float>());
	TEST_ASSERT(!floatVali2->hasMin());
	TEST_ASSERT(!floatVali2->hasMax());
	TEST_ASSERT(floatVali2->getMin() == -std::numeric_limits<float>::max());
	TEST_ASSERT(floatVali2->getMax() == std::numeric_limits<float>::max());
	TEST_ASSERT(floatVali2->getStep() == Teuchos::EnhancedNumberValidator<float>::floatDefaultStep);
	TEST_ASSERT(floatVali2->getPrecision() == Teuchos::EnhancedNumberValidator<float>::floatDefaultPrecision);
	floatList->set("Float Parameter", (float)5.0, "float parameter", floatVali);
	TEST_NOTHROW(floatList->validateParameters(*floatList));
	TEST_THROW(floatList->set("Float Parameter", (float)11.0), Teuchos::Exceptions::InvalidParameterValue);
	TEST_THROW(floatList->set("Int Parameter", 5, "int parameter", floatVali), Teuchos::Exceptions::InvalidParameterType);

	/*
	 * Testing Double Validator.
	 */
	Teuchos::RCP<Teuchos::ParameterList> doubleList = Teuchos::rcp(new Teuchos::ParameterList("Double List"));
	Teuchos::RCP<Teuchos::EnhancedNumberValidator<double> > doubleVali = Teuchos::rcp(new Teuchos::EnhancedNumberValidator<double>(0,10.0,4.0,6));
	TEST_ASSERT(doubleVali->getMin() == 0.0);
	TEST_ASSERT(doubleVali->getMax() == 10.0);
	TEST_ASSERT(doubleVali->getStep() == 4.0);
	TEST_ASSERT(doubleVali->getPrecision() == 6);
	TEST_ASSERT(doubleVali->hasMin());
	TEST_ASSERT(doubleVali->hasMax());
	Teuchos::RCP<Teuchos::EnhancedNumberValidator<double> > doubleVali2 = Teuchos::rcp(new Teuchos::EnhancedNumberValidator<double>());
	TEST_ASSERT(!doubleVali2->hasMin());
	TEST_ASSERT(!doubleVali2->hasMax());
	TEST_ASSERT(doubleVali2->getMin() == -std::numeric_limits<double>::max());
	TEST_ASSERT(doubleVali2->getMax() == std::numeric_limits<double>::max());
	TEST_ASSERT(doubleVali2->getStep() == Teuchos::EnhancedNumberValidator<double>::doubleDefaultStep);
	TEST_ASSERT(doubleVali2->getPrecision() == Teuchos::EnhancedNumberValidator<double>::doubleDefaultPrecision);
	doubleList->set("Double Parameter", (double)5.0, "double parameter", doubleVali);
	TEST_NOTHROW(doubleList->validateParameters(*doubleList));
	TEST_THROW(doubleList->set("Double Parameter", (double)11.0), Teuchos::Exceptions::InvalidParameterValue);
	TEST_THROW(doubleList->set("Int Parameter", 5, "int parameter", doubleVali), Teuchos::Exceptions::InvalidParameterType);
	return (success ? 0:1);
}

/*
 * Testing StringValidator.
 */
int testStringVali(Teuchos::FancyOStream &out){
	bool success = true;
	Teuchos::RCP<Teuchos::ParameterList> stringList = Teuchos::rcp(new Teuchos::ParameterList("String List"));
	Teuchos::Array<std::string> stringVals = Teuchos::tuple<std::string>("str1", "str2", "str3");
	Teuchos::RCP<Teuchos::StringValidator> stringVali = Teuchos::rcp(new Teuchos::StringValidator(stringVals));
	Teuchos::RCP<const Teuchos::Array<std::string> > valiVals = stringVali->validStringValues();
	bool local_success = true;
	for(int i=0;i<stringVals.size(); i++){
		TEST_ARRAY_ELE_EQUALITY(*valiVals, i, stringVals[i]);
	}
	if(!local_success){
		success = false;
	}
	TEST_NOTHROW(stringList->set("String param1", "str1", "a string parameter", stringVali));
	TEST_THROW(stringList->set("String param2", "not in list", "a string parameter", stringVali), Teuchos::Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("int param", 5, "a int parameter", stringVali), Teuchos::Exceptions::InvalidParameterType);

	return (success ? 0:1);
}

/*
 * Testing FileNameValidator.
 */
int testFileNameVali(Teuchos::FancyOStream &out){
	bool success = true;
	Teuchos::RCP<Teuchos::ParameterList> fileNameList = Teuchos::rcp(new Teuchos::ParameterList("Filename List"));
	Teuchos::RCP<Teuchos::FileNameValidator> fileNameVali = Teuchos::rcp(new Teuchos::FileNameValidator(true));
	TEST_ASSERT(fileNameVali->fileMustExist());
	fileNameVali->setFileMustExist(false);
	TEST_ASSERT(!fileNameVali->fileMustExist());
	TEST_NOTHROW(fileNameList->set("File name param", "../path", "file name parameter", fileNameVali));
	TEST_THROW(fileNameList->set("int param", 5, "int parameter", fileNameVali), Teuchos::Exceptions::InvalidParameterType);
	fileNameVali->setFileMustExist(true);
	TEST_NOTHROW(fileNameList->set("file name param", "testFile.txt", "a file name", fileNameVali));
	TEST_THROW(fileNameList->set("file name param", "doesntexist.txt", "a file name", fileNameVali), Teuchos::Exceptions::InvalidParameterValue);

	return (success ? 0:1);
}

/*
 * Testing Array Validators.
 */
int testArrayValis(Teuchos::FancyOStream &out){
	bool success = true;
	/*
	 * Testing StringArrayValidator.
	 */
	Teuchos::RCP<Teuchos::ParameterList> stringList = Teuchos::rcp(new Teuchos::ParameterList("String List"));
	Teuchos::Array<std::string> stringVals = Teuchos::tuple<std::string>("str1", "str2", "str3");
	Teuchos::RCP<Teuchos::StringValidator> stringVali = Teuchos::rcp(new Teuchos::StringValidator(stringVals));
	Teuchos::RCP<Teuchos::ArrayStringValidator> stringArrayVali = Teuchos::rcp(new Teuchos::ArrayStringValidator(stringVali));
	TEST_ASSERT(stringVali.get() == stringArrayVali->getPrototype().get());
	Teuchos::Array<std::string> stringArray = Teuchos::tuple<std::string>("str2","str3","str1","str3","str2");
	TEST_NOTHROW(stringList->set("String Array Param", stringArray, "string array parameter", stringArrayVali));
	Teuchos::Array<std::string> badStringArray = Teuchos::tuple<std::string>("not valid","str3","str1","str3","str2");
	TEST_THROW(stringList->set("String Array Param", badStringArray, "string array parameter", stringArrayVali), Teuchos::Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("Int param", 5, "int parameter", stringArrayVali), Teuchos::Exceptions::InvalidParameterType);
	Teuchos::Array<long> longArray = Teuchos::tuple<long>((long)5,(long)5,(long)3);
	TEST_THROW(stringList->set("Long array param", longArray, "long array parameter", stringArrayVali), Teuchos::Exceptions::InvalidParameterType);

	
	/*
	 * Testing Int ArrayValidator.
	 */
	Teuchos::RCP<Teuchos::ParameterList> intList = Teuchos::rcp(new Teuchos::ParameterList("Int List"));
	Teuchos::RCP<Teuchos::EnhancedNumberValidator<int> > intVali = Teuchos::rcp(new Teuchos::EnhancedNumberValidator<int>(0, 10));
	Teuchos::RCP<Teuchos::ArrayNumberValidator<int> > intArrayVali = Teuchos::rcp(new Teuchos::ArrayNumberValidator<int>(intVali));
	TEST_ASSERT(intVali.get() == intArrayVali->getPrototype().get());
	Teuchos::Array<int> intArray = Teuchos::tuple<int>(1,4,2,5);
	TEST_NOTHROW(intList->set("int array param", intArray, "int array parameter", intArrayVali));
	Teuchos::Array<int> intBadArray = Teuchos::tuple<int>(11,4,2,5);
	TEST_THROW(intList->set("int bad array param", intBadArray, "int bad array parameter", intArrayVali), Teuchos::Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("Long array param", longArray, "long array parameter", intArrayVali), Teuchos::Exceptions::InvalidParameterType);

	/*
	 * Testing Short ArrayValidator.
	 */
	Teuchos::RCP<Teuchos::ParameterList> shortList = Teuchos::rcp(new Teuchos::ParameterList("Short List"));
	Teuchos::RCP<Teuchos::EnhancedNumberValidator<short> > shortVali = Teuchos::rcp(new Teuchos::EnhancedNumberValidator<short>(0, 10));
	Teuchos::RCP<Teuchos::ArrayNumberValidator<short> > shortArrayVali = Teuchos::rcp(new Teuchos::ArrayNumberValidator<short>(shortVali));
	TEST_ASSERT(shortVali.get() == shortArrayVali->getPrototype().get());
	Teuchos::Array<short> shortArray = Teuchos::tuple<short>(1,4,2,5);
	TEST_NOTHROW(shortList->set("short array param", shortArray, "short array parameter", shortArrayVali));
	Teuchos::Array<short> shortBadArray = Teuchos::tuple<short>(11,4,2,5);
	TEST_THROW(shortList->set("short bad array param", shortBadArray, "short bad array parameter", shortArrayVali), Teuchos::Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("Long array param", longArray, "long array parameter", shortArrayVali), Teuchos::Exceptions::InvalidParameterType);

	/*
	 * Testing Float ArrayValidator.
	 */
	Teuchos::RCP<Teuchos::ParameterList> floatList = Teuchos::rcp(new Teuchos::ParameterList("Float List"));
	Teuchos::RCP<Teuchos::EnhancedNumberValidator<float> > floatVali = Teuchos::rcp(new Teuchos::EnhancedNumberValidator<float>(0.0, 10.0));
	Teuchos::RCP<Teuchos::ArrayNumberValidator<float> > floatArrayVali = Teuchos::rcp(new Teuchos::ArrayNumberValidator<float>(floatVali));
	TEST_ASSERT(floatVali.get() == floatArrayVali->getPrototype().get());
	Teuchos::Array<float> floatArray = Teuchos::tuple<float>(1.0,4.0,2.0,5.0);
	TEST_NOTHROW(floatList->set("float array param", floatArray, "float array parameter", floatArrayVali));
	Teuchos::Array<float> floatBadArray = Teuchos::tuple<float>(11.0,4.0,2.0,5.0);
	TEST_THROW(floatList->set("float bad array param", floatBadArray, "float bad array parameter", floatArrayVali), Teuchos::Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("Long array param", longArray, "long array parameter", floatArrayVali), Teuchos::Exceptions::InvalidParameterType);

	/*
	 * Testing Double ArrayValidator.
	 */
	Teuchos::RCP<Teuchos::ParameterList> doubleList = Teuchos::rcp(new Teuchos::ParameterList("Double List"));
	Teuchos::RCP<Teuchos::EnhancedNumberValidator<double> > doubleVali = Teuchos::rcp(new Teuchos::EnhancedNumberValidator<double>(0.0, 10.0));
	Teuchos::RCP<Teuchos::ArrayNumberValidator<double> > doubleArrayVali = Teuchos::rcp(new Teuchos::ArrayNumberValidator<double>(doubleVali));
	TEST_ASSERT(doubleVali.get() == doubleArrayVali->getPrototype().get());
	Teuchos::Array<double> doubleArray = Teuchos::tuple<double>(1.0,4.0,2.0,5.0);
	TEST_NOTHROW(doubleList->set("double array param", doubleArray, "double array parameter", doubleArrayVali));
	Teuchos::Array<double> doubleBadArray = Teuchos::tuple<double>(11.0,4.0,2.0,5.0);
	TEST_THROW(doubleList->set("double bad array param", doubleBadArray, "double bad array parameter", doubleArrayVali), Teuchos::Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("Long array param", longArray, "long array parameter", doubleArrayVali), Teuchos::Exceptions::InvalidParameterType);

	/*
	 * Testing FileName ArrayValidator.
	 */
	Teuchos::RCP<Teuchos::ParameterList> fileNameList = Teuchos::rcp(new Teuchos::ParameterList("Filename List"));
	Teuchos::RCP<Teuchos::FileNameValidator> fileNameVali = Teuchos::rcp(new Teuchos::FileNameValidator(true));
	Teuchos::RCP<Teuchos::ArrayFileNameValidator> arrayFileNameVali = Teuchos::rcp(new Teuchos::ArrayFileNameValidator(fileNameVali));
	TEST_ASSERT(arrayFileNameVali->getPrototype().get() == fileNameVali.get());
	Teuchos::Array<std::string> fileNameArray = Teuchos::tuple<std::string>("testFile.txt", "testFile2.txt", "testFile3.txt");
	Teuchos::Array<std::string> fileNameBadArray = Teuchos::tuple<std::string>("doesnexist.txt", "testFile2.txt", "testFile3.txt");
	TEST_NOTHROW(fileNameList->set("File name array", fileNameArray, "file name array parameter", arrayFileNameVali));
	TEST_THROW(fileNameList->set("Bad File name array", fileNameBadArray, "bad file name array parameter", arrayFileNameVali), Teuchos::Exceptions::InvalidParameterValue);
	TEST_THROW(stringList->set("Long array param", longArray, "long array parameter", arrayFileNameVali), Teuchos::Exceptions::InvalidParameterType);


	return (success ? 0:1);
}

int main(int argc, char* argv[]){
	bool success = true;
	Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
	if(testNumberValis(*out) == 1){
		success = false;
	}
	if(testStringVali(*out) == 1){
		success = false;
	}
	if(testFileNameVali(*out) == 1){
		success = false;
	}
	if(testArrayValis(*out) == 1){
		success = false;
	}
	return (success ? 0:1);
}

