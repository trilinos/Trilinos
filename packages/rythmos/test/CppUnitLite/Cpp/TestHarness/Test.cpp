

#include "Test.h"
#include "TestRegistry.h"
#include "TestResult.h"



Test::Test (const std::string& testName) : name (testName) 
{
	TestRegistry::addTest (this);
}


void Test::run (TestResult& result) 
{
	runTest (result);
	result.testWasRun();
}

