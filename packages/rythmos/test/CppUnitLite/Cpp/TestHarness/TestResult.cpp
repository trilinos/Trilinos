
#include "TestResult.h"
#include "Failure.h"

#include <iostream>



void TestResult::testWasRun()
{
	testCount++;
}

void TestResult::startTests () 
{
}

void TestResult::addFailure (Failure failure) 
{
	using namespace std;

	cout << failure << endl;
	failureCount++;
}

void TestResult::endTests () 
{
	using namespace std;

	cout << testCount << " tests run" << endl;
	if (failureCount > 0)
		cout << "There were " << failureCount << " failures" << endl;
	else
		cout << "There were no test failures" << endl;
}
