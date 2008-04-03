
#ifndef TESTREGISTRY_H
#define TESTREGISTRY_H

// TestRegistry is a primitive singleton which collects all of the
// tests in a system and allows them to be executed.  To see how
// tests are added to the TestRegistry look at the Test.h file


#include <vector>

class Test;
class TestResult;

class TestRegistry
{
public:
	static void addTest (Test *test);
	static void runAllTests (TestResult& result);

private:

	static TestRegistry&	instance ();
	void					add (Test *test);
	void					run (TestResult& result);

	std::vector<Test *>		tests;

};



#endif