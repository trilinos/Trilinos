

#ifndef TESTRESULT_H
#define TESTRESULT_H

// TestResult is a simple little class which accepts notifications from
// tests.  In its current form, it takes these notifications and prints
// them on the standard output.
//
// If you need different behavior, you can alter it in place or subclass 
// TestResult to provide it.  The TestResult hierarchy can be as simple
// or complex as you require to support your application.


class Failure;

class TestResult
{
public:
		 TestResult () : failureCount (0), testCount(0) {}
	void testWasRun ();
	void startTests ();
	void addFailure (Failure failure);
	void endTests ();

private:
	int failureCount;
	int testCount;
};

#endif
