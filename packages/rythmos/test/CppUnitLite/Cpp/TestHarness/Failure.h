

#ifndef FAILURE_H
#define FAILURE_H

// Failure records the circumstances surrounding a test failure.  Using C++
// macros were are able to record the name of the file where the failure 
// occurred, the line number, and the text of the condition which provoked
// the failure.

 
#include <string>

class Failure
{
public:
	Failure (std::string theCondition, std::string theTestName, std::string theFileName, long theLineNumber) 
		: condition (theCondition), testName (theTestName), fileName (theFileName), lineNumber (theLineNumber)
	{
	}

	std::string condition;
	std::string testName;
	std::string fileName;
	long lineNumber;
};


inline std::ostream& operator<< (std::ostream& stream, Failure& failure)
{
	stream 
		<< "Failure: \"" << failure.condition.c_str () << "\" " 
		<< "line " << failure.lineNumber << " in "
		<< failure.fileName.c_str ()
		<< std::endl;

	return stream;
}

#endif