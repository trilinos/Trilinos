//@HEADER

/*!
  \file NOX_TestUtils.H
  \brief Header for the library that reads NOX input parameter
  information.

  This header defines the routine NOX::parseTextInputFile,
  used to read NOX parameter input file information. It also defines
  a number of required string parsing routines.

  These routines are based on a similar code in APPSPACK. Functionality
  for vector-typed input parameters has been commented out in this
  file, but may be reinsated in the future.
*/

#ifndef NOX_TEST_UTILS_H
#define NOX_TEST_UTILS_H

#include "NOX_Common.H"
#include "NOX_Parameter_List.H"


namespace NOX
{

/*!  Get the next quoted string on the given line, starting at
  position pos. 

  \param line - Line of text from which to read

  \param pos - On input, the starting position in the line. On output,
  the next position after the quoted string (which may be
  std::string::npos). If there is any sort of error, this is set to
  std::string::npos upon return.

  \param value - On output, filled in with the quoted string (without
  the quotes). This is an empty string if no quoted string is found.

  \retval Returns true if the quoted string is successfully found,
  false otherwise.
  
*/
bool NOX::getNextQuotedString(const string& line, string::size_type& pos, string& value);

/*!  Get the next string on the given line, starting at
  position pos. 

  \param line - Line of text from which to read

  \param pos - On input, the starting position in the line. On output,
  the next position after the string (which may be
  std::string::npos). If there is any sort of error, this is set to
  std::string::npos upon return.

  \param value - On output, filled in with the next string (i.e., the
  next contguous block of non-space characters). This is an empty
  string if no string is found.

  \retval Returns true if the string is successfully found,
  false otherwise.
  
*/
bool NOX::getNextString(const string& line, string::size_type& pos, string& value);

/*!  Get the next string on the given line, starting at
  position pos, and convert it to a double. 

  \param line - Line of text from which to read

  \param pos - On input, the starting position in the line. On output,
  the next position after the string (which may be
  std::string::npos). If there is any sort of error in reading the
  next string, this is set to std::string::npos upon return.

  \param value - On output, filled in with the double value constained
  in the next string (i.e., the next contguous block of non-space
  characters).
  
  \retval Returns true if the next string contains a double,
  false otherwise.
  
*/
bool NOX::getNextDouble(const string& line, string::size_type& pos, double& value);


/*!  Get the next string on the given line, starting at
  position pos, and convert it to a int. 

  \param line - Line of text from which to read

  \param pos - On input, the starting position in the line. On output,
  the next position after the string (which may be
  std::string::npos). If there is any sort of error in reading the
  next string, this is set to std::string::npos upon return.

  \param value - On output, filled in with the int value constained
  in the next string (i.e., the next contguous block of non-space
  characters).
  
  \retval Returns true if the next string contains a int,
  false otherwise.
  
*/
bool NOX::getNextInt(const string& line, string::size_type& pos, int& value);


/*!
  \brief Parse a NOX input file and store the data in the given parameter list.

  \param filename - The file name.  See \ref pageExecutables_step2 for details on
  the structure of the input file.
  
  \param params - The parameter list that is to be filled in by this
  function

  \return Returns false if there are any problems parsing the input
  file, true otherwise.


*/


bool NOX::parseTextInputFile(const string filename, NOX::Parameter::List& params);

//! Process a single line of the text input file.
bool NOX::processTextInputFileLine(const string& line, NOX::Parameter::List& params, 
			      NOX::Parameter::List*& subPtr);

} // namespace NOX

#endif
