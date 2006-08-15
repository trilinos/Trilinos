//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
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
  file, but may be reinstated in the future.
*/

#ifndef NOX_TEST_UTILS_H
#define NOX_TEST_UTILS_H

#include "NOX_Common.H"
#include "Teuchos_ParameterList.hpp"


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
bool getNextQuotedString(const string& line, string::size_type& pos, string& value);

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
bool getNextString(const string& line, string::size_type& pos, string& value);

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
bool getNextDouble(const string& line, string::size_type& pos, double& value);


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
bool getNextInt(const string& line, string::size_type& pos, int& value);


/*!
  \brief Parse a NOX input file and store the data in the given parameter list.

  \param filename - The file name.  See \ref pageExecutables_step2 for details on
  the structure of the input file.
  
  \param params - The parameter list that is to be filled in by this
  function

  \return Returns false if there are any problems parsing the input
  file, true otherwise.


*/


bool parseTextInputFile(const string filename, Teuchos::ParameterList& params);

//! Process a single line of the text input file.
bool processTextInputFileLine(const string& line, Teuchos::ParameterList& params, 
			      Teuchos::ParameterList*& subPtr);

} // namespace NOX

#endif
