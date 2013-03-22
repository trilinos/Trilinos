//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
  a number of required std::string parsing routines.

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

/*!  Get the next quoted std::string on the given line, starting at
  position pos. 

  \param line - Line of text from which to read

  \param pos - On input, the starting position in the line. On output,
  the next position after the quoted std::string (which may be
  std::string::npos). If there is any sort of error, this is set to
  std::string::npos upon return.

  \param value - On output, filled in with the quoted std::string (without
  the quotes). This is an empty std::string if no quoted std::string is found.

  \retval Returns true if the quoted std::string is successfully found,
  false otherwise.
  
*/
bool getNextQuotedString(const std::string& line, std::string::size_type& pos, std::string& value);

/*!  Get the next std::string on the given line, starting at
  position pos. 

  \param line - Line of text from which to read

  \param pos - On input, the starting position in the line. On output,
  the next position after the std::string (which may be
  std::string::npos). If there is any sort of error, this is set to
  std::string::npos upon return.

  \param value - On output, filled in with the next std::string (i.e., the
  next contguous block of non-space characters). This is an empty
  std::string if no std::string is found.

  \retval Returns true if the std::string is successfully found,
  false otherwise.
  
*/
bool getNextString(const std::string& line, std::string::size_type& pos, std::string& value);

/*!  Get the next std::string on the given line, starting at
  position pos, and convert it to a double. 

  \param line - Line of text from which to read

  \param pos - On input, the starting position in the line. On output,
  the next position after the std::string (which may be
  std::string::npos). If there is any sort of error in reading the
  next std::string, this is set to std::string::npos upon return.

  \param value - On output, filled in with the double value constained
  in the next std::string (i.e., the next contguous block of non-space
  characters).
  
  \retval Returns true if the next std::string contains a double,
  false otherwise.
  
*/
bool getNextDouble(const std::string& line, std::string::size_type& pos, double& value);


/*!  Get the next std::string on the given line, starting at
  position pos, and convert it to a int. 

  \param line - Line of text from which to read

  \param pos - On input, the starting position in the line. On output,
  the next position after the std::string (which may be
  std::string::npos). If there is any sort of error in reading the
  next std::string, this is set to std::string::npos upon return.

  \param value - On output, filled in with the int value constained
  in the next std::string (i.e., the next contguous block of non-space
  characters).
  
  \retval Returns true if the next std::string contains a int,
  false otherwise.
  
*/
bool getNextInt(const std::string& line, std::string::size_type& pos, int& value);


/*!
  \brief Parse a NOX input file and store the data in the given parameter list.

  \param filename - The file name.  See \ref pageExecutables_step2 for details on
  the structure of the input file.
  
  \param params - The parameter list that is to be filled in by this
  function

  \return Returns false if there are any problems parsing the input
  file, true otherwise.


*/


bool parseTextInputFile(const std::string filename, Teuchos::ParameterList& params);

//! Process a single line of the text input file.
bool processTextInputFileLine(const std::string& line, Teuchos::ParameterList& params, 
			      Teuchos::ParameterList*& subPtr);

} // namespace NOX

#endif
