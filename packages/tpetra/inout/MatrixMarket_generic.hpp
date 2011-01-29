//@HEADER
// ************************************************************************
// 
//               Tpetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
//@HEADER

#ifndef __MatrixMarket_generic_hpp
#define __MatrixMarket_generic_hpp

#include "Teuchos_ScalarTraits.hpp"
#ifdef HAVE_TEUCHOS_COMPLEX
#  include <complex>
#endif // HAVE_TEUCHOS_COMPLEX
#include <sstream>
#include <stdexcept>

namespace MatrixMarket {

  /// Maximum number of ASCII characters per line allowed in a Matrix
  /// Market file.
  int maxLineLength();

  //! True if line is a comment line, false otherwise.
  bool 
  checkCommentLine (const std::string& line, 
		    size_t& start, 
		    size_t& end, 
		    const size_t lineNumber,
		    const bool tolerant);

  template<class Ordinal>
  bool
  readPatternData (std::istream& istr, 
		   Ordinal& rowIndex, 
		   Ordinal& colIndex, 
		   const size_t lineNumber,
		   const bool tolerant)
  {
    Ordinal __rowIndex, __colIndex;

    if (istr.eof() || istr.fail())
      {
	if (tolerant)
	  return false;
	else
	  {
	    std::ostringstream os;
	    os << "Unable to read any data from line " << lineNumber << " of input";
	    throw std::invalid_argument(os.str());
	  }
      }
    istr >> __rowIndex;
    if (istr.fail())
      {
	if (tolerant)
	  return false;
	else
	  {
	    std::ostringstream os;
	    os << "Failed to get row index from line " << lineNumber << " of input";
	    throw std::invalid_argument(os.str());
	  }
      }
    else if (istr.eof())
      {
	if (tolerant)
	  return false;
	else
	  {
	    std::ostringstream os;
	    os << "No more data after row index on line " << lineNumber << " of input";
	    throw std::invalid_argument(os.str());
	  }
      }
    istr >> __colIndex;
    if (istr.fail())
      {
	if (tolerant)
	  return false;
	else
	  {
	    std::ostringstream os;
	    os << "Failed to get column index from line " << lineNumber << " of input";
	    throw std::invalid_argument(os.str());
	  }
      }
    rowIndex = __rowIndex;
    colIndex = __colIndex;
    return true;
  }

  template<class Ordinal, class Real>
  bool
  readRealData (std::istream& istr, 
		Ordinal& rowIndex, 
		Ordinal& colIndex, 
		Real& realValue, 
		const size_t lineNumber,
		const bool tolerant)
  {
    Real __realValue;
    if (! readPatternData (istr, rowIndex, colIndex, lineNumber, tolerant))
      {
	if (tolerant) 
	  return false;
	else
	  {
	    std::ostringstream os;
	    os << "Failed to read pattern data from line " << lineNumber << " of input";
	    throw std::invalid_argument(os.str());
	  }
      }
    if (istr.eof())
      {
	if (tolerant)
	  return false;
	else
	  {
	    std::ostringstream os;
	    os << "No more data after pattern data on line " << lineNumber << " of input";
	    throw std::invalid_argument(os.str());
	  }
      }
    istr >> __realValue;
    if (istr.fail())
      {
	if (tolerant)
	  return false;
	else
	  {
	    std::ostringstream os;
	    os << "Failed to get real value from line " << lineNumber << " of input";
	    throw std::invalid_argument(os.str());
	  }
      }
    realValue = __realValue;
    return true;
  }

#ifdef HAVE_TEUCHOS_COMPLEX
  template<class Ordinal, class Complex>
  bool
  readComplexData (std::istream& istr, 
		   Ordinal& rowIndex, 
		   Ordinal& colIndex, 
		   Complex& value,
		   const size_t lineNumber,
		   const bool tolerant)
  {
    typedef typename Teuchos::ScalarTraits<Complex>::magnitudeType Real;

    Real realValue, imagValue;
    if (! readRealData (istr, rowIndex, colIndex, realValue, lineNumber, tolerant))
      {
	if (tolerant) 
	  return false;
	else
	  {
	    std::ostringstream os;
	    os << "Failed to read pattern data and/or real value from line " << lineNumber << " of input";
	    throw std::invalid_argument(os.str());
	  }
      }
    if (istr.eof())
      {
	if (tolerant)
	  return false;
	else
	  {
	    std::ostringstream os;
	    os << "No more data after real value on line " << lineNumber << " of input";
	    throw std::invalid_argument(os.str());
	  }
      }
    istr >> imagValue;
    if (istr.fail())
      {
	if (tolerant)
	  return false;
	else
	  {
	    std::ostringstream os;
	    os << "Failed to get imaginary value from line " << lineNumber << " of input";
	    throw std::invalid_argument(os.str());
	  }
      }
    value = Complex (realValue, imagValue);
    return true;
  }
#endif // HAVE_TEUCHOS_COMPLEX

  template<class Ordinal>
  bool
  readPatternLine (const std::string& line, 
		   Ordinal& rowIndex, 
		   Ordinal& colIndex, 
		   const size_t lineNumber,
		   const bool tolerant)
  {
    size_t start, end;
    if (checkCommentLine (line, start, end, lineNumber, tolerant)) 
      return false; // It's a comment line
    // If it's an empty line, checkCommentLine() will throw an
    // exception if non-tolerant parsing is being performed, so
    // we need only return false otherwise.
    if (end == 0) 
      {
	if (tolerant)
	  throw std::logic_error("Should never get here! checkCommentLine() "
				 "is supposed to catch empty lines.");
	else
	  return false;
      }
    // The part of the line that contains data
    std::istringstream istr (line.substr (start, end));
    return readPatternLineData (istr, rowIndex, colIndex, lineNumber, tolerant);
  }

  template<class Ordinal, class Real>
  bool
  readRealLine (const std::string& line, 
		Ordinal& rowIndex, 
		Ordinal& colIndex, 
		Real& realValue, 
		const size_t lineNumber,
		const bool tolerant)
  {
    size_t start, end;
    if (checkCommentLine (line, start, end, lineNumber, tolerant)) 
      return false; // It's a comment line
    // If it's an empty line, checkCommentLine() will throw an
    // exception if non-tolerant parsing is being performed, so
    // we need only return false otherwise.
    if (end == 0) 
      {
	if (tolerant)
	  throw std::logic_error("Should never get here! checkCommentLine() "
				 "is supposed to catch empty lines.");
	else
	  return false;
      }
    // The part of the line that contains data
    std::istringstream istr (line.substr (start, end));
    return readRealLineData (istr, rowIndex, colIndex, realValue, lineNumber, tolerant);
  }

#ifdef HAVE_TEUCHOS_COMPLEX
  template<class Ordinal, class Complex>
  bool
  readComplexLine (const std::string& line, 
		   Ordinal& rowIndex, 
		   Ordinal& colIndex,
		   Complex& value,
		   const size_t lineNumber,
		   const bool tolerant)
  {
    size_t start, end;
    if (checkCommentLine (line, start, end, lineNumber, tolerant)) 
      return false; // It's a comment line
    // If it's an empty line, checkCommentLine() will throw an
    // exception if non-tolerant parsing is being performed, so
    // we need only return false otherwise.
    if (end == 0) 
      {
	if (tolerant)
	  throw std::logic_error("Should never get here! checkCommentLine() "
				 "is supposed to catch empty lines.");
	else
	  return false;
      }
    // The part of the line that contains data
    std::istringstream istr (line.substr (start, end));
    return readComplexLineData (istr, rowIndex, colIndex, value, lineNumber, tolerant);
  }
#endif // HAVE_TEUCHOS_COMPLEX

  template<class Ordinal, class PatternCallback>
  std::pair<bool, std::vector<size_t> >
  readPatternCoordinateData (std::istream& in, 
			     PatternCallback add,
			     const size_t startingLineNumber, 
			     const bool tolerant)
  {
    std::string line;
    size_t lineNumber = startingLineNumber;
    bool anySucceeded = false;
    bool allSucceeded = true;
    std::vector<size_t> badLineNumbers; 
    size_t validDataLines = 0;
    while (in.getline(line)) // ???
      {
	Ordinal rowIndex, colIndex;
	const bool localSuccess = readPatternLine (line, rowIndex, colIndex, lineNumber++, tolerant);
	anySucceeded = anySucceeded || localSuccess;
	allSucceeded = allSucceeded && localSuccess;
	if (! localSuccess)
	  badLineNumbers.push_back (lineNumber);
	else
	  {
	    add (rowIndex, colIndex);
	    validDataLines++;
	  }
      }
    if (lineNumber == startingLineNumber)
      anySucceeded = true; // Trivially true
    
    return std::make_pair (allSucceeded, badLineNumbers);
  }

  template<class Ordinal, class Real, class RealCallback>
  std::pair<bool, std::vector<size_t> >
  readRealCoordinateData (std::istream& in, 
			  RealCallback add,
			  const size_t startingLineNumber, 
			  const bool tolerant)
  {
    std::string line;
    size_t lineNumber = startingLineNumber;
    bool anySucceeded = false;
    bool allSucceeded = true;
    std::vector<size_t> badLineNumbers; 
    size_t validDataLines = 0;
    while (in.getline(line)) // ???
      {
	Ordinal rowIndex, colIndex;
	Real value;
	const bool localSuccess = readRealLine (line, rowIndex, colIndex, value, lineNumber++, tolerant);
	anySucceeded = anySucceeded || localSuccess;
	allSucceeded = allSucceeded && localSuccess;
	if (! localSuccess)
	  badLineNumbers.push_back (lineNumber);
	else
	  {
	    add (rowIndex, colIndex, value);
	    validDataLines++;
	  }
      }
    if (lineNumber == startingLineNumber)
      anySucceeded = true; // Trivially true
    
    return std::make_pair (allSucceeded, badLineNumbers);
  }

#ifdef HAVE_TEUCHOS_COMPLEX
  template<class Ordinal, class Complex, class ComplexCallback>
  std::pair<bool, std::vector<size_t> >
  readComplexCoordinateData (std::istream& in, 
			     ComplexCallback add,
			     const size_t startingLineNumber, 
			     const bool tolerant)
  {
    std::string line;
    size_t lineNumber = startingLineNumber;
    bool anySucceeded = false;
    bool allSucceeded = true;
    std::vector<size_t> badLineNumbers; 
    size_t validDataLines = 0;
    while (in.getline(line)) // ???
      {
	Ordinal rowIndex, colIndex;
	Complex complexValue;
	const bool localSuccess = readComplexLine (line, rowIndex, colIndex, complexValue, lineNumber++, tolerant);
	anySucceeded = anySucceeded || localSuccess;
	allSucceeded = allSucceeded && localSuccess;
	if (! localSuccess)
	  badLineNumbers.push_back (lineNumber);
	else
	  {
	    add (rowIndex, colIndex, complexValue);
	    validDataLines++;
	  }
      }
    if (lineNumber == startingLineNumber)
      anySucceeded = true; // Trivially true
    
    return std::make_pair (allSucceeded, badLineNumbers);
  }
#endif // HAVE_TEUCHOS_COMPLEX

  template<class Ordinal>
  bool
  readCoordinateDimensions (std::istream& istr, 
			    Ordinal& numRows, 
			    Ordinal& numCols,
			    Ordinal& numNonzeros,
			    const size_t lineNumber,
			    const bool tolerant)
  {
    Ordinal __numRows, __numCols, __numNonzeros;

    if (istr.eof() || istr.fail())
      {
	std::ostringstream os;
	os << "Unable to read any data from line " << lineNumber << " of input";
	throw std::invalid_argument(os.str());
      }
    istr >> __numRows;
    if (istr.fail())
      {
	std::ostringstream os;
	os << "Failed to get number of rows from line " << lineNumber << " of input";
	throw std::invalid_argument(os.str());
      }
    else if (istr.eof())
      {
	std::ostringstream os;
	os << "No more data after number of rows on line " << lineNumber << " of input";
	throw std::invalid_argument(os.str());
      }
    istr >> __numCols;
    if (istr.fail())
      {
	std::ostringstream os;
	os << "Failed to get number of columns from line " << lineNumber << " of input";
	throw std::invalid_argument(os.str());
      }
    else if (istr.eof())
      {
	std::ostringstream os;
	os << "No more data after number of columns on line " << lineNumber << " of input";
	throw std::invalid_argument(os.str());
      }
    istr >> __numNonzeros;
    if (istr.fail())
      {
	std::ostringstream os;
	os << "Failed to get number of nonzeros from line " << lineNumber << " of input";
	throw std::invalid_argument(os.str());
      }
    numRows = __numRows;
    numCols = __numCols;
    numNonzeros = __numNonzeros;
    return true;
  }

} // namespace MatrixMarket

#endif // __MatrixMarket_generic_hpp
