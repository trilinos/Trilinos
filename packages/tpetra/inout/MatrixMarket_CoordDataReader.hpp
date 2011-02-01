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

#ifndef __MatrixMarket_CoordDataReader_hpp
#define __MatrixMarket_CoordDataReader_hpp

#include "MatrixMarket_generic.hpp"

namespace MatrixMarket {

  template<class Callback, class Ordinal, class Scalar, bool isComplex = Teuchos::ScalarTraits<Scalar>::isComplex>
  class CoordDataReaderBase {
  private:
    Callback& adder_;

  public:
    CoordDataReaderBase (Callback& adder) : adder_ (adder) {}

    /// Read in the data from a single line of the file.
    ///
    /// This method has different implementation, depending on whether
    /// Scalar is complex or not.
    virtual bool
    readLine (const std::string& theLine, 
	      Ordinal& rowIndex, 
	      Ordinal& colIndex, 
	      Scalar& value, 
	      const size_t lineNumber,
	      const bool tolerant) = 0;

    std::pair<bool, std::vector<size_t> >
    read (std::istream& in, 
	  const size_t startingLineNumber, 
	  const bool tolerant,
	  const bool debug = false)
    {
      using std::cerr;
      using std::endl;
      typedef Teuchos::ScalarTraits<Scalar> STS;

      if (isComplex != STS::isComplex)
	throw std::logic_error("Should never get here!");
      else if (! in)
	throw std::invalid_argument("Input stream is invalid");

      std::string line;
      size_t lineNumber = startingLineNumber;
      bool allSucceeded = true;
      std::vector<size_t> badLineNumbers; 
      size_t validDataLines = 0;
      while (getline (in, line))
	{
	  size_t start, size;
	  if (checkCommentLine (line, start, size, lineNumber, tolerant))
	    {
	      // if (debug)
	      // 	cerr << "Comment line: " << lineNumber << endl;
	      lineNumber++;
	      continue; // it's a comment line
	    }
	  const std::string theLine = line.substr (start, size);
	  // if (debug)
	  //   cerr << "Possible data line " << lineNumber << ": " << line << endl;

	  Ordinal rowIndex, colIndex;
	  Scalar value;
	  const bool localSuccess = readLine (theLine, rowIndex, colIndex, 
					      value, lineNumber, tolerant);
	  lineNumber++;
	  allSucceeded = allSucceeded && localSuccess;
	  if (! localSuccess)
	    badLineNumbers.push_back (lineNumber);
	  else
	    {
	      adder_ (rowIndex, colIndex, value);
	      validDataLines++;
	    }
	}
      return std::make_pair (allSucceeded, badLineNumbers);
    }

    bool
    readDimensions (std::istream& in, 
		    Ordinal& numRows, 
		    Ordinal& numCols,
		    Ordinal& numNonzeros,
		    const size_t lineNumber,
		    const bool tolerant = false)
    {
      Ordinal __numRows, __numCols, __numNonzeros;
      std::string line;

      if (! getline(in, line))
	{
	  if (tolerant)
	    return false;
	  std::ostringstream os;
	  os << "Failed to read line " << lineNumber << " from input stream"
	    "; it should contains coordinate matrix dimensions.";
	  throw std::invalid_argument (os.str());
	}
      // Read in <numRows> <numCols> <numNonzeros> from input line
      std::istringstream istr (line);

      if (istr.eof() || istr.fail())
	{
	  if (tolerant)
	    return false;
	  std::ostringstream os;
	  os << "Unable to read any data from line " << lineNumber 
	     << " of input; it should contain coordinate matrix dimensions.";
	  throw std::invalid_argument(os.str());
	}
      istr >> __numRows;
      if (istr.fail())
	{
	  if (tolerant)
	    return false;
	  std::ostringstream os;
	  os << "Failed to get number of rows from line " << lineNumber << " of input";
	  throw std::invalid_argument(os.str());
	}
      else if (istr.eof())
	{
	  if (tolerant)
	    return false;
	  std::ostringstream os;
	  os << "No more data after number of rows on line " << lineNumber << " of input";
	  throw std::invalid_argument(os.str());
	}
      istr >> __numCols;
      if (istr.fail())
	{
	  if (tolerant)
	    return false;
	  std::ostringstream os;
	  os << "Failed to get number of columns from line " << lineNumber << " of input";
	  throw std::invalid_argument(os.str());
	}
      else if (istr.eof())
	{
	  if (tolerant)
	    return false;
	  std::ostringstream os;
	  os << "No more data after number of columns on line " << lineNumber << " of input";
	  throw std::invalid_argument(os.str());
	}
      istr >> __numNonzeros;
      if (istr.fail())
	{
	  if (tolerant)
	    return false;
	  std::ostringstream os;
	  os << "Failed to get number of nonzeros from line " << lineNumber << " of input";
	  throw std::invalid_argument(os.str());
	}
      // It would be nice to validate the read-in data further.  The
      // only thing we can do now is test if it's negative.  However,
      // we don't know syntactically whether Ordinal is a signed or
      // unsigned type.
      numRows = __numRows;
      numCols = __numCols;
      numNonzeros = __numNonzeros;
      return true;
    }
  };


  template<class Callback, class Ordinal, class Scalar, bool isComplex = Teuchos::ScalarTraits<Scalar>::isComplex>
  class CoordDataReader :
    public CoordDataReaderBase<Callback, Ordinal, Scalar, isComplex>
  {
  public:
    CoordDataReader (Callback& adder);

    bool
    readLine (const std::string& theLine, 
	      Ordinal& rowIndex, 
	      Ordinal& colIndex, 
	      Scalar& value, 
	      const size_t lineNumber,
	      const bool tolerant);
  };

#ifdef HAVE_TEUCHOS_COMPLEX
  template<class Callback, class Ordinal, class Scalar>
  class CoordDataReader<Callback, Ordinal, Scalar, true> : 
    public CoordDataReaderBase<Callback, Ordinal, Scalar, true>
  {
  public:
    CoordDataReader (Callback& adder) :
      CoordDataReaderBase<Callback, Ordinal, Scalar, true>(adder) 
    {}

    bool
    readLine (const std::string& theLine, 
	      Ordinal& rowIndex, 
	      Ordinal& colIndex, 
	      Scalar& value, 
	      const size_t lineNumber,
	      const bool tolerant)
    {
      typedef Teuchos::ScalarTraits<Scalar> STS;
      typedef typename STS::magnitudeType Real;
      Real realPart, imagPart;
      const bool localSuccess = 
	readComplexLine (theLine, rowIndex, colIndex, realPart, imagPart,
			 lineNumber, tolerant);
      if (localSuccess)
	{
	  // Assume that assignment from std::complex<Real> to
	  // Scalar (which itself is complex-valued) is valid.
	  // We have to do this, since the C++ compiler may not
	  // be smart enough to assume here (when it
	  // instantiates the templates) that Scalar is an
	  // std::complex<Real> -- even though it has to be, if
	  // STS::isComplex is true (which as of 31 Jan 2011,
	  // only holds for std::complex<T>).
	  const std::complex<Real> __value (realPart, imagPart);
	  value = __value;
	}
      return localSuccess;
    }
  };
#endif // HAVE_TEUCHOS_COMPLEX

  template<class Callback, class Ordinal, class Scalar>
  class CoordDataReader<Callback, Ordinal, Scalar, false> : 
    public CoordDataReaderBase<Callback, Ordinal, Scalar, false>
  {
  public:
    CoordDataReader (Callback& adder) :
      CoordDataReaderBase<Callback, Ordinal, Scalar, false>(adder) 
    {}

    bool
    readLine (const std::string& theLine, 
	      Ordinal& rowIndex, 
	      Ordinal& colIndex, 
	      Scalar& value, 
	      const size_t lineNumber,
	      const bool tolerant)
    {
      return readRealLine (theLine, rowIndex, colIndex, value, 
			   lineNumber, tolerant);
    }
  };

} // namespace MatrixMarket
#endif // __MatrixMarket_CoordDataReader_hpp
