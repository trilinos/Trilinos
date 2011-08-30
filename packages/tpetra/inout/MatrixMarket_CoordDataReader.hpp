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

namespace Tpetra {
  namespace MatrixMarket {

    /// \class CoordDataReaderBase 
    /// \brief Coordinate-format sparse matrix data reader.
    ///
    /// This class is an interface and partial implementation of
    /// reading coordinate-format sparse matrix data from a Matrix
    /// Market file.  There are two concrete implementations: one for
    /// real-valued data, and the other for complex-valued data.
    template<class Callback, class Ordinal, class Scalar, bool isComplex = ScalarTraits<Scalar>::isComplex>
    class CoordDataReaderBase {
    private:
      //! Closure that knows how to add entries to the sparse matrix.
      RCP<Callback> adder_;

    public:
      /// \brief Constructor with "adder" argument.
      ///
      /// This is the favored way to construct an instance of this
      /// type.  Only use the no-argument constructor if you have a
      /// "chicken-and-egg" problem (where the adder needs the matrix
      /// dimensions).
      /// 
      /// \param adder [in/out] Closure whose operator() adds an entry
      ///   to the sparse matrix on each invocation.
      CoordDataReaderBase (const RCP<Callback>& adder) : 
	adder_ (adder) {}

      /// \brief No-argument constructor.
      ///
      /// We offer this option in case the adder's constructor needs
      /// the matrix dimensions, so that it's necessary to call
      /// readDimensions() first before constructing the adder.  You
      /// should call setAdder() with a non-null argument before
      /// calling read() or readLine().
      CoordDataReaderBase () : adder_ (null) {}

      //! Virtual destructor for safety and happy compilers.
      virtual ~CoordDataReaderBase () {}

      /// \brief Set the Adder object.
      /// 
      /// Please don't call this after calling read() or readLine().
      /// The right time to call this is right after calling the
      /// no-argument constructor, if it's not possible to supply an
      /// Adder object before calling readDimensions().
      void setAdder (const RCP<Callback>& adder) {
	adder_ = adder;
      }

      /// \brief Read in the data from a single line of the input stream.
      ///
      /// This method has a different implementation, depending on
      /// whether Scalar is complex or not.
      virtual bool
      readLine (const std::string& theLine, 
		Ordinal& rowIndex, 
		Ordinal& colIndex, 
		Scalar& value, 
		const size_t lineNumber,
		const bool tolerant) = 0;

      /// \brief Read in all the data from the given input stream.
      ///
      /// \param in [in/out] The input stream from which to read
      ///
      /// \param startingLineNumber [in] The line number of the file
      ///   from which we begin reading.  (This is used for
      ///   informative error output, if an error is detected in the
      ///   file.)
      ///
      /// \param tolerant [in] If true, parse tolerantly.  The
      ///   resulting read-in data may be incorrect, but the parser
      ///   won't throw an exception if it can just let bad data
      ///   through and continue.
      ///
      /// \param debug [in] If true, print verbose debugging output.
      ///
      std::pair<bool, std::vector<size_t> >
      read (std::istream& in, 
	    const size_t startingLineNumber, 
	    const bool tolerant,
	    const bool debug = false)
      {
	using std::cerr;
	using std::endl;
	typedef ScalarTraits<Scalar> STS;

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
		(*adder_) (rowIndex, colIndex, value);
		validDataLines++;
	      }
	  }
	return std::make_pair (allSucceeded, badLineNumbers);
      }

      /// \brief Read (numRows, numCols, numNonzeros).
      ///
      /// Read one line from the given input stream, and parse it into
      /// the number of rows, the number of columns, and the number of
      /// nonzeros in the sparse matrix.  We assume that those data
      /// are in whitespace-delimited format and are read from a
      /// Matrix Market - format file.
      ///
      /// \param in [in/out] The input stream from which to attempt to
      ///   read one line.
      /// 
      /// \param lineNumber The starting line number from which we
      ///   begin reading from the input stream.  Used for diagnostic
      ///   error output.
      ///
      /// \param tolerant [in] Whether to read "tolerantly" (setting
      ///   defaults and returning whether we were successful) or
      ///   "intolerantly" (throwing an exception on any deviation
      ///   from the expected format).
      ///
      /// \return ((numRows, numCols, numNonzeros), success).  In
      ///   tolerant mode, success may be false, meaning that the
      ///   read-in triple may not be valid.
      ///
      std::pair<Tuple<Ordinal, 3>, bool>
      readDimensions (std::istream& in, 
		      size_t& lineNumber,
		      const bool tolerant = false)
      {
	Tuple<Ordinal, 3> dims;
	// Fill in (numRows, numCols, numNonzeros) with reasonable
	// defaults.  If we don't succeed in reading all the data
	// from the current line of the input stream, the remaining
	// values not read will be these default values.
	dims[0] = 0;
	dims[1] = 0;
	dims[2] = 0;

	// Keep reading lines from the input stream until we find a
	// non-comment line, or until we run out of lines.  The latter
	// is an error, since every "coordinate" format Matrix Market
	// file must have a dimensions line after the banner (even if
	// the matrix has zero rows or columns, or zero entries).
	std::string line;
	bool commentLine = true;
	while (commentLine)
	  {
	    // Is it even valid to read from the input stream?
	    if (in.eof() || in.fail())
	      {
		if (tolerant)
		  return std::make_pair (dims, false);
		else
		  {
		    std::ostringstream os;
		    os << "Unable to get coordinate dimensions line (at all) "
		      "from (line " << lineNumber << ") of input stream; the "
		      "input stream claims that it is at \"end-of-file\" or has "
		      "an otherwise \"fail\"ed state.";
		    throw std::invalid_argument(os.str());
		  }
	      }
	    // Try to get the next line from the input stream.
	    if (getline(in, line))
	      lineNumber++; // We did actually read a line
	    else
	      {
		if (tolerant)
		  return std::make_pair (dims, false);
		else
		  {
		    std::ostringstream os;
		    os << "Failed to read coordinate dimensions line (at all) "
		      "from (line " << lineNumber << " from input stream.  The "
		      "line should contain the coordinate matrix dimensions in "
		       << " the form \"<numRows> <numCols> <numNonzeros>\".";
		    throw std::invalid_argument (os.str());
		  }
	      }
	    // Is the current line a comment line?  Ignore start and
	    // size; they are only useful for reading the actual matrix
	    // entries.  (We could use them here as an optimization, but
	    // we've chosen not to.)
	    size_t start = 0, size = 0;
	    commentLine = checkCommentLine (line, start, size, 
					    lineNumber, tolerant);
	  }
	//
	// Read in <numRows> <numCols> <numNonzeros> from input line
	//
	std::istringstream istr (line);
	// Does line contain anything at all?  Can we safely read from
	// the input stream wrapping the line?
	if (istr.eof() || istr.fail())
	  {
	    if (tolerant)
	      return std::make_pair (dims, false);
	    std::ostringstream os;
	    os << "Unable to read any data from line " << lineNumber 
	       << " of input; the line should contain the coordinate matrix "
	       << "dimensions \"<numRows> <numCols> <numNonzeros>\".";
	    throw std::invalid_argument(os.str());
	  }
	// Read in <numRows>
	{
	  Ordinal theNumRows = 0;
	  istr >> theNumRows;
	  if (istr.fail())
	    {
	      if (tolerant)
		return std::make_pair (dims, false);
	      std::ostringstream os;
	      os << "Failed to get number of rows from line " << lineNumber 
		 << " of input; the line should contain the coordinate matrix "
		 << " dimensions \"<numRows> <numCols> <numNonzeros>\".";
	      throw std::invalid_argument(os.str());
	    }
	  else // Capture the validly read result before checking for eof.
	    dims[0] = theNumRows;
	}
	// There should be two more things to read.
	if (istr.eof())
	  {
	    if (tolerant)
	      return std::make_pair (dims, false);
	    std::ostringstream os;
	    os << "No more data after number of rows on line " << lineNumber
	       << " of input; the line should contain the coordinate matrix "
	       << " dimensions \"<numRows> <numCols> <numNonzeros>\".";
	    throw std::invalid_argument(os.str());
	  }
	// Read in <numCols>
	{
	  Ordinal theNumCols = 0;
	  istr >> theNumCols;
	  if (istr.fail())
	    {
	      if (tolerant)
		return std::make_pair (dims, false);
	      std::ostringstream os;
	      os << "Failed to get number of columns from line " << lineNumber 
		 << " of input; the line should contain the coordinate matrix "
		 << " dimensions \"<numRows> <numCols> <numNonzeros>\".";
	      throw std::invalid_argument(os.str());
	    }
	  else // Capture the validly read result before checking for eof.
	    dims[1] = theNumCols;
	}
	// There should be one more thing to read.
	if (istr.eof())
	  {
	    if (tolerant)
	      return std::make_pair (dims, false);
	    std::ostringstream os;
	    os << "No more data after number of columns on line " << lineNumber
	       << " of input; the line should contain the coordinate matrix "
	       << " dimensions \"<numRows> <numCols> <numNonzeros>\".";
	    throw std::invalid_argument(os.str());
	  }
	// Read in <numNonzeros>
	{
	  Ordinal theNumNonzeros = 0;
	  istr >> theNumNonzeros;
	  if (istr.fail())
	    {
	      if (tolerant)
		return std::make_pair (dims, false);
	      std::ostringstream os;
	      os << "Failed to get number of (structural) nonzeros from line " 
		 << lineNumber 
		 << " of input; the line should contain the coordinate matrix "
		 << " dimensions \"<numRows> <numCols> <numNonzeros>\".";
	      throw std::invalid_argument(os.str());
	    }
	  else // Capture the validly read result
	    dims[2] = theNumNonzeros;
	}
	// It would be nice to validate the read-in data further.  The
	// only thing we can do now is test if it's negative.  However,
	// we don't know syntactically whether Ordinal is a signed or
	// unsigned type, so we shouldn't even test for negativity.
	return std::make_pair (dims, true);
      }
    };


    /// \class CoordDataReader
    /// \brief Coordinate-format sparse matrix data reader.
    ///
    /// This class completes the implementation of
    /// CoordDataReaderBase.  There are two concrete specializations:
    /// one for real-valued data, and the other for complex-valued
    /// data.
    template<class Callback, class Ordinal, class Scalar, bool isComplex = ScalarTraits<Scalar>::isComplex>
    class CoordDataReader :
      public CoordDataReaderBase<Callback, Ordinal, Scalar, isComplex>
    {
    public:
      CoordDataReader (const RCP<Callback>& adder);
      CoordDataReader ();

      //! Virtual destructor for safety and happy compilers.
      virtual ~CoordDataReader();

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
      CoordDataReader (const RCP<Callback>& adder) :
	CoordDataReaderBase<Callback, Ordinal, Scalar, true> (adder)
      {}

      CoordDataReader() : 
	CoordDataReaderBase<Callback, Ordinal, Scalar, true> (null)
      {}

      //! Virtual destructor for safety and happy compilers.
      virtual ~CoordDataReader() {};

      bool
      readLine (const std::string& theLine, 
		Ordinal& rowIndex, 
		Ordinal& colIndex, 
		Scalar& value, 
		const size_t lineNumber,
		const bool tolerant)
      {
	typedef ScalarTraits<Scalar> STS;
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
	    const std::complex<Real> theValue (realPart, imagPart);
	    value = theValue;
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
      CoordDataReader (const RCP<Callback>& adder) :
	CoordDataReaderBase<Callback, Ordinal, Scalar, false> (adder)
      {}

      CoordDataReader() : 
	CoordDataReaderBase<Callback, Ordinal, Scalar, false> (null)
      {}

      //! Virtual destructor for safety and happy compilers.
      virtual ~CoordDataReader() {};

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
} // namespace Tpetra

#endif // __MatrixMarket_CoordDataReader_hpp
