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

#ifndef __MatrixMarket_raw_hpp
#define __MatrixMarket_raw_hpp

#include "MatrixMarket_util.hpp"
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>

namespace MatrixMarket {
  namespace Raw {

    /// \class Element
    /// \author Mark Hoemmen
    /// \brief One structural nonzero of a sparse matrix
    ///
    template<class Scalar, class Ordinal>
    class Element {
    public:
      Element (const Ordinal i, const Ordinal j, const Scalar& Aij) :
	rowIndex_ (i), colIndex_ (j), value_ (Aij) {}

      //! Ignore the nonzero value for comparisons.
      bool operator== (const Element& rhs) {
	return rowIndex_ == rhs.rowIndex_ && colIndex_ == rhs.colIndex_;
      }

      //! Ignore the nonzero value for comparisons.
      bool operator!= (const Element& rhs) {
	return ! (*this == rhs);
      }

      //! Lex order first by row index, then by column index.
      bool operator< (const Element& rhs) const {
	if (rowIndex_ < rhs.rowIndex_)
	  return true;
	else if (rowIndex_ > rhs.rowIndex_)
	  return false;
	else { // equal
	  return colIndex_ < rhs.colIndex_;
	}
      }

      void merge (const Element& rhs, const bool replace=false) {
	if (*this != rhs)
	  throw std::logic_error("Can only merge elements at the same "
				 "location in the sparse matrix");
	else if (replace)
	  value_ = rhs.value_;
	else
	  value_ += rhs.value_;
      }

      Ordinal rowIndex() const { return rowIndex_; }
      Ordinal colIndex() const { return colIndex_; }
      Scalar value() const { return value_; }

    private:
      Ordinal rowIndex_, colIndex_;
      Scalar value_;
    };

    //! Print out an Element to the given output stream
    std::ostream& 
    operator<< (std::ostream& out, const Element<Scalar, Ordinal>& elt) {
      out << elt.rowIndex() << " " << elt.colIndex() << " " << elt.value();
    }

    template<class Scalar, class Ordinal>
    class Adder {
    public:
      typedef Ordinal index_type;
      typedef Scalar value_type;
      typedef Element<Scalar, Ordinal> element_type;

      Adder () : numRows_(0), numCols_(0), numNonzeros_(0) {}

      //! Add an element to the sparse matrix at location (i,j) (one-based indexing).
      void operator() (const Ordinal i, const Ordinal j, const Scalar& Aij) {
	// i and j are 1-based
	elts_.push_back (Element (i-1, j-1, Aij));
	// Keep track of the rightmost column containing a nonzero,
	// and the bottommost row containing a nonzero.  This gives us
	// a lower bound for the dimensions of the matrix, and a check
	// for the reported dimensions of the matrix in the Matrix
	// Market file.
	numRows_ = std::max(numRows_, i);
	numCols_ = std::max(numCols_, j);
	numNonzeros_++;
      }

      /// \brief Print the sparse matrix data.  
      ///
      /// We always print the data sorted.  You may also merge
      /// duplicate entries if you prefer.
      /// 
      /// \param out [out] Output stream to which to print
      /// \param doMerge [in] Whether to merge entries before printing
      /// \param replace [in] If merging, whether to replace duplicate
      ///   entries; otherwise their values are added together.
      void print (std::ostream& out, const bool doMerge, const bool replace=false) {
	if (doMerge)
	  merge (replace);
	else
	  std::sort (elts_.begin(), elts_.end());
	// Print out the results, delimited by newlines.
	std::copy (elts_.begin(), elts_.end(), std::ostream_iterator (out, "\n"));
      }

      /// \brief Merge duplicate elements 
      ///
      /// Merge duplicate elements of the sparse matrix, where
      /// "duplicate" means at the same (i,j) location in the sparse
      /// matrix.  Resize the array of elements to fit just the
      /// "unique" (meaning "nonduplicate") elements.
      ///
      /// \return (# unique elements, # removed elements)
      std::pair<std::vector<element_type>::size_type, std::vector<element_type>::size_type>
      merge (const bool replace=false) 
      {
	typedef std::vector<element_type>::iterator iter_type;
	typedef std::vector<element_type>::size_type size_type;

	// Start with a sorted container.  It may be sorted already,
	// but we just do the extra work.
	std::sort (elts_.begin(), elts_.end());

	// Walk through the array in place, merging duplicates and
	// pushing unique elements up to the front of the array.  We
	// can't use std::unique for this because it doesn't let us
	// merge duplicate elements; it only removes them from the
	// sequence.
	size_type numUnique = 0;
	iter_type cur = elts_.begin();
	if (cur == elts_.end())
	  // There are no elements to merge
	  return std::make_pair (numUnique, size_type(0));
	else {
	  iter_type next = cur;
	  ++next; // There is one unique element
	  ++numUnique;
	  while (next != elts_.end()) {
	    if (*cur == *next) {
	      // Merge in the duplicated element *next
	      cur->merge (*next, replace);
	    } else {
	      // *cur is already a unique element.  Move over one to
	      // *make space for the new unique element.
	      ++cur; 
	      *cur = *next; // Add the new unique element
	      ++numUnique;
	    }
	    // Look at the "next" not-yet-considered element
	    ++next;
	  }
	  // Remember how many elements we removed before resizing.
	  const size_type numRemoved = elts_.size() - numUnique;
	  elts_.resize (numUnique);
	  return std::make_pair (numUnique, numRemoved);
	}
      }

    private:
      Ordinal numRows_, numCols_, numNonzeros_;
      std::vector<element_type> elts_;
    };


    template<class Scalar, class Ordinal>
    class Reader {
    public:
      static void
      readFile (const std::string& filename,
		const bool tolerant=false)
      {
	std::ifstream in (filename.c_str());
	return read (in, tolerant);
      }
      
      static void
      read (const std::istream& in,	
	    const bool tolerant=false)
      {
	typedef Teuchos::ScalarTraits<Scalar> STS;

	const Ordinal numRows, numCols, numNonzeros;
	bool success = true;
	std::string line;
	if (! in.getline(line))
	  throw std::invalid_argument ("Failed to get first (banner) line");

	Banner banner (line, tolerant);
	if (matrixType() != "coordinate")
	  throw std::invalid_argument ("Matrix Market input file must contain a "
				       "\"coordinate\"-format sparse matrix in "
				       "order to create a sparse matrix object "
				       "from it.");
	else if (! STS::isComplex && banner.dataType() == "complex")
	  throw std::invalid_argument ("Matrix Market file contains complex-"
				       "valued data, but your chosen Scalar "
				       "type is real.");
	else if (banner.dataType() != "real" && banner.dataType() != "complex")
	  throw std::invalid_argument ("Only real or complex data types (no "
				       "pattern or integer matrices) are "
				       "currently supported");
	// The rest of the file starts at line 2, after the banner line.
	size_t lineNumber = 2;
	// Read in the dimensions of the sparse matrix and the number
	// of nonzeros.
	success = readCoordinateDimensions (in, numRows, numCols, numNonzeros, 
					    lineNumber, tolerant);
	TEST_FOR_EXCEPTION(! success, std::invalid_argument,
			   "Error reading Matrix Market sparse matrix "
			   "file: failed to read coordinate dimensions.");
	// Read the sparse matrix entries
	std::pair<bool, std::vector<size_t> > results;	    
	typedef Adder<Scalar, Ordinal> raw_adder_type;
	typedef SymmetrizingAdder<raw_adder_type> adder_type;
	raw_adder_type rawAdder;
	adder_type adder (rawAdder, banner.symmType());

	if (banner.dataType() == "real")
	  results = readRealCoordinateData (in, adder, lineNumber, tolerant);
	else if (banner.dataType() == "complex")
	  {
	    TEST_FOR_EXCEPTION(! STS::isComplex, std::invalid_argument,
			       "The Matrix Market sparse matrix file contains "
			       "complex-valued data, but you are trying to read"
			       " the data into a sparse matrix of real values.");
	    results = readComplexCoordinateData (in, adder, lineNumber, tolerant);
	  }
	else
	  throw std::logic_error ("Should never get here!");

	// In tolerant mode, report any bad line number(s)
	if (tolerant && ! results.first)
	  reportBadness (std::cerr, results);

	// We're done reading in the sparse matrix.  Now print out the
	// nonzero entry/ies.
	rawAdder.print (std::cout);
	std::cout << std::endl;
      }

      static void 
      reportBadness (std::ostream& out, 
		     const std::pair<bool, std::vector<size_t> >& results) 
      {
	using std::endl;
	const size_t numErrors = results.second.size();
	const size_t maxNumErrorsToReport = 100;
	out << numErrors << " errors when reading Matrix Market sparse matrix file." << endl;
	if (numErrors > maxNumErrorsToReport)
	  out << "-- We do not report individual errors when there "
	    "are more than " << maxNumErrorsToReport << ".";
	else if (numErrors == 1)
	  out << "Error on line " << results.second[0] << endl;
	else if (numErrors > 1)
	  {
	    out << "Errors on lines {";
	    for (size_t k = 0; k < numErrors-1; ++k)
	      out << results.second[k] << ", ";
	    out << results.second[numErrors-1] << "}" << endl;
	  }
      }
    };



  } // namespace Raw
} // namespace MatrixMarket


#endif // __MatrixMarket_raw_hpp
