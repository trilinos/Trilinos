// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

#ifndef __Teuchos_MatrixMarket_Raw_Writer_hpp
#define __Teuchos_MatrixMarket_Raw_Writer_hpp

#include "Teuchos_MatrixMarket_SetScientific.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_ParameterList.hpp"
#include <fstream>
#include <iostream>

namespace Teuchos {
  namespace MatrixMarket {
    namespace Raw {
      /// \class Writer
      /// \brief Write a sparse matrix from raw CSR (compressed sparse
      ///   row) storage to a Matrix Market file.
      ///
      /// \tparam ScalarType The type of entries of the sparse matrix.
      /// \tparam OrdinalType The type of indices of the sparse matrix.
      ///
      /// This class is useful for testing local sparse kernels.  It
      /// should only be called by one MPI process at a time and is
      /// not aware of parallel communication.  Use
      /// Tpetra::MatrixMarket::Writer if you want to write a
      /// Tpetra::CrsMatrix to a Matrix Market file.
      template<class ScalarType, class OrdinalType>
      class Writer {
      public:
        /// \brief Write the sparse matrix to the given file.
        ///
        /// The input arrays rowptr, colind, values together form the
        /// common three-arrays representation of compressed sparse
        /// row (CSR) storage.
        ///
        /// \param filename [in] Name of the Matrix Market file to
        ///   which to write the sparse matrix.
        /// \param rowptr [in] Array of numRows+1 offsets, where
        ///   numRows is the number of rows in the sparse matrix.  For
        ///   row i (zero-based indexing), the entries of that row are
        ///   in indices rowptr[i] .. rowptr[i+1]-1 of colind and
        ///   values.
        /// \param colind [in] Column indices of the matrix.  Same
        ///   number of entries as values.  colind[k] is the column
        ///   index of values[k].
        /// \param values [in] Values stored in the matrix.
        /// \param numRows [in] Number of rows in the sparse matrix.
        ///   This is redundant, because rowptr.size() == numRows on
        ///   output.
        /// \param numCols [in] Number of columns in the sparse matrix.
	void
        writeFile (const std::string& filename,
		   const ArrayView<const OrdinalType>& rowptr,
		   const ArrayView<const OrdinalType>& colind,
		   const ArrayView<const ScalarType>& values,
		   const OrdinalType numRows,
		   const OrdinalType numCols)
        {
          std::ofstream out (filename.c_str ());
          TEUCHOS_TEST_FOR_EXCEPTION(! out, std::runtime_error,
            "Failed to open file \"" << filename << "\" for writing.");
          write (out, rowptr, colind, values, numRows, numCols);
        }

        /// \brief Write the sparse matrix to the given output stream.
        ///
        /// The input arrays rowptr, colind, values together form the
        /// common three-arrays representation of compressed sparse
        /// row (CSR) storage.
        ///
        /// \param out [in/out] Output stream to which to write the 
        ///   sparse matrix.
        /// \param rowptr [in] Array of numRows+1 offsets, where
        ///   numRows is the number of rows in the sparse matrix.  For
        ///   row i (zero-based indexing), the entries of that row are
        ///   in indices rowptr[i] .. rowptr[i+1]-1 of colind and
        ///   values.
        /// \param colind [in] Column indices of the matrix.  Same
        ///   number of entries as values.  colind[k] is the column
        ///   index of values[k].
        /// \param values [in] Values stored in the matrix.
        /// \param numRows [in] Number of rows in the sparse matrix.
        ///   This is redundant, because rowptr.size() == numRows on
        ///   output.
        /// \param numCols [in] Number of columns in the sparse matrix.
        ///
        /// \return If parsing tolerantly: false if the Matrix Market
        ///   file has any syntax errors, else true.  If not parsing
        ///   tolerantly, this should always return true.
	void
        write (std::ostream& out,
	       const ArrayView<const OrdinalType>& rowptr,
	       const ArrayView<const OrdinalType>& colind,
	       const ArrayView<const ScalarType>& values,
	       const OrdinalType numRows,
	       const OrdinalType numCols)
        {
          using std::endl;
          typedef ScalarTraits<ScalarType> STS;
	  typedef typename ArrayView<const OrdinalType>::size_type size_type;

	  // Make the output stream write floating-point numbers in
	  // scientific notation.  It will politely put the output
	  // stream back to its state on input, when this scope
	  // terminates.
	  Teuchos::MatrixMarket::details::SetScientific<ScalarType> sci (out);

	  // Data type string for ScalarType.
	  std::string dataType;
	  if (STS::isComplex) {
	    dataType = "complex"; 
	  } else if (STS::isOrdinal) {
	    dataType = "integer";
	  } else {
	    dataType = "real";
	  }

          // Print the Matrix Market banner line.  We assume
          // nonsymmetric storage ("general").
          out << "%%MatrixMarket matrix coordinate " << dataType << " general"
	      << endl;

          // // Print comments (the matrix name and / or description).
          // if (matrixName != "") {
          //   printAsComment (out, matrixName);
	  // }
          // if (matrixDescription != "") {
          //   printAsComment (out, matrixDescription);
	  // }

          // Write the dimensions of the sparse matrix: (# rows, #
          // columns, # matrix entries (counting duplicates as
          // separate entries)).  
	  out << numRows << " " << numCols << " " << rowptr[numRows] << endl;

	  for (size_type i = 0; i < numRows; ++i) {
	    for (OrdinalType k = rowptr[i]; k < rowptr[i+1]; ++k) {
	      const OrdinalType j = colind[k];
	      const ScalarType& A_ij = values[k];

	      // Matrix Market files use 1-based row and column indices.
	      out << (i+1) << " " << (j+1) << " ";
	      if (STS::isComplex) {
		out << STS::real (A_ij) << " " << STS::imag (A_ij);
	      } else {
		out << A_ij;
	      }
	      out << endl;
	    }
	  }
        }
      }; // end of class Writer
    } // namespace Raw
  } // namespace MatrixMarket
} // namespace Teuchos

#endif // __Teuchos_MatrixMarket_Raw_Writer_hpp
