// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teuchos_MatrixMarket_Raw_Writer_hpp
#define __Teuchos_MatrixMarket_Raw_Writer_hpp

#include <Teuchos_SetScientific.hpp>
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
	  Teuchos::SetScientific<ScalarType> sci (out);

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
