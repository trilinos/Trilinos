//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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

#include "Kokkos_Mkl_MatrixDescriptor.hpp"


namespace Kokkos {
  namespace Mkl {
    void MatrixDescriptor::print (Teuchos::FancyOStream& out) const
    {
      using std::endl;
      TEUCHOS_TEST_FOR_EXCEPTION(descr_.size() < 6, std::invalid_argument,
        "MKL documentation requires that the sparse matrix descriptor array "
        "have at least 6 elements.  MKL only uses the first 4 currently.");

      const char s0 = descr_[0];
      std::string structure;
      if (s0 == 'G') {
        structure = "General";
      } else if (s0 == 'S') {
        structure = "Symmetric";
      } else if (s0 == 'H') {
        structure = "Hermitian";
      } else if (s0 == 'T') {
        structure = "Triangular";
      } else if (s0 == 'A') {
        structure = "Antisymmetric";
      } else if (s0 == 'D') {
        structure = "Diagonal";
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
          "First element of the matrix descriptor array is invalid.  "
          "Valid values are 'G', 'S', 'H', 'T', 'A', and 'D', but you "
          "provided '" << s0 << "'.");
      }

      // Symmetric, Hermitian, and Triangular matrices all use
      // triangular storage, in the sense that MKL ignores all entries
      // in the lower resp. upper triangle of the matrix in these
      // cases.
      const char s1 = descr_[1];
      std::string upperLower;
      if (s1 == 'L') {
        upperLower = "Lower";
      } else if (s1 == 'U') {
        upperLower = "Upper";
      } else {
        upperLower = "Neither";
      }

      const char s2 = descr_[2];
      const bool unitDiag = (s2 == 'U');

      // 'C' means zero-based indexing (as in C); 'F' means one-based
      // indexing (as in Fortran).
      const char s3 = descr_[3];
      int indexBase = 0;
      if (s3 == 'C') {
        indexBase = 0;
      } else if (s3 == 'F') {
        indexBase = 1;
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(
          s3 != 'C' && s3 != 'F',
          std::invalid_argument,
          "Fourth element of matrix descriptor has invalid index base "
          "specification '" << s3 << "'.  Valid values are 'C' (for zero-based, "
          "C-style indices) and 'F' (for one-based, Fortran-style indices).");
      }

      // YAML (Yet Another Markup Language) output.
      // We defer all output until end of routine, for exception safety.
      out << "MKL sparse matrix descriptor:" << endl;
      Teuchos::OSTab tab (Teuchos::rcpFromRef (out));
      out << "Structure: \"" << structure << "\"" << endl
          << "UpperLower: \"" << upperLower << "\"" << endl
          << "UnitDiagonal: " << (unitDiag ? "true" : "false") << endl
          << "IndexBase: " << indexBase << endl;
    }

    MatrixDescriptor::
    MatrixDescriptor (const bool symmetric,
                      const bool hermitian,
                      const bool triangular,
                      const bool skewSymmetric,
                      const bool diagonal,
                      const Teuchos::EUplo uplo,
                      const Teuchos::EDiag diag,
                      const int indexBase) :
      descr_ ("XXXXXX") // initialize with correct length but invalid values
    {
      fill (descr_, symmetric, hermitian, triangular, skewSymmetric,
            diagonal, uplo, diag, indexBase);
    }

    MatrixDescriptor::MatrixDescriptor () :
      descr_ ("XXXXXX") // initialize with correct length but invalid values
    {
      const bool symmetric = false;
      const bool hermitian = false;
      const bool triangular = false;
      const bool skewSymmetric = false;
      const bool diagonal = false;
      const Teuchos::EUplo uplo = Teuchos::UNDEF_TRI;
      const Teuchos::EDiag diag = Teuchos::NON_UNIT_DIAG;
      // We use one-based indices by default, so that MKL's mat-vec
      // routines that take multiple vectors accept the vectors in
      // column-major order.  With zero-based indices, those routines
      // require the vectors in row-major order.
      const int indexBase = 1;
      fill (descr_, symmetric, hermitian, triangular, skewSymmetric,
            diagonal, uplo, diag, indexBase);
    }

    void
    MatrixDescriptor::fill (std::string& descr,
                            const bool symmetric,
                            const bool hermitian,
                            const bool triangular,
                            const bool skewSymmetric,
                            const bool diagonal,
                            const Teuchos::EUplo uplo,
                            const Teuchos::EDiag diag,
                            const int indexBase)
    {
      // MKL docs require 6 elements, of which only the first 4 are
      // used currently.
      if (descr.size() < 6) {
        descr.resize (6);
        // Just put something here that's different, as a marker.
        // MKL doesn't use it, so we can put whatever we want.
        descr[4] = 'X';
        descr[5] = 'X';
      }

      char s0 = 'G';
      // Diagonal test goes first, since a diagonal matrix may also be
      // symmetric, Hermitian, or skew-symmetric, depending on the
      // diagonal's values.  The user is responsible for verifying
      // symmetry, Hermitian-ness, or skew-symmetry of a diagonal
      // matrix.
      if (diagonal) {
        s0 = 'D';
      } else if (symmetric) {
        s0 = 'S';
      } else if (hermitian) {
        s0 = 'H';
      } else if (uplo != Teuchos::UNDEF_TRI) {
        s0 = 'T';
      } else if (skewSymmetric) {
        s0 = 'A'; // A for antisymmetric
      } else {
        s0 = 'G'; // G for general
      }

      // If the matrix is neither upper nor lower triangular, we put
      // 'N' (for "Neither").  Note that MKL's routines dealing with
      // triangular, (skew-)symmetric, or Hermitian matrices only read
      // the elements of the matrix in the lower resp. upper triangle,
      // even if storage includes both lower and upper triangles.
      char s1 = 'N';
      if (uplo == Teuchos::LOWER_TRI) {
        s1 = 'L';
      } else if (uplo == Teuchos::UPPER_TRI) {
        s1 = 'U';
      } else {
        s1 = 'N';
      }

      // Whether or not the matrix has unit diagonal.
      char s2 = 'N';
      if (diag == Teuchos::NON_UNIT_DIAG) {
        s2 = 'N';
      } else if (diag == Teuchos::UNIT_DIAG) {
        s2 = 'U';
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
          "Invalid Teuchos::EDiag enum value " << diag << ".  Valid values are "
          "Teuchos::NON_UNIT_DIAG and Teuchos::UNIT_DIAG.");
      }

      // 'C' means zero-based indexing (as in C); 'F' means one-based
      // indexing (as in Fortran).
      char s3 = 'C';
      if (indexBase == 0) {
        s3 = 'C';
      } else if (indexBase == 1) {
        s3 = 'F';
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(indexBase != 0 && indexBase != 1,
          std::invalid_argument, "Invalid index base " << indexBase
          << ".  Valid values are 0 (for zero-based indexing) and 1 "
          "(for one-based indexing).");
      }

      // Now that we've validated the input, commit changes.
      descr[0] = s0;
      descr[1] = s1;
      descr[2] = s2;
      descr[3] = s3;
    }

    int
    MatrixDescriptor::getIndexBase () const
    {
      if (descr_[3] == 'C') {
        return 0;
      } else if (descr_[3] == 'F') {
        return 1;
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::invalid_argument, "Invalid MKL matrix descriptor: fourth "
          "entry of the character array must be either 'C' or 'F', but is '"
          << descr_[3] << "' instead.");
      }
    }

    Teuchos::EDiag
    MatrixDescriptor::getDiag () const
    {
      const char s2 = descr_[2];
      if (s2 == 'U') {
        return Teuchos::UNIT_DIAG;
      } else {
        return Teuchos::NON_UNIT_DIAG;
      }
    }

    Teuchos::EUplo
    MatrixDescriptor::getUplo () const
    {
      const char s1 = descr_[1];
      if (s1 == 'L') {
        return Teuchos::LOWER_TRI;
      } else if (s1 == 'U') {
        return Teuchos::UPPER_TRI;
      } else {
        return Teuchos::UNDEF_TRI;
      }
    }

  } // namespace Mkl
} // namespace Kokkos

