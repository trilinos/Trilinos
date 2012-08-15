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

#ifndef __Kokkos_Mkl_MatrixDescriptor_hpp
#define __Kokkos_Mkl_MatrixDescriptor_hpp

#include <Teuchos_BLAS_types.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <string>

namespace Kokkos {
  namespace Mkl {
    /// \class MatrixDescriptor
    /// \brief Wrapper for MKL's sparse matrix descriptor array.
    ///
    /// \note This class is not meant for end users.  We use it in the
    ///   implementation of Kokkos::MklSparseOps.
    ///
    /// MKL's sparse kernels use an array of characters to describe
    /// the storage format of a sparse matrix.  This class wraps the
    /// matrix descriptor array and provides methods for filling in
    /// the array, validation, and printing.
    class MatrixDescriptor {
    public:
      /// \brief Construct a matrix descriptor.
      ///
      /// The properties symmetric, Hermitian, skew-symmetric, and
      /// triangular below all refer to the storage of the matrix as
      /// well as to its mathematical properties.  In all these cases,
      /// MKL will only read the lower or upper triangle of the
      /// matrix, depending on the value of the uplo parameter.  If
      /// you specify one of these properties, you must set uplo to
      /// either Teuchos::LOWER_TRI or Teuchos::UPPER_TRI.
      ///
      /// \param symmetric [in] Whether the matrix \f$A\f$ is symmetric,
      ///   i.e., \f$A = A^T\f$.
      /// \param hermitian [in] Whether the matrix \f$A\f$ is Hermitian,
      ///   i.e., \f$A = A^H\f$.  If Scalar is complex, then symmetric
      ///   means something different than Hermitian.  Otherwise, they
      ///   mean the same thing.
      /// \param triangular [in] Whether the matrix \f$A\f$ is
      ///   triangular.  MKL interprets this as meaning that it should
      ///   ignore any entries not in the part of the matrix specified
      ///   by the uplo argument.
      /// \param skewSymmetric [in] Whether the matrix \f$A\f$ is
      ///   skew-symmetric (a.k.a. antisymmetric), i.e., \f$A = -A^T\f$.
      /// \param diagonal [in] Whether the matrix is diagonal.
      /// \param uplo [in] If the matrix is triangular, then this refers
      ///   to whether the matrix is lower (Teuchos::LOWER_TRI) or upper
      ///   (Teuchos::UPPER_TRI) triangular.  If the matrix is
      ///   symmetric, Hermitian, or skew-symmetric, then this refers to
      ///   whether MKL should read only the lower or the upper triangle
      ///   of the matrix.  Otherwise, this parameter is ignored.  You
      ///   may set uplo=Teuchos::UNDEF_TRI in that case.
      /// \param diag [in] If diag=Teuchos::UNIT_DIAG, then the diagonal
      ///   entries of the matrix are not stored explicitly, but
      ///   implicitly all have the value one.  Otherwise, if
      ///   diag=Teuchos::NON_UNIT_DIAG, any diagonal entries of the
      ///   matrix must be stored explicitly (or they will be considered
      ///   to be zero).  This is ignored if the matrix is
      ///   skew-symmetric, or if it is "general" (neither symmetric,
      ///   nor Hermitiian, nor triangular, nor diagonal).
      /// \param indexBase [in] 0 for zero-based indices, 1 for
      ///   one-based indices.
      MatrixDescriptor (const bool symmetric,
                        const bool hermitian,
                        const bool triangular,
                        const bool skewSymmetric,
                        const bool diagonal,
                        const Teuchos::EUplo uplo,
                        const Teuchos::EDiag diag,
                        const int indexBase);

      //! Default constructor (fills in the descriptor with default values).
      MatrixDescriptor ();

      /// \brief Get a raw pointer to the matrix descriptor array.
      ///
      /// The array has at least 6 entries, even though MKL currently
      /// only uses the first 4 entries.
      ///
      /// \warning This pointer is only valid during the lifetime of
      ///   the MatrixDescriptor object.
      const char* const getRawPtr () const {
        return &descr_[0];
      }

      //! Get whether or not the matrix has an implicitly stored unit diagonal.
      Teuchos::EDiag getDiag () const;

      //! Get whether or not the matrix is lower or upper triangular (or neither).
      Teuchos::EUplo getUplo () const;

      //! Return the index base specified by the given matrix descriptor.
      int getIndexBase () const;

      /// \brief Print the matrix descriptor in human-readable YAML format.
      ///
      /// \param out [out] Output stream to which to print.  We
      ///   guarantee not to write to this stream or otherwise have
      ///   externally visible side effects unless the matrix descriptor
      ///   is valid.
      /// \param descr [in] The MKL sparse matrix descriptor array, as
      ///   would be created by fillMatrixDescriptor().
      ///
      /// \note We use FancyOStream instead of plain std::ostream because
      ///   FancyOStream lets us set tabs to indent each line of output.
      ///   This is nice for formatted human-readable YAML output.
      void print (Teuchos::FancyOStream& out) const;

    private:
      /// \brief Fill the given matrix descriptor array.
      ///
      /// \param descr [in] The matrix descriptor character array, as
      ///   a string.  The array must be at least 6 characters long,
      ///   although MKL currently only uses the first 4 characters.
      ///   The string will be resized if necessary.
      ///
      /// All other arguments are the same as those of the constructor.
      ///
      /// \note This method is meant to be called only by the constructor.
      static void
      fill (std::string& descr,
            const bool symmetric,
            const bool hermitian,
            const bool triangular,
            const bool skewSymmetric,
            const bool diagonal,
            const Teuchos::EUplo uplo,
            const Teuchos::EDiag diag,
            const int indexBase);

      //! The matrix descriptor character array, as a string.
      std::string descr_;
    };
  } // namespace Mkl
} // namespace Kokkos

#endif // __Kokkos_Mkl_MatrixDescriptor_hpp
