//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
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

#ifndef __Belos_ProjectedLeastSquaresSolver_hpp
#define __Belos_ProjectedLeastSquaresSolver_hpp

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

/// \file BelosProjectedLeastSquaresSolver.hpp
/// \brief Methods for solving GMRES' projected least-squares problem.
/// \author Mark Hoemmen

namespace Belos {

  /// \namespace details
  /// \brief Namespace containing implementation details of Belos solvers.
  /// \author Mark Hoemmen
  ///
  /// \warning Belos users should not use anything in this namespace.
  ///   They should not even assume that the namespace will continue to
  ///   exist between releases.  The namespace's name itself or anything
  ///   it contains may change at any time.
  namespace details {

    // Anonymous namespace restricts contents to file scope.
    namespace {
      /// \fn printMatrix
      /// \brief Print A, a dense matrix, in Matlab-readable ASCII format.
      template<class Scalar>
      void
      printMatrix (std::ostream& out,
                   const std::string& name,
                   const Teuchos::SerialDenseMatrix<int, Scalar>& A)
      {
        using std::endl;

        const int numRows = A.numRows();
        const int numCols = A.numCols();

        out << name << " = " << endl << '[';
        if (numCols == 1) {
          // Compact form for column vectors; valid Matlab.
          for (int i = 0; i < numRows; ++i) {
            out << A(i,0);
            if (i < numRows-1) {
              out << "; ";
            }
          }
        } else {
          for (int i = 0; i < numRows; ++i) {
            for (int j = 0; j < numCols; ++j) {
              out << A(i,j);
              if (j < numCols-1) {
                out << ", ";
              } else if (i < numRows-1) {
                out << ';' << endl;
              }
            }
          }
        }
        out << ']' << endl;
      }

      /// \fn print
      /// \brief Print A, a dense matrix, in Matlab-readable ASCII format.
      template<class Scalar>
      void
      print (std::ostream& out,
             const Teuchos::SerialDenseMatrix<int, Scalar>& A,
             const std::string& linePrefix)
      {
        using std::endl;

        const int numRows = A.numRows();
        const int numCols = A.numCols();

        out << linePrefix << '[';
        if (numCols == 1) {
          // Compact form for column vectors; valid Matlab.
          for (int i = 0; i < numRows; ++i) {
            out << A(i,0);
            if (i < numRows-1) {
              out << "; ";
            }
          }
        } else {
          for (int i = 0; i < numRows; ++i) {
            for (int j = 0; j < numCols; ++j) {
              if (numRows > 1) {
                out << linePrefix << "  ";
              }
              out << A(i,j);
              if (j < numCols-1) {
                out << ", ";
              } else if (i < numRows-1) {
                out << ';' << endl;
              }
            }
          }
        }
        out << linePrefix << ']' << endl;
      }
    } // namespace (anonymous)

    /// \class ProjectedLeastSquaresProblem
    /// \brief "Container" for the GMRES projected least-squares problem.
    /// \author Mark Hoemmen
    ///
    /// See the documentation of \c ProjectedLeastSquaresSolver to
    /// learn how to use this class in a GMRES-type iterative solver.
    template<class Scalar>
    class ProjectedLeastSquaresProblem {
    public:
      /// \typedef scalar_type
      /// \brief The type of the entries in the projected least-squares problem.
      typedef Scalar scalar_type;
      /// \typedef magnitude_type
      /// \brief The type of the magnitude of \c scalar_type values.
      typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

      /// \brief The upper Hessenberg matrix from GMRES.
      ///
      /// This matrix's number of rows is one more than its number of
      /// columns.  The updating methods never modify H; they just
      /// copy out the relevant data into R.  This allows GMRES
      /// implementations to implement features like backtracking
      /// (throwing away iterations).
      Teuchos::SerialDenseMatrix<int,Scalar> H;

      /// \brief Upper triangular factor from the QR factorization of H.
      ///
      /// R is a matrix with the same dimensions as H.  It is used for
      /// computing and storing the incrementally computed upper
      /// triangular factor from the QR factorization of H.  R must
      /// have the same dimensions as H (the number of rows is one
      /// more than the number of columns).
      ///
      /// H[0:k, 0:k-1] (inclusive zero-based index ranges) is the
      /// upper Hessenberg matrix for the first k iterations of GMRES
      /// (where k = 0, 1, 2, ...).
      Teuchos::SerialDenseMatrix<int,Scalar> R;

      /// \brief Current solution of the projected least-squares problem.
      ///
      /// The vector (matrix with one column) y has the same number of
      /// rows as H.  It is used to store the solution of the
      /// projected least-squares problem at each step.  The vector
      /// should have one more entry than necessary for the solution,
      /// because of the way we solve the least-squares problem.
      /// (Most of the methods require copying the right-hand side
      /// vector into y, and the right-hand side has one more entry
      /// than the solution.)
      Teuchos::SerialDenseMatrix<int,Scalar> y;

      /// \brief Current right-hand side of the projected least-squares problem.
      ///
      /// The vector (one-column matrix) z has the same number of rows
      /// as H.  It stores the current right-hand side of the
      /// projected least-squares problem.  The z vector may be
      /// updated either progressively (if a Givens rotation method is
      /// used) or all at once (if an LAPACK factorization method is
      /// used).
      Teuchos::SerialDenseMatrix<int,Scalar> z;

      /// \brief Array of cosines from the computed Givens rotations.
      ///
      /// This array is only filled in if a Givens rotation method is
      /// used for updating the least-squares problem.
      ///
      /// In the complex-arithmetic case, it is possible to compute
      /// entirely real-valued cosines.  However, we choose to
      /// represent them as the original scalar type, rather than the
      /// magnitude type.
      Teuchos::Array<Scalar> theCosines;

      /// \brief Array of sines from the computed Givens rotations.
      ///
      /// This array is only filled in if a Givens rotation method is
      /// used for updating the least-squares problem.
      Teuchos::Array<Scalar> theSines;

      /// \brief Constructor.
      ///
      /// Reserve space for a projected least-squares problem of dimension
      /// at most (maxNumIterations+1) by maxNumIterations.  "Iterations"
      /// refers to GMRES iterations.  We assume that after the first
      /// iteration (<i>not</i> counting the computation of the initial
      /// residual as an iteration), the projected least-squares problem
      /// has dimension 2 by 1.
      ProjectedLeastSquaresProblem (const int maxNumIterations) :
        H (maxNumIterations+1, maxNumIterations),
        R (maxNumIterations+1, maxNumIterations),
        y (maxNumIterations+1, 1),
        z (maxNumIterations+1, 1),
        theCosines (maxNumIterations+1),
        theSines (maxNumIterations+1)
      {}

      /// \brief Reset the projected least-squares problem.
      ///
      /// "Reset" means that the right-hand side of the projected
      /// least-squares problem is restored to \f$\beta e_1\f$.  The
      /// value \f$\beta\f$ is typically the norm of the initial or
      /// restarted residual vector.
      ///
      /// This method does not reallocate or resize any of the
      /// matrices or vectors.  Since this class keeps the original
      /// upper Hessenberg matrix, the application may choose to
      /// recompute its QR factorization up to the desired point and
      /// apply the resulting Givens rotations to the right-hand size
      /// z resulting from this reset operation.  Alternately, the
      /// application may simply overwrite the upper Hessenberg
      /// matrix's entries with new data.
      ///
      /// \param beta [in] The initial residual norm of the
      ///   (non-projected) linear system \f$Ax=b\f$.
      void
      reset (const typename Teuchos::ScalarTraits<Scalar>::magnitudeType beta)
      {
        typedef Teuchos::ScalarTraits<Scalar> STS;

        // Zero out the right-hand side of the least-squares problem.
        z.putScalar (STS::zero());

        // Promote the initial residual norm from a magnitude type to
        // a scalar type, so we can assign it to the first entry of z.
        const Scalar initialResidualNorm (beta);
        z(0,0) = initialResidualNorm;
      }

      /// \brief (Re)allocate and reset the projected least-squares problem.
      ///
      /// "(Re)allocate" means to (re)size H, R, y, and z to their
      /// appropriate maximum dimensions, given the maximum number of
      /// iterations that GMRES may execute.  "Reset" means to do what
      /// the \c reset() method does.  Reallocation happens first,
      /// then reset.
      ///
      /// \param beta [in] The initial residual norm of the
      ///   (non-projected) linear system \f$Ax=b\f$.
      ///
      /// \param maxNumIterations [in] The maximum number of
      ///   iterations that GMRES may execute.
      void
      reallocateAndReset (const typename Teuchos::ScalarTraits<Scalar>::magnitudeType beta,
                          const int maxNumIterations)
      {
        typedef Teuchos::ScalarTraits<Scalar> STS;
        typedef Teuchos::ScalarTraits<magnitude_type> STM;

        TEUCHOS_TEST_FOR_EXCEPTION(beta < STM::zero(), std::invalid_argument,
                           "ProjectedLeastSquaresProblem::reset: initial "
                           "residual beta = " << beta << " < 0.");
        TEUCHOS_TEST_FOR_EXCEPTION(maxNumIterations <= 0, std::invalid_argument,
                           "ProjectedLeastSquaresProblem::reset: maximum number "
                           "of iterations " << maxNumIterations << " <= 0.");

        if (H.numRows() < maxNumIterations+1 || H.numCols() < maxNumIterations) {
          const int errcode = H.reshape (maxNumIterations+1, maxNumIterations);
          TEUCHOS_TEST_FOR_EXCEPTION(errcode != 0, std::runtime_error,
                             "Failed to reshape H into a " << (maxNumIterations+1)
                             << " x " << maxNumIterations << " matrix.");
        }
        (void) H.putScalar (STS::zero());

        if (R.numRows() < maxNumIterations+1 || R.numCols() < maxNumIterations) {
          const int errcode = R.reshape (maxNumIterations+1, maxNumIterations);
          TEUCHOS_TEST_FOR_EXCEPTION(errcode != 0, std::runtime_error,
                             "Failed to reshape R into a " << (maxNumIterations+1)
                             << " x " << maxNumIterations << " matrix.");
        }
        (void) R.putScalar (STS::zero());

        if (y.numRows() < maxNumIterations+1 || y.numCols() < 1) {
          const int errcode = y.reshape (maxNumIterations+1, 1);
          TEUCHOS_TEST_FOR_EXCEPTION(errcode != 0, std::runtime_error,
                             "Failed to reshape y into a " << (maxNumIterations+1)
                             << " x " << 1 << " matrix.");
        }
        (void) y.putScalar (STS::zero());

        if (z.numRows() < maxNumIterations+1 || z.numCols() < 1) {
          const int errcode = z.reshape (maxNumIterations+1, 1);
          TEUCHOS_TEST_FOR_EXCEPTION(errcode != 0, std::runtime_error,
                             "Failed to reshape z into a " << (maxNumIterations+1)
                             << " x " << 1 << " matrix.");
        }
        reset (beta);
      }

    };


    /// \class LocalDenseMatrixOps
    /// \brief Low-level operations on non-distributed dense matrices.
    /// \author Mark Hoemmen
    ///
    /// This class provides a convenient wrapper around some BLAS
    /// operations, operating on non-distributed (hence "local") dense
    /// matrices.
    template<class Scalar>
    class LocalDenseMatrixOps {
    public:
      /// \typedef scalar_type
      /// \brief The template parameter of this class.
      typedef Scalar scalar_type;
      /// \typedef magnitude_type
      /// \brief The type of the magnitude of a \c scalar_type value.
      typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
      /// \typedef mat_type
      /// \brief The type of a dense matrix (or vector) of \c scalar_type.
      typedef Teuchos::SerialDenseMatrix<int,Scalar> mat_type;

    private:
      typedef Teuchos::ScalarTraits<scalar_type> STS;
      typedef Teuchos::ScalarTraits<magnitude_type> STM;
      typedef Teuchos::BLAS<int, scalar_type> blas_type;
      typedef Teuchos::LAPACK<int, scalar_type> lapack_type;

    public:
      //! A_star := (conjugate) transpose of A.
      void
      conjugateTranspose (mat_type& A_star, const mat_type& A) const
      {
        for (int i = 0; i < A.numRows(); ++i) {
          for (int j = 0; j < A.numCols(); ++j) {
            A_star(j,i) = STS::conjugate (A(i,j));
          }
        }
      }

      //! L := (conjugate) transpose of R (upper triangular).
      void
      conjugateTransposeOfUpperTriangular (mat_type& L, const mat_type& R) const
      {
        const int N = R.numCols();

        for (int j = 0; j < N; ++j) {
          for (int i = 0; i <= j; ++i) {
            L(j,i) = STS::conjugate (R(i,j));
          }
        }
      }

      //! Zero out everything below the diagonal of A.
      void
      zeroOutStrictLowerTriangle (mat_type& A) const
      {
        const int N = std::min (A.numRows(), A.numCols());

        for (int j = 0; j < N; ++j) {
          for (int i = j+1; i < A.numRows(); ++i) {
            A(i,j) = STS::zero();
          }
        }
      }

      /// \brief A -> [A_11, A_21, A_12, A_22].
      ///
      /// The first four arguments are the output arguments.  They are
      /// views of their respective submatrices of A.
      ///
      /// \note SerialDenseMatrix's operator= and copy constructor
      ///   both always copy the matrix deeply, so I can't make the
      ///   output arguments "mat_type&".
      void
      partition (Teuchos::RCP<mat_type>& A_11,
                 Teuchos::RCP<mat_type>& A_21,
                 Teuchos::RCP<mat_type>& A_12,
                 Teuchos::RCP<mat_type>& A_22,
                 mat_type& A,
                 const int numRows1,
                 const int numRows2,
                 const int numCols1,
                 const int numCols2)
      {
        using Teuchos::rcp;
        using Teuchos::View;

        A_11 = rcp (new mat_type (View, A, numRows1, numCols1, 0, 0));
        A_21 = rcp (new mat_type (View, A, numRows2, numCols1, numRows1, 0));
        A_12 = rcp (new mat_type (View, A, numRows1, numCols2, 0, numCols1));
        A_22 = rcp (new mat_type (View, A, numRows2, numCols2, numRows1, numCols1));
      }

      //! A := alpha * A.
      void
      matScale (mat_type& A, const scalar_type& alpha) const
      {
        const int LDA = A.stride();
        const int numRows = A.numRows();
        const int numCols = A.numCols();

        if (numRows == 0 || numCols == 0) {
          return;
        } else {
          for (int j = 0; j < numCols; ++j) {
            scalar_type* const A_j = &A(0,j);

            for (int i = 0; i < numRows; ++i) {
              A_j[i] *= alpha;
            }
          }
        }
      }

      /// \brief Y := Y + alpha * X.
      ///
      /// "AXPY" stands for "alpha times X plus y," and is the
      /// traditional abbreviation for this operation.
      void
      axpy (mat_type& Y,
            const scalar_type& alpha,
            const mat_type& X) const
      {
        const int numRows = Y.numRows();
        const int numCols = Y.numCols();

        TEUCHOS_TEST_FOR_EXCEPTION(numRows != X.numRows() || numCols != X.numCols(),
                           std::invalid_argument, "Dimensions of X and Y don't "
                           "match.  X is " << X.numRows() << " x " << X.numCols()
                           << ", and Y is " << numRows << " x " << numCols << ".");
        for (int j = 0; j < numCols; ++j) {
          for (int i = 0; i < numRows; ++i) {
            Y(i,j) += alpha * X(i,j);
          }
        }
      }

      //! A := A + B.
      void
      matAdd (mat_type& A, const mat_type& B) const
      {
        const int LDA = A.stride();
        const int LDB = B.stride();
        const int numRows = A.numRows();
        const int numCols = A.numCols();

        TEUCHOS_TEST_FOR_EXCEPTION(B.numRows() != numRows || B.numCols() != numCols,
                           std::invalid_argument,
                           "matAdd: The input matrices A and B have "
                           "incompatible dimensions.  A is " << numRows
                           << " x " << numCols << ", but B is " << B.numRows()
                           << " x " << B.numCols() << ".");
        if (numRows == 0 || numCols == 0) {
          return;
        } else {
          for (int j = 0; j < numCols; ++j) {
            scalar_type* const A_j = &A(0,j);
            const scalar_type* const B_j = &B(0,j);

            for (int i = 0; i < numRows; ++i) {
              A_j[i] += B_j[i];
            }
          }
        }
      }

      //! A := A - B.
      void
      matSub (mat_type& A, const mat_type& B) const
      {
        const int LDA = A.stride();
        const int LDB = B.stride();
        const int numRows = A.numRows();
        const int numCols = A.numCols();

        TEUCHOS_TEST_FOR_EXCEPTION(B.numRows() != numRows || B.numCols() != numCols,
                           std::invalid_argument,
                           "matSub: The input matrices A and B have "
                           "incompatible dimensions.  A is " << numRows
                           << " x " << numCols << ", but B is " << B.numRows()
                           << " x " << B.numCols() << ".");
        if (numRows == 0 || numCols == 0) {
          return;
        } else {
          for (int j = 0; j < numCols; ++j) {
            scalar_type* const A_j = &A(0,j);
            const scalar_type* const B_j = &B(0,j);

            for (int i = 0; i < numRows; ++i) {
              A_j[i] -= B_j[i];
            }
          }
        }
      }

      /// \brief In Matlab notation: B = B / R, where R is upper triangular.
      ///
      /// This method only looks at the upper left R.numCols() by
      /// R.numCols() part of R.
      void
      rightUpperTriSolve (mat_type& B,
                          const mat_type& R) const
      {
        TEUCHOS_TEST_FOR_EXCEPTION(B.numCols() != R.numRows(),
                           std::invalid_argument,
                           "rightUpperTriSolve: R and B have incompatible "
                           "dimensions.  B has " << B.numCols() << " columns, "
                           "but R has " << R.numRows() << " rows.");
        blas_type blas;
        blas.TRSM (Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI,
                   Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG,
                   R.numCols(), B.numCols(),
                   STS::one(), R.values(), R.stride(),
                   B.values(), B.stride());
      }

      /// \brief C := beta*C + alpha*A*B.
      ///
      /// This method is a thin wrapper around the BLAS' _GEMM
      /// routine.  The matrix C is NOT allowed to alias the matrices
      /// A or B.  This method makes no effort to check for aliasing.
      void
      matMatMult (const scalar_type& beta,
                  mat_type& C,
                  const scalar_type& alpha,
                  const mat_type& A,
                  const mat_type& B) const
      {
        using Teuchos::NO_TRANS;

        TEUCHOS_TEST_FOR_EXCEPTION(A.numCols() != B.numRows(),
                           std::invalid_argument,
                           "matMatMult: The input matrices A and B have "
                           "incompatible dimensions.  A is " << A.numRows()
                           << " x " << A.numCols() << ", but B is "
                           << B.numRows() << " x " << B.numCols() << ".");
        TEUCHOS_TEST_FOR_EXCEPTION(A.numRows() != C.numRows(),
                           std::invalid_argument,
                           "matMatMult: The input matrix A and the output "
                           "matrix C have incompatible dimensions.  A has "
                           << A.numRows() << " rows, but C has " << C.numRows()
                           << " rows.");
        TEUCHOS_TEST_FOR_EXCEPTION(B.numCols() != C.numCols(),
                           std::invalid_argument,
                           "matMatMult: The input matrix B and the output "
                           "matrix C have incompatible dimensions.  B has "
                           << B.numCols() << " columns, but C has "
                           << C.numCols() << " columns.");
        blas_type blas;
        blas.GEMM (NO_TRANS, NO_TRANS, C.numRows(), C.numCols(), A.numCols(),
                   alpha, A.values(), A.stride(), B.values(), B.stride(),
                   beta, C.values(), C.stride());
      }

      /// \brief Return the number of Inf or NaN entries in the matrix A.
      ///
      /// \param A [in] The matrix to check.
      ///
      /// \param upperTriangular [in] If true, only check the upper
      ///   triangle / trapezoid of A.  Otherwise, check all entries
      ///   of A.
      int
      infNaNCount (const mat_type& A, const bool upperTriangular=false) const
      {
        int count = 0;
        for (int j = 0; j < A.numCols(); ++j) {
          if (upperTriangular) {
            for (int i = 0; i <= j && i < A.numRows(); ++i) {
              if (STS::isnaninf (A(i,j))) {
                ++count;
              }
            }
          } else {
            for (int i = 0; i < A.numRows(); ++i) {
              if (STS::isnaninf (A(i,j))) {
                ++count;
              }
            }
          }
        }
        return count;
      }

      /// \brief Is the matrix A upper triangular / trapezoidal?
      ///
      /// \return (is upper triangular?, (squared Frobenius norm of
      ///   strict lower triangle, squared Frobenius norm of the whole
      ///   matrix))
      std::pair<bool, std::pair<magnitude_type, magnitude_type> >
      isUpperTriangular (const mat_type& A) const
      {
        magnitude_type lowerTri = STM::zero();
        magnitude_type upperTri = STM::zero();
        int count = 0;

        for (int j = 0; j < A.numCols(); ++j) {
          // Compute the Frobenius norm of the upper triangle /
          // trapezoid of A.  The second clause of the loop upper
          // bound is for matrices with fewer rows than columns.
          for (int i = 0; i <= j && i < A.numRows(); ++i) {
            const magnitude_type A_ij_mag = STS::magnitude (A(i,j));
            upperTri += A_ij_mag * A_ij_mag;
          }
          // Scan the strict lower triangle / trapezoid of A.
          for (int i = j+1; i < A.numRows(); ++i) {
            const magnitude_type A_ij_mag = STS::magnitude (A(i,j));
            lowerTri += A_ij_mag * A_ij_mag;
            if (A_ij_mag != STM::zero()) {
              ++count;
            }
          }
        }
        return std::make_pair (count == 0, std::make_pair (lowerTri, upperTri));
      }


      /// \brief Is the matrix A upper Hessenberg?
      ///
      /// \return (is upper Hessenberg?, (squared Frobenius norm of
      ///   the part of A that should be zero if A is upper
      ///   Hessenberg, squared Frobenius norm of the whole matrix))
      std::pair<bool, std::pair<magnitude_type, magnitude_type> >
      isUpperHessenberg (const mat_type& A) const
      {
        magnitude_type lower = STM::zero();
        magnitude_type upper = STM::zero();
        int count = 0;

        for (int j = 0; j < A.numCols(); ++j) {
          // Compute the Frobenius norm of the upper Hessenberg part
          // of A.  The second clause of the loop upper bound is for
          // matrices with fewer rows than columns.
          for (int i = 0; i <= j+1 && i < A.numRows(); ++i) {
            const magnitude_type A_ij_mag = STS::magnitude (A(i,j));
            upper += A_ij_mag * A_ij_mag;
          }
          // Scan the strict lower part of A.
          for (int i = j+2; i < A.numRows(); ++i) {
            const magnitude_type A_ij_mag = STS::magnitude (A(i,j));
            lower += A_ij_mag * A_ij_mag;
            if (A_ij_mag != STM::zero()) {
              ++count;
            }
          }
        }
        return std::make_pair (count == 0, std::make_pair (lower, upper));
      }

      /// \brief Throw an exception if A is not upper triangular / trapezoidal.
      ///
      /// \param A [in] The matrix to test.
      ///
      /// \param matrixName [in] Name of the matrix.  Used only to
      ///   make the exception message more informative.
      void
      ensureUpperTriangular (const mat_type& A,
                             const char* const matrixName) const
      {
        std::pair<bool, std::pair<magnitude_type, magnitude_type> > result =
          isUpperTriangular (A);

        TEUCHOS_TEST_FOR_EXCEPTION(! result.first, std::invalid_argument,
                           "The " << A.numRows() << " x " << A.numCols()
                           << " matrix " << matrixName << " is not upper "
                           "triangular.  ||tril(A)||_F = "
                           << result.second.first << " and ||A||_F = "
                           << result.second.second << ".");
      }

      /// \brief Throw an exception if A is not (strictly) upper Hessenberg.
      ///
      /// \param A [in] The matrix to test.
      ///
      /// \param matrixName [in] Name of the matrix.  Used only to
      ///   make the exception message more informative.
      void
      ensureUpperHessenberg (const mat_type& A,
                             const char* const matrixName) const
      {
        std::pair<bool, std::pair<magnitude_type, magnitude_type> > result =
          isUpperHessenberg (A);

        TEUCHOS_TEST_FOR_EXCEPTION(! result.first, std::invalid_argument,
                           "The " << A.numRows() << " x " << A.numCols()
                           << " matrix " << matrixName << " is not upper "
                           "triangular.  ||tril(A(2:end, :))||_F = "
                           << result.second.first << " and ||A||_F = "
                           << result.second.second << ".");
      }

      /// \brief Throw an exception if A is not "approximately" upper Hessenberg.
      ///
      /// "Approximately" in this case means that the Frobenius norm
      /// of the part of A that should be zero, divided by the
      /// Frobenius norm of all of A, is less than or equal to the
      /// given relative tolerance.
      ///
      /// \param A [in] The matrix to test.
      ///
      /// \param matrixName [in] Name of the matrix.  Used only to
      ///   make the exception message more informative.
      ///
      /// \param relativeTolerance [in] Amount by which we allow the
      ///   norm of the part of A that should be zero to deviate from
      ///   zero, relative to the norm of A.
      void
      ensureUpperHessenberg (const mat_type& A,
                             const char* const matrixName,
                             const magnitude_type relativeTolerance) const
      {
        std::pair<bool, std::pair<magnitude_type, magnitude_type> > result =
          isUpperHessenberg (A);

        if (result.first) {
          // Mollified relative departure from upper Hessenberg.
          const magnitude_type err = (result.second.second == STM::zero() ?
                                      result.second.first :
                                      result.second.first / result.second.second);
          TEUCHOS_TEST_FOR_EXCEPTION(err > relativeTolerance, std::invalid_argument,
                             "The " << A.numRows() << " x " << A.numCols()
                             << " matrix " << matrixName << " is not upper "
                             "triangular.  ||tril(A(2:end, :))||_F "
                             << (result.second.second == STM::zero() ? "" : " / ||A||_F")
                             << " = " << err << " > " << relativeTolerance << ".");
        }
      }

      /// \brief Ensure that the matrix A is at least minNumRows by
      ///   minNumCols.
      ///
      /// If A has fewer rows than minNumRows, or fewer columns than
      /// minNumCols, this method throws an informative exception.
      ///
      /// \param A [in] The matrix whose dimensions to check.
      /// \param matrixName [in] Name of the matrix; used to make the
      ///   exception message more informative.
      /// \param minNumRows [in] Minimum number of rows allowed in A.
      /// \param minNumCols [in] Minimum number of columns allowed in A.
      void
      ensureMinimumDimensions (const mat_type& A,
                               const char* const matrixName,
                               const int minNumRows,
                               const int minNumCols) const
      {
        TEUCHOS_TEST_FOR_EXCEPTION(A.numRows() < minNumRows || A.numCols() < minNumCols,
                           std::invalid_argument,
                           "The matrix " << matrixName << " is " << A.numRows()
                           << " x " << A.numCols() << ", and therefore does not "
                           "satisfy the minimum dimensions " << minNumRows
                           << " x " << minNumCols << ".");
      }

      /// \brief Ensure that the matrix A is exactly numRows by numCols.
      ///
      /// If A has a different number of rows than numRows, or a
      /// different number of columns than numCols, this method throws
      /// an informative exception.
      ///
      /// \param A [in] The matrix whose dimensions to check.
      /// \param matrixName [in] Name of the matrix; used to make the
      ///   exception message more informative.
      /// \param numRows [in] Number of rows that A must have.
      /// \param numCols [in] Number of columns that A must have.
      void
      ensureEqualDimensions (const mat_type& A,
                             const char* const matrixName,
                             const int numRows,
                             const int numCols) const
      {
        TEUCHOS_TEST_FOR_EXCEPTION(A.numRows() != numRows || A.numCols() != numCols,
                           std::invalid_argument,
                           "The matrix " << matrixName << " is supposed to be "
                           << numRows << " x " << numCols << ", but is "
                           << A.numRows() << " x " << A.numCols() << " instead.");
      }

    };

    /// \enum ERobustness
    /// \brief Robustness level of projected least-squares solver operations.
    /// \author Mark Hoemmen
    ///
    /// "Robustness" refers in particular to dense triangular solves.
    /// ROBUSTNESS_NONE means use the BLAS' _TRSM, which may result in
    /// inaccurate or invalid (e.g., NaN or Inf) results if the upper
    /// triangular matrix R is singular.  ROBUSTNESS_LOTS means use an
    /// SVD-based least-squares solver for upper triangular solves.
    /// This will work even if R is singular, but is much more
    /// expensive than _TRSM (in fact, it's at least N times more
    /// expensive, where N is the number of columns in R).
    /// ROBUSTNESS_SOME means some algorithmic point in between those
    /// two extremes.
    ///
    /// The projected least-squares solver may use this enum to
    /// control the robustness of other operations besides dense
    /// triangular solves.
    enum ERobustness {
      ROBUSTNESS_NONE,
      ROBUSTNESS_SOME,
      ROBUSTNESS_LOTS,
      ROBUSTNESS_INVALID
    };

    //! Convert the given ERobustness enum value to a string.
    inline std::string
    robustnessEnumToString (const ERobustness x)
    {
      const char* strings[] = {"None", "Some", "Lots"};
      TEUCHOS_TEST_FOR_EXCEPTION(x < ROBUSTNESS_NONE || x >= ROBUSTNESS_INVALID,
                         std::invalid_argument,
                         "Invalid enum value " << x << ".");
      return std::string (strings[x]);
    }

    //! Convert the given robustness string value to an ERobustness enum.
    ERobustness
    inline robustnessStringToEnum (const std::string& x)
    {
      const char* strings[] = {"None", "Some", "Lots"};
      for (int r = 0; r < static_cast<int> (ROBUSTNESS_INVALID); ++r) {
        if (x == strings[r]) {
          return static_cast<ERobustness> (r);
        }
      }
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
                         "Invalid robustness string " << x << ".");
    }

    /// \brief Make a ParameterList validator for ERobustness.
    ///
    /// Use this validator when you are setting up a default
    /// ParameterList for a solver that accepts a ERobustness
    /// parameter.  This will let users supply a string that turns
    /// into an enum value.
    ///
    /// The validator converts strings ("None", "Some", "Lots") to the
    /// corresponding ERobustness enum values.  It also includes
    /// documentation for each value.
    inline Teuchos::RCP<Teuchos::ParameterEntryValidator>
    robustnessValidator ()
    {
      using Teuchos::stringToIntegralParameterEntryValidator;

      Teuchos::Array<std::string> strs (3);
      strs[0] = robustnessEnumToString (ROBUSTNESS_NONE);
      strs[1] = robustnessEnumToString (ROBUSTNESS_SOME);
      strs[2] = robustnessEnumToString (ROBUSTNESS_LOTS);
      Teuchos::Array<std::string> docs (3);
      docs[0] = "Use the BLAS' triangular solve.  This may result in Inf or "
        "NaN output if the triangular matrix is rank deficient.";
      docs[1] = "Robustness somewhere between \"None\" and \"Lots\".";
      docs[2] = "Solve the triangular system in a least-squares sense, using "
        "an SVD-based algorithm.  This will always succeed, though the "
        "solution may not make sense for GMRES.";
      Teuchos::Array<ERobustness> ints (3);
      ints[0] = ROBUSTNESS_NONE;
      ints[1] = ROBUSTNESS_SOME;
      ints[2] = ROBUSTNESS_LOTS;
      const std::string pname ("Robustness of Projected Least-Squares Solve");

      return stringToIntegralParameterEntryValidator<ERobustness> (strs, docs,
                                                                   ints, pname);
    }

    /// \class ProjectedLeastSquaresSolver
    /// \brief Methods for solving GMRES' projected least-squares problem.
    /// \author Mark Hoemmen
    ///
    /// \tparam Scalar The type of the matrix and vector entries in the
    ///   least-squares problem.
    ///
    /// Expected use of this class:
    /// -# Use a \c ProjectedLeastSquaresProblem struct instance to
    ///    store the projected problem in your GMRES solver.
    /// -# Instantiate a ProjectedLeastSquaresSolver:
    ///    \code
    ///    ProjectedLeastSquaresSolver<Scalar> solver;
    ///    \endcode
    /// -# Update the current column(s) of the QR factorization
    ///    of GMRES' upper Hessenberg matrix:
    ///    \code
    ///    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
    ///    // Standard GMRES: update one column at a time.
    ///    MT resNorm = solver.updateColumn (problem, endCol);
    ///    \endcode
    ///    or
    ///    \code
    ///    // Some GMRES variants update multiple columns at a time.
    ///    // You can also use this feature to recompute the upper
    ///    // Hessenberg's QR factorization from scratch.
    ///    MT resNorm = solver.updateColumns (problem, startCol, endCol);
    ///    \endcode
    /// -# Solve for the current GMRES solution update coefficients:
    ///    \code
    ///    solver.solve (problem, endCol);
    ///    \endcode
    ///
    /// You can defer Step 4 as long as you want.  Step 4 must always
    /// follow Step 3.
    ///
    /// Purposes of this class:
    /// -# Isolate and factor out BLAS and LAPACK dependencies.  This
    ///    makes it easier to write custom replacements for routines
    ///    for which no Scalar specialization is available.
    /// -# Encapsulate common functionality of many GMRES-like
    ///    solvers.  Avoid duplicated code and simplify debugging,
    ///    testing, and implementation of new GMRES-like solvers.
    /// -# Provide an option for more robust implementations of
    ///    solvers for the projected least-squares problem.
    ///
    /// "Robust" here means regularizing the least-squares solve, so
    /// that the solution is well-defined even if the problem is
    /// ill-conditioned.  Many distributed-memory iterative solvers,
    /// including those in Belos, currently solve the projected
    /// least-squares problem redundantly on different processes.  If
    /// those processes are heterogeneous or implement the BLAS and
    /// LAPACK themselves in parallel (via multithreading, for
    /// example), then different BLAS or LAPACK calls on different
    /// processes may result in different answers.  The answers may be
    /// significantly different if the projected problem is singular
    /// or ill-conditioned.  This is bad because GMRES variants use
    /// the projected problem's solution as the coefficients for the
    /// solution update.  The solution update coefficients must be
    /// (almost nearly) the same on all processes.  Regularizing the
    /// projected problem is one way to ensure that different
    /// processes compute (almost) the same solution.
    template<class Scalar>
    class ProjectedLeastSquaresSolver {
    public:
      /// \typedef scalar_type
      /// \brief The template parameter of this class.
      ///
      /// The type of the matrix and vector entries in the projected
      /// least-squares problem to solve, and the type of the
      /// resulting GMRES solution update coefficients.
      typedef Scalar scalar_type;
      /// \typedef magnitude_type
      /// \brief The type of the magnitude of a \c scalar_type value.
      ///
      /// If scalar_type is complex-valued, then magnitude_type is real.
      typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
      /// \typedef mat_type
      /// \brief The type of a dense matrix (or vector) of \c scalar_type.
      typedef Teuchos::SerialDenseMatrix<int,Scalar> mat_type;

    private:
      typedef Teuchos::ScalarTraits<scalar_type> STS;
      typedef Teuchos::ScalarTraits<magnitude_type> STM;
      typedef Teuchos::BLAS<int, scalar_type> blas_type;
      typedef Teuchos::LAPACK<int, scalar_type> lapack_type;

    public:
      /// \brief Constructor.
      ///
      /// \param warnStream [out] Stream to which to output warnings.
      ///   Set to a Teuchos::oblackholestream if you don't want to
      ///   display warnings.
      ///
      /// \param defaultRobustness [in] Default robustness level for
      ///   operations like triangular solves.  For example, at a low
      ///   robustness level, triangular solves might fail if the
      ///   triangular matrix is singular or ill-conditioned in a
      ///   particular way.  At a high robustness level, triangular
      ///   solves will always succeed, but may only be solved in a
      ///   least-squares sense.
      ProjectedLeastSquaresSolver (std::ostream& warnStream,
                                   const ERobustness defaultRobustness=ROBUSTNESS_NONE) :
        warn_ (warnStream),
        defaultRobustness_ (defaultRobustness)
      {}

      /// \brief Update column curCol of the projected least-squares problem.
      ///
      /// The upper Hessenberg matrix H is read but not touched.  The
      /// R factor, the cosines and sines, and the right-hand side z
      /// are updated.  This method does <i>not</i> compute the
      /// solution of the least-squares problem; call \c solve() for
      /// that.
      ///
      /// \param problem [in/out] The projected least-squares problem.
      /// \param curCol [in] Zero-based index of the current column to update.
      ///
      /// \return 2-norm of the absolute residual of the projected
      ///   least-squares problem.
      magnitude_type
      updateColumn (ProjectedLeastSquaresProblem<Scalar>& problem,
                    const int curCol)
      {
        return updateColumnGivens (problem.H, problem.R, problem.y, problem.z,
                                   problem.theCosines, problem.theSines, curCol);
      }

      /// \brief Update columns [startCol,endCol] of the projected least-squares problem.
      ///
      /// The upper Hessenberg matrix H is read but not touched.  The
      /// R factor, the cosines and sines, and the right-hand side z
      /// are updated.  This method does <i>not</i> compute the
      /// solution of the least-squares problem; call \c solve() for
      /// that.
      ///
      /// \param problem [in/out] The projected least-squares problem.
      /// \param startCol [in] Zero-based index of the first column to update.
      /// \param endCol [in] Zero-based index of the last column (inclusive) to update.
      ///
      /// \return 2-norm of the absolute residual of the projected
      ///   least-squares problem.
      magnitude_type
      updateColumns (ProjectedLeastSquaresProblem<Scalar>& problem,
                     const int startCol,
                     const int endCol)
      {
        return updateColumnsGivens (problem.H, problem.R, problem.y, problem.z,
                                    problem.theCosines, problem.theSines,
                                    startCol, endCol);
      }

      /// \brief Solve the projected least-squares problem.
      ///
      /// Call this method only after calling \c updateColumn() or \c
      /// updateColumns().  If you call \c updateColumn(), use the
      /// same column index when calling this method.  If you call \c
      /// updateColumns(), use the endCol argument as the column index
      /// for calling this method.
      ///
      /// \param problem [in/out] The projected least-squares problem.
      ///
      /// \param curCol [in] Zero-based index of the most recently
      ///   updated column of the least-squares problem.
      void
      solve (ProjectedLeastSquaresProblem<Scalar>& problem,
             const int curCol)
      {
        solveGivens (problem.y, problem.R, problem.z, curCol);
      }

      /// \brief Solve the given square upper triangular linear system(s).
      ///
      /// This method does the same thing as the five-argument
      /// overload, but uses the default robustness level.
      std::pair<int, bool>
      solveUpperTriangularSystem (Teuchos::ESide side,
                                  mat_type& X,
                                  const mat_type& R,
                                  const mat_type& B)
      {
        return solveUpperTriangularSystem (side, X, R, B, defaultRobustness_);
      }

      /// \brief Solve the given square upper triangular linear system(s).
      ///
      /// \param side [in] Whether to solve with the triangular matrix
      ///   R on the left side or right side of the matrix X of
      ///   unknowns.  side = Teuchos::LEFT_SIDE or
      ///   Teuchos::RIGHT_SIDE.
      ///
      /// \param X [out] The matrix of unknowns.  (We consider a
      ///   matrix with one column a vector.)
      ///
      /// \param R [in] The square upper triangular matrix.  If R has
      ///   more rows than columns, this method uses the number of
      ///   columns of R as the dimension of the square matrix, and
      ///   ignores the "extra" rows of R.  If R has more columns than
      ///   rows, this method throws an exception.
      ///
      /// \param B [in] The matrix of right-hand sides.  B must have
      ///   at least as many columns as X, otherwise this method will
      ///   throw an exception.  If B has more columns than X, those
      ///   columns will be ignored.
      ///
      /// \param robustness [in] The robustness level to use for the
      ///   solve.  ROBUSTNESS_NONE is fastest but does not attempt to
      ///   detect rank deficiency in R.  ROBUSTNESS_LOTS is slowest,
      ///   but can "solve" the problem (in a least-squares
      ///   minimum-norm-solution sense) even if R is singular.
      ///   ROBUSTNESS_SOME may do something in between.
      ///
      /// Depending on the specified robustness level, this method may
      /// also attempt to compute the numerical rank and detect rank
      /// deficiency in R.  The first return value is the computed
      /// numerical rank, and the second return value is a Boolean
      /// indicating whether any rank deficiency was detected.  If the
      /// method does not attempt to compute the numerical rank or
      /// detect rank deficiency, it will return (N, false), where N
      /// is the number of columns in R.
      ///
      /// \return (detectedRank, foundRankDeficiency).
      std::pair<int, bool>
      solveUpperTriangularSystem (Teuchos::ESide side,
                                  mat_type& X,
                                  const mat_type& R,
                                  const mat_type& B,
                                  const ERobustness robustness)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(X.numRows() != B.numRows(), std::invalid_argument,
                           "The output X and right-hand side B have different "
                           "numbers of rows.  X has " << X.numRows() << " rows"
                           ", and B has " << B.numRows() << " rows.");
        // If B has more columns than X, we ignore the remaining
        // columns of B when solving the upper triangular system.  If
        // B has _fewer_ columns than X, we can't solve for all the
        // columns of X, so we throw an exception.
        TEUCHOS_TEST_FOR_EXCEPTION(X.numCols() > B.numCols(), std::invalid_argument,
                           "The output X has more columns than the "
                           "right-hand side B.  X has " << X.numCols()
                           << " columns and B has " << B.numCols()
                           << " columns.");
        // See above explaining the number of columns in B_view.
        mat_type B_view (Teuchos::View, B, B.numRows(), X.numCols());

        // Both the BLAS' _TRSM and LAPACK's _LATRS overwrite the
        // right-hand side with the solution, so first copy B_view
        // into X.
        X.assign (B_view);

        // Solve the upper triangular system.
        return solveUpperTriangularSystemInPlace (side, X, R, robustness);
      }

      /// \brief Solve square upper triangular linear system(s) in place.
      ///
      /// This method does the same thing as its four-argument
      /// overload, but uses the default robustness level.
      std::pair<int, bool>
      solveUpperTriangularSystemInPlace (Teuchos::ESide side,
                                         mat_type& X,
                                         const mat_type& R)
      {
        return solveUpperTriangularSystemInPlace (side, X, R, defaultRobustness_);
      }

      /// \brief Solve square upper triangular linear system(s) in place.
      ///
      /// This is the "in-place" version of \c
      /// solveUpperTriangularSystem().  The difference between that
      /// method and this one is that this method assumes that the
      /// right-hand side(s) is/are already stored in X on input.  It
      /// then overwrites X with the solution.
      std::pair<int, bool>
      solveUpperTriangularSystemInPlace (Teuchos::ESide side,
                                         mat_type& X,
                                         const mat_type& R,
                                         const ERobustness robustness)
      {
        using Teuchos::Array;
        using Teuchos::Copy;
        using Teuchos::LEFT_SIDE;
        using Teuchos::RIGHT_SIDE;
        LocalDenseMatrixOps<Scalar> ops;

        const int M = R.numRows();
        const int N = R.numCols();
        TEUCHOS_TEST_FOR_EXCEPTION(M < N, std::invalid_argument,
                           "The input matrix R has fewer columns than rows.  "
                           "R is " << M << " x " << N << ".");
        // Ignore any additional rows of R by working with a square view.
        mat_type R_view (Teuchos::View, R, N, N);

        if (side == LEFT_SIDE) {
          TEUCHOS_TEST_FOR_EXCEPTION(X.numRows() < N, std::invalid_argument,
                             "The input/output matrix X has only "
                             << X.numRows() << " rows, but needs at least "
                             << N << " rows to match the matrix for a "
                             "left-side solve R \\ X.");
        } else if (side == RIGHT_SIDE) {
          TEUCHOS_TEST_FOR_EXCEPTION(X.numCols() < N, std::invalid_argument,
                             "The input/output matrix X has only "
                             << X.numCols() << " columns, but needs at least "
                             << N << " columns to match the matrix for a "
                             "right-side solve X / R.");
        }
        TEUCHOS_TEST_FOR_EXCEPTION(robustness < ROBUSTNESS_NONE ||
                           robustness >= ROBUSTNESS_INVALID,
                           std::invalid_argument,
                           "Invalid robustness value " << robustness << ".");

        // In robust mode, scan the matrix and right-hand side(s) for
        // Infs and NaNs.  Only look at the upper triangle of the
        // matrix.
        if (robustness > ROBUSTNESS_NONE) {
          int count = ops.infNaNCount (R_view, true);
          TEUCHOS_TEST_FOR_EXCEPTION(count > 0, std::runtime_error,
                             "There " << (count != 1 ? "are" : "is")
                             << " " << count << " Inf or NaN entr"
                             << (count != 1 ? "ies" : "y")
                             << " in the upper triangle of R.");
          count = ops.infNaNCount (X, false);
          TEUCHOS_TEST_FOR_EXCEPTION(count > 0, std::runtime_error,
                             "There " << (count != 1 ? "are" : "is")
                             << " " << count << " Inf or NaN entr"
                             << (count != 1 ? "ies" : "y") << " in the "
                             "right-hand side(s) X.");
        }

        // Pair of values to return from this method.
        int rank = N;
        bool foundRankDeficiency = false;

        // Solve for X.
        blas_type blas;

        if (robustness == ROBUSTNESS_NONE) {
          // Fast triangular solve using the BLAS' _TRSM.  This does
          // no checking for rank deficiency.
          blas.TRSM(side, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, X.numRows(), X.numCols(),
                    STS::one(), R.values(), R.stride(),
                    X.values(), X.stride());
        } else if (robustness < ROBUSTNESS_INVALID) {
          // Save a copy of X, since X contains the right-hand side on
          // input.
          mat_type B (Copy, X, X.numRows(), X.numCols());

          // Fast triangular solve using the BLAS' _TRSM.  This does
          // no checking for rank deficiency.
          blas.TRSM(side, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                    Teuchos::NON_UNIT_DIAG, X.numRows(), X.numCols(),
                    STS::one(), R.values(), R.stride(),
                    X.values(), X.stride());

          // Check for Infs or NaNs in X.  If there are any, then
          // assume that TRSM failed, and use a more robust algorithm.
          if (ops.infNaNCount (X, false) != 0) {

            warn_ << "Upper triangular solve: Found Infs and/or NaNs in the "
              "solution after using the fast algorithm.  Retrying using a more "
              "robust algorithm." << std::endl;

            // Restore X from the copy.
            X.assign (B);

            // Find the minimum-norm solution to the least-squares
            // problem $\min_x \|RX - B\|_2$, using the singular value
            // decomposition (SVD).
            LocalDenseMatrixOps<Scalar> ops;
            if (side == LEFT_SIDE) {
              // _GELSS overwrites its matrix input, so make a copy.
              mat_type R_copy (Teuchos::Copy, R_view, N, N);

              // Zero out the lower triangle of R_copy, since the
              // mat_type constructor copies all the entries, not just
              // the upper triangle.  _GELSS will read all the entries
              // of the input matrix.
              ops.zeroOutStrictLowerTriangle (R_copy);

              // Solve the least-squares problem.
              rank = solveLeastSquaresUsingSVD (R_copy, X);
            } else {
              // If solving with R on the right-hand side, the interface
              // requires that instead of solving $\min \|XR - B\|_2$,
              // we have to solve $\min \|R^* X^* - B^*\|_2$.  We
              // compute (conjugate) transposes in newly allocated
              // temporary matrices X_star resp. R_star.  (B is already
              // in X and _GELSS overwrites its input vector X with the
              // solution.)
              mat_type X_star (X.numCols(), X.numRows());
              ops.conjugateTranspose (X_star, X);
              mat_type R_star (N, N); // Filled with zeros automatically.
              ops.conjugateTransposeOfUpperTriangular (R_star, R);

              // Solve the least-squares problem.
              rank = solveLeastSquaresUsingSVD (R_star, X_star);

              // Copy the transpose of X_star back into X.
              ops.conjugateTranspose (X, X_star);
            }
            if (rank < N) {
              foundRankDeficiency = true;
            }
          }
        } else {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                             "Should never get here!  Invalid robustness value "
                             << robustness << ".  Please report this bug to the "
                             "Belos developers.");
        }
        return std::make_pair (rank, foundRankDeficiency);
      }


    public:
      /// \brief Test Givens rotations.
      ///
      /// This routine tests both computing Givens rotations (via \c
      /// computeGivensRotation()) and applying them.
      ///
      /// \param out [out] Stream to which to write test output.
      ///
      /// \return true if the test succeeded, else false.
      bool
      testGivensRotations (std::ostream& out)
      {
        using std::endl;

        out << "Testing Givens rotations:" << endl;
        Scalar x = STS::random();
        Scalar y = STS::random();
        out << "  x = " << x << ", y = " << y << endl;

        Scalar theCosine, theSine, result;
        blas_type blas;
        computeGivensRotation (x, y, theCosine, theSine, result);
        out << "-- After computing rotation:" << endl;
        out << "---- cos,sin = " << theCosine << "," << theSine << endl;
        out << "---- x = " << x << ", y = " << y
            << ", result = " << result << endl;

        blas.ROT (1, &x, 1, &y, 1, &theCosine, &theSine);
        out << "-- After applying rotation:" << endl;
        out << "---- cos,sin = " << theCosine << "," << theSine << endl;
        out << "---- x = " << x << ", y = " << y << endl;

        // Allow only a tiny bit of wiggle room for zeroing-out of y.
        if (STS::magnitude(y) > 2*STS::eps())
          return false;
        else
          return true;
      }

      /// \brief Test update and solve using Givens rotations.
      ///
      /// \param out [out] Output stream to which to print results.
      ///
      /// \param numCols [in] Number of columns in the projected
      ///   least-squares problem to test.  (The number of rows is one
      ///   plus the number of columns.)
      ///
      /// \param testBlockGivens [in] Whether to test the "block"
      ///   (i.e., panel) version of the Givens rotations update.
      ///
      /// \param extraVerbose [in] Whether to print extra-verbose
      ///   output (e.g., the test problem and results).
      ///
      /// \return Whether the test succeeded, meaning that none of the
      ///   solves reported failure and the least-squares solution
      ///   error was within the expected bound.
      ///
      /// Test updating and solving the least-squares problem using
      /// Givens rotations by comparison against LAPACK's
      /// least-squares solver.  First generate a random least-squares
      /// problem that looks like it comes from GMRES.  The matrix is
      /// upper Hessenberg, and the right-hand side starts out with
      /// the first entry being nonzero with nonnegative real part and
      /// zero imaginary part, and all the other entries being zero.
      /// Then compare the results of \c updateColumnGivens() (applied
      /// to each column in turn) and solveGivens() (applied at the
      /// end) with the results of \c solveLapack() (applied at the
      /// end).  Print out the results to the given output stream.
      bool
      testUpdateColumn (std::ostream& out,
                        const int numCols,
                        const bool testBlockGivens=false,
                        const bool extraVerbose=false)
      {
        using Teuchos::Array;
        using std::endl;

        TEUCHOS_TEST_FOR_EXCEPTION(numCols <= 0, std::invalid_argument,
                           "numCols = " << numCols << " <= 0.");
        const int numRows = numCols + 1;

        mat_type H (numRows, numCols);
        mat_type z (numRows, 1);

        mat_type R_givens (numRows, numCols);
        mat_type y_givens (numRows, 1);
        mat_type z_givens (numRows, 1);
        Array<Scalar> theCosines (numCols);
        Array<Scalar> theSines (numCols);

        mat_type R_blockGivens (numRows, numCols);
        mat_type y_blockGivens (numRows, 1);
        mat_type z_blockGivens (numRows, 1);
        Array<Scalar> blockCosines (numCols);
        Array<Scalar> blockSines (numCols);
        const int panelWidth = std::min (3, numCols);

        mat_type R_lapack (numRows, numCols);
        mat_type y_lapack (numRows, 1);
        mat_type z_lapack (numRows, 1);

        // Make a random least-squares problem.
        makeRandomProblem (H, z);
        if (extraVerbose) {
          printMatrix<Scalar> (out, "H", H);
          printMatrix<Scalar> (out, "z", z);
        }

        // Set up the right-hand side copies for each of the methods.
        // Each method is free to overwrite its given right-hand side.
        z_givens.assign (z);
        if (testBlockGivens) {
          z_blockGivens.assign (z);
        }
        z_lapack.assign (z);

        //
        // Imitate how one would update the least-squares problem in a
        // typical GMRES implementation, for each updating method.
        //
        // Update using Givens rotations, one at a time.
        magnitude_type residualNormGivens = STM::zero();
        for (int curCol = 0; curCol < numCols; ++curCol) {
          residualNormGivens = updateColumnGivens (H, R_givens, y_givens, z_givens,
                                                   theCosines, theSines, curCol);
        }
        solveGivens (y_givens, R_givens, z_givens, numCols-1);

        // Update using the "panel left-looking" Givens approach, with
        // the given panel width.
        magnitude_type residualNormBlockGivens = STM::zero();
        if (testBlockGivens) {
          const bool testBlocksAtATime = true;
          if (testBlocksAtATime) {
            // Blocks of columns at a time.
            for (int startCol = 0; startCol < numCols; startCol += panelWidth) {
              int endCol = std::min (startCol + panelWidth - 1, numCols - 1);
              residualNormBlockGivens =
                updateColumnsGivens (H, R_blockGivens, y_blockGivens, z_blockGivens,
                                     blockCosines, blockSines, startCol, endCol);
            }
          } else {
            // One column at a time.  This is good as a sanity check
            // to make sure updateColumnsGivens() with a single column
            // does the same thing as updateColumnGivens().
            for (int startCol = 0; startCol < numCols; ++startCol) {
              residualNormBlockGivens =
                updateColumnsGivens (H, R_blockGivens, y_blockGivens, z_blockGivens,
                                     blockCosines, blockSines, startCol, startCol);
            }
          }
          // The panel version of Givens should compute the same
          // cosines and sines as the non-panel version, and should
          // update the right-hand side z in the same way.  Thus, we
          // should be able to use the same triangular solver.
          solveGivens (y_blockGivens, R_blockGivens, z_blockGivens, numCols-1);
        }

        // Solve using LAPACK's least-squares solver.
        const magnitude_type residualNormLapack =
          solveLapack (H, R_lapack, y_lapack, z_lapack, numCols-1);

        // Compute the condition number of the least-squares problem.
        // This requires a residual, so use the residual from the
        // LAPACK method.  All that the method needs for an accurate
        // residual norm is forward stability.
        const magnitude_type leastSquaresCondNum =
          leastSquaresConditionNumber (H, z, residualNormLapack);

        // Compute the relative least-squares solution error for both
        // Givens methods.  We assume that the LAPACK solution is
        // "exact" and compare against the Givens rotations solution.
        // This is taking liberties with the definition of condition
        // number, but it's the best we can do, since we don't know
        // the exact solution and don't have an extended-precision
        // solver.

        // The solution lives only in y[0 .. numCols-1].
        mat_type y_givens_view (Teuchos::View, y_givens, numCols, 1);
        mat_type y_blockGivens_view (Teuchos::View, y_blockGivens, numCols, 1);
        mat_type y_lapack_view (Teuchos::View, y_lapack, numCols, 1);

        const magnitude_type givensSolutionError =
          solutionError (y_givens_view, y_lapack_view);
        const magnitude_type blockGivensSolutionError = testBlockGivens ?
          solutionError (y_blockGivens_view, y_lapack_view) :
          STM::zero();

        // If printing out the matrices, copy out the upper triangular
        // factors for printing.  (Both methods are free to leave data
        // below the lower triangle.)
        if (extraVerbose) {
          mat_type R_factorFromGivens (numCols, numCols);
          mat_type R_factorFromBlockGivens (numCols, numCols);
          mat_type R_factorFromLapack (numCols, numCols);

          for (int j = 0; j < numCols; ++j) {
            for (int i = 0; i <= j; ++i) {
              R_factorFromGivens(i,j) = R_givens(i,j);
              if (testBlockGivens) {
                R_factorFromBlockGivens(i,j) = R_blockGivens(i,j);
              }
              R_factorFromLapack(i,j) = R_lapack(i,j);
            }
          }

          printMatrix<Scalar> (out, "R_givens", R_factorFromGivens);
          printMatrix<Scalar> (out, "y_givens", y_givens_view);
          printMatrix<Scalar> (out, "z_givens", z_givens);

          if (testBlockGivens) {
            printMatrix<Scalar> (out, "R_blockGivens", R_factorFromBlockGivens);
            printMatrix<Scalar> (out, "y_blockGivens", y_blockGivens_view);
            printMatrix<Scalar> (out, "z_blockGivens", z_blockGivens);
          }

          printMatrix<Scalar> (out, "R_lapack", R_factorFromLapack);
          printMatrix<Scalar> (out, "y_lapack", y_lapack_view);
          printMatrix<Scalar> (out, "z_lapack", z_lapack);
        }

        // Compute the (Frobenius) norm of the original matrix H.
        const magnitude_type H_norm = H.normFrobenius();

        out << "||H||_F = " << H_norm << endl;

        out << "||H y_givens - z||_2 / ||H||_F = "
            << leastSquaresResidualNorm (H, y_givens_view, z) / H_norm << endl;
        if (testBlockGivens) {
          out << "||H y_blockGivens - z||_2 / ||H||_F = "
              << leastSquaresResidualNorm (H, y_blockGivens_view, z) / H_norm << endl;
        }
        out << "||H y_lapack - z||_2 / ||H||_F = "
            << leastSquaresResidualNorm (H, y_lapack_view, z) / H_norm << endl;

        out << "||y_givens - y_lapack||_2 / ||y_lapack||_2 = "
            << givensSolutionError << endl;
        if (testBlockGivens) {
          out << "||y_blockGivens - y_lapack||_2 / ||y_lapack||_2 = "
              << blockGivensSolutionError << endl;
        }

        out << "Least-squares condition number = "
            << leastSquaresCondNum << endl;

        // Now for the controversial part of the test: judging whether
        // we succeeded.  This includes the problem's condition
        // number, which is a measure of the maximum perturbation in
        // the solution for which we can still say that the solution
        // is valid.  We include a little wiggle room by including a
        // factor proportional to the square root of the number of
        // floating-point operations that influence the last entry
        // (the conventional Wilkinsonian heuristic), times 10 for
        // good measure.
        //
        // (The square root looks like it has something to do with an
        // average-case probabilistic argument, but doesn't really.
        // What's an "average problem"?)
        const magnitude_type wiggleFactor =
          10 * STM::squareroot( numRows*numCols );
        const magnitude_type solutionErrorBoundFactor =
          wiggleFactor * leastSquaresCondNum;
        const magnitude_type solutionErrorBound =
          solutionErrorBoundFactor * STS::eps();
        out << "Solution error bound: " << solutionErrorBoundFactor
            << " * eps = " << solutionErrorBound << endl;

        // Remember that NaN is not greater than, not less than, and
        // not equal to any other number, including itself.  Some
        // compilers will rudely optimize away the "x != x" test.
        if (STM::isnaninf (solutionErrorBound)) {
          // Hm, the solution error bound is Inf or NaN.  This
          // probably means that the test problem was generated
          // incorrectly.  We should return false in this case.
          return false;
        } else { // solution error bound is finite.
          if (STM::isnaninf (givensSolutionError)) {
            return false;
          } else if (givensSolutionError > solutionErrorBound) {
            return false;
          } else if (testBlockGivens) {
            if (STM::isnaninf (blockGivensSolutionError)) {
              return false;
            } else if (blockGivensSolutionError > solutionErrorBound) {
              return false;
            } else { // Givens and Block Givens tests succeeded.
              return true;
            }
          } else { // Not testing block Givens; Givens test succeeded.
            return true;
          }
        }
      }

      /// \brief Test upper triangular solves.
      ///
      /// \param out [out] Stream to which this method writes output.
      /// \param testProblemSize [in] Dimension (number of rows and
      ///   columns) of the upper triangular matrices tested by this
      ///   method.
      /// \param robustness [in] Robustness level of the upper
      ///   triangular solves.
      /// \param verbose [in] Whether to print more verbose output.
      ///
      /// \return True if test passed, else false.
      ///
      /// \warning Before calling this routine, seed the random number
      ///   generator via Teuchos::ScalarTraits<Scalar>::seedrandom
      ///   (seed).
      bool
      testTriangularSolves (std::ostream& out,
                            const int testProblemSize,
                            const ERobustness robustness,
                            const bool verbose=false)
      {
        using Teuchos::LEFT_SIDE;
        using Teuchos::RIGHT_SIDE;
        using std::endl;
        typedef Teuchos::SerialDenseMatrix<int, scalar_type> mat_type;

        Teuchos::oblackholestream blackHole;
        std::ostream& verboseOut = verbose ? out : blackHole;

        verboseOut << "Testing upper triangular solves" << endl;
        //
        // Construct an upper triangular linear system to solve.
        //
        verboseOut << "-- Generating test matrix" << endl;
        const int N = testProblemSize;
        mat_type R (N, N);
        // Fill the upper triangle of R with random numbers.
        for (int j = 0; j < N; ++j) {
          for (int i = 0; i <= j; ++i) {
            R(i,j) = STS::random ();
          }
        }
        // Compute the Frobenius norm of R for later use.
        const magnitude_type R_norm = R.normFrobenius ();
        // Fill the right-hand side B with random numbers.
        mat_type B (N, 1);
        B.random ();
        // Compute the Frobenius norm of B for later use.
        const magnitude_type B_norm = B.normFrobenius ();

        // Save a copy of the original upper triangular system.
        mat_type R_copy (Teuchos::Copy, R, N, N);
        mat_type B_copy (Teuchos::Copy, B, N, 1);

        // Solution vector.
        mat_type X (N, 1);

        // Solve RX = B.
        verboseOut << "-- Solving RX=B" << endl;
        // We're ignoring the return values for now.
        (void) solveUpperTriangularSystem (LEFT_SIDE, X, R, B, robustness);
        // Test the residual error.
        mat_type Resid (N, 1);
        Resid.assign (B_copy);
        Belos::details::LocalDenseMatrixOps<scalar_type> ops;
        ops.matMatMult (STS::one(), Resid, -STS::one(), R_copy, X);
        verboseOut << "---- ||R*X - B||_F = " << Resid.normFrobenius() << endl;
        verboseOut << "---- ||R||_F ||X||_F + ||B||_F = "
                   << (R_norm * X.normFrobenius() + B_norm)
                   << endl;

        // Restore R and B.
        R.assign (R_copy);
        B.assign (B_copy);

        //
        // Set up a right-side test problem: YR = B^*.
        //
        mat_type Y (1, N);
        mat_type B_star (1, N);
        ops.conjugateTranspose (B_star, B);
        mat_type B_star_copy (1, N);
        B_star_copy.assign (B_star);
        // Solve YR = B^*.
        verboseOut << "-- Solving YR=B^*" << endl;
        // We're ignoring the return values for now.
        (void) solveUpperTriangularSystem (RIGHT_SIDE, Y, R, B_star, robustness);
        // Test the residual error.
        mat_type Resid2 (1, N);
        Resid2.assign (B_star_copy);
        ops.matMatMult (STS::one(), Resid2, -STS::one(), Y, R_copy);
        verboseOut << "---- ||Y*R - B^*||_F = " << Resid2.normFrobenius() << endl;
        verboseOut << "---- ||Y||_F ||R||_F + ||B^*||_F = "
                   << (Y.normFrobenius() * R_norm + B_norm)
                   << endl;

        // FIXME (mfh 14 Oct 2011) The test always "passes" for now;
        // you have to inspect the output in order to see whether it
        // succeeded.  We really should fix the above to use the
        // infinity-norm bounds in Higham's book for triangular
        // solves.  That would automate the test.
        return true;
      }

    private:
      //! Stream to which to output warnings.
      std::ostream& warn_;

      /// \brief Default robustness level, for things like triangular solves.
      ///
      /// \note Many methods, including the triangular solve methods,
      ///   let you override the default robustness level.
      ERobustness defaultRobustness_;

    private:
      /// \brief Solve \f$\min_X \|AX - B\|_2\f$ using the SVD.
      ///
      /// \param X [in/out] On input: the right-hand side(s) B.  On
      ///   output: the solution to the least-squares problem.
      ///
      /// \param A [in/out] On input: the matrix in the least-squares
      ///   problem.  On output, A will be overwritten with
      ///   factorization information.
      ///
      /// \return Numerical rank of A, computed using Scalar's machine
      ///   precision as the rank tolerance.
      int
      solveLeastSquaresUsingSVD (mat_type& A, mat_type& X)
      {
        using Teuchos::Array;
        LocalDenseMatrixOps<Scalar> ops;

        if (defaultRobustness_ > ROBUSTNESS_SOME) {
          int count = ops.infNaNCount (A);
          TEUCHOS_TEST_FOR_EXCEPTION(count != 0, std::invalid_argument,
                             "solveLeastSquaresUsingSVD: The input matrix A "
                             "contains " << count << "Inf and/or NaN entr"
                             << (count != 1 ? "ies" : "y") << ".");
          count = ops.infNaNCount (X);
          TEUCHOS_TEST_FOR_EXCEPTION(count != 0, std::invalid_argument,
                             "solveLeastSquaresUsingSVD: The input matrix X "
                             "contains " << count << "Inf and/or NaN entr"
                             << (count != 1 ? "ies" : "y") << ".");
        }
        const int N = std::min (A.numRows(), A.numCols());
        lapack_type lapack;

        // Rank of A; to be computed by _GELSS and returned.
        int rank = N;

        // Use Scalar's machine precision for the rank tolerance,
        // not magnitude_type's machine precision.
        const magnitude_type rankTolerance = STS::eps();

        // Array of singular values.
        Array<magnitude_type> singularValues (N);

        // Extra workspace.  This is only used by _GELSS if Scalar is
        // complex.  Teuchos::LAPACK presents a unified interface to
        // _GELSS that always includes the RWORK argument, even though
        // LAPACK's SGELSS and DGELSS don't have the RWORK argument.
        // We always allocate at least one entry so that &rwork[0]
        // makes sense.
        Array<magnitude_type> rwork (1);
        if (STS::isComplex) {
          rwork.resize (std::max (1, 5 * N));
        }
        //
        // Workspace query
        //
        Scalar lworkScalar = STS::one(); // To be set by workspace query
        int info = 0;
        lapack.GELSS (A.numRows(), A.numCols(), X.numCols(),
                      A.values(), A.stride(), X.values(), X.stride(),
                      &singularValues[0], rankTolerance, &rank,
                      &lworkScalar, -1, &rwork[0], &info);
        TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error,
                           "_GELSS workspace query returned INFO = "
                           << info << " != 0.");
        const int lwork = static_cast<int> (STS::real (lworkScalar));
        TEUCHOS_TEST_FOR_EXCEPTION(lwork < 0, std::logic_error,
                           "_GELSS workspace query returned LWORK = "
                           << lwork << " < 0.");
        // Allocate workspace.  Size > 0 means &work[0] makes sense.
        Array<Scalar> work (std::max (1, lwork));
        // Solve the least-squares problem.
        lapack.GELSS (A.numRows(), A.numCols(), X.numCols(),
                      A.values(), A.stride(), X.values(), X.stride(),
                      &singularValues[0], rankTolerance, &rank,
                      &work[0], lwork, &rwork[0], &info);
        TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
                           "_GELSS returned INFO = " << info << " != 0.");
        return rank;
      }

      /// \brief Solve the projected least-squares problem, assuming
      ///   Givens rotations updates.
      ///
      /// Call this method after invoking either \c updateColumnGivens()
      /// with the same curCol, or \c updateColumnsGivens() with curCol =
      /// endCol.
      void
      solveGivens (mat_type& y,
                   mat_type& R,
                   const mat_type& z,
                   const int curCol)
      {
        const int numRows = curCol + 2;

        // Now that we have the updated R factor of H, and the updated
        // right-hand side z, solve the least-squares problem by
        // solving the upper triangular linear system Ry=z for y.
        const mat_type R_view (Teuchos::View, R, numRows-1, numRows-1);
        const mat_type z_view (Teuchos::View, z, numRows-1, z.numCols());
        mat_type y_view (Teuchos::View, y, numRows-1, y.numCols());

        (void) solveUpperTriangularSystem (Teuchos::LEFT_SIDE, y_view,
                                           R_view, z_view, defaultRobustness_);
      }

      //! Make a random projected least-squares problem.
      void
      makeRandomProblem (mat_type& H, mat_type& z)
      {
        // In GMRES, z always starts out with only the first entry
        // being nonzero.  That entry always has nonnegative real part
        // and zero imaginary part, since it is the initial residual
        // norm.
        H.random ();
        // Zero out the entries below the subdiagonal of H, so that it
        // is upper Hessenberg.
        for (int j = 0; j < H.numCols(); ++j) {
          for (int i = j+2; i < H.numRows(); ++i) {
            H(i,j) = STS::zero();
          }
        }
        // Initialize z, the right-hand side of the least-squares
        // problem.  Make the first entry of z nonzero.
        {
          // It's still possible that a random number will come up
          // zero after 1000 trials, but unlikely.  Nevertheless, it's
          // still important not to allow an infinite loop, for
          // example if the pseudorandom number generator is broken
          // and always returns zero.
          const int numTrials = 1000;
          magnitude_type z_init = STM::zero();
          for (int trial = 0; trial < numTrials && z_init == STM::zero(); ++trial) {
            z_init = STM::random();
          }
          TEUCHOS_TEST_FOR_EXCEPTION(z_init == STM::zero(), std::runtime_error,
                             "After " << numTrials << " trial"
                             << (numTrials != 1 ? "s" : "")
                             << ", we were unable to generate a nonzero pseudo"
                             "random real number.  This most likely indicates a "
                             "broken pseudorandom number generator.");
          const magnitude_type z_first = (z_init < 0) ? -z_init : z_init;

          // NOTE I'm assuming here that "scalar_type = magnitude_type"
          // assignments make sense.
          z(0,0) = z_first;
        }
      }

      /// \brief Compute the Givens rotation corresponding to [x; y].
      ///
      /// The result of applying the rotation is [result; 0].
      void
      computeGivensRotation (const Scalar& x,
                             const Scalar& y,
                             Scalar& theCosine,
                             Scalar& theSine,
                             Scalar& result)
      {
        // _LARTG, an LAPACK aux routine, is slower but more accurate
        // than the BLAS' _ROTG.
        const bool useLartg = false;

        if (useLartg) {
          lapack_type lapack;
          // _LARTG doesn't clobber its input arguments x and y.
          lapack.LARTG (x, y, &theCosine, &theSine, &result);
        } else {
          // _ROTG clobbers its first two arguments.  x is overwritten
          // with the result of applying the Givens rotation: [x; y] ->
          // [x (on output); 0].  y is overwritten with the "fast"
          // Givens transform (see Golub and Van Loan, 3rd ed.).
          Scalar x_temp = x;
          Scalar y_temp = y;
          blas_type blas;
          blas.ROTG (&x_temp, &y_temp, &theCosine, &theSine);
          result = x_temp;
        }
      }

      //! Compute the singular values of A.  Store them in the sigmas array.
      void
      singularValues (const mat_type& A,
                      Teuchos::ArrayView<magnitude_type> sigmas)
      {
        using Teuchos::Array;
        using Teuchos::ArrayView;

        const int numRows = A.numRows();
        const int numCols = A.numCols();
        TEUCHOS_TEST_FOR_EXCEPTION(sigmas.size() < std::min (numRows, numCols),
                           std::invalid_argument,
                           "The sigmas array is only of length " << sigmas.size()
                           << ", but must be of length at least "
                           << std::min (numRows, numCols)
                           << " in order to hold all the singular values of the "
                           "matrix A.");

        // Compute the condition number of the matrix A, using a singular
        // value decomposition (SVD).  LAPACK's SVD routine overwrites the
        // input matrix, so make a copy.
        mat_type A_copy (numRows, numCols);
        A_copy.assign (A);

        // Workspace query.
        lapack_type lapack;
        int info = 0;
        Scalar lworkScalar = STS::zero();
        Array<magnitude_type> rwork (std::max (std::min (numRows, numCols) - 1, 1));
        lapack.GESVD ('N', 'N', numRows, numCols,
                      A_copy.values(), A_copy.stride(), &sigmas[0],
                      (Scalar*) NULL, 1, (Scalar*) NULL, 1,
                      &lworkScalar, -1, &rwork[0], &info);

        TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error,
                           "LAPACK _GESVD workspace query failed with INFO = "
                           << info << ".");
        const int lwork = static_cast<int> (STS::real (lworkScalar));
        TEUCHOS_TEST_FOR_EXCEPTION(lwork < 0, std::logic_error,
                           "LAPACK _GESVD workspace query returned LWORK = "
                           << lwork << " < 0.");
        // Make sure that the workspace array always has positive
        // length, so that &work[0] makes sense.
        Teuchos::Array<Scalar> work (std::max (1, lwork));

        // Compute the singular values of A.
        lapack.GESVD ('N', 'N', numRows, numCols,
                      A_copy.values(), A_copy.stride(), &sigmas[0],
                      (Scalar*) NULL, 1, (Scalar*) NULL, 1,
                      &work[0], lwork, &rwork[0], &info);
        TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error,
                           "LAPACK _GESVD failed with INFO = " << info << ".");
      }

      /// \brief The (largest, smallest) singular values of the given matrix.
      ///
      /// We use these for computing the 2-norm condition number of the
      /// matrix A.  We separate out the singular values rather than
      /// returning their quotient, so that you can see the value of the
      /// largest singular value, even if the smallest singular value is
      /// zero.
      std::pair<magnitude_type, magnitude_type>
      extremeSingularValues (const mat_type& A)
      {
        using Teuchos::Array;

        const int numRows = A.numRows();
        const int numCols = A.numCols();

        Array<magnitude_type> sigmas (std::min (numRows, numCols));
        singularValues (A, sigmas);
        return std::make_pair (sigmas[0], sigmas[std::min(numRows, numCols) - 1]);
      }

      /// \brief Normwise 2-norm condition number of the least-squares problem.
      ///
      /// \param A [in] The matrix in the least-squares problem.
      /// \param b [in] The right-hand side in the least-squares problem.
      ///
      /// \param residualNorm [in] Residual norm from solving the
      ///   least-squares problem using a known good method.
      ///
      /// For details on the condition number formula, see Section 3.3 of
      /// J. W. Demmel, "Applied Numerical Linear Algebra," SIAM Press.
      magnitude_type
      leastSquaresConditionNumber (const mat_type& A,
                                   const mat_type& b,
                                   const magnitude_type& residualNorm)
      {
        // Extreme singular values of A.
        const std::pair<magnitude_type, magnitude_type> sigmaMaxMin =
          extremeSingularValues (A);

        // Our solvers currently assume that H has full rank.  If the
        // test matrix doesn't have full rank, we stop right away.
        TEUCHOS_TEST_FOR_EXCEPTION(sigmaMaxMin.second == STM::zero(), std::runtime_error,
                           "The test matrix is rank deficient; LAPACK's _GESVD "
                           "routine reports that its smallest singular value is "
                           "zero.");
        // 2-norm condition number of A.  We checked above that the
        // denominator is nonzero.
        const magnitude_type A_cond = sigmaMaxMin.first / sigmaMaxMin.second;

        // "Theta" in the variable names below refers to the angle between
        // the vectors b and A*x, where x is the computed solution.  It
        // measures whether the residual norm is large (near ||b||) or
        // small (near 0).
        const magnitude_type sinTheta = residualNorm / b.normFrobenius();

        // \sin^2 \theta + \cos^2 \theta = 1.
        //
        // The range of sine is [-1,1], so squaring it won't overflow.
        // We still have to check whether sinTheta > 1, though.  This
        // is impossible in exact arithmetic, assuming that the
        // least-squares solver worked (b-A*0 = b and x minimizes
        // ||b-A*x||_2, so ||b-A*0||_2 >= ||b-A*x||_2).  However, it
        // might just be possible in floating-point arithmetic.  We're
        // just looking for an estimate, so if sinTheta > 1, we cap it
        // at 1.
        const magnitude_type cosTheta = (sinTheta > STM::one()) ?
          STM::zero() : STM::squareroot (1 - sinTheta * sinTheta);

        // This may result in Inf, if cosTheta is zero.  That's OK; in
        // that case, the condition number of the (full-rank)
        // least-squares problem is rightfully infinite.
        const magnitude_type tanTheta = sinTheta / cosTheta;

        // Condition number for the full-rank least-squares problem.
        return 2 * A_cond / cosTheta + tanTheta * A_cond * A_cond;
      }

      //! \f$\| b - A x \|_2\f$ (Frobenius norm if b has more than one column).
      magnitude_type
      leastSquaresResidualNorm (const mat_type& A,
                                const mat_type& x,
                                const mat_type& b)
      {
        mat_type r (b.numRows(), b.numCols());

        // r := b - A*x
        r.assign (b);
        LocalDenseMatrixOps<Scalar> ops;
        ops.matMatMult (STS::one(), r, -STS::one(), A, x);
        return r.normFrobenius ();
      }

      /// \brief ||x_approx - x_exact||_2 // ||x_exact||_2.
      ///
      /// Use the Frobenius norm if more than one column.
      /// Don't scale if ||x_exact|| == 0.
      magnitude_type
      solutionError (const mat_type& x_approx,
                     const mat_type& x_exact)
      {
        const int numRows = x_exact.numRows();
        const int numCols = x_exact.numCols();

        mat_type x_diff (numRows, numCols);
        for (int j = 0; j < numCols; ++j) {
          for (int i = 0; i < numRows; ++i) {
            x_diff(i,j) = x_exact(i,j) - x_approx(i,j);
          }
        }
        const magnitude_type scalingFactor = x_exact.normFrobenius();

        // If x_exact has zero norm, just use the absolute difference.
        return x_diff.normFrobenius() /
          (scalingFactor == STM::zero() ? STM::one() : scalingFactor);
      }

      /// \brief Update current column using Givens rotations.
      ///
      /// Update the current column of the projected least-squares
      /// problem, using Givens rotations.  This updates the QR
      /// factorization of the upper Hessenberg matrix H.  The
      /// resulting R factor is stored in the matrix R.  The Q factor
      /// is stored implicitly in the list of cosines and sines,
      /// representing the Givens rotations applied to the problem.
      /// These Givens rotations are also applied to the right-hand
      /// side z.
      ///
      /// \return The residual of the resulting least-squares problem,
      ///   assuming that the upper triangular system Ry=z can be
      ///   solved exactly (with zero residual).  (This may not be the
      ///   case if R is singular and the system Ry=z is
      ///   inconsistent.)
      ///
      /// \param H [in] The upper Hessenberg matrix from GMRES.  We only
      ///   view H( 1:curCol+2, 1:curCol+1 ).  It's copied and not
      ///   overwritten, so as not to disturb any backtracking or other
      ///   features of GMRES.
      ///
      /// \param R [in/out] Matrix with the same dimensions as H, used for
      ///   storing the incrementally computed upper triangular factor from
      ///   the QR factorization of H.
      ///
      /// \param y [out] Vector (one-column matrix) with the same number
      ///   of rows as H.  On output: the solution of the projected
      ///   least-squares problem.  The vector should have one more entry
      ///   than necessary for the solution, because of the way we solve
      ///   the least-squares problem.
      ///
      /// \param z [in/out] Vector (one-column matrix) with the same
      ///   number of rows as H.  On input: the current right-hand side of
      ///   the projected least-squares problem.  On output: the updated
      ///   right-hand side.
      ///
      /// \param theCosines [in/out] On input: Array of cosines from the
      ///   previously computed Givens rotations.  On output: the same
      ///   cosines, with the new Givens rotation's cosine appended.
      ///
      /// \param theSines [in/out] On input: Array of sines from the
      ///   previously computed Givens rotations.  On output: the same
      ///   sines, with the new Givens rotation's sine appended.
      ///
      /// \param curCol [in] Index of the current column to update.
      magnitude_type
      updateColumnGivens (const mat_type& H,
                          mat_type& R,
                          mat_type& y,
                          mat_type& z,
                          Teuchos::ArrayView<scalar_type> theCosines,
                          Teuchos::ArrayView<scalar_type> theSines,
                          const int curCol)
      {
        using std::cerr;
        using std::endl;

        const int numRows = curCol + 2; // curCol is zero-based
        const int LDR = R.stride();
        const bool extraDebug = false;

        if (extraDebug) {
          cerr << "updateColumnGivens: curCol = " << curCol << endl;
        }

        // View of H( 1:curCol+1, curCol ) (in Matlab notation, if
        // curCol were a one-based index, as it would be in Matlab but
        // is not here).
        const mat_type H_col (Teuchos::View, H, numRows, 1, 0, curCol);

        // View of R( 1:curCol+1, curCol ) (again, in Matlab notation,
        // if curCol were a one-based index).
        mat_type R_col (Teuchos::View, R, numRows, 1, 0, curCol);

        // 1. Copy the current column from H into R, where it will be
        //    modified.
        R_col.assign (H_col);

        if (extraDebug) {
          printMatrix<Scalar> (cerr, "R_col before ", R_col);
        }

        // 2. Apply all the previous Givens rotations, if any, to the
        //    current column of the matrix.
        blas_type blas;
        for (int j = 0; j < curCol; ++j) {
          // BLAS::ROT really should take "const Scalar*" for these
          // arguments, but it wants a "Scalar*" instead, alas.
          Scalar theCosine = theCosines[j];
          Scalar theSine = theSines[j];

          if (extraDebug) {
            cerr << "  j = " << j << ": (cos,sin) = "
                 << theCosines[j] << "," << theSines[j] << endl;
          }
          blas.ROT (1, &R_col(j,0), LDR, &R_col(j+1,0), LDR,
                    &theCosine, &theSine);
        }
        if (extraDebug && curCol > 0) {
          printMatrix<Scalar> (cerr, "R_col after applying previous "
                               "Givens rotations", R_col);
        }

        // 3. Calculate new Givens rotation for R(curCol, curCol),
        //    R(curCol+1, curCol).
        Scalar theCosine, theSine, result;
        computeGivensRotation (R_col(curCol,0), R_col(curCol+1,0),
                               theCosine, theSine, result);
        theCosines[curCol] = theCosine;
        theSines[curCol] = theSine;

        if (extraDebug) {
          cerr << "  New cos,sin = " << theCosine << "," << theSine << endl;
        }

        // 4. _Apply_ the new Givens rotation.  We don't need to
        //    invoke _ROT here, because computeGivensRotation()
        //    already gives us the result: [x; y] -> [result; 0].
        R_col(curCol, 0) = result;
        R_col(curCol+1, 0) = STS::zero();

        if (extraDebug) {
          printMatrix<Scalar> (cerr, "R_col after applying current "
                               "Givens rotation", R_col);
        }

        // 5. Apply the resulting Givens rotation to z (the right-hand
        //    side of the projected least-squares problem).
        //
        // We prefer overgeneralization to undergeneralization by assuming
        // here that z may have more than one column.
        const int LDZ = z.stride();
        blas.ROT (z.numCols(),
                  &z(curCol,0), LDZ, &z(curCol+1,0), LDZ,
                  &theCosine, &theSine);

        if (extraDebug) {
          //mat_type R_after (Teuchos::View, R, numRows, numRows-1);
          //printMatrix<Scalar> (cerr, "R_after", R_after);
          //mat_type z_after (Teuchos::View, z, numRows, z.numCols());
          printMatrix<Scalar> (cerr, "z_after", z);
        }

        // The last entry of z is the nonzero part of the residual of the
        // least-squares problem.  Its magnitude gives the residual 2-norm
        // of the least-squares problem.
        return STS::magnitude( z(numRows-1, 0) );
      }

      /// \brief Solve the least-squares problem using LAPACK.
      ///
      /// In particular, this method uses LAPACK's _GELS dense
      /// least-squares solver.  _GELS assumes that the input matrix
      /// has full rank, so this method might fail.  This method is
      /// inefficient compared to the Givens rotations approach, but
      /// useful for testing the latter.
      ///
      /// \param H [in] The upper Hessenberg matrix from GMRES.  We only
      ///   view H( 1:curCol+2, 1:curCol+1 ).  It's copied and not
      ///   overwritten, so as not to disturb any backtracking or other
      ///   features of GMRES.
      ///
      /// \param R [out] Matrix with the same dimensions as H.  On
      ///   output: the upper triangle of R contains the R factor in
      ///   the QR factorization of H.
      ///
      /// \param y [out] Vector (one-column matrix) with the same
      ///   number of rows as H.  On output: the solution of the
      ///   projected least-squares problem.  The vector should have
      ///   one more entry than necessary for the solution, because of
      ///   the way we solve the least-squares problem.
      ///
      /// \param z [in] Vector (one-column matrix) with the same
      ///   number of rows as H.  On input: the current right-hand
      ///   side of the projected least-squares problem.
      ///
      /// \param curCol [in] Index of the current column of H.
      ///
      /// \return The 2-norm of the residual of the least-squares
      ///   problem, assuming that it can be solved using LAPACK's
      ///   _GELS least-squares solver.  (This may not be the case if
      ///   H is singular.)
      magnitude_type
      solveLapack (const mat_type& H,
                   mat_type& R,
                   mat_type& y,
                   mat_type& z,
                   const int curCol)
      {
        const int numRows = curCol + 2;
        const int numCols = curCol + 1;
        const int LDR = R.stride();

        // Copy H( 1:curCol+1, 1:curCol ) into R( 1:curCol+1, 1:curCol ).
        const mat_type H_view (Teuchos::View, H, numRows, numCols);
        mat_type R_view (Teuchos::View, R, numRows, numCols);
        R_view.assign (H_view);

        // The LAPACK least-squares solver overwrites the right-hand side
        // vector with the solution, so first copy z into y.
        mat_type y_view (Teuchos::View, y, numRows, y.numCols());
        mat_type z_view (Teuchos::View, z, numRows, y.numCols());
        y_view.assign (z_view);

        // Workspace query for the least-squares routine.
        int info = 0;
        Scalar lworkScalar = STS::zero();
        lapack_type lapack;
        lapack.GELS ('N', numRows, numCols, y_view.numCols(),
                     NULL, LDR, NULL, y_view.stride(),
                     &lworkScalar, -1, &info);
        TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error,
                           "LAPACK _GELS workspace query failed with INFO = "
                           << info << ", for a " << numRows << " x " << numCols
                           << " matrix with " << y_view.numCols()
                           << " right hand side"
                           << ((y_view.numCols() != 1) ? "s" : "") << ".");
        TEUCHOS_TEST_FOR_EXCEPTION(STS::real(lworkScalar) < STM::zero(),
                           std::logic_error,
                           "LAPACK _GELS workspace query returned an LWORK with "
                           "negative real part: LWORK = " << lworkScalar
                           << ".  That should never happen.  Please report this "
                           "to the Belos developers.");
        TEUCHOS_TEST_FOR_EXCEPTION(STS::isComplex && STS::imag(lworkScalar) != STM::zero(),
                           std::logic_error,
                           "LAPACK _GELS workspace query returned an LWORK with "
                           "nonzero imaginary part: LWORK = " << lworkScalar
                           << ".  That should never happen.  Please report this "
                           "to the Belos developers.");
        // Cast workspace from Scalar to int.  Scalar may be complex,
        // hence the request for the real part.  Don't ask for the
        // magnitude, since computing the magnitude may overflow due
        // to squaring and square root to int.  Hopefully LAPACK
        // doesn't ever overflow int this way.
        const int lwork = std::max (1, static_cast<int> (STS::real (lworkScalar)));

        // Allocate workspace for solving the least-squares problem.
        Teuchos::Array<Scalar> work (lwork);

        // Solve the least-squares problem.  The ?: operator prevents
        // accessing the first element of the work array, if it has
        // length zero.
        lapack.GELS ('N', numRows, numCols, y_view.numCols(),
                     R_view.values(), R_view.stride(),
                     y_view.values(), y_view.stride(),
                     (lwork > 0 ? &work[0] : (Scalar*) NULL),
                     lwork, &info);

        TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error,
                           "Solving projected least-squares problem with LAPACK "
                           "_GELS failed with INFO = " << info << ", for a "
                           << numRows << " x " << numCols << " matrix with "
                           << y_view.numCols() << " right hand side"
                           << (y_view.numCols() != 1 ? "s" : "") << ".");
        // Extract the projected least-squares problem's residual error.
        // It's the magnitude of the last entry of y_view on output from
        // LAPACK's least-squares solver.
        return STS::magnitude( y_view(numRows-1, 0) );
      }

      /// \brief Update columns [startCol,endCol] of the projected least-squares problem.
      ///
      /// This method implements a "left-looking panel QR
      /// factorization" of the upper Hessenberg matrix in the
      /// projected least-squares problem.  It's "left-looking"
      /// because we don't updating anything to the right of columns
      /// [startCol, endCol], which is the "panel."
      ///
      /// \return 2-norm of the absolute residual of the projected
      ///   least-squares problem.
      magnitude_type
      updateColumnsGivens (const mat_type& H,
                           mat_type& R,
                           mat_type& y,
                           mat_type& z,
                           Teuchos::ArrayView<scalar_type> theCosines,
                           Teuchos::ArrayView<scalar_type> theSines,
                           const int startCol,
                           const int endCol)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(startCol > endCol, std::invalid_argument,
                           "updateColumnGivens: startCol = " << startCol
                           << " > endCol = " << endCol << ".");
        magnitude_type lastResult = STM::zero();
        // [startCol, endCol] is an inclusive range.
        for (int curCol = startCol; curCol <= endCol; ++curCol) {
          lastResult = updateColumnGivens (H, R, y, z, theCosines, theSines, curCol);
        }
        return lastResult;
      }

      /// \brief Update columns [startCol,endCol] of the projected least-squares problem.
      ///
      /// This method implements a "left-looking panel QR factorization"
      /// of the upper Hessenberg matrix in the projected least-squares
      /// problem.  It's "left-looking" because we don't updating anything
      /// to the right of columns [startCol, endCol], which is the
      /// "panel."
      ///
      /// \return 2-norm of the absolute residual of the projected
      ///   least-squares problem.
      ///
      /// \warning This method doesn't work!
      magnitude_type
      updateColumnsGivensBlock (const mat_type& H,
                                mat_type& R,
                                mat_type& y,
                                mat_type& z,
                                Teuchos::ArrayView<scalar_type> theCosines,
                                Teuchos::ArrayView<scalar_type> theSines,
                                const int startCol,
                                const int endCol)
      {
        const int numRows = endCol + 2;
        const int numColsToUpdate = endCol - startCol + 1;
        const int LDR = R.stride();

        // 1. Copy columns [startCol, endCol] from H into R, where they
        //    will be modified.
        {
          const mat_type H_view (Teuchos::View, H, numRows, numColsToUpdate, 0, startCol);
          mat_type R_view (Teuchos::View, R, numRows, numColsToUpdate, 0, startCol);
          R_view.assign (H_view);
        }

        // 2. Apply all the previous Givens rotations, if any, to
        //    columns [startCol, endCol] of the matrix.  (Remember
        //    that we're using a left-looking QR factorization
        //    approach; we haven't yet touched those columns.)
        blas_type blas;
        for (int j = 0; j < startCol; ++j) {
          blas.ROT (numColsToUpdate,
                    &R(j, startCol), LDR, &R(j+1, startCol), LDR,
                    &theCosines[j], &theSines[j]);
        }

        // 3. Update each column in turn of columns [startCol, endCol].
        for (int curCol = startCol; curCol < endCol; ++curCol) {
          // a. Apply the Givens rotations computed in previous
          //    iterations of this loop to the current column of R.
          for (int j = startCol; j < curCol; ++j) {
            blas.ROT (1, &R(j, curCol), LDR, &R(j+1, curCol), LDR,
                      &theCosines[j], &theSines[j]);
          }
          // b. Calculate new Givens rotation for R(curCol, curCol),
          //    R(curCol+1, curCol).
          Scalar theCosine, theSine, result;
          computeGivensRotation (R(curCol, curCol), R(curCol+1, curCol),
                                 theCosine, theSine, result);
          theCosines[curCol] = theCosine;
          theSines[curCol] = theSine;

          // c. _Apply_ the new Givens rotation.  We don't need to
          //    invoke _ROT here, because computeGivensRotation()
          //    already gives us the result: [x; y] -> [result; 0].
          R(curCol+1, curCol) = result;
          R(curCol+1, curCol) = STS::zero();

          // d. Apply the resulting Givens rotation to z (the right-hand
          //    side of the projected least-squares problem).
          //
          // We prefer overgeneralization to undergeneralization by
          // assuming here that z may have more than one column.
          const int LDZ = z.stride();
          blas.ROT (z.numCols(),
                    &z(curCol,0), LDZ, &z(curCol+1,0), LDZ,
                    &theCosine, &theSine);
        }

        // The last entry of z is the nonzero part of the residual of the
        // least-squares problem.  Its magnitude gives the residual 2-norm
        // of the least-squares problem.
        return STS::magnitude( z(numRows-1, 0) );
      }
    }; // class ProjectedLeastSquaresSolver
  } // namespace details
} // namespace Belos

#endif // __Belos_ProjectedLeastSquaresSolver_hpp
