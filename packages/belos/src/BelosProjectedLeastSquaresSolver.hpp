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
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

/// \file BelosProjectedLeastSquaresSolver.hpp 
/// \brief Methods for solving GMRES' projected least-squares problem.
/// \author Mark Hoemmen

namespace Belos {

/// \namespace details
/// \brief Namespace containing implementation details of Belos solvers.
/// 
/// \warning Belos users should not use anything in this namespace.
///   They should not even assume that the namespace will continue to
///   exist between releases.  The namespace's name itself or anything
///   it contains may change at any time.
namespace details {

/// \struct ProjectedLeastSquaresProblem
/// \brief "Container" for the data representing the projected least-squares problem.
/// \author Mark Hoemmen
///
template<class Scalar>
struct ProjectedLeastSquaresProblem {
  Teuchos::SerialDenseMatrix<int,Scalar> H;
  Teuchos::SerialDenseMatrix<int,Scalar> R;
  Teuchos::SerialDenseMatrix<int,Scalar> y;
  Teuchos::SerialDenseMatrix<int,Scalar> z;
  Teuchos::Array<Scalar> theCosines;
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
    theCosines (maxNumIterations),
    theSines (maxNumIterations)
  {}
};


/// \class ProjectedLeastSquaresSolver
/// \brief Methods for solving GMRES' projected least-squares problem.
/// \author Mark Hoemmen
///
/// \tparam Scalar The type of the matrix and vector entries in the
///   least-squares problem.
///
template<class Scalar>
class ProjectedLeastSquaresSolver {
public:
  typedef Scalar scalar_type;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
  typedef Teuchos::SerialDenseMatrix<int,Scalar> mat_type;

private:
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;
  typedef Teuchos::BLAS<int, scalar_type> blas_type;
  typedef Teuchos::LAPACK<int, scalar_type> lapack_type;

  /// \brief The (largest, smallest) singular values of the given matrix.
  ///
  /// We use these for computing the 2-norm condition number of the
  /// matrix A.  We separate out the singular values rather than
  /// returning their quotient, so that you can see the value of the
  /// largest singular value, even if the smallest singular value is
  /// zero.
  std::pair<magnitude_type, magnitude_type>
  extremeSingularValues (const mat_type& A) const
  {
    const int numRows = A.numRows();
    const int numCols = A.numCols();

    // Compute the condition number of the matrix A, using a singular
    // value decomposition (SVD).  LAPACK's SVD routine overwrites the
    // input matrix, so make a copy.
    mat_type A_copy (numRows, numCols);
    A_copy.assign (A);

    // Workspace query.
    int info = 0;
    Scalar lworkScalar = STS::zero();
    lapack.GESVD ('N', 'N', numRows, numCols, NULL, A_copy.stride(), NULL,
		  NULL, 1, NULL, 1, &lworkScalar, -1, &info);
    TEST_FOR_EXCEPTION(info != 0, std::logic_error,
		       "LAPACK _GESVD workspace query failed with INFO = " 
		       << info << ".");
    const int lwork = static_cast<int> (STS::real (lworkScalar));
    TEST_FOR_EXCEPTION(lwork < 0, std::logic_error,
		       "LAPACK _GESVD workspace query returned LWORK = " 
		       << lwork << " < 0.");
    Teuchos::Array<Scalar> work (lwork);

    // Compute the singular values of A.
    Teuchos::Array<magnitude_type> singularValues (numCols);
    lapack.GESVD ('N', 'N', numRows, numCols, A_copy.values(), A_copy.stride(),
		  &singularValues[0], NULL, 1, NULL, 1, 
		  (lwork > 0 ? &work[0] : NULL), lwork, &info);
    TEST_FOR_EXCEPTION(info != 0, std::logic_error,
		       "LAPACK _GESVD failed with INFO = " << info << ".");
    return std::make_pair (singularValues[0], singularValues[numCols-1]);
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
			       const magnitude_type& residualNorm) const
  {
    // Extreme singular values of A.
    const std::pair<magnitude_type, magnitude_type> sigmaMaxMin = 
      extremeSingularValues (A);

    // Our solvers currently assume that H has full rank.  If the
    // test matrix doesn't have full rank, we stop right away.
    TEST_FOR_EXCEPTION(sigmaMaxMin.second == STM::zero(), std::runtime_error,
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
    const magnitude_type cosTheta = STM::squareroot (1 - sinTheta * sinTheta);

    // This may result in Inf, if cosTheta is zero.  That's OK; in
    // that case, the condition number of the (full-rank)
    // least-squares problem is rightfully infinite.
    const magnitude_type tanTheta = sinTheta / cosTheta;

    // Condition number for the full-rank least-squares problem.
    return 2 * A_cond / cosTheta + tanTheta * A_cond * A_cond;
  }

public:
  /// \brief Update column curCol of the projected least-squares problem.
  ///
  /// \param problem [in/out] The projected least-squares problem.
  /// \param curCol [in] Zero-based index of the current column to update.
  ///
  /// \return 2-norm of the absolute residual of the projected
  ///   least-squares problem.
  magnitude_type
  updateColumn (ProjectedLeastSquaresProblem<Scalar>& problem,
		const int curCol) const
  {
    return updateColumnGivens (problem.H, problem.R, problem.y, problem.z,
			       problem.theCosines, problem.theSines, curCol);
  }

  /// \brief Test \c updateColumnGivens() against \c updateColumnLapack().
  ///
  /// \param out [out] Output stream to which to print results.
  ///
  /// \param numCols [in] Number of columns in the projected
  ///   least-squares problem to test.  (The number of rows is one
  ///   plus the number of columns.)
  ///
  /// \return Whether the least-squares solution error is within the
  ///   expected bound.
  ///
  /// We test updating the least-squares problem by generating a
  /// random least-squares problem that looks like it comes from
  /// GMRES.  The matrix is upper Hessenberg, and the right-hand side
  /// starts out with the first entry being nonzero with nonnegative
  /// real part and zero imaginary part, and all the other entries
  /// being zero.  Then, we compare the results of \c
  /// updateColumnGivens() (applied to each column in turn) with the
  /// results of \c updateColumnLapack() (again, applied to each
  /// column in turn).  We print out the results to the given output
  /// stream.
  bool
  testUpdateColumn (std::ostream& out, const int numCols) const
  {
    using Teuchos::Array;

    TEST_FOR_EXCEPTION(numCols <= 0, std::invalid_argument,
		       "numCols = " << numCols << " <= 0.");
    const int numRows = numCols + 1;

    mat_type H (numRows, numCols);
    mat_type z (numRows, numCols);

    mat_type R_givens (numRows, numCols);
    mat_type y_givens (numRows, 1);
    mat_type z_givens (numRows, 1);
    Array<Scalar> theCosines (numCols);
    Array<Scalar> theSines (numCols);

    mat_type R_lapack (numRows, numCols);
    mat_type y_lapack (numRows, 1);
    mat_type z_lapack (numRows, 1);

    // Make a random least-squares problem.  In GMRES, z always starts
    // out with only the first entry being nonzero.  That entry always
    // has nonnegative real part and zero imaginary part, since it is
    // the initial residual norm.
    H.random ();
    {
      const magnitude_type z_init = STM::random ();
      const magnitude_type z_first = (z_init < 0) ? -z_init : z_init;

      // NOTE I'm assuming here that "scalar_type = magnitude_type"
      // assignments make sense.
      z(0,0) = z_first;
      z_givens(0,0) = z_first;
      z_lapack(0,0) = z_first;
    }

    //
    // Imitate how one would update the least-squares problem in a
    // typical GMRES implementation, for each updating method.
    //
    // Update using Givens rotations.
    magnitude_type residualNormGivens = STM::zero();
    for (int curCol = 0; curCol < numCols; ++curCol) {
      residualNormGivens = updateColumnGivens (H, R_givens, y_givens, z_givens, 
					       theCosines, theSines, curCol);
    }
    // Update using LAPACK's least-squares solver.
    magnitude_type residualNormLapack = STM::zero();
    for (int curCol = 0; curCol < numCols; ++curCol) {
      residualNormLapack = updateColumnLapack (H, R_lapack, y_lapack, z_lapack, 
					       curCol);
    }

    // Compute the condition number of the least-squares problem,
    // using the residual from the LAPACK method.
    const magnitude_type leastSquaresCondNum = 
      leastSquaresConditionNumber (H, z, residualNormLapack);

    // Compute the relative least-squares solution error.  We assume
    // that the LAPACK solution is "exact" and compare against the
    // Givens rotations solution.  This is taking liberties with the
    // definition of condition number, but it's the best we can do,
    // since we don't know the exact solution.
    magnitude_type leastSquaresSolutionError = STM::zero();
    {
      mat_type y (numCols, 1);
      for (int i = 0; i < numCols; ++i) {
	y(i,0) = y_givens(i,0) - y_lapack(i,0);
      }
      const magnitude_type scalingFactor = y_lapack.normFrobenius();

      // If y_lapack has zero norm, just use the absolute difference.
      leastSquaresSolutionError = y.normFrobenius() / 
	(scalingFactor == STM::zero() ? STM::one() : scalingFactor);
    }    

    // Compute the normwise difference between the two computed R factors.
    magnitude_type R_factorDeviation = STM::zero();
    {
      // Copy out the R factors.  R_lapack, if not also R_givens,
      // contains some additional stuff below the upper triangle, so
      // we copy out the upper triangle (including the diagonal).
      mat_type R_factorFromGivens (numCols, numCols);
      mat_type R_factorFromLapack (numCols, numCols);
      for (int j = 0; j < numCols; ++j) {
	for (int i = 0; i <= j; ++i) {
	  R_factorFromGivens(i,j) = R_givens(i,j);
	  R_factorFromLapack(i,j) = R_lapack(i,j);
	}
      }
      mat_type R_diff (numCols, numCols);
      for (int j = 0; j < numCols; ++j) {
	for (int i = 0; i < numCols; ++i) {
	  R_diff(i,j) = R_factorFromGivens(i,j) - R_factorFromLapack(i,j);
	}
      } 
      R_factorDeviation = R_diff.normFrobenius();
    }

    // Compute the (Frobenius) norm of the original matrix H.
    const magnitude_type H_norm = H.normFrobenius();

    out << "||R_givens - R_lapack||_F = " << R_factorDeviation << endl;
    out << "||H||_F = " << H_norm << endl;
    out << "||y_givens - y_lapack||_2 / ||y_lapack||_2 = " << leastSquaresSolutionError << endl;
    out << "Least-squares condition number = " << leastSquaresCondNum << endl;

    // Now for the controversial part of the test: judging whether we
    // succeeded.  This includes the problem's condition number, which
    // is a measure of the maximum perturbation in the solution for
    // which we can still say that the solution is valid.  We include
    // a little wiggle room, computing by the conventional
    // Wilkinsonian square-root root, times 10 for good measure.
    // (This doesn't really have anything to do with probability.
    // What's an "average problem"?  However, you can explain it to
    // yourself with a standard deviation argument.)
    const magnitude_type solutionErrorBound = 
      10 * STM::squareroot( numRows*numCols ) * leastSquaresCondNum;
    out << "Solution error bound: 10*sqrt(numRows*numCols)*cond = " 
	<< solutionErrorBound << endl;

    if (leastSquaresSolutionError > solutionErrorBound)
      return false;
    else
      return true;
  }

private:
  /// \brief Update current column using Givens rotations.
  ///
  /// Update the current column of the projected least-squares
  /// problem, using Givens rotations.
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
		      const int curCol) const
  {
    // View of H( 1:curCol+1, curCol ) (in Matlab notation).
    const mat_type H_col (Teuchos::View, H, curCol+2, 1, 0, curCol);

    // View of R( 1:curCol+1, curCol ).
    mat_type R_col (Teuchos::View, R, curCol+2, 1, 0, curCol);

    // Copy the current column of the upper Hessenberg matrix H into
    // the corresponding column of R.
    R_col.assign (H_col);

    // Apply the previous Givens rotations, if any, to the current
    // column.  (This happens in R, not in H; we keep around the
    // original H for various reasons, e.g., for safe backtracking.)
    blas_type blas;
    for (int j = 0; j < curCol; ++j) {
      blas.ROT (1, &R_col(j,0), R_col.stride(), &R_col(j+1,0), R_col.stride());
    }

    // Calculate new Givens rotation for [R(j,j), R(j+1, j)], 
    // where j = curCol.
    magnitude_type theCosine;
    scalar_type theSine;
    blas.ROTG (&R_col(curCol,0), &R_col(curCol+1,0), &theCosine, &theSine);
    theCosines_[curCol] = theCosine;
    theSines_[curCol] = theSine;

    // Clear the subdiagonal element R(j+1,j)
    R_col(curCol+1,0) = STS::zero();

    // Apply the resulting Givens rotation to z (the right-hand side
    // of the projected least-squares problem).
    blas.ROT (1, &z(curCol,0), z.stride(), &z(curCol+1,0), z.stride(),
	      &theCosine, &theSine);

    // Now that we have the updated R factor of H, and the updated
    // right-hand side z, solve the least-squares problem by solving
    // the upper triangular linear system Ry=z for y.
    // 
    // TRSM() overwrites the right-hand side vector with the solution,
    // so first copy z (the right-hand side vector) into y.
    {
      const mat_type z_view (Teuchos::View, z, curCol+1, 1);
      mat_type y_view (Teuchos::View, y, curCol+1, 1);
      y_view.assign (z_view);
    }
    // Solve Ry = z.  (Only use the first curCol+1 entries of z.)
    const scalar_type one = STS::one();
    blas.TRSM(Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
	      Teuchos::NON_UNIT_DIAG, curCol+1, curCol+1,
	      one, R->values(), R->stride(), y->values(), y->stride());

    // The last entry of z is the nonzero part of the residual of the
    // least-squares problem.  Its magnitude gives the residual 2-norm
    // of the least-squares problem.
    return STS::magnitude( z(curCol+1, 0) );
  }

  /// \brief Update current column using LAPACK least-squares solver.
  ///
  /// Update the current column of the projected least-squares
  /// problem, using LAPACK's _GELS least-squares solver.  This is
  /// inefficient, but useful for testing \c updateColumnGivens().
  /// See that method's documentation for an explanation of the
  /// arguments.
  magnitude_type
  updateColumnLapack (const mat_type& H,
		      mat_type& R,
		      mat_type& y,
		      mat_type& z,
		      const int curCol) const
  {
    // Copy H( 1:curCol+1, 1:curCol ) into R( 1:curCol+1, 1:curCol ).
    const mat_type H_view (Teuchos::View, H, curCol+2, curCol+1, 0, 0);
    const mat_type R_view (Teuchos::View, R, curCol+2, curCol+1, 0, 0);
    R_view.assign (H_view);

    // The LAPACK least-squares solver overwrites the right-hand side
    // vector with the solution, so first copy z into y.
    mat_type y_view (Teuchos::View, y, curCol+2, 1);
    mat_type z_view (Teuchos::View, z, curCol+2, 1);
    y_view.assign (z_view);

    // Workspace query for the least-squares routine.
    int info = 0;
    Scalar lworkScalar = STS::zero();
    lapack.GELS ('N', curCol+2, curCol+1, 1, 
		 NULL, R_view.stride(), NULL, y_view.stride(), 
		 &lworkScalar, -1, &info);
    
    TEST_FOR_EXCEPTION(info != 0, std::logic_error,
		       "LAPACK _GELS workspace query failed with INFO = " 
		       << info << ", for a " << (curCol+2) << " x " 
		       << (curCol+1) << " matrix with 1 right-hand side.");
    TEST_FOR_EXCEPTION(STS::isComplex && STS::imag(lworkScalar) != STM::zero(),
		       std::logic_error,
		       "LAPACK _GELS workspace query returned an LWORK with "
		       "nonzero imaginary part.  That's totally weird and "
		       "should never happen.  It could be that we Belos "
		       "developers made a programming error; please report "
		       "this to us.");

    // Cast workspace from Scalar (which may be complex, hence the
    // request for the real part -- don't call magnitude() since it
    // may overflow due to squaring and square root) to int.
    // Hopefully LAPACK doesn't ever overflow int this way.
    const int lwork = static_cast<int> (STS::real (lworkScalar));

    TEST_FOR_EXCEPTION(lwork < 0, std::logic_error,
		       "LAPACK _GELS workspace query returned LWORK = " 
		       << lwork << " < 0.  This should never happen.");

    Teuchos::Array<Scalar> work (lwork, STS::zero()); // Workspace

    // Solve the least-squares problem.  The ?: operator prevents
    // accessing the first element of the work array, if it has length
    // zero.
    lapack.GELS ('N', curCol+2, curCol+1, 1, 
		 R_view.values(), R_view.stride(),
		 y_view.values(), y_view.stride(),
		 (lwork > 0 ? &work[0] : (Scalar*) NULL), 
		 lwork, &info);
    TEST_FOR_EXCEPTION(info != 0, std::logic_error,
		       "Solving projected least-squares problem with LAPACK "
		       "_GELS failed with INFO = " << info << ", for a " 
		       << (curCol+2) << " x " << (curCol+1) 
		       << " matrix with 1 right-hand side.");

    // Extract the projected least-squares problem's residual error.
    // It's the magnitude of the last entry of y_view on output from
    // LAPACK's least-squares solver.  
    return STS::magnitude( y_view(curCol+2, 0) );
  }



}; // class ProjectedLeastSquaresSolver

} // namespace Belos

#endif // __Belos_ProjectedLeastSquaresSolver_hpp
