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
    } // namespace (anonymous)

    /// \struct ProjectedLeastSquaresProblem
    /// \brief "Container" for the data representing the projected least-squares problem.
    /// \author Mark Hoemmen
    ///
    template<class Scalar>
    struct ProjectedLeastSquaresProblem {
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
      Teuchos::SerialDenseMatrix<int,Scalar> R;

      /// \brief Current solution of the projected least-squares problem.
      ///
      /// The vector (matrix with one column) y has the same number of
      /// rows as H.  It is used to store the solution of the
      /// projected least-squares problem at each step.  The vector
      /// should have one more entry than necessary for the solution,
      /// because of the way we solve the least-squares problem.
      Teuchos::SerialDenseMatrix<int,Scalar> y;

      /// \brief Current right-hand side of the projected least-squares problem.
      ///
      /// The vector (one-column matrix) z has the same number of rows
      /// as H.  It stores the current right-hand side of the
      /// projected least-squares problem, which may be updated
      /// either progressively (if a Givens rotation method is used)
      /// or all at once (if an LAPACK factorization method is
      /// used).
      Teuchos::SerialDenseMatrix<int,Scalar> z;

      /// \brief Array of cosines from the computed Givens rotations.
      ///
      /// This array is only filled in if a Givens rotation method is
      /// used for updating the least-squares problem.
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
    };
    
    /// \class ProjectedLeastSquaresSolver
    /// \brief Methods for solving GMRES' projected least-squares problem.
    /// \author Mark Hoemmen
    ///
    /// \tparam Scalar The type of the matrix and vector entries in the
    ///   least-squares problem.
    ///
    /// Purposes of this class:
    /// 1. Isolate and factor out BLAS and LAPACK dependencies.  This
    ///    makes it easier to write custom replacements for routines
    ///    for which no Scalar specialization is available.
    /// 2. Encapsulate common functionality of many GMRES-like
    ///    solvers.  Avoid duplicated code and simplify debugging,
    ///    testing, and implementation of new GMRES-like solvers.
    /// 3. Provide more robust implementations of solvers for the
    ///    projected least-squares problem.
    ///
    /// "Robust" here means regularizing the least-squares solve, so
    /// that the solution is well-defined even if the problem is
    /// ill-conditioned.  Many distributed-memory iterative solvers,
    /// including Belos, currently solve the projected least-squares
    /// problem redundantly on different processes.  If those
    /// processes are heterogeneous or implement the BLAS and LAPACK
    /// themselves in parallel (via multithreading, for example), then
    /// different BLAS or LAPACK calls on different processes may
    /// result in different answers.  The answers may be significantly
    /// different if the projected problem is singular or
    /// ill-conditioned.  This is bad because iterative methods assume
    /// that the projected problem has the same solution on all
    /// processes.  Regularizing the projected problem is one way to
    /// ensure that different processes compute (almost) the same
    /// solution.
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

    public:
      //! Constructor.
      ProjectedLeastSquaresSolver () {}

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

      /// \brief Update columns [startCol,endCol] of the projected least-squares problem.
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
		     const int endCol) const
      {
	return updateColumnsGivens (problem.H, problem.R, problem.y, problem.z,
				    problem.theCosines, problem.theSines, 
				    startCol, endCol);
      }

      /// \brief Test Givens rotations.
      ///
      /// This routine tests both computing Givens rotations (via \c
      /// computeGivensRotation()) and applying them.
      ///
      /// \param out [out] Stream to which to write test output.
      ///
      void
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
      }

      /// \brief Test \c updateColumnGivens() against \c updateColumnLapack().
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
      /// \param verbose [in] Whether to print verbose output (e.g.,
      ///   the test problem and results).
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
      testUpdateColumn (std::ostream& out, 
			const int numCols,
			const bool testBlockGivens=false,
			const bool verbose=false) const
      {
	using Teuchos::Array;
	using std::endl;
	
	TEST_FOR_EXCEPTION(numCols <= 0, std::invalid_argument,
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
	if (verbose) {
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
	// Update using the "panel left-looking" Givens approach, with
	// the given panel width.
	magnitude_type residualNormBlockGivens = STM::zero();
	if (testBlockGivens) {
	  for (int startCol = 0; startCol < numCols; startCol += panelWidth) {
	    int endCol = std::min (startCol + panelWidth - 1, numCols - 1);
	    residualNormBlockGivens = 
	      updateColumnsGivens (H, R_blockGivens, y_blockGivens, z_blockGivens, 
				   blockCosines, blockSines, startCol, endCol);
	  }
	}
	// Update using LAPACK's least-squares solver.
	magnitude_type residualNormLapack = STM::zero();
	for (int curCol = 0; curCol < numCols; ++curCol) {
	  residualNormLapack = updateColumnLapack (H, R_lapack, y_lapack, z_lapack, 
						   curCol);
	}

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
	if (verbose) {
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

	if (givensSolutionError > solutionErrorBound)
	  return false;
	else if (testBlockGivens && blockGivensSolutionError > solutionErrorBound)
	  return false;
	else
	  return true;
      }

    private:
      //! Make a random projected least-squares problem.
      void
      makeRandomProblem (mat_type& H, mat_type& z) const
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
	  TEST_FOR_EXCEPTION(z_init == STM::zero(), std::runtime_error,
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
			     Scalar& result) const
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
		      Teuchos::ArrayView<magnitude_type> sigmas) const
      {
	using Teuchos::Array;
	using Teuchos::ArrayView;

	const int numRows = A.numRows();
	const int numCols = A.numCols();
	TEST_FOR_EXCEPTION(sigmas.size() < std::min (numRows, numCols),
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

	TEST_FOR_EXCEPTION(info != 0, std::logic_error,
			   "LAPACK _GESVD workspace query failed with INFO = " 
			   << info << ".");
	const int lwork = static_cast<int> (STS::real (lworkScalar));
	TEST_FOR_EXCEPTION(lwork < 0, std::logic_error,
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
	TEST_FOR_EXCEPTION(info != 0, std::logic_error,
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
      extremeSingularValues (const mat_type& A) const
      {
	using Teuchos::Array;

	const int numRows = A.numRows();
	const int numCols = A.numCols();

	Array<magnitude_type> sigmas (std::min (numRows, numCols));
	singularValues (A, sigmas);
	return std::make_pair (sigmas[0], sigmas[std::min(numRows, numCols) - 1]);
      }

      
      /// \brief Solve the square upper triangular linear system \f$Rx = b\f$.
      ///
      /// This method uses the number of columns of R as the dimension
      /// of the linear system, so R may have more rows than columns;
      /// we just won't use the "extra" rows in the solve.
      ///
      /// \return (detectedRank, foundRankDeficiency).
      std::pair<int, bool>
      solveUpperTriangularSystem (mat_type& x,
				  const mat_type& R,
				  const mat_type& b,
				  const int robustness=0) const
      {
	using Teuchos::Array;

	const int M = R.numRows();
	const int N = R.numCols();
	const int LDR = R.stride();
	const int NRHS = x.numCols();
	const int LDX = x.stride();

	// If b has more columns than x, we ignore the remaining
	// columns of b when solving the upper triangular system.  If
	// b has _fewer_ columns than x, we can't solve for all the
	// columns of x, so we throw an exception.
	TEST_FOR_EXCEPTION(NRHS > b.numCols(), std::invalid_argument,
			   "The solution vector x has more columns than the "
			   "right-hand side vector b.  x has " << x.numCols() 
			   << " columns and b has " << b.numCols() 
			   << " columns.");
	TEST_FOR_EXCEPTION(b.numRows() < N, std::invalid_argument,
			   "The right-hand side vector b has only " 
			   << b.numRows() << " rows, but needs at least " 
			   << N << " rows to match the matrix.");
	TEST_FOR_EXCEPTION(x.numRows() < N, std::invalid_argument,
			   "The solution vector x has only " << x.numRows() 
			   << " rows, but needs at least " << N 
			   << " rows to match the matrix.");
	TEST_FOR_EXCEPTION(M < N, std::invalid_argument,
			   "R is " << M << " x " << N << ", but "
			   "solveUpperTriangularSystem needs R to have at "
			   "least as many rows as columns.");
	TEST_FOR_EXCEPTION(M < N, std::invalid_argument,
			   "R is " << M << " x " << N << ", but "
			   "solveUpperTriangularSystem needs R to have at "
			   "least as many rows as columns.");
	TEST_FOR_EXCEPTION(robustness < 0 || robustness > 2, 
			   std::invalid_argument,
			   "Invalid robustness value " << robustness << ".");
	lapack_type lapack;
	blas_type blas;
	int detectedRank = N;
	bool foundRankDeficiency = false;

	// Both the BLAS' _TRSM and LAPACK's _LATRS overwrite the
	// right-hand side with the solution, so first copy b into x.
	if (x.numCols() == b.numCols()) {
	  x.assign (b);
	} else {
	  mat_type b_view (Teuchos::View, b, b.numRows(), NRHS);
	  x.assign (b_view);
	}

	// Solve Rx = b.
	if (robustness == 0) 
	  { // Fast BLAS triangular solve.  No rank checks.
	    blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
		      Teuchos::NON_UNIT_DIAG, N, NRHS, STS::one(), R.values(),
		      LDR, x.values(), LDX);
	  }
	else if (robustness == 1) 
	  { // Robust triangular solve using LAPACK's LATRS routine.
	    // Rudimentary rank detection, using diagonal entries of R
	    // and the norms of the off-diagonal entries of each
	    // column as computed by LATRS.
	    Array<magnitude_type> cnorm (N);
	    magnitude_type scaleFactor = STM::one();
	    char NORMIN = 'N';
	    int info = 0;

	    for (int j = 0; j < x.numCols(); ++j) {
	      // _LATRS only solves from the left side, and only one
	      // column at a time.
	      lapack.LATRS ('U', 'N', 'N', NORMIN, N, R.values(), LDR,
			    &x(0,j), &scaleFactor, &cnorm[0], &info);
	      TEST_FOR_EXCEPTION(info != 0, std::logic_error,
				 "LAPACK's _LATRS routine returned INFO = " 
				 << info << " != 0.");
	      // LATRS computes cnorm if NORMIN='N'.  We don't need to
	      // compute them again, so tell LATRS to reuse cnorm the
	      // next time around.
	      NORMIN = 'Y'; 

	      if (scaleFactor == STM::zero()) {
		// LATRS doesn't tell us the actual rank, just that the
		// matrix is either rank-deficient or badly scaled.
		foundRankDeficiency = true;
	      }
	    }
	    // However, _LATRS _does_ return the 1-norms of the
	    // off-diagonal parts of the columns of R (in the cnorm
	    // array).  We can at least use this to detect zero
	    // columns of R.
	    int rank = N;
	    for (int j = 0; j < N; ++j) {
	      if (R(j,j) == STS::zero() && (j == 0 || cnorm[j] == STM::zero())) {
		--rank;
	      }
	    }
	    if (rank < N) {
	      foundRankDeficiency = true;
	    }
	    detectedRank = rank;
	  } 
	else if (robustness == 2) 
	  { // Find the minimum-norm solution to the least-squares
	    // problem $\min_x \|Rx-b\|_2$, using the singular value
	    // decomposition (SVD).
	    //
	    // _GELSS overwrites its matrix input, so make a copy.
	    mat_type R_copy (Teuchos::Copy, R, N, N);
	    Array<magnitude_type> singularValues (N);
	    int rank = N; // to be set by _GELSS

	    // Use Scalar's machine precision for the rank tolerance,
	    // not magnitude_type's machine precision.
	    const magnitude_type rankTolerance = STS::eps();

	    // Extra workspace.  This is only used if Scalar is
	    // complex, by CGELSS or ZGELSS.  Teuchos::LAPACK presents
	    // a unified interface to _GELSS that always includes the
	    // RWORK argument, even though SGELSS and DGELSS don't
	    // have the RWORK argument.  We always allocate at least
	    // one entry so that &rwork[0] makes sense.
	    Array<magnitude_type> rwork (1);
	    if (STS::isComplex) {
	      rwork.resize (std::max (1, 5 * N));
	    }

	    //
	    // Workspace query
	    //
	    Scalar lworkScalar = STS::one(); // To be set by workspace query
	    int info = 0;
	    lapack.GELSS (N, N, NRHS, R_copy.values(), LDR, x.values(), LDX,
			  &singularValues[0], rankTolerance, &rank, 
			  &lworkScalar, -1, &rwork[0], &info);
	    TEST_FOR_EXCEPTION(info != 0, std::logic_error,
			       "_GELSS workspace query returned INFO = " 
			       << info << " != 0.");
	    const int lwork = static_cast<int> (STS::real (lworkScalar));
	    TEST_FOR_EXCEPTION(lwork < 0, std::logic_error,
			       "_GELSS workspace query returned LWORK = " 
			       << lwork << " < 0.");
	    // Allocate workspace.  Size > 0 means &work[0] makes sense.
	    Array<Scalar> work (std::max (1, lwork));
	    // Solve the least-squares problem.
	    lapack.GELSS (N, N, NRHS, R_copy.values(), LDR, x.values(), LDX,
			  &singularValues[0], rankTolerance, &rank,
			  &lworkScalar, -1, &rwork[0], &info);
	    TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
			       "_GELSS returned INFO = " << info << " != 0.");
	    if (rank < N) {
	      foundRankDeficiency = true;
	    }
	    detectedRank = rank;
	  }
	else 
	  {
	    TEST_FOR_EXCEPTION(true, std::logic_error, 
			       "Should never get here!  Invalid robustness value " 
			       << robustness << ".  Please report this bug to the "
			       "Belos developers.");
	  }

	return std::make_pair (detectedRank, foundRankDeficiency);
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
      
      //! \f$\| b - A x \|_2\f$ (Frobenius norm if b has more than one column).
      magnitude_type
      leastSquaresResidualNorm (const mat_type& A,
				const mat_type& x,
				const mat_type& b) const
      {
	mat_type r (b.numRows(), b.numCols());

	// r := b - A*x
	r.assign (b);
	blas_type blas;
	blas.GEMM (Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
		   b.numRows(), b.numCols(), A.numCols(),
		   -STS::one(), A.values(), A.stride(), x.values(), x.stride(),
		   STS::one(), r.values(), r.stride());
	return r.normFrobenius ();
      }

      /// \brief ||x_approx - x_exact||_2 // ||x_exact||_2.
      ///
      /// Use the Frobenius norm if more than one column.
      /// Don't scale if ||x_exact|| == 0.
      magnitude_type
      solutionError (const mat_type& x_approx,
		     const mat_type& x_exact) const
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
      /// Return the residual of the resulting least-squares problem,
      /// assuming that the upper triangular system Ry=z can be solved
      /// exactly (with zero residual).  (This may not be the case if
      /// R is singular and the system Ry=z is inconsistent.)
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

	// 6. Now that we have the updated R factor of H, and the
	//    updated right-hand side z, solve the least-squares
	//    problem by solving the upper triangular linear system
	//    Ry=z for y.
	{
	  const mat_type R_view (Teuchos::View, R, numRows-1, numRows-1);
	  const mat_type z_view (Teuchos::View, z, numRows-1, z.numCols());
	  mat_type y_view (Teuchos::View, y, numRows-1, y.numCols());

	  (void) solveUpperTriangularSystem (y_view, R_view, z_view);
	}

#if 0
	// 
	// The BLAS' _TRSM overwrites the right-hand side with the
	// solution, so first copy z into y.  We need to keep z, since
	// each call to one of this class' update routines updates it.
	{
	  TEST_FOR_EXCEPTION(y.numCols() > z.numCols(), std::logic_error,
			     "y.numCols() = " << y.numCols() 
			     << " > z.numCols() = " << z.numCols() << ".");
	  // If z has more columns than y, we ignore the remaining
	  // columns of z when solving the upper triangular system.
	  // They will get updated correctly, but they just won't be
	  // counted in the solution.  Of course, z should always have
	  // been initialized with the same number of columns as y,
	  // but we relax this unnecessary restriction here.
	  const mat_type z_view (Teuchos::View, z, numRows, y.numCols());
	  mat_type y_view (Teuchos::View, y, numRows, y.numCols());
	  y_view.assign (z_view);

	  // Solve Ry = z.
	  blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
		    Teuchos::NON_UNIT_DIAG, numRows - 1, y.numCols(), 
		    STS::one(), R.values(), LDR, y.values(), y.stride());
	}
#endif // 0

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

      /// \brief Update current column using LAPACK least-squares solver.
      ///
      /// Update the current column of the projected least-squares
      /// problem, using LAPACK's _GELS least-squares solver.  This is
      /// inefficient, but useful for testing \c updateColumnGivens().
      /// See that method's documentation for an explanation of the
      /// arguments.
      magnitude_type
      updateColumnLapack (const mat_type& H, mat_type& R, mat_type& y, 
			  mat_type& z, const int curCol) const
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
	TEST_FOR_EXCEPTION(info != 0, std::logic_error,
			   "LAPACK _GELS workspace query failed with INFO = " 
			   << info << ", for a " << numRows << " x " << numCols
			   << " matrix with " << y_view.numCols() 
			   << " right hand side"
			   << ((y_view.numCols() != 1) ? "s" : "") << ".");
	TEST_FOR_EXCEPTION(STS::real(lworkScalar) < STM::zero(),
			   std::logic_error,
			   "LAPACK _GELS workspace query returned an LWORK with "
			   "negative real part: LWORK = " << lworkScalar 
			   << ".  That should never happen.  Please report this "
			   "to the Belos developers.");
	TEST_FOR_EXCEPTION(STS::isComplex && STS::imag(lworkScalar) != STM::zero(),
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

	TEST_FOR_EXCEPTION(info != 0, std::logic_error,
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
      /// This method implements a "left-looking panel QR factorization"
      /// of the upper Hessenberg matrix in the projected least-squares
      /// problem.  It's "left-looking" because we don't updating anything
      /// to the right of columns [startCol, endCol], which is the
      /// "panel."
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
			   const int endCol) const
      {
	const int numRows = endCol + 2;
	const int numCols = endCol - startCol + 1;
	const int LDR = R.stride();

	// 1. Copy columns [startCol, endCol] from H into R, where they
	//    will be modified.
	{
	  const mat_type H_view (Teuchos::View, H, numRows, numCols, 0, startCol);
	  mat_type R_view (Teuchos::View, R, numRows, numCols, 0, startCol);
	  R_view.assign (H_view);
	}

	// 2. Apply all the previous Givens rotations, if any, to columns
	//    [startCol, endCol] of the matrix.
	blas_type blas;
	for (int j = 0; j < startCol; ++j) {
	  blas.ROT (numCols - j, 
		    &R(j, startCol), LDR, &R(j+1, startCol), LDR, 
		    &theCosines[j], &theSines[j]);
	}

	// 3. Update each column in turn of columns [startCol, endCol].
	for (int curCol = startCol; curCol < endCol; ++curCol) {
	  // a. Apply the previous Givens rotations computed in this loop
	  //    (Step 3 of the current invocation of this method) to the
	  //    current column of R.
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

	// 4. Now that we have the updated R factor of H, and the updated
	//    right-hand side z, solve the least-squares problem by
	//    solving the upper triangular linear system Ry=z for y.
	{
	  const mat_type R_view (Teuchos::View, R, numRows-1, numRows-1);
	  const mat_type z_view (Teuchos::View, z, numRows-1, z.numCols());
	  mat_type y_view (Teuchos::View, y, numRows-1, y.numCols());

	  (void) solveUpperTriangularSystem (y_view, R_view, z_view);
	}

#if 0
	// 
	// The BLAS' _TRSM overwrites the right-hand side with the
	// solution, so first copy z into y.  We need to keep z, since
	// each call to one of this class' update routines updates it.
	{
	  // If z has more columns than y, we ignore the remaining columns
	  // of z when solving the upper triangular system.  They will get
	  // updated correctly, but they just won't be counted in the
	  // solution.  Of course, z should always have been initialized
	  // with the same number of columns as y, but we relax this
	  // unnecessary restriction here.
	  const mat_type z_view (Teuchos::View, z, numRows, y.numCols());
	  mat_type y_view (Teuchos::View, y, numRows, y.numCols());
	  y_view.assign (z_view);
	}
	// Solve Ry = z.
	blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
		  Teuchos::NON_UNIT_DIAG, numRows - 1, y.numCols(),
		  STS::one(), R.values(), LDR, y.values(), y.stride());
#endif // 0

	// The last entry of z is the nonzero part of the residual of the
	// least-squares problem.  Its magnitude gives the residual 2-norm
	// of the least-squares problem.
	return STS::magnitude( z(numRows-1, 0) );
      }
    }; // class ProjectedLeastSquaresSolver
  } // namespace details
} // namespace Belos

#endif // __Belos_ProjectedLeastSquaresSolver_hpp
