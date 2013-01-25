// ***********************************************************************
// 
//      Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************


/*! \file Ifpack2_UnitTestChebyshev2.cpp
\brief A convergence test for Ifpack2::Chebyshev.
\author Mark Hoemmen

This test compares Ifpack2's implementation of Chebyshev iteration
(which as of 25 Jan 2013 is a direct imitation of Ifpack's
implementation) against the following implementations:
1. A textbook version of the algorithm
2. A textbook implementation of CG

All three implementations use left diagonal scaling.  "Textbook" means
"Templates for the Solution of Linear Systems," 2nd edition.  We
include CG just to give us a rough measure of how fast the methods
"should" converge.

The test exercises all three algorithms with a 1-D Poisson equation.
We know the eigenvalues of the matrix exactly as a function of its
dimensions (see Chapter 6 of "Applied Numerical Linear Algebra," James
Demmel, SIAM), so we can give perfect eigenvalue bounds.  We also
experiment with variations of max and min eigenvalue estimates and the
eigenvalue ratio parameter.

This test has the following command-line arguments:
- numIters: The number of iterations of Chebyshev or CG.
- localNumRows: The number of rows of the matrix on each process.
- numEigIters: The number of iterations of eigenvalue analysis (the
  power method) to find the max eigenvalue of the matrix.

The textbook implementation of Chebyshev converges faster if the
eigenvalue bounds are good, but it is much more sensitive than
Ifpack2::Chebyshev to an incorrect upper bound on the eigenvalues.
This gives me confidence that Ifpack2's version is correct.

There is also dead code for imitating ML's Chebyshev implementation
(ML_Cheby(), in packages/ml/src/Smoother/ml_smoother.c) in
Ifpack2::Details::Chebyshev.  I couldn't get it to converge, and
didn't want to waste time.  ML uses Ifpack::Chebyshev for the
top-level smoother if you give it an Epetra matrix, so Ifpack's
implementation (which Ifpack2 imitates) been tested in the field.
*/

#include <Ifpack2_ConfigDefs.hpp>
#include <Ifpack2_Chebyshev.hpp>
#include <Ifpack2_UnitTestHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_UnitTestRepository.hpp>

// mfh 24 Jan 2013: Hack to make the linker stop complaining that it
// can't find the Ifpack2::Details::Chebyshev methods.
#ifdef HAVE_IFPACK2_EXPLICIT_INSTANTIATION
#include "Ifpack2_Details_Chebyshev_def.hpp"
#endif // HAVE_IFPACK2_EXPLICIT_INSTANTIATION

namespace {

/// \class CG
/// \brief Method of conjugate gradients, with left-scaling preconditioning.
/// \tparam ScalarType The type of entries in the matrix and vectors.
/// \tparam MV Specialization of Tpetra::MultiVector.
/// \tparam MAT Corresponding specialization of Tpetra::CrsMatrix.
///
/// This class requires that the matrix A be symmetric
/// (resp. Hermitian) positive definite.
template<class ScalarType, class MV, class MAT>
class CG {
public:
  typedef ScalarType ST;
  typedef Teuchos::ScalarTraits<ScalarType> STS;
  typedef typename STS::magnitudeType MT;
  typedef Tpetra::Vector<typename MV::scalar_type,
			 typename MV::local_ordinal_type,
			 typename MV::global_ordinal_type,
			 typename MV::node_type> V;
  /// Constructor.
  ///
  /// \param A [in] The matrix A in the linear system to solve.
  ///   A must be symmetric (resp. Hermitian) positive definite.
  CG (Teuchos::RCP<const MAT> A) : 
    A_ (A), 
    D_ (getDiagonal (*A)),
    numIters_ (1)
  {}

  /// Constructor with parameters.
  ///
  /// \param A [in] The matrix A in the linear system to solve.
  ///   A must be symmetric (resp. Hermitian) positive definite.
  /// \param params [in/out] On input: the parameters.  On output:
  ///   filled with the current parameter settings.
  CG (Teuchos::RCP<const MAT> A, Teuchos::ParameterList& params) : 
    A_ (A), 
    D_ (getDiagonal (*A)),
    numIters_ (1)
  {
    setParameters (params);
  }

  /// \brief Set (or reset) parameters.
  ///
  /// This method accepts the following parameters:
  /// - "relaxation: sweeps" (\c int): numIters, the number of iterations.
  ///
  /// \pre numIters >= 0
  void setParameters (Teuchos::ParameterList& plist) {
    int numIters = numIters_;
    if (plist.isParameter ("chebyshev: degree")) {
      numIters = plist.get<int> ("chebyshev: degree");
    } else if (plist.isParameter ("CG: iterations")) {
      numIters = plist.get<int> ("CG: iterations");
    } else {
      numIters = plist.get ("relaxation: sweeps", numIters);
    }
    numIters_ = numIters;
  }

  /// Solve Ax=b for x with Chebyshev iteration, using diagonal left preconditioning.
  ///
  /// \param b [in] Right-hand side(s) in the linear system to solve.
  /// \param x [in] Initial guess(es) for the linear system to solve.
  ///
  /// \return Max (over all columns) absolute residual 2-norm after iterating.
  MT apply (const MV& b, MV& x) {
    return leftScaledCG (*A_, b, x, numIters_, *D_);
  }

private:
  Teuchos::RCP<const MAT> A_;
  Teuchos::RCP<const V> D_;
  int numIters_;

  //! r = b - A * x
  static void 
  computeResidual (MV& r, const MV& b, const MAT& A, const MV& x, 
		   const Teuchos::ETransp mode = Teuchos::NO_TRANS) 
  {
    r = b;
    A.apply (x, r, mode, -STS::one(), STS::one());
  }

  //! z = D_inv * r, = D \ r.
  static void solve (MV& z, const V& D_inv, const MV& r) {
    z.elementWiseMultiply (STS::one(), D_inv, r, STS::zero());
  }

  //! Get a copy of the diagonal of the matrix, as a row Map vector.
  static Teuchos::RCP<V> getDiagonal (const MAT& A) {
    Teuchos::RCP<V> D (new V (A.getGraph ()->getRowMap ()));
    A.getLocalDiagCopy (*D);

    typedef Kokkos::MultiVector<ST, typename MAT::node_type> KMV;
    KMV& localDiag = D->getLocalMVNonConst ();
    typedef Kokkos::DefaultArithmetic<KMV> KMVT;
    KMVT::ReciprocalThreshold (localDiag, STS::eps ());

    return D;
  }

  /// Solve Ax=b for x with CG, using diagonal left preconditioning.
  ///
  /// \pre A must be real-valued and symmetric positive definite.
  /// \pre iterNum >= 0.
  /// \pre All entries of D_inv are positive.
  ///
  /// \param A [in] The matrix A in the linear system to solve.
  /// \param b [in] Right-hand side(s) in the linear system to solve.
  /// \param x [in] Initial guess(es) for the linear system to solve.
  /// \param numIters [in] Number of iterations.
  /// \param D_inv [in] Vector of diagonal entries of A.  It must have
  ///   the same distribution as b.
  ///
  /// \return Max (over all columns) absolute residual 2-norm after iterating.
  static MT
  leftScaledCG (const MAT& A,
		const MV& B,
		MV& X,
		const int numIters,
		const V& D_inv)
  {
    Teuchos::Array<MT> norms (B.getNumVectors ());
    for (size_t j = 0; j < B.getNumVectors (); ++j) {
      Teuchos::RCP<const V> b_j = B.getVector (j);
      Teuchos::RCP<V> x_j = X.getVectorNonConst (j);
      norms[j] = oneVecLeftScaledCG (A, *b_j, *x_j, numIters, D_inv);
    }
    return *std::max_element (norms.begin (), norms.end ());
  }

  /// Solve Ax=b for x with CG, using diagonal left preconditioning.
  ///
  /// \pre A must be real-valued and symmetric positive definite.
  /// \pre iterNum >= 0.
  /// \pre All entries of D_inv are positive.
  ///
  /// \param A [in] The matrix A in the linear system to solve.
  /// \param b [in] Right-hand side(s) in the linear system to solve.
  /// \param x [in] Initial guess(es) for the linear system to solve.
  /// \param numIters [in] Number of iterations.
  /// \param D_inv [in] Vector of diagonal entries of A.  It must have
  ///   the same distribution as b.
  ///
  /// \return Max (over all columns) absolute residual 2-norm after iterating.
  static MT
  oneVecLeftScaledCG (const MAT& A,
		      const V& b,
		      V& x,
		      const int numIters,
		      const V& D_inv)
  {
    const ST one = STS::one ();
    V r (b.getMap ());
    V p (b.getMap ());
    V q (b.getMap ());
    V z (b.getMap ());

    ST alpha, beta, rho, rho_prev;
    computeResidual (r, b, A, x); // r = b - A*x
    for (int i = 0; i < numIters; ++i) {
      solve (z, D_inv, r); // z = D_inv * r, that is, D \ r.
      rho = r.dot (z); // rho = r^T z; not sure if the order is right for complex arithmetic.
      if (i == 0) {
	p = z;
      } else {
	beta = rho / rho_prev;
	p.update (one, z, beta); // p = z + beta*p
      }
      A.apply (p, q);
      const ST p_dot_q = p.dot (q); // p_dot_q = p^T q; not sure if the order is right for complex arithmetic.
      alpha = rho / p_dot_q;
      x.update (+alpha, p, one); // x = x + alpha*p
      r.update (-alpha, q, one); // r = r - alpha*q
      rho_prev = rho;
    }

    const bool computeResidualNorm = true;
    if (computeResidualNorm) {
      computeResidual (r, b, A, x);
      return r.norm2 ();
    }
  }
};

//////////////////////////////////////////////////////////////////////
// Command-line arguments
//////////////////////////////////////////////////////////////////////

// They have long names so I don't confuse them with the shorter-named
// actual options in the body of the test.
int numberOfIterations = 50;
int numberOfEigenanalysisIterations = 15;
int localNumberOfRows = 10000;

} // namespace (anonymous) 


TEUCHOS_STATIC_SETUP()
{
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP ();
  clp.setOption ("numIters", &numberOfIterations, 
		 "Number of Chebyshev iterations");
  clp.setOption ("numEigIters", &numberOfEigenanalysisIterations, 
		 "Number of iterations of eigenvalue analysis (e.g., power method)");
  clp.setOption ("localNumRows", &localNumberOfRows, 
		 "Number of rows per process in the sparse matrix.");
}


TEUCHOS_UNIT_TEST(Ifpack2Chebyshev, Convergence)
{
  // We are now in a class method declared by the above macro, and
  // that method has these input arguments:
  // Teuchos::FancyOStream& out, bool& success

  using Tpetra::global_size_t;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::Comm;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
  using std::cout;
  using std::endl;

  // Typedefs for basic Tpetra template parameters.
  typedef double ST;
  typedef int LO;
  //typedef long GO;
  typedef int GO;
  //typedef Kokkos::SerialNode NT;
  typedef Kokkos::DefaultNode::DefaultNodeType NT;

  // Convenience typedefs.
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Tpetra::CrsMatrix<ST, LO, GO, NT> crs_matrix_type;
  typedef Tpetra::MultiVector<ST, LO, GO, NT> MV;
  typedef Tpetra::Vector<ST, LO, GO, NT> V;
  typedef Ifpack2::Chebyshev<crs_matrix_type> prec_type;
  typedef Teuchos::ScalarTraits<ST> STS;
  typedef STS::magnitudeType MT;

  const ST zero = STS::zero ();
  const ST one = STS::one ();
  const ST two = one + one;

  // Prepare arguments for creating the Map.
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm ();
  RCP<NT> node;
  {
    ParameterList junk;
    node = rcp (new NT (junk));
  }
  const size_t localNumRows = as<size_t> (localNumberOfRows);
  const global_size_t globalNumRows = localNumRows * comm->getSize ();
  const GO indexBase = 0;
  const Tpetra::LocalGlobal lg = Tpetra::GloballyDistributed;

  // Create the row Map of the matrix, and the matrix's other Maps.
  RCP<const map_type> rowMap (new map_type (globalNumRows, indexBase, comm, lg, node));
  RCP<const map_type> rangeMap = rowMap;
  RCP<const map_type> domainMap = rowMap;

  // Create the matrix, with static profile.
  RCP<crs_matrix_type> A (new crs_matrix_type (rowMap, 3, Tpetra::StaticProfile));

  // Fill the matrix.
  Array<GO> cols (3);
  Array<ST> vals (3);
  for (GO globalRow = rowMap->getMinGlobalIndex (); 
       globalRow <= rowMap->getMaxGlobalIndex (); ++globalRow) {
    size_t numEntries = 3;
    if (globalRow == rowMap->getMinAllGlobalIndex ()) {
      numEntries = 2;
      cols[0] = globalRow;
      cols[1] = globalRow+1;
      vals[0] = two;
      vals[1] = -one;
    }
    else if (globalRow == rowMap->getMaxAllGlobalIndex ()) {
      numEntries = 2;
      cols[0] = globalRow-1;
      cols[1] = globalRow;
      vals[0] = -one;
      vals[1] = two;
    }
    else {
      numEntries = 3;
      cols[0] = globalRow-1;
      cols[1] = globalRow;
      cols[2] = globalRow+1;
      vals[0] = -one;
      vals[1] = two;
      vals[2] = -one;
    }
    ArrayView<const GO> colsView = cols.view (0, numEntries);
    ArrayView<const ST> valsView = vals.view (0, numEntries);
    A->insertGlobalValues (globalRow, colsView, valsView);
  }
  A->fillComplete (domainMap, rangeMap);

  // See James Demmel, "Applied Numerical Linear Algebra," SIAM,
  // pp. 267-8.  The eigenvalue approximations apply to the N x N
  // matrix with typical row (-1, 2, -1), representing the
  // discretization of the 1-D Poisson equation with Dirichlet
  // boundary conditions.
  const ST pi = acos (-1.0);
  const ST N = as<ST> (globalNumRows);
  const ST lambdaMax = two * (one - cos ((pi*N) / (N+one)));
  const ST lambdaMin = (pi / N) * (pi / N);
  ST eigRatio = lambdaMax / lambdaMin;
  const int numIters = numberOfIterations;
  int numEigIters = numberOfEigenanalysisIterations;

  // Set up the linear system to solve.
  V x_exact (domainMap), x (domainMap), b (rangeMap);
  x_exact.randomize ();
  A->apply (x_exact, b);
  x.putScalar (zero);

  V r (rangeMap); // for storing residual vector(s).
  Array<MT> norms (b.getNumVectors ());

  // Compute max initial absolute residual 2-norm.
  r = b;
  A->apply (x, r, Teuchos::NO_TRANS, -one, one);
  r.norm2 (norms ());
  const MT maxInitResNorm = *std::max_element (norms.begin (), norms.end ());

  // Make an output stream that only prints on Proc 0.
  // The usual 'out' stream in unit tests only prints if the test failed.
  Teuchos::oblackholestream blackHole;
  std::ostream& os2 = (comm->getRank () == 0) ? cout : blackHole;

  os2 << std::scientific;
  os2 << endl
      << "numIters: " << numIters << endl
      << "localNumRows: " << localNumRows << endl
      << "globalNumRows: " << globalNumRows << endl
      << "lambdaMin: " << lambdaMin << endl
      << "lambdaMax: " << lambdaMax << endl
      << "eigRatio: " << eigRatio << endl
      << "Initial residual norm: " << maxInitResNorm << endl 
      << endl;

  Teuchos::ParameterList params;
  // Set parameters for the various Chebyshev implementations.  The
  // above Chebyshev class understands many of the same parameters as
  // Ifpack2, Ifpack, and ML.  For this first pass, we only set the
  // max eigenvalue.  Below, we'll experiment with also setting the
  // min eigenvalue and the min / max eigenvalue ratio.
  params.set ("chebyshev: eigenvalue max iterations", numEigIters);
  params.set ("chebyshev: degree", numIters);
  params.set ("chebyshev: max eigenvalue", lambdaMax);

  // Create the operators: Ifpack2, textbook Chebyshev, and custom CG.
  prec_type ifpack2Cheby (A);
  Ifpack2::Details::Chebyshev<ST, MV, crs_matrix_type> myCheby (A);
  CG<ST, MV, crs_matrix_type> cg (A);

  // Residual 2-norms for comparison.
  MT maxResNormIfpack2, maxResNormTextbook, maxResNormCg;

  ////////////////////////////////////////////////////////////////////
  // Test 1: set lambdaMax exactly, use default values of eigRatio and
  // lambdaMin.  Run each version of Chebyshev and compare results.
  ////////////////////////////////////////////////////////////////////

  // Run Ifpack2's version of Chebyshev.
  ifpack2Cheby.setParameters (params);
  ifpack2Cheby.initialize ();
  ifpack2Cheby.compute ();
  ifpack2Cheby.apply (b, x);
  r = b;
  A->apply (x, r, Teuchos::NO_TRANS, -one, one);
  r.norm2 (norms ());
  maxResNormIfpack2 = *std::max_element (norms.begin (), norms.end ());

  // Run our custom version of Chebyshev.
  x.putScalar (zero); // Reset the initial guess(es).
  params.set ("chebyshev: textbook algorithm", true);
  myCheby.setParameters (params);
  myCheby.compute ();
  maxResNormTextbook = myCheby.apply (b, x);

  // Run CG, just to compare.
  x.putScalar (zero); // Reset the initial guess(es).
  cg.setParameters (params);
  maxResNormCg = cg.apply (b, x);

  os2 << "Results with lambdaMax = " << lambdaMax 
      << ", default lambdaMin and eigRatio:" << endl
      << "- Ifpack2::Chebyshev:         " << maxResNormIfpack2 / maxInitResNorm << endl
      << "- Textbook Chebyshev:         " << maxResNormTextbook / maxInitResNorm << endl
      << "- CG:                         " << maxResNormCg / maxInitResNorm << endl;

  ////////////////////////////////////////////////////////////////////
  // Test 2: set lambdaMax and lambdaMin exactly, and set eigRatio =
  // lambdaMax / lambdaMin.
  ////////////////////////////////////////////////////////////////////

  // Reset parameters.
  params.set ("chebyshev: textbook algorithm", false);
  params.set ("chebyshev: min eigenvalue", lambdaMin);
  params.set ("chebyshev: ratio eigenvalue", eigRatio);

  // Run Ifpack2's version of Chebyshev.
  ifpack2Cheby.setParameters (params);
  ifpack2Cheby.initialize ();
  ifpack2Cheby.compute ();
  ifpack2Cheby.apply (b, x);
  r = b;
  A->apply (x, r, Teuchos::NO_TRANS, -one, one);
  r.norm2 (norms ());
  maxResNormIfpack2 = *std::max_element (norms.begin (), norms.end ());

  // Run our custom version of Chebyshev.
  x.putScalar (zero); // Reset the initial guess(es).
  params.set ("chebyshev: textbook algorithm", true);
  myCheby.setParameters (params);
  myCheby.compute ();
  maxResNormTextbook = myCheby.apply (b, x);

  // Run CG, just to compare.
  x.putScalar (zero); // Reset the initial guess(es).
  cg.setParameters (params);
  maxResNormCg = cg.apply (b, x);

  os2 << "Results with lambdaMax = " << lambdaMax 
      << ", lambdaMin = " << lambdaMin << ", eigRatio = " << eigRatio << endl
      << "- Ifpack2::Chebyshev:         " << maxResNormIfpack2 / maxInitResNorm << endl
      << "- Textbook Chebyshev:         " << maxResNormTextbook / maxInitResNorm << endl
      << "- CG:                         " << maxResNormCg / maxInitResNorm << endl;

  // Reset parameters.
  params.set ("chebyshev: textbook algorithm", false);
  // ParameterList is NOT a delta.  That is, if we remove these
  // parameters from the list, setParameters() will use default
  // values, rather than letting the current settings remain.
  params.remove ("chebyshev: min eigenvalue", false);
  params.remove ("chebyshev: ratio eigenvalue", false);

  ////////////////////////////////////////////////////////////////////
  // Test 3: set lambdaMax exactly, and set eigRatio = 20 (ML's
  // default).
  ////////////////////////////////////////////////////////////////////

  // Set new parameter values.
  eigRatio = Teuchos::as<ST> (20);
  params.set ("chebyshev: ratio eigenvalue", eigRatio);

  //
  // Run each version of Chebyshev and compare their results.
  //
  // Run Ifpack2's version of Chebyshev.
  ifpack2Cheby.setParameters (params);
  ifpack2Cheby.initialize ();
  ifpack2Cheby.compute ();
  ifpack2Cheby.apply (b, x);
  r = b;
  A->apply (x, r, Teuchos::NO_TRANS, -one, one);
  r.norm2 (norms ());
  maxResNormIfpack2 = *std::max_element (norms.begin (), norms.end ());

  // Run our custom version of Chebyshev.
  x.putScalar (zero); // Reset the initial guess(es).
  params.set ("chebyshev: textbook algorithm", true);
  myCheby.setParameters (params);
  myCheby.compute ();
  maxResNormTextbook = myCheby.apply (b, x);

  // Run CG, just to compare.
  x.putScalar (zero); // Reset the initial guess(es).
  cg.setParameters (params);
  maxResNormCg = cg.apply (b, x);

  os2 << "Results with lambdaMax = " << lambdaMax 
      << ", default lambdaMin, eigRatio = " << eigRatio << endl
      << "- Ifpack2::Chebyshev:         " << maxResNormIfpack2 / maxInitResNorm << endl
      << "- Textbook Chebyshev:         " << maxResNormTextbook / maxInitResNorm << endl
      << "- CG:                         " << maxResNormCg / maxInitResNorm << endl;

  // Reset parameters.
  params.set ("chebyshev: textbook algorithm", false);

  ////////////////////////////////////////////////////////////////////
  // Test 4: set lambdaMax exactly, and set eigRatio = 30 (Ifpack's
  // default).
  ////////////////////////////////////////////////////////////////////
  eigRatio = Teuchos::as<ST> (30);
  params.set ("chebyshev: ratio eigenvalue", eigRatio);

  //
  // Run each version of Chebyshev and compare their results.
  //
  // Run Ifpack2's version of Chebyshev.
  ifpack2Cheby.setParameters (params);
  ifpack2Cheby.initialize ();
  ifpack2Cheby.compute ();
  ifpack2Cheby.apply (b, x);
  r = b;
  A->apply (x, r, Teuchos::NO_TRANS, -one, one);
  r.norm2 (norms ());
  maxResNormIfpack2 = *std::max_element (norms.begin (), norms.end ());

  // Run our custom version of Chebyshev.
  x.putScalar (zero); // Reset the initial guess(es).
  params.set ("chebyshev: textbook algorithm", true);
  myCheby.setParameters (params);
  myCheby.compute ();
  maxResNormTextbook = myCheby.apply (b, x);

  // Run CG, just to compare.
  x.putScalar (zero); // Reset the initial guess(es).
  cg.setParameters (params);
  maxResNormCg = cg.apply (b, x);

  os2 << "Results with lambdaMax = " << lambdaMax 
      << ", default lambdaMin, eigRatio = " << eigRatio << endl
      << "- Ifpack2::Chebyshev:         " << maxResNormIfpack2 / maxInitResNorm << endl
      << "- Textbook Chebyshev:         " << maxResNormTextbook / maxInitResNorm << endl
      << "- CG:                         " << maxResNormCg / maxInitResNorm << endl;

  // Reset parameters to their original values.
  params.set ("chebyshev: textbook algorithm", false);
  params.remove ("chebyshev: ratio eigenvalue", false);
  eigRatio = lambdaMax / lambdaMin;

  ////////////////////////////////////////////////////////////////////
  // Test 5: Clear lambdaMax, lambdaMin, and eigRatio.  Let the
  // smoother do eigenanalysis to estimate lambdaMax.
  ////////////////////////////////////////////////////////////////////

  params.remove ("chebyshev: min eigenvalue", false);
  params.remove ("chebyshev: max eigenvalue", false);
  params.remove ("chebyshev: ratio eigenvalue", false);

  // Run Ifpack2's version of Chebyshev.
  ifpack2Cheby.setParameters (params);
  ifpack2Cheby.initialize ();
  ifpack2Cheby.compute ();
  ifpack2Cheby.apply (b, x);
  r = b;
  A->apply (x, r, Teuchos::NO_TRANS, -one, one);
  r.norm2 (norms ());
  maxResNormIfpack2 = *std::max_element (norms.begin (), norms.end ());

  // Run our custom version of Chebyshev.
  x.putScalar (zero); // Reset the initial guess(es).
  params.set ("chebyshev: textbook algorithm", true);
  myCheby.setParameters (params);
  myCheby.compute ();
  maxResNormTextbook = myCheby.apply (b, x);

  // Run CG, just to compare.
  x.putScalar (zero); // Reset the initial guess(es).
  cg.setParameters (params);
  maxResNormCg = cg.apply (b, x);

  os2 << "Results with default lambdaMax, lambdaMin, and eigRatio, "
    "with numEigIters = " << numEigIters << ":" << endl
      << "- Ifpack2::Chebyshev:         " << maxResNormIfpack2 / maxInitResNorm << endl
      << "- Textbook Chebyshev:         " << maxResNormTextbook / maxInitResNorm << endl
      << "- CG:                         " << maxResNormCg / maxInitResNorm << endl;

  // Print the computed max and min eigenvalues, and other details.
  os2 << endl;
  myCheby.print (os2);

  // Reset parameters.
  params.set ("chebyshev: textbook algorithm", false);

  ////////////////////////////////////////////////////////////////////
  // Test 6: Clear lambdaMax, lambdaMin, and eigRatio.  Let the
  // smoother do eigenanalysis to estimate lambdaMax, with more
  // iterations.
  ////////////////////////////////////////////////////////////////////

  params.remove ("chebyshev: min eigenvalue", false);
  params.remove ("chebyshev: max eigenvalue", false);
  params.remove ("chebyshev: ratio eigenvalue", false);
  numEigIters = 2 * numEigIters;
  params.set ("chebyshev: eigenvalue max iterations", numEigIters);

  // Run Ifpack2's version of Chebyshev.
  ifpack2Cheby.setParameters (params);
  ifpack2Cheby.initialize ();
  ifpack2Cheby.compute ();
  ifpack2Cheby.apply (b, x);
  r = b;
  A->apply (x, r, Teuchos::NO_TRANS, -one, one);
  r.norm2 (norms ());
  maxResNormIfpack2 = *std::max_element (norms.begin (), norms.end ());

  // Run our custom version of Chebyshev.
  x.putScalar (zero); // Reset the initial guess(es).
  params.set ("chebyshev: textbook algorithm", true);
  myCheby.setParameters (params);
  myCheby.compute ();
  maxResNormTextbook = myCheby.apply (b, x);

  // Run CG, just to compare.
  x.putScalar (zero); // Reset the initial guess(es).
  cg.setParameters (params);
  maxResNormCg = cg.apply (b, x);

  os2 << "Results with default lambdaMax, lambdaMin, and eigRatio, "
    "with numEigIters = " << numEigIters << ":" << endl
      << "- Ifpack2::Chebyshev:         " << maxResNormIfpack2 / maxInitResNorm << endl
      << "- Textbook Chebyshev:         " << maxResNormTextbook / maxInitResNorm << endl
      << "- CG:                         " << maxResNormCg / maxInitResNorm << endl;


  // For this case, if there are enough eigenanalysis iterations,
  // Ifpack2 should do quite a bit better than the textbook version of
  // the algorithm.  We'll be generous and say that it does "no worse"
  // than the textbook version.  We give "wiggle room" of four digits
  // for defining "no worse than."
  if (numEigIters >= 15) {
    const MT tol = Teuchos::as<MT> (1.0e-4);
    // Avoid division by zero when computing relative accuracy.
    const MT relDiff = maxResNormTextbook == zero ? 
      STS::magnitude (maxResNormIfpack2 - maxResNormTextbook) :
      STS::magnitude (maxResNormIfpack2 - maxResNormTextbook) / maxResNormTextbook;
    TEUCHOS_TEST_FOR_EXCEPTION(
      maxResNormIfpack2 > maxResNormTextbook && relDiff > tol, 
      std::runtime_error, 
      "After " << numIters << " iterations of Chebyshev, with lambdaMax = " 
      << lambdaMax << " and default lambdaMin and eigRatio, Ifpack2::Chebyshev "
      "does quite a bit worse than the textbook version of the algorithm.  The "
      "former has a max relative residual norm of " << maxResNormIfpack2 << ", "
      "and the latter of " << maxResNormTextbook << ".");
  }

  // Print the computed max and min eigenvalues, and other details.
  os2 << endl;
  myCheby.print (os2);
}


