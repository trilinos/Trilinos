// ***********************************************************************
// 
//      Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
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

\brief A more thorough test for Ifpack2::Chebyshev.
*/

#include <Ifpack2_ConfigDefs.hpp>
#include <Ifpack2_Chebyshev.hpp>
#include <Ifpack2_UnitTestHelpers.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_UnitTestRepository.hpp>
#include <cmath>

namespace {
/// \class Chebyshev
/// \brief Left-scaled Chebyshev iteration.
/// \tparam ScalarType The type of entries in the matrix and vectors.
/// \tparam MV Specialization of Tpetra::MultiVector.
/// \tparam MAT Corresponding specialization of Tpetra::CrsMatrix.
///
/// This class requires that the matrix A be real valued and symmetric
/// positive definite.  If users could provide the ellipse parameters
/// ("d" and "c" in the literature, where d is the real-valued center
/// of the ellipse, and d-c and d+c the two foci), the code would work
/// fine with nonsymmetric A, as long as the eigenvalues of A can be
/// bounded in an ellipse that is entirely to the right of the origin.
template<class ScalarType, class MV, class MAT>
class Chebyshev {
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
  ///   A must be real-valued and symmetric positive definite.
  Chebyshev (Teuchos::RCP<const MAT> A) : 
    A_ (A), 
    D_ (getDiagonal (*A)),
    lambdaMax_ (Teuchos::as<ST> (4)),
    lambdaMin_ (Teuchos::as<ST> (4) / Teuchos::as<ST> (30)),
    numIters_ (1)
  {}

  /// Constructor with parameters.
  ///
  /// \param A [in] The matrix A in the linear system to solve.
  ///   A must be real-valued and symmetric positive definite.
  /// \param params [in/out] On input: the parameters.  On output:
  ///   filled with the current parameter settings.
  Chebyshev (Teuchos::RCP<const MAT> A, Teuchos::ParameterList& params) : 
    A_ (A), 
    D_ (getDiagonal (*A)),
    lambdaMax_ (Teuchos::as<ST> (4)),
    lambdaMin_ (Teuchos::as<ST> (4) / Teuchos::as<ST> (30)),
    numIters_ (1)
  {
    setParameters (params);
  }

  /// \brief Set (or reset) parameters.
  ///
  /// This method accepts the following parameters:
  /// - "chebyshev: max eigenvalue" (\c ScalarType): lambdaMax, a
  ///   lower bound of real part of bounding ellipse of eigenvalues of
  ///   the matrix A.
  /// - "chebyshev: min eigenvalue" (\c ScalarType): lambdaMin, an
  ///   upper bound of real part of bounding ellipse of eigenvalues of
  ///   the matrix A.
  /// - "chebyshev: degree" (\c int): numIters, the number of iterations.
  ///
  /// \pre 0 < lambdaMin <= lambdaMax
  /// \pre lambdaMin and lambdaMax are real
  /// \pre numIters >= 0
  void setParameters (Teuchos::ParameterList& plist) {
    lambdaMax_ = plist.get ("chebyshev: max eigenvalue", lambdaMax_);
    lambdaMin_ = plist.get ("chebyshev: min eigenvalue", lambdaMin_);
    numIters_ = plist.get ("chebyshev: degree", numIters_);
  }

  /// Solve Ax=b for x with Chebyshev iteration, using diagonal left preconditioning.
  ///
  /// \pre numIters >= 0.
  ///
  /// \param b [in] Right-hand side(s) in the linear system to solve.
  /// \param x [in] Initial guess(es) for the linear system to solve.
  ///
  /// \return Max (over all columns) absolute residual 2-norm after iterating.
  MT apply (const MV& b, MV& x) {
    const ST eigRatio = Teuchos::as<ST> (30);
    const ST lambdaMax = lambdaMax_ * Teuchos::as<ST> (1.1);
    const ST lambdaMin = lambdaMax / eigRatio;

    return leftScaledChebyshevImpl (*A_, b, x, numIters_, lambdaMax, lambdaMin, *D_);
  }

private:
  Teuchos::RCP<const MAT> A_;
  Teuchos::RCP<const V> D_;
  ST lambdaMax_;
  ST lambdaMin_;
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

  /// Solve Ax=b for x with Chebyshev iteration, using diagonal left preconditioning.
  ///
  /// \pre A must be real-valued and symmetric positive definite.
  /// \pre iterNum >= 0.
  /// \pre 0 < lMin <= lMax
  /// \pre All entries of D_inv are positive.
  ///
  /// \param A [in] The matrix A in the linear system to solve.
  /// \param b [in] Right-hand side(s) in the linear system to solve.
  /// \param x [in] Initial guess(es) for the linear system to solve.
  /// \param iterNum [in] Number of Chebyshev iterations.
  /// \param lMax [in] Estimate of max eigenvalue of A.
  /// \param lMin [in] Estimate of min eigenvalue of A.
  /// \param D_inv [in] Vector of diagonal entries of A.  It must have
  ///   the same distribution as b.
  ///
  /// \return Max (over all columns) absolute residual 2-norm after iterating.
  static MT
  leftScaledChebyshevImpl (const MAT& A,
			   const MV& b,
			   MV& x,
			   const int iterNum,
			   const ST lMax,
			   const ST lMin,
			   const V& D_inv)
  {
    const ST one = Teuchos::as<ST> (1);
    const ST two = Teuchos::as<ST> (2);
    const ST d = (lMax + lMin) / two;
    const ST c = (lMax - lMin) / two;

    MV r (b.getMap (), b.getNumVectors (), false);
    MV p (b.getMap (), b.getNumVectors (), false);
    MV z (b.getMap (), b.getNumVectors (), false);
    ST alpha, beta;
    for (int i = 0; i < iterNum; ++i) {
      computeResidual (r, b, A, x); // r = b - A*x
      solve (z, D_inv, r); // z = D_inv * r, that is, D \ r.
      if (i == 0) {
	p = z;
	alpha = two / d;
      } else {
	//beta = (c * alpha / two)^2;
	//const ST sqrtBeta = c * alpha / two;
	//beta = sqrtBeta * sqrtBeta;
	beta = alpha * (c/two) * (c/two);
	alpha = one / (d - beta);
	p.update (one, z, beta); // p = z + beta*p
      }
      x.update (alpha, p, one); // x = x + alpha*p
      // If we compute the residual here, we could either do r = b -
      // A*x, or r = r - alpha*A*p.  Since we choose the former, we
      // can move the computeResidual call to the top of the loop.
    }

    const bool computeResidualNorm = true;
    if (computeResidualNorm) {
      computeResidual (r, b, A, x);
      Teuchos::Array<MT> norms (b.getNumVectors ());
      r.norm2 (norms ());
      return *std::max_element (norms.begin (), norms.end ());
    }
  }
};

} // namespace (anonymous) 

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
  typedef typename STS::magnitudeType MT;

  const ST zero = STS::zero ();
  const ST one = STS::one ();
  const ST two = one + one;

  // Any real-valued, symmetric positive definite matrix is spectrally
  // equivalent to a diagonal matrix.  Thus, it suffices to test
  // Chebyshev with a diagonal matrix.

  // Prepare arguments for creating the Map.
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm ();
  RCP<NT> node;
  {
    ParameterList junk;
    node = rcp (new NT (junk));
  }
  const size_t localNumRows = 10;
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
  for (GO globalRow = rowMap->getMinGlobalIndex (); globalRow <= rowMap->getMaxGlobalIndex (); ++globalRow) {
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
  const ST lambdaMin = (pi / globalNumRows) * (pi / globalNumRows);
  const ST lambdaMax = two * (one - cos ((pi*N) / (N+one)));
  const ST eigRatio = lambdaMax / lambdaMin;
  const int numIters = 50;

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

  // Create the Ifpack2::Chebyshev iteration operator.
  prec_type M1 (A);
  Teuchos::ParameterList params;
  params.set ("chebyshev: degree", numIters);
  params.set("chebyshev: min eigenvalue", lambdaMin);
  params.set("chebyshev: max eigenvalue", lambdaMax);
  params.set("chebyshev: ratio eigenvalue", eigRatio);
  M1.setParameters (params);
  M1.initialize ();
  M1.compute ();

  // Create the comparison Chebyshev iteration operator.
  Chebyshev<ST, MV, crs_matrix_type> M2 (A, params);

  //
  // Run each version of Chebyshev and compare their results.
  //
  // Run Ifpack2's version of Chebyshev.
  M1.apply (b, x);
  r = b;
  A->apply (x, r, Teuchos::NO_TRANS, -one, one);
  r.norm2 (norms ());
  const MT maxResNormM1 = *std::max_element (norms.begin (), norms.end ());

  // Run our custom version of Chebyshev.
  x.putScalar (zero); // Reset the initial guess(es).
  const MT maxResNormM2 = M2.apply (b, x);

  using std::cerr;
  cerr << "Initial residual norm: " << maxInitResNorm << endl
       << "Number of iterations:  " << numIters << endl
       << "Ifpack2::Chebyshev:    " << maxResNormM1 / maxInitResNorm << endl
       << "Chebyshev:             " << maxResNormM2 / maxInitResNorm << endl;
}


