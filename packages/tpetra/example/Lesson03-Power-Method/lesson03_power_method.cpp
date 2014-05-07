//@HEADER
// ************************************************************************
//
//               Tpetra: Linear Algebra Services Package
//                 Copyright (2009) Sandia Corporation
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

/*!
\example lesson03_power_method.cpp
\brief Example showing how to fill and compute with a Tpetra sparse matrix

\ref Tpetra_Lesson03 explains this example in detail.
*/

#include <Teuchos_Array.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"

// Implementation of the power method for estimating the eigenvalue of
// maximum magnitude of a matrix.  This function returns the
// eigenvalue estimate.
//
// TpetraOperatorType: the type of the Tpetra::Operator specialization
//   used to represent the sparse matrix or operator A.
//
// A: The sparse matrix or operator, as a Tpetra::Operator.
// niters: Maximum number of iterations of the power method.
// tolerance: If the 2-norm of the residual A*x-lambda*x (for the
//   current eigenvalue estimate lambda) is less than this, stop
//   iterating.  The complicated expression for the type ensures that
//   if the type of entries in the matrix A (scalar_type) is complex,
//   then we'll be using a real-valued type ("magnitude") for the
//   tolerance.  (You can't compare complex numbers using less than,
//   so you can't test for convergence using a complex number.)
// out: output stream to which to print the current status of the
//   power method.
template <class TpetraOperatorType>
typename TpetraOperatorType::scalar_type
powerMethod (const TpetraOperatorType& A,
             const size_t niters,
             const typename Teuchos::ScalarTraits<typename TpetraOperatorType::scalar_type>::magnitudeType tolerance,
             std::ostream& out);

int
main (int argc, char *argv[])
{
  // global_size_t: Tpetra defines this unsigned integer type big
  // enough to hold any global dimension or amount of data.
  using Tpetra::global_size_t;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::ArrayRCP;
  using Teuchos::arcp;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
  using std::cerr;
  using std::cout;
  using std::endl;

  // Template parameters for Tpetra objects.  Tpetra::CrsMatrix has
  // five template parameters, but all but the first of them have
  // defaults.
  typedef double scalar_type;
  typedef int local_ordinal_type;
  typedef long global_ordinal_type;

  // Convenience typedef.
  typedef Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackhole);
  RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  const size_t myRank = comm->getRank();
  //const size_t numProcs = comm->getSize();

  // Make an output stream (for verbose output) that only prints on
  // Proc 0 of the communicator.
  Teuchos::oblackholestream blackHole;
  std::ostream& out = (myRank == 0) ? std::cout : blackHole;

  // Print the current version of Tpetra.
  out << Tpetra::version() << endl << endl;

  // The number of rows and columns in the matrix.
  const global_size_t numGlobalElements = 50;

  // Construct a Map that puts approximately the same number of
  // equations on each processor.
  const global_ordinal_type indexBase = 0;
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type> map_type;
  RCP<const map_type> map = rcp (new map_type (numGlobalElements, indexBase, comm));

  // Get the list of global indices that this process owns.  In this
  // example, this is unnecessary, because we know that we created a
  // contiguous Map (see above).  (Thus, we really only need the min
  // and max global index on this process.)  However, in general, we
  // don't know what global indices the Map owns, so if we plan to add
  // entries into the sparse matrix using global indices, we have to
  // get the list of global indices this process owns.
  const size_t numMyElements = map->getNodeNumElements();
  ArrayView<const global_ordinal_type> myGlobalElements = map->getNodeElementList ();

  out << endl << "Creating the sparse matrix" << endl;

  // Create a Tpetra sparse matrix whose rows have distribution given by the Map.
  typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type> mat_type;
  RCP<mat_type> A = rcp (new mat_type (map, 0));

  // Fill the sparse matrix, one row at a time.
  const scalar_type two    = static_cast<scalar_type>( 2.0);
  const scalar_type negOne = static_cast<scalar_type>(-1.0);
  for (size_t i = 0; i < numMyElements; ++i) {
    // A(0, 0:1) = [2, -1]
    if (myGlobalElements[i] == 0) {
      A->insertGlobalValues (myGlobalElements[i],
                             tuple<global_ordinal_type> (myGlobalElements[i], myGlobalElements[i]+1),
                             tuple<scalar_type> (two, negOne));
    }
    // A(N-1, N-2:N-1) = [-1, 2]
    else if (static_cast<global_size_t> (myGlobalElements[i]) == numGlobalElements - 1) {
      A->insertGlobalValues (myGlobalElements[i],
                             tuple<global_ordinal_type> (myGlobalElements[i]-1, myGlobalElements[i]),
                             tuple<scalar_type> (negOne, two));
    }
    // A(i, i-1:i+1) = [-1, 2, -1]
    else {
      A->insertGlobalValues (myGlobalElements[i],
                             tuple<global_ordinal_type> (myGlobalElements[i]-1, myGlobalElements[i], myGlobalElements[i]+1),
                             tuple<scalar_type> (negOne, two, negOne));
    }
  }

  // Tell the sparse matrix that we are done adding entries to it.
  A->fillComplete ();

  // Number of iterations
  const size_t niters = static_cast<size_t> (500);
  // Desired (absolute) residual tolerance
  const magnitude_type tolerance = 1.0e-2;

  // Run the power method and report the result.
  scalar_type lambda = powerMethod<mat_type> (*A, niters, tolerance, out);
  out << endl << "Estimated max eigenvalue: " << lambda << endl;

  //
  // Now we're going to change values in the sparse matrix and run the
  // power method again.  In Tpetra, if fillComplete() has been
  // called, you have to call resumeFill() before you may change the
  // matrix (either its values or its structure).
  //

  //
  // Increase diagonal dominance
  //
  out << endl << "Increasing magnitude of A(0,0), solving again" << endl;

  // Must call resumeFill() before changing the matrix, even its values.
  A->resumeFill();

  if (A->getRowMap()->isNodeGlobalElement (0)) {
    // Get a copy of the row with with global index 0.  Modify the
    // diagonal entry of that row.  Submit the modified values to the
    // matrix.
    const global_ordinal_type idOfFirstRow = 0;
    size_t numEntriesInRow = A->getNumEntriesInGlobalRow (idOfFirstRow);
    Array<scalar_type>         rowvals (numEntriesInRow);
    Array<global_ordinal_type> rowinds (numEntriesInRow);

    // Fill rowvals and rowinds with the values resp. (global) column
    // indices of the sparse matrix entries owned by the calling
    // process.
    //
    // Note that it's legal (though we don't exercise it in this
    // example) for the row Map of the sparse matrix not to be one to
    // one.  This means that more than one process might own entries
    // in the first row.  In general, multiple processes might own the
    // (0,0) entry, so that the global A(0,0) value is really the sum
    // of all processes' values for that entry.  However, scaling the
    // entry by a constant factor distributes across that sum, so it's
    // OK to do so.
    //
    // The parentheses after rowinds and rowvalues indicate "a view of
    // the Array's data."  Array::operator() returns an ArrayView.
    A->getGlobalRowCopy (idOfFirstRow, rowinds (), rowvals (), numEntriesInRow);
    for (size_t i = 0; i < numEntriesInRow; i++) {
      if (rowinds[i] == idOfFirstRow) {
        // We have found the diagonal entry; modify it.
        rowvals[i] *= 10.0;
      }
    }
    // "Replace global values" means modify the values, but not the
    // structure of the sparse matrix.  If the specified columns
    // aren't already populated in this row on this process, then this
    // method throws an exception.  If you want to modify the
    // structure (by adding new entries), you'll need to call
    // insertGlobalValues().
    A->replaceGlobalValues (idOfFirstRow, rowinds (), rowvals ());
  }

  // Call fillComplete() again to signal that we are done changing the
  // matrix.
  A->fillComplete ();

  // Run the power method again.
  lambda = powerMethod<mat_type> (*A, niters, tolerance, out);
  out << endl << "Estimated max eigenvalue: " << lambda << endl;

  // This tells the Trilinos test framework that the test passed.
  if (myRank == 0) {
    std::cout << "End Result: TEST PASSED" << std::endl;
  }

  return 0;
}


template <class TpetraOperatorType>
typename TpetraOperatorType::scalar_type
powerMethod (const TpetraOperatorType& A,
             const size_t niters,
             const typename Teuchos::ScalarTraits<typename TpetraOperatorType::scalar_type>::magnitudeType tolerance,
             std::ostream& out)
{
  using std::endl;
  typedef typename TpetraOperatorType::scalar_type scalar_type;
  typedef typename TpetraOperatorType::local_ordinal_type local_ordinal_type;
  typedef typename TpetraOperatorType::global_ordinal_type global_ordinal_type;
  typedef typename TpetraOperatorType::node_type node_type;

  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef typename STS::magnitudeType magnitude_type;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;

  // The type of a Tpetra vector with the same template parameters as those of A.
  typedef Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> vec_type;

  // Create three vectors for iterating the power method.  Since the
  // power method computes z = A*q, q should be in the domain of A and
  // z should be in the range.  (Obviously the power method requires
  // that the domain and the range are equal, but it's a good idea to
  // get into the habit of thinking whether a particular vector
  // "belongs" in the domain or range of the matrix.)  The residual
  // vector "resid" is of course in the range of A.
  vec_type q (A.getDomainMap ());
  vec_type z (A.getRangeMap ());
  vec_type resid (A.getRangeMap ());

  // Fill the iteration vector z with random numbers to start.  Don't
  // have grand expectations about the quality of our pseudorandom
  // number generator; this is usually good enough for eigensolvers.
  z.randomize();

  // lambda: the current approximation of the eigenvalue of maximum magnitude.
  // normz: the 2-norm of the current iteration vector z.
  // residual: the 2-norm of the current residual vector "resid"
  scalar_type lambda = STS::zero();
  magnitude_type normz = STM::zero();
  magnitude_type residual = STM::zero();

  const scalar_type one  = STS::one();
  const scalar_type zero = STS::zero();

  // How often to report progress in the power method.  Reporting
  // progress requires computing a residual which can be expensive.
  // However, if you don't compute the residual often enough, you
  // might keep iterating even after you've converged.
  const size_t reportFrequency = 10;

  // Do the power method, until the method has converged or the
  // maximum iteration count has been reached.
  for (size_t iter = 0; iter < niters; ++iter) {
    normz = z.norm2 (); // Compute the 2-norm of z
    q.scale (one / normz, z); // q := z / normz
    A.apply (q, z);        // z := A * q
    lambda = q.dot (z);     // Approx. max eigenvalue

    // Compute and report the residual norm every reportFrequency
    // iterations, or if we've reached the maximum iteration count.
    if (iter % reportFrequency == 0 || iter + 1 == niters) {
      resid.update (one, z, -lambda, q, zero); // z := A*q - lambda*q
      residual = resid.norm2(); // 2-norm of the residual vector
      out << "Iteration " << iter << ":" << endl
          << "- lambda = " << lambda << endl
          << "- ||A*q - lambda*q||_2 = " << residual << endl;
    }
    if (residual < tolerance) {
      out << "Converged after " << iter << " iterations" << endl;
      break;
    } else if (iter+1 == niters) {
      out << "Failed to converge after " << niters << " iterations" << endl;
      break;
    }
  }
  return lambda;
}
