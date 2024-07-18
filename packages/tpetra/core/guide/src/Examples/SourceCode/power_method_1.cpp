// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_Map.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
// Power method for estimating the eigenvalue of maximum magnitude of
// a matrix.  This function returns the eigenvalue estimate.
//
// We don't intend for you to write your own eigensolvers; the Anasazi
// package provides them.  You should instead see this class as a
// surrogate for a Tpetra interface to a Trilinos package.
//
// TpetraOperatorType: the type of the Tpetra::Operator specialization
//   used to represent the sparse matrix or operator A.
//
// Tpetra::Operator implements a function from one
// Tpetra::(Multi)Vector to another Tpetra::(Multi)Vector.
// Tpetra::CrsMatrix implements Tpetra::Operator; its apply() method
// computes a sparse matrix-(multi)vector multiply.  It's typical for
// numerical algorithms that use Tpetra objects to be templated on the
// type of the Tpetra::Operator specialization.  We do so here, and
// thus demonstrate how you can use the public typedefs in Tpetra
// classes to write generic code.
//
// One could use a templated function here instead of a templated
// class with a static (class) method.  I prefer the class approach
// because one can lift typedefs out of the function into the class.
// It tends to makes the function declaration easier to read.
template <class TpetraOperatorType>
class PowerMethod {
public:
  typedef typename TpetraOperatorType::scalar_type scalar_type;
  typedef typename TpetraOperatorType::local_ordinal_type local_ordinal_type;
  typedef typename TpetraOperatorType::global_ordinal_type global_ordinal_type;
  typedef typename TpetraOperatorType::node_type node_type;
  // The type of a Tpetra vector with the same template parameters as
  // those of TpetraOperatorType.
  typedef Tpetra::Vector<scalar_type, local_ordinal_type,
                         global_ordinal_type, node_type> vec_type;
  // The type of the norm of the above Tpetra::Vector specialization.
  typedef typename vec_type::mag_type magnitude_type;
  // Run the power method and return the eigenvalue estimate.
  //
  // Input arguments:
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
  static scalar_type
  run (const TpetraOperatorType& A,
       const int niters,
       const magnitude_type tolerance,
       std::ostream& out)
  {
    using std::endl;
    typedef Teuchos::ScalarTraits<scalar_type> STS;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;
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
    // Fill the iteration vector z with random numbers to start.
    // Don't have grand expectations about the quality of our
    // pseudorandom number generator, but it is usually good enough
    // for eigensolvers.
    z.randomize ();
    // lambda: Current approximation of the eigenvalue of maximum magnitude.
    // normz: 2-norm of the current iteration vector z.
    // residual: 2-norm of the current residual vector 'resid'.
    //
    // Teuchos::ScalarTraits defines what zero and one means for any
    // type.  Most number types T know how to turn a 0 or a 1 (int)
    // into a T.  I have encountered some number types in C++ that do
    // not.  These tend to be extended-precision types that define
    // number operators and know how to convert from a float or
    // double, but don't have conversion operators for int.  Thus,
    // using Teuchos::ScalarTraits makes this code maximally general.
    scalar_type lambda = STS::zero ();
    magnitude_type normz = STM::zero ();
    magnitude_type residual = STM::zero ();
    const scalar_type one  = STS::one ();
    const scalar_type zero = STS::zero ();
    // How often to report progress in the power method.  Reporting
    // progress requires computing a residual, which can be expensive.
    // However, if you don't compute the residual often enough, you
    // might keep iterating even after you've converged.
    const int reportFrequency = 10;
    // Do the power method, until the method has converged or the
    // maximum iteration count has been reached.
    for (int iter = 0; iter < niters; ++iter) {
      normz = z.norm2 ();       // Compute the 2-norm of z
      q.scale (one / normz, z); // q := z / normz
      A.apply (q, z);           // z := A * q
      lambda = q.dot (z);       // Approx. max eigenvalue
      // Compute and report the residual norm every reportFrequency
      // iterations, or if we've reached the maximum iteration count.
      if (iter % reportFrequency == 0 || iter + 1 == niters) {
        resid.update (one, z, -lambda, q, zero); // z := A*q - lambda*q
        residual = resid.norm2 (); // 2-norm of the residual vector
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
};
int
main (int argc, char *argv[])
{
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
  typedef Tpetra::Map<> map_type;
  typedef Tpetra::Vector<>::scalar_type scalar_type;
  typedef Tpetra::Vector<>::local_ordinal_type local_ordinal_type;
  typedef Tpetra::Vector<>::global_ordinal_type global_ordinal_type;
  typedef Tpetra::Vector<>::mag_type magnitude_type;
  typedef Tpetra::CrsMatrix<> crs_matrix_type;
  typedef Tpetra::CrsMatrix<>::nonconst_global_inds_host_view_type nonconst_global_inds_host_view_type;
  typedef Tpetra::CrsMatrix<>::nonconst_values_host_view_type nonconst_values_host_view_type;

  Teuchos::oblackholestream blackhole;

  Tpetra::ScopeGuard tpetraScope(&argc, &argv);
  {
    // Get a communicator corresponding to MPI_COMM_WORLD
    Teuchos::RCP<const Teuchos::Comm<int>> comm = Tpetra::getDefaultComm();

    const size_t myRank = comm->getRank();
    // Make an output stream (for verbose output) that only prints on
    // Proc 0 of the communicator.
    Teuchos::oblackholestream blackHole;
    std::ostream& out = (myRank == 0) ? std::cout : blackHole;

    // The number of rows and columns in the matrix.
    const Tpetra::global_size_t numGblIndices = 50;

    // Construct a Map that puts approximately the same number of
    // equations on each processor.
    const global_ordinal_type indexBase = 0;
    RCP<const map_type> map = rcp (new map_type (numGblIndices, indexBase, comm));
    const size_t numMyElements = map->getLocalNumElements ();

    // If you like, you may get the list of global indices that the
    // calling process owns.  This is unnecessary if you don't mind
    // converting local indices to global indices.
    //
    // ArrayView<const global_ordinal_type> myGlobalElements =
    //   map->getLocalElementList ();
    out << endl << "Creating the sparse matrix" << endl;

    // Create a Tpetra sparse matrix whose rows have distribution given by the Map.
    RCP<crs_matrix_type> A (new crs_matrix_type (map, 3));

    // Fill the sparse matrix, one row at a time.
    const scalar_type two = static_cast<scalar_type> (2.0);
    const scalar_type negOne = static_cast<scalar_type> (-1.0);
    for (local_ordinal_type lclRow = 0;
         lclRow < static_cast<local_ordinal_type> (numMyElements);
         ++lclRow) {
      const global_ordinal_type gblRow = map->getGlobalElement (lclRow);
      // A(0, 0:1) = [2, -1]
      if (gblRow == 0) {
        A->insertGlobalValues(gblRow,
                              tuple<global_ordinal_type> (gblRow, gblRow + 1),
                              tuple<scalar_type> (two, negOne));
      }
      // A(N-1, N-2:N-1) = [-1, 2]
      else if (static_cast<Tpetra::global_size_t>(gblRow) == numGblIndices - 1) {
        A->insertGlobalValues(gblRow,
                              tuple<global_ordinal_type> (gblRow - 1, gblRow),
                              tuple<scalar_type> (negOne, two));
      }
      // A(i, i-1:i+1) = [-1, 2, -1]
      else {
        A->insertGlobalValues(gblRow,
                              tuple<global_ordinal_type> (gblRow - 1, gblRow, gblRow + 1),
                              tuple<scalar_type> (negOne, two, negOne));
      }
    }

    // Tell the sparse matrix that we are done adding entries to it.
    A->fillComplete ();

    // Number of iterations
    const int niters = 500;

    // Desired (absolute) residual tolerance
    const magnitude_type tolerance = 1.0e-2;
    // Run the power method and report the result.
    scalar_type lambda = PowerMethod<crs_matrix_type>::run (*A, niters, tolerance, out);
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
    A->resumeFill ();
    if (A->getRowMap ()->isNodeGlobalElement (0)) {
      // Get a copy of the row with with global index 0.  Modify the
      // diagonal entry of that row.  Submit the modified values to the
      // matrix.
      const global_ordinal_type idOfFirstRow = 0;
      size_t numEntriesInRow = A->getNumEntriesInGlobalRow (idOfFirstRow);
      nonconst_global_inds_host_view_type rowinds ("indices",numEntriesInRow);
      nonconst_values_host_view_type rowvals ("vals",numEntriesInRow);
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
      A->getGlobalRowCopy(idOfFirstRow, rowinds, rowvals, numEntriesInRow);
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
      A->replaceGlobalValues (idOfFirstRow, rowinds, rowvals);
    }
    // Call fillComplete() again to signal that we are done changing the
    // matrix.
    A->fillComplete ();
    // Run the power method again.
    lambda = PowerMethod<crs_matrix_type>::run (*A, niters, tolerance, out);
    out << endl << "Estimated max eigenvalue: " << lambda << endl;

    // This tells the Trilinos test framework that the test passed.
    if (myRank == 0) {
      cout << "End Result: TEST PASSED" << endl;
    }
    // ScopeGuard's destructor calls MPI_Finalize, if its constructor
    // called MPI_Init.  Likewise, it calls Kokkos::finalize, if its
    // constructor called Kokkos::initialize.
  }
  return 0;
}
