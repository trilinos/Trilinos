#include <Tpetra_Core.hpp>
#include <Tpetra_Version.hpp>

#include "Solve.hpp"

int
main (int argc, char **argv)
{
  typedef Ifpack2::Test::ST ST;
  typedef Ifpack2::Test::GO GO;
  typedef Ifpack2::Test::STS STS;
  typedef Ifpack2::Test::map_type map_type;
  typedef Ifpack2::Test::multivector_type multivector_type;
  typedef Ifpack2::Test::sparse_mat_type sparse_mat_type;

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

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  // The number of rows and columns in the matrix.
  const global_size_t numGlobalElements = 50;

  // Create the matrix's row Map.
  RCP<const map_type> map (new map_type (numGlobalElements, 0, comm));
  const size_t numMyElements = map->getNodeNumElements ();
  ArrayView<const GO> myGlobalElements = map->getNodeElementList ();

  // Create a Tpetra sparse matrix with the given row Map.
  sparse_mat_type A (map, 0);

  // Fill the sparse matrix, one row at a time.
  const ST two    = static_cast<ST> ( 2.0);
  const ST negOne = static_cast<ST> (-1.0);

  for (size_t i = 0; i < numMyElements; ++i) {
    // A(0, 0:1) = [2, -1]
    if (myGlobalElements[i] == 0) {
      A.insertGlobalValues (myGlobalElements[i],
                            tuple<GO> (myGlobalElements[i], myGlobalElements[i]+1),
                            tuple<ST> (two, negOne));
    }
    // A(N-1, N-2:N-1) = [-1, 2]
    else if (static_cast<global_size_t> (myGlobalElements[i]) == numGlobalElements - 1) {
      A.insertGlobalValues (myGlobalElements[i],
                            tuple<GO> (myGlobalElements[i]-1, myGlobalElements[i]),
                            tuple<ST> (negOne, two));
    }
    // A(i, i-1:i+1) = [-1, 2, -1]
    else {
      A.insertGlobalValues (myGlobalElements[i],
                            tuple<GO> (myGlobalElements[i]-1, myGlobalElements[i], myGlobalElements[i]+1),
                            tuple<ST> (negOne, two, negOne));
    }
  }

  // Tell the sparse matrix that we are done adding entries to it.
  A.fillComplete ();

  // Create the right-hand side and initial guess of the linear system to solve.
  multivector_type b (map, 1);
  b.putScalar (STS::one ());
  multivector_type x (map, 1, true);

  // Attempt to solve the linear system, with the following parameters:
  //
  //   - ILU(k) with k = 1
  //   - Additive Schwarz overlap level: 1
  //   - GMRES restart length: 100
  //   - Maximum number of GMRES iterations (over all restarts): 1000
  //   - GMRES relative residual tolerance: 1.0e-8
  //   - Do not reorder

  using Teuchos::reduceAll;
  using Teuchos::REDUCE_MIN;
  using Teuchos::outArg;
  int lclSuccess = 1;
  int gblSuccess = 1;
  std::ostringstream errStrm;
  try {
    Ifpack2::Test::solve (A, b, x, 1, 1, 100, 1000, 1.0e-8, false, "RILUK");
  }
  catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << e.what ();
  }
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  if (gblSuccess != 1) {
    // We assume that it's OK for MPI processes other than Proc 0 in
    // MPI_COMM_WORLD to print to stderr.  That's generally true when
    // we run tests.
    cerr << "Belos solve with Ifpack2 preconditioner threw an exception "
      "on one or more processes!" << endl;
    for (int r = 0; r < numProcs; ++r) {
      if (r == myRank) {
        std::ostringstream os;
        os << "Process " << myRank << ": " << errStrm.str () << endl;
        cerr << os.str ();
      }
      comm->barrier (); // wait for output to finish
      comm->barrier ();
      comm->barrier ();
    }
  }

  if (comm->getRank () == 0) {
    if (gblSuccess == 1) {
      cout << "End Result: TEST PASSED" << endl;
    }
    else {
      cout << "End Result: TEST FAILED" << endl;
    }
  }
  return EXIT_SUCCESS;
}
