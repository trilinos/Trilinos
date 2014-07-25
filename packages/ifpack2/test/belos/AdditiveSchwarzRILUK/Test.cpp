#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Version.hpp>

#include "Solve.hpp"

int
main (int argc, char **argv)
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

  Teuchos::oblackholestream blackHole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackHole);
  RCP<const Teuchos::Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();

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

  // Solve the linear system.
  Ifpack2::Test::solve (A, b, x, 1, 1, 100, 1000, 1.0e-8, false, "RILUK");

  if (comm->getRank () == 0) {
    cout << "End Result: TEST PASSED" << endl;
  }
  return EXIT_SUCCESS;
}
