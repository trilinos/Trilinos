// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include <Tpetra_Core.hpp>
#include <Tpetra_Version.hpp>
#include <MatrixMarket_Tpetra.hpp>

#include "TypedefsAndIncludes.hpp"
#include "Solve.hpp"

int main(int argc, char** argv) {

  typedef Tpetra::MatrixMarket::Reader<Ifpack2::Test::sparse_mat_type> reader_type;

  // global_size_t: Tpetra defines this unsigned integer type big
  // enough to hold any global dimension or amount of data.
  using std::cerr;
  using std::cout;
  using std::endl;
  using Teuchos::arcp;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
  using Tpetra::global_size_t;

  Tpetra::ScopeGuard tpetraScope(&argc, &argv);
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const int myRank                    = comm->getRank();
  const int numProcs                  = comm->getSize();

  Teuchos::RCP<Ifpack2::Test::sparse_mat_type> A   = reader_type::readSparseFile("opamp_mtx1.mm", comm);
  Teuchos::RCP<Ifpack2::Test::sparse_mat_type> A2  = reader_type::readSparseFile("opamp_mtx2.mm", comm);
  RCP<const Ifpack2::Test::map_type> map = A->getRangeMap();
  Teuchos::RCP<Ifpack2::Test::multivector_type> b  = reader_type::readDenseFile("opamp_rhs1.mm", comm, map); 
  Teuchos::RCP<Ifpack2::Test::multivector_type> b2 = reader_type::readDenseFile("opamp_rhs2.mm", comm, map); 

  // Attempt to solve the linear system, with the following parameters:
  //   - Amesos2 - KLU subdomain solver
  //   - GMRES restart length: 100
  //   - Maximum number of GMRES iterations (over all restarts): 1000
  //   - GMRES relative residual tolerance: 1.0e-8
  //   - Do not reorder

  using Teuchos::outArg;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  int lclSuccess = 1;
  int gblSuccess = 1;
  int overlapLevel = ((numProcs==1)? 0 : 1);
  std::ostringstream errStrm;
  bool expectOneIteration = ((numProcs==1)? true : false);
  bool success = false;
  try {
    success = Ifpack2::Test::solveTwice(*A, *A2, *b, *b2, overlapLevel, 100, 1000,
                                        1.0e-8, false, "AMESOS2", 
                                        expectOneIteration);
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << e.what();
  }
  reduceAll<int, int>(*comm, REDUCE_MIN, lclSuccess, outArg(gblSuccess));
  if (gblSuccess != 1) {
    // We assume that it's OK for MPI processes other than Proc 0 in
    // MPI_COMM_WORLD to print to stderr.  That's generally true when
    // we run tests.
    cerr << "Belos solve with Ifpack2 Amesos2 preconditioner threw an exception "
            "on one or more processes!"
         << endl;
    for (int r = 0; r < numProcs; ++r) {
      if (r == myRank) {
        std::ostringstream os;
        os << "Process " << myRank << ": " << errStrm.str() << endl;
        cerr << os.str();
      }
      comm->barrier();  // wait for output to finish
      comm->barrier();
      comm->barrier();
    }
  }

  // If the test did not succeed, then global success is false (0);
  if (!success && gblSuccess)
    gblSuccess = 0;

  if (comm->getRank() == 0) {
    if (gblSuccess == 1) {
      cout << "End Result: TEST PASSED" << endl;
    } else {
      cout << "End Result: TEST FAILED" << endl;
    }
  }
  return EXIT_SUCCESS;
}
