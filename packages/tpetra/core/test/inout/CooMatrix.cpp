// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Details_CooMatrix.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#ifdef HAVE_TPETRACORE_MPI
#  include "Tpetra_Details_extractMpiCommFromTeuchos.hpp"
#  include "Teuchos_DefaultMpiComm.hpp"
#endif // HAVE_TPETRACORE_MPI

namespace { // (anonymous)

using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::rcp;
using std::endl;
typedef double SC;
typedef Tpetra::DistObject<char>::local_ordinal_type LO;
typedef Tpetra::DistObject<char>::global_ordinal_type GO;
typedef Tpetra::global_size_t GST;
typedef Tpetra::Export<> export_type;
typedef Tpetra::Map<> map_type;

void
testCooMatrix (bool& success,
               Teuchos::FancyOStream& out,
               const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
{
  using Tpetra::Details::CooMatrix;
  using Teuchos::outArg;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  int lclSuccess = 1;
  int gblSuccess = 0; // output argument

  out << "Test CooMatrix" << endl;
  Teuchos::OSTab tab1 (out);

  TEST_ASSERT( comm->getSize () >= 2 );
  if (comm->getSize () < 2) {
    out << "This test needs at least 2 MPI processes!" << endl;
    return;
  }

  out << "CooMatrix default constructor" << endl;
  CooMatrix<SC, LO, GO> A_in;
  TEST_ASSERT( A_in.getMap ().is_null () );
  TEST_EQUALITY( A_in.getLclNumEntries (), static_cast<std::size_t> (0) );

  out << "Add entries locally to CooMatrix" << endl;
  const int myRank = comm->getRank ();
  if (myRank == 0) {
    A_in.sumIntoGlobalValues ({666, 31, 31, 31}, {11, 6, 5, 6}, {-1.0, 1.0, 2.0, 111.0}, 4);
    TEST_EQUALITY( A_in.getLclNumEntries (), static_cast<std::size_t> (3) );
  }
  else if (myRank == 1) {
    A_in.sumIntoGlobalValues ({418, 31}, {11, 5}, {11.0, 5.0}, 2);
    TEST_EQUALITY( A_in.getLclNumEntries (), static_cast<std::size_t> (2) );
  }

  lclSuccess = success ? 1 : 0;
  gblSuccess = 0; // output argument
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
  if (gblSuccess != 1) {
    out << "A_in not in a consistent state before fillComplete, "
      "so don't bother continuing the test." << endl;
    return;
  }

  out << "Call fillComplete on CooMatrix" << endl;
  A_in.fillComplete (comm);
  TEST_ASSERT( ! A_in.getMap ().is_null () );

  out << "Create output Map" << endl;
  RCP<const map_type> outMap;
  const GO indexBase = 31; // the smallest global index in the Map
  const GST numGblInds = 3;
  if (myRank == 0) {
    const GO myGblInds[] = {418, 666};
    const LO numLclInds = 2;
    outMap = rcp (new map_type (numGblInds, myGblInds, numLclInds, indexBase, comm));
  }
  else if (myRank == 1) {
    const GO myGblInds[] = {31};
    const LO numLclInds = 1;
    outMap = rcp (new map_type (numGblInds, myGblInds, numLclInds, indexBase, comm));
  }
  else {
    const GO* myGblInds = NULL;
    const LO numLclInds = 0;
    outMap = rcp (new map_type (numGblInds, myGblInds, numLclInds, indexBase, comm));
  }

  out << "Create output CooMatrix" << endl;
  CooMatrix<SC, LO, GO> A_out (outMap);
  TEST_EQUALITY( A_out.getLclNumEntries (), static_cast<std::size_t> (0) );
  TEST_ASSERT( ! A_out.getMap ().is_null () );
  const bool outMapsSame = outMap->isSameAs (* (A_out.getMap ()));
  TEST_ASSERT( outMapsSame );
  const bool outMapIsOneToOne = outMap->isOneToOne ();
  TEST_ASSERT( outMapIsOneToOne );

  out << "Create Export object" << endl;
  export_type exporter (A_in.getMap (), A_out.getMap ());

  out << "Call doExport on CooMatrix" << endl;
  A_out.doExport (A_in, exporter, Tpetra::ADD);

  out << "Test global success" << endl;
  lclSuccess = (success && ! A_out.localError ()) ? 1 : 0;
  gblSuccess = 0; // output argument
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));

  TEST_EQUALITY( gblSuccess, 1 );
  if (gblSuccess != 1) {
    std::ostringstream os;
    os << "Process " << myRank << ": " << A_out.errorMessages () << endl;
    Tpetra::Details::gathervPrint (out, os.str (), *comm);
  }
}

TEUCHOS_UNIT_TEST( CooMatrix, doubleIntLongLong )
{
  using Tpetra::TestingUtilities::getDefaultComm;
  RCP<const Comm<int> > comm = getDefaultComm ();
  TEST_ASSERT( ! comm.is_null () );
  if (comm.is_null ()) {
    return;
  }

  // Throw away map, just to make sure that Kokkos is initialized
  RCP<const map_type> throwaway_map;
  throwaway_map = rcp(new map_type(static_cast<GST>(0),
                                   static_cast<GO>(0),
                                   comm));

#ifdef HAVE_TPETRACORE_MPI
  // Set the MPI error handler so that errors return, instead of
  // immediately causing MPI_Abort.  This will help us catch any bugs
  // with how we use, e.g., MPI_Pack and MPI_Unpack.
  {
    using Teuchos::MpiComm;
    using Teuchos::rcp_const_cast;
    using Teuchos::rcp_dynamic_cast;

    constexpr bool throwOnFail = true;
    auto mpiComm = rcp_dynamic_cast<const MpiComm<int> > (comm, throwOnFail);
    // We have to cast away const to call setErrorHandler.
    auto mpiCommNonConst = rcp_const_cast<MpiComm<int> > (mpiComm);
    auto errHandler =
      rcp (new Teuchos::OpaqueWrapper<MPI_Errhandler> (MPI_ERRORS_RETURN));
    mpiCommNonConst->setErrorHandler (errHandler);
  }
#endif // HAVE_TPETRACORE_MPI


  testCooMatrix (success, out, comm);
}

} // namespace (anonymous)
