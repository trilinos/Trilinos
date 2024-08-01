// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_Comm.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include <memory>
#include <sstream>

namespace {

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Teuchos::OrdinalTraits;
  using Teuchos::Comm;
  using Teuchos::null;
  using Tpetra::Import;
  using Tpetra::Map;
  using Tpetra::Vector;
  using Tpetra::global_size_t;

  bool testMpi = true;
  double errorTolSlack = 1e+1;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST( MultiVector, BlockedViews_7234 )
  {
    // Test that replicates Trilinos issue #7234 
    // https://github.com/trilinos/Trilinos/issues/7234
    // On CUDA platforms, subviews of Kokkos::DualView did not perform properly
    // when pointing to the end of a DualView, as the shared Vector does
    // below.  See Kokkos issues #2981 and #2979 for more details
    // https://github.com/kokkos/kokkos/issues/2981
    // https://github.com/kokkos/kokkos/issues/2979
    using Teuchos::RCP;
    using Teuchos::rcp;
    using namespace Tpetra;
    using namespace Kokkos;
    using GST = Tpetra::global_size_t;
    using LO = int;
    using GO = Tpetra::Map<>::global_ordinal_type;
    using SC = Tpetra::Vector<>::scalar_type;
    using NodeType = Tpetra::Map<>::node_type;
    using TpetraMap = Tpetra::Map<LO,GO,NodeType>;
    using TpetraImport = Tpetra::Import<LO,GO,NodeType>;
    using TpetraExport = Tpetra::Export<LO,GO,NodeType>;
    using TpetraVec = Tpetra::Vector<SC,LO,GO,NodeType>;

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
    const int rank = comm->getRank();
    const int size = comm->getSize();

    // Build the array: GIDs are BLOCKED - owned first, then shared
    Teuchos::Array<GO> owned_gids;
    Teuchos::Array<GO> shared_gids;
    Teuchos::Array<GO> owned_and_shared_gids;
    owned_gids.push_back(2*rank);
    owned_gids.push_back(2*rank+1);
    owned_and_shared_gids.push_back(2*rank);
    owned_and_shared_gids.push_back(2*rank+1);
    if (rank > 0) {
      shared_gids.push_back(owned_gids[0]-1);
      owned_and_shared_gids.push_back(shared_gids[0]);
    }

    out.setShowProcRank(true);
    out.setOutputToRootOnly(rank);
    out << "owned_map\n";
    for (int i=0; i < owned_gids.size(); ++i)
      out << i << " " << owned_gids[i] << std::endl;
    out << "owned_and_shared_map\n";
    for (int i=0; i < owned_and_shared_gids.size(); ++i)
      out << i << " " << owned_and_shared_gids[i] << std::endl;
    out.setOutputToRootOnly(0);

    RCP<TpetraMap> owned_map = rcp(new TpetraMap(Teuchos::OrdinalTraits<GST>::invalid(),owned_gids,GO(0),comm));
    RCP<TpetraMap> shared_map = rcp(new TpetraMap(Teuchos::OrdinalTraits<GST>::invalid(),shared_gids,GO(0),comm));
    RCP<TpetraMap> owned_and_shared_map = rcp(new TpetraMap(Teuchos::OrdinalTraits<GST>::invalid(),owned_and_shared_gids,GO(0),comm));

    bool zeroOut = true;
    RCP<TpetraVec> owned_and_shared = 
                   rcp(new TpetraVec(owned_and_shared_map,zeroOut));
    RCP<TpetraExport> exporter = rcp(new TpetraExport(shared_map,owned_map));
    RCP<TpetraImport> importer = rcp(new TpetraImport(owned_map,shared_map));

    // If we block the owned_and_shared, we can make the code much
    // more efficient with subviews.
    RCP<TpetraVec> owned = 
                  owned_and_shared->offsetViewNonConst(owned_map,0);
    RCP<TpetraVec> shared = 
                  owned_and_shared->offsetViewNonConst(shared_map,
                                                       owned->getLocalLength());

    owned->putScalar(1.0);
    shared->doImport(*owned, *importer, Tpetra::REPLACE);

    const double tol = Teuchos::ScalarTraits<double>::eps() * 100.0;
    {
      auto host_shared = shared->getLocalViewHost(Tpetra::Access::ReadOnly);
      if (rank > 0) {
        TEST_FLOATING_EQUALITY(host_shared(0,0),1.0,tol);
      }

      auto host_owned_and_shared = 
                owned_and_shared->getLocalViewHost(Tpetra::Access::ReadOnly);
      TEST_FLOATING_EQUALITY(host_owned_and_shared(0,0),1.0,tol);
      TEST_FLOATING_EQUALITY(host_owned_and_shared(1,0),1.0,tol);
      if (rank > 0) {
        TEST_FLOATING_EQUALITY(host_owned_and_shared(2,0),1.0,tol);
      }
    }

    // now test the export
    owned->doExport(*shared, *exporter, Tpetra::ADD);

    {
      auto host_owned = owned->getLocalViewHost(Tpetra::Access::ReadOnly);
      auto host_owned_and_shared = 
                owned_and_shared->getLocalViewHost(Tpetra::Access::ReadOnly);
      TEST_FLOATING_EQUALITY(host_owned(0,0),1.0,tol);
      // check owned entries only
      TEST_FLOATING_EQUALITY(host_owned_and_shared(0,0),1.0,tol); 
      if (rank != size-1) {
        TEST_FLOATING_EQUALITY(host_owned(1,0),2.0,tol);
        // check owned entries only
        TEST_FLOATING_EQUALITY(host_owned_and_shared(1,0),2.0,tol); 
      }
      else {
        TEST_FLOATING_EQUALITY(host_owned(1,0),1.0,tol);
        // check owned entries only
        TEST_FLOATING_EQUALITY(host_owned_and_shared(1,0),1.0,tol); 
      }
    }
  }
}
