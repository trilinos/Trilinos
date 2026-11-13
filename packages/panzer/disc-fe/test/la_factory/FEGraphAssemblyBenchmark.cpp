// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Time building the FE graph of a TpetraLinearObjFactory on a structured hex mesh:
//   - "device":    TpetraLinearObjFactory::getFEGraph(), which assembles the graph on
//                  device (panzer::buildFEGraphGeometry + Tpetra::assembleFECrsGraphWithMaps);
//   - "insertion": the previous implementation, inserting each element's couplings into
//                  a Tpetra::FECrsGraph on host between beginAssembly() and endAssembly().
// Both use the same, already built, maps.  With --compare (the default), also check
// that the two graphs are identical.
//
// The fields are a velocity (three components, HGRAD_HEX_C<order>) and a pressure
// (HGRAD_HEX_C1), like a Taylor-Hood discretization.

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Tpetra_Core.hpp"
#include "Kokkos_Core.hpp"

#include "Panzer_TpetraLinearObjFactory.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_FEGraphGeometry.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"

#include "CartesianConnManager.hpp"
#include "FEGraphTestUtils.hpp"

#include <algorithm>
#include <iostream>
#include <limits>

using namespace panzer::fe_graph_test;
using LOFType = panzer::TpetraLinearObjFactory<panzer::Traits, double, LO, GO>;

namespace {
template <typename Intrepid2Type>
RCP<const panzer::FieldPattern> buildFieldPattern() {
  RCP<Intrepid2::Basis<PHX::Device, double, double>> basis = rcp(new Intrepid2Type);
  return rcp(new panzer::Intrepid2FieldPattern(basis));
}
}  // namespace

int main(int argc, char* argv[]) {
  Tpetra::ScopeGuard scope(&argc, &argv);
  {
    auto comm         = rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    const int myRank  = comm->getRank();
    const int numProc = comm->getSize();

    int nx = 40, ny = 40, nz = 40;
    int px = numProc, py = 1, pz = 1;
    int order = 2, repetitions = 3;
    bool compare = true;
    Teuchos::CommandLineProcessor clp;
    clp.setOption("nx", &nx, "Number of elements in the x direction.");
    clp.setOption("ny", &ny, "Number of elements in the y direction.");
    clp.setOption("nz", &nz, "Number of elements in the z direction.");
    clp.setOption("px", &px, "Number of processes in the x direction.");
    clp.setOption("py", &py, "Number of processes in the y direction.");
    clp.setOption("pz", &pz, "Number of processes in the z direction.");
    clp.setOption("order", &order, "Polynomial order of the velocity (1 or 2).");
    clp.setOption("repetitions", &repetitions, "Number of times to build each graph (the minimum time is reported).");
    clp.setOption("compare", "no-compare", &compare, "Check that the two graphs are identical.");
    if (clp.parse(argc, argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) return -1;
    if (px * py * pz != numProc || (order != 1 && order != 2)) {
      if (myRank == 0) std::cerr << "Need px*py*pz == number of processes, and order 1 or 2." << std::endl;
      return -1;
    }

    // The mesh and the DOFs.
    using CCM = panzer::unit_test::CartesianConnManager;
    RCP<CCM> connManager = rcp(new CCM);
    connManager->initialize(*comm, GO(nx), GO(ny), GO(nz), px, py, pz, 1, 1, 1);
    RCP<const panzer::FieldPattern> velocity =
        order == 2 ? buildFieldPattern<Intrepid2::Basis_HGRAD_HEX_C2_FEM<PHX::Device, double, double>>()
                   : buildFieldPattern<Intrepid2::Basis_HGRAD_HEX_C1_FEM<PHX::Device, double, double>>();
    RCP<const panzer::FieldPattern> pressure =
        buildFieldPattern<Intrepid2::Basis_HGRAD_HEX_C1_FEM<PHX::Device, double, double>>();
    RCP<panzer::DOFManager> dofManager = rcp(new panzer::DOFManager);
    dofManager->setConnManager(connManager, MPI_COMM_WORLD);
    dofManager->addField("UX", velocity);
    dofManager->addField("UY", velocity);
    dofManager->addField("UZ", velocity);
    dofManager->addField("P", pressure);
    dofManager->buildGlobalUnknowns();

    // Time a function: the maximum over processes, the minimum over repetitions.
    auto timeIt = [&](auto&& f) {
      double best = std::numeric_limits<double>::max();
      for (int r = 0; r < repetitions; ++r) {
        Kokkos::fence();
        comm->barrier();
        const double start = MPI_Wtime();
        f();
        Kokkos::fence();
        double elapsed = MPI_Wtime() - start, maxElapsed = 0;
        Teuchos::reduceAll<int, double>(*comm, Teuchos::REDUCE_MAX, elapsed, Teuchos::outArg(maxElapsed));
        best = std::min(best, maxElapsed);
      }
      return best;
    };

    // The factory's maps (built once, outside the timings).
    RCP<LOFType> factory = rcp(new LOFType(comm, dofManager, true));
    auto ownedRowMap = factory->getMap(), ownedPlusSharedRowMap = factory->getGhostedMap();
    auto ownedDomainMap = factory->getColMap(), ownedPlusSharedDomainMap = factory->getGhostedColMap();
    factory->getGhostedImport();

    RCP<FEGraph> deviceGraph, insertionGraph;
    const double deviceTime = timeIt([&] {
      // What getFEGraph() does (it caches its graph, so it can't be timed repeatedly).
      deviceGraph = Tpetra::assembleFECrsGraphWithMaps(
          panzer::buildFEGraphGeometry<LO, GO, NT>(*dofManager, *dofManager, *dofManager, connManager),
          ownedRowMap, ownedPlusSharedRowMap, Teuchos::null, ownedPlusSharedDomainMap, factory->getGhostedImport(),
          ownedDomainMap);
    });
    const double insertionTime = timeIt([&] {
      insertionGraph = insertionFEGraph(*dofManager, *dofManager, *dofManager, connManager, ownedRowMap,
                                        ownedPlusSharedRowMap, ownedPlusSharedDomainMap, ownedDomainMap);
    });

    bool ok = true;
    if (compare) {
      ok = identicalFEGraphs(*insertionGraph, *deviceGraph, *comm, std::cerr);
      // The factory's own graph, too.
      ok = identicalFEGraphs(*insertionGraph, *factory->getFEGraph(), *comm, std::cerr) && ok;
    }

    if (myRank == 0) {
      std::cout << "FE graph assembly: " << nx << "x" << ny << "x" << nz << " hexes, " << numProc
                << " processes, velocity order " << order << ", "
                << dofManager->getElementBlockGIDCount("eblock-0_0_0") << " DOFs per element, "
                << ownedRowMap->getGlobalNumElements() << " rows, " << deviceGraph->getGlobalNumEntries()
                << " entries" << std::endl;
      std::cout << "  insertion (host): " << insertionTime << " s" << std::endl;
      std::cout << "  device:           " << deviceTime << " s" << std::endl;
      std::cout << "  speedup:          " << insertionTime / deviceTime << std::endl;
      if (compare) std::cout << (ok ? "Graphs are identical." : "ERROR: graphs differ.") << std::endl;
      std::cout << (ok ? "End Result: TEST PASSED" : "End Result: TEST FAILED") << std::endl;
    }
    if (!ok) return 1;
  }
  return 0;
}
