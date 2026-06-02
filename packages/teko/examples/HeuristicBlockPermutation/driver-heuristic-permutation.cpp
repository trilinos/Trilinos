// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Teuchos includes
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

// Tpetra includes
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Core.hpp"
#include "MatrixMarket_Tpetra.hpp"

// Teko-Package includes
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_BlockedTpetraOperator.hpp"
#include "Teko_TpetraInverseFactoryOperator.hpp"
#include "Teko_HeuristicBlockPermutation.hpp"

#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraVectorSpace.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

// In a regualar user-app, one can essentially always include this assuming MueLu is enabled.
// Here, we use this ifdef since this particular driver is used as a test inside Teko,
// leading to circular library dependencies if we were to add a dependency to MueLu here.
#ifdef HAVE_MUELU_STRATIMIKOS
#include <Stratimikos_MueLuHelpers.hpp>
#endif

// Belos includes
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosSolverFactory_Tpetra.hpp"
#include "BelosSolverManager.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosThyraAdapter.hpp"  // Requires Stratimikos...

#include <iostream>

// for simplicity
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;

// Generates proc-by-proc block gid lists
template <class LO, class GO, class NO>
std::vector<std::vector<GO>> read_block_gids(std::string partitionFile,
                                             RCP<const Tpetra::Map<LO, GO, NO>> &rowMap) {
  using V   = Tpetra::Vector<LO, LO, GO, NO>;
  using CRS = Tpetra::CrsMatrix<LO, LO, GO, NO>;

  RCP<V> pfile =
      Tpetra::MatrixMarket::Reader<CRS>::readVectorFile(partitionFile, rowMap->getComm(), rowMap);

  int num_blocks = 1 + pfile->normInf();

  if (!rowMap->getComm()->getRank())
    std::cout << "Reading partition file: Found " << num_blocks << " blocks" << std::endl;

  std::vector<std::vector<GO>> block_gids(num_blocks);
  auto vv = pfile->get1dView();
  for (LO i = 0; i < (LO)vv.size(); i++) {
    LO block_id = vv[i];
    block_gids[block_id].push_back(rowMap->getGlobalElement(i));
  }

  return block_gids;
}

template <class SC, class LO, class GO, class NO>
void ReadSplittingFromDisk(const std::string &partitionFile,
                           const RCP<Tpetra::CrsMatrix<SC, LO, GO, NO>> &crsMat,
                           RCP<Teko::TpetraHelpers::BlockedTpetraOperator> &A,
                           std::vector<std::vector<GO>> &block_gids) {
  // The format for the partition file is a number of lines of the form:
  // proc_id block_id gid0, gid1,...
  //
  // Blocks do not need to be uniquely owned by a processor.

  //  Each rank is going to read the file, one at a time (to avoid hammering on the disk)
  auto comm                                 = crsMat->getRowMap()->getComm();
  RCP<const Tpetra::Map<LO, GO, NO>> rowmap = crsMat->getRowMap();
  block_gids                                = read_block_gids<LO, GO, NO>(partitionFile, rowmap);

  RCP<Teko::TpetraHelpers::BlockedTpetraOperator> rA =
      rcp(new Teko::TpetraHelpers::BlockedTpetraOperator(block_gids, crsMat));

  A = rA;
}

template <class SC, class LO, class GO, class NO>
Teuchos::ParameterList augment_muelu_parameters_with_coordinates(
    const Teuchos::ParameterList &initial_params,
    RCP<Tpetra::MultiVector<SC, LO, GO, NO>> &coords) {
  Teuchos::ParameterList params_with_coords(initial_params);
  if (!coords) return params_with_coords;

  auto thyra_coords = Thyra::createMultiVector(coords);
  for (auto &&param : params_with_coords) {
    if (!params_with_coords.isSublist(param.key)) continue;

    auto &sublist = params_with_coords.sublist(param.key, true);
    if (sublist.isParameter("Type") && sublist.get<std::string>("Type") == "MueLu") {
      auto &muelu_sublist = sublist.sublist("user data");
      muelu_sublist.set("Coordinates", thyra_coords);
    }
  }

  return params_with_coords;
}

// Optional ParameterList to generate_heuristic_permutation, but provided here nevertheless for
// demonstration purposes
Teuchos::ParameterList default_heuristic_settings() {
  Teuchos::ParameterList tekoHeuristicSettings;
  tekoHeuristicSettings.set("Heuristic Method", "Greedy Block Merging Heuristic");
  tekoHeuristicSettings.set("Block Inverse Type", "Block Gauss-Seidel");
  tekoHeuristicSettings.set("Use Block Upper Triangle", true);

  const double maxWalltime = 0.01;
  tekoHeuristicSettings.set("Max Heuristic Walltime", maxWalltime);

  const double targetNormLoss = 1e-3;
  tekoHeuristicSettings.set("Target Norm Loss", targetNormLoss);
  return tekoHeuristicSettings;
}

template <class SC, class LO, class GO, class NO>
int solve_tpetra(RCP<Tpetra::CrsMatrix<SC, LO, GO, NO>> &crsMat,
                 RCP<Tpetra::Vector<SC, LO, GO, NO>> &b,
                 RCP<Tpetra::MultiVector<SC, LO, GO, NO>> &coords,
                 const std::string &partitionFile) {
  typedef Tpetra::MultiVector<SC, LO, GO, NO> MV;
  typedef Tpetra::Operator<SC, LO, GO, NO> OP;
  auto comm = crsMat->getRowMap()->getComm();

  // tell Stratimikos => Teko about MueLu
  RCP<Stratimikos::DefaultLinearSolverBuilder> linearSolverBuilder =
      Teuchos::rcp(new Stratimikos::DefaultLinearSolverBuilder);

// In a regualar user-app, one can essentially always include the `Stratimikos::enableMueLu` call
// below. Here, we use ifdef this call since this particular driver is used as a test inside Teko,
// leading to circular library dependencies if we were to add a dependency to MueLu here.
#ifdef HAVE_MUELU_STRATIMIKOS
  Stratimikos::enableMueLu<SC, LO, GO, NO>(*linearSolverBuilder);
#endif

  /////////////////////////////////////////////////////////
  // Build the Thyra operators
  /////////////////////////////////////////////////////////
  RCP<Teko::TpetraHelpers::BlockedTpetraOperator> A;
  RCP<Tpetra::Vector<SC, LO, GO, NO>> x = rcp(new Tpetra::Vector<SC, LO, GO, NO>(b->getMap()));

  std::vector<std::vector<GO>> dof_block_gids;
  ReadSplittingFromDisk<SC, LO, GO, NO>(partitionFile, crsMat, A, dof_block_gids);
  x->putScalar(Teuchos::ScalarTraits<SC>::zero());

  RCP<Thyra::MultiVectorBase<SC>> xt = Thyra::createVector(x);
  RCP<Thyra::MultiVectorBase<SC>> bt = Thyra::createVector(b);

  /////////////////////////////////////////////////////////
  // Build the preconditioner
  /////////////////////////////////////////////////////////

  // Optional
  Teuchos::ParameterList tekoHeuristicSettings = default_heuristic_settings();
  auto [permutation, score] =
      Teko::generate_heuristic_permutation(A, Teuchos::rcpFromRef(tekoHeuristicSettings));

  auto block_gids = Teko::construct_block_gids_from_permutation(permutation, dof_block_gids);

  // Reconstruct blocked representation with updated block gids from heuristic permutation
  A = rcp(new Teko::TpetraHelpers::BlockedTpetraOperator(block_gids, crsMat));

  const std::string tekoInverseName = "MyTekoPreconditioner";

  // build an InverseLibrary
  RCP<Teuchos::ParameterList> xmlList =
      Teko::generate_parameters_from_permutation(permutation, tekoInverseName);

  // Add coordinates if we need to (user app provides this information)
  Teuchos::ParameterList paramsWithCoords =
      augment_muelu_parameters_with_coordinates(*xmlList, coords);

  RCP<Teko::InverseLibrary> invLib =
      Teko::InverseLibrary::buildFromParameterList(paramsWithCoords, linearSolverBuilder);
  RCP<Teko::InverseFactory> inverse = invLib->getInverseFactory(tekoInverseName);
  auto preconditioner = Teuchos::rcp(new Teko::TpetraHelpers::InverseFactoryOperator(inverse));
  preconditioner->initInverse();

  auto blocked_operator =
      Teuchos::rcp_dynamic_cast<const Tpetra::Operator<Teko::ST, Teko::LO, Teko::GO, Teko::NT>>(A);
  preconditioner->rebuildInverseOperator(blocked_operator);

  // Setup the Belos solver
  /////////////////////////////////////////////////////////

  Teuchos::ParameterList belosList;
  belosList.set("Num Blocks", 200);              // Maximum number of blocks in Krylov factorization
  belosList.set("Block Size", 1);                // Blocksize to be used by iterative solver
  belosList.set("Maximum Iterations", 200);      // Maximum number of iterations allowed
  belosList.set("Maximum Restarts", 1);          // Maximum number of restarts allowed
  belosList.set("Convergence Tolerance", 1e-8);  // Relative convergence tolerance requested
  belosList.set(
      "Verbosity",
      33);  // Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails );
  belosList.set("Output Frequency", 1);
  belosList.set("Output Style", 1);
  belosList.set("Flexible Gmres", true);  // NOTE: required due to sub-iterative schemes

  Teuchos::RCP<const OP> tpetraOp = crsMat;
  auto problem                    = rcp(new Belos::LinearProblem<double, MV, OP>(tpetraOp, x, b));
  problem->setRightPrec(preconditioner);
  problem->setProblem();

  RCP<Belos::SolverManager<double, MV, OP>> solver =
      rcp(new Belos::BlockGmresSolMgr<double, MV, OP>(problem, rcpFromRef(belosList)));

  //
  // Perform solve
  //
  Belos::ReturnType ret = solver->solve();

  if (ret != Belos::Converged) {
    std::cout << std::endl << "ERROR:  Belos did not converge!" << std::endl;
    return -1;
  } else {
    std::cout << std::endl << " Belos converged" << std::endl;
  }

  return 0;
}

int main(int argc, char *argv[]) {
  typedef double SC;

  typedef Tpetra::Map<> TP_Map;
  typedef Tpetra::Vector<SC> TP_Vec;
  typedef Tpetra::MultiVector<SC> TP_MV;
  typedef Tpetra::CrsMatrix<SC> TP_Crs;
  // typedef Tpetra::Operator<SC>        TP_Op;

  typedef TP_Vec::local_ordinal_type LO;
  typedef TP_Vec::global_ordinal_type GO;
  typedef TP_Vec::node_type NO;

  // calls MPI_Init and MPI_Finalize
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Parse CLI options
  Teuchos::CommandLineProcessor clp(false);
  std::string rhsFile;
  clp.setOption("rhs", &rhsFile, "rhs data file");

  std::string rowMapFile;
  clp.setOption("rowmap", &rowMapFile, "map data file");
  std::string matrixFile = "../data/nsjac.mm";
  clp.setOption("matrix", &matrixFile, "matrix data file");

  std::string partitionFile = "";
  clp.setOption("partition", &partitionFile, "partition file which defines the blocks");

  std::string coordsFile;
  clp.setOption("coords", &coordsFile, "coords data file");

  std::string coordsMapFile;
  clp.setOption("coordsmap", &coordsMapFile, "coords data file");

  clp.recogniseAllOptions(true);
  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
  }

  /////////////////////////////////////////////////////////
  // Build the Tpetra matrices and vectors
  /////////////////////////////////////////////////////////

  // read in the CRS matrix
  RCP<TP_Crs> crsMat;
  RCP<TP_Vec> rhs;
  RCP<TP_MV> coords;
  if (rowMapFile == "") {
    crsMat =
        Tpetra::MatrixMarket::Reader<TP_Crs>::readSparseFile(matrixFile, Tpetra::getDefaultComm());
  } else {
    RCP<const TP_Map> rowMap =
        Tpetra::MatrixMarket::Reader<TP_Crs>::readMapFile(rowMapFile, Tpetra::getDefaultComm());
    RCP<const TP_Map> colMap;
    crsMat = Tpetra::MatrixMarket::Reader<TP_Crs>::readSparseFile(matrixFile, rowMap, colMap,
                                                                  rowMap, rowMap);
  }

  auto rowMap = crsMat->getRowMap();

  if (rhsFile != "") {
    rhs = Tpetra::MatrixMarket::Reader<TP_Crs>::readVectorFile(rhsFile, rowMap->getComm(), rowMap);
  }

  if (coordsFile != "") {
    RCP<const TP_Map> coordsMap = rowMap;
    if (coordsMapFile != "")
      coordsMap =
          Tpetra::MatrixMarket::Reader<TP_Crs>::readMapFile(coordsMapFile, rowMap->getComm());
    coords = Tpetra::MatrixMarket::Reader<TP_Crs>::readDenseFile(coordsFile, coordsMap->getComm(),
                                                                 coordsMap);
  }

  // Sanity checks
  if (rhs.is_null()) throw std::runtime_error("rhs is null");
  if (crsMat.is_null()) throw std::runtime_error("crsMat is null");

  auto comm = crsMat->getRowMap()->getComm();
  int rv    = solve_tpetra<SC, LO, GO, NO>(crsMat, rhs, coords, partitionFile);

  return rv;
}