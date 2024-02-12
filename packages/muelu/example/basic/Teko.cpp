/*
// @HEADER
//
// ***********************************************************************
//
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

*/

// Teuchos includes /*@ \label{lned:being-includes} @*/
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

// Tpetra includes
#include "mpi.h"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Core.hpp"
#include "MatrixMarket_Tpetra.hpp"

// Teko-Package includes
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_BlockedTpetraOperator.hpp"
#include "Teko_LSCPreconditionerFactory.hpp"
#include "Teko_InvLSCStrategy.hpp"
#include "Teko_SIMPLEPreconditionerFactory.hpp"
#include "Teko_TpetraInverseFactoryOperator.hpp"

#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraVectorSpace.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Stratimikos_MueLuHelpers.hpp"

// Belos includes
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosThyraAdapter.hpp"  // Requires Stratimikos...

#include <iostream> /*@ \label{lned:end-includes} @*/

// for simplicity
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;

RCP<Teuchos::ParameterList> buildLibPL();

template <class SC, class LO, class GO, class NO>
void GenerateDefault_2x2_Splitting(const RCP<Tpetra::CrsMatrix<SC, LO, GO, NO> > &crsMat,
                                   Teko::LinearOp &A,
                                   Teko::MultiVector &x,
                                   Teko::MultiVector &b) {
  typedef Tpetra::Vector<SC, LO, GO, NO> TP_Vec;
  typedef Tpetra::CrsMatrix<SC, LO, GO, NO> TP_Crs;
  typedef Tpetra::Operator<SC, LO, GO, NO> TP_Op;
  RCP<TP_Crs> zeroCrsMat = rcp(new TP_Crs(*crsMat, Teuchos::Copy));
  zeroCrsMat->setAllToScalar(0.0);

  RCP<TP_Op> Mat     = crsMat;
  RCP<TP_Op> zeroMat = zeroCrsMat;

  // Allocate some right handside vectors
  RCP<TP_Vec> x0_tp = rcp(new TP_Vec(Mat->getDomainMap()));
  RCP<TP_Vec> x1_tp = rcp(new TP_Vec(Mat->getDomainMap()));
  RCP<TP_Vec> b0_tp = rcp(new TP_Vec(Mat->getRangeMap()));
  RCP<TP_Vec> b1_tp = rcp(new TP_Vec(Mat->getRangeMap()));
  b0_tp->randomize();
  b1_tp->randomize();

  RCP<const Thyra::TpetraVectorSpace<SC, LO, GO, NO> > domain = Thyra::tpetraVectorSpace<SC>(Mat->getDomainMap());
  RCP<const Thyra::TpetraVectorSpace<SC, LO, GO, NO> > range  = Thyra::tpetraVectorSpace<SC>(Mat->getRangeMap());

  /////////////////////////////////////////////////////////
  // Build Teko compatible matrices and vectors
  /////////////////////////////////////////////////////////

  // convert them to teko compatible sub vectors
  Teko::MultiVector x0_th = Thyra::tpetraVector(domain, x0_tp);
  Teko::MultiVector x1_th = Thyra::tpetraVector(domain, x1_tp);
  Teko::MultiVector b0_th = Thyra::tpetraVector(range, b0_tp);
  Teko::MultiVector b1_th = Thyra::tpetraVector(range, b1_tp);
  std::vector<Teko::MultiVector> x_vec;
  x_vec.push_back(x0_th);
  x_vec.push_back(x1_th);
  std::vector<Teko::MultiVector> b_vec;
  b_vec.push_back(b0_th);
  b_vec.push_back(b1_th);

  x = Teko::buildBlockedMultiVector(x_vec);  // these will be used in the Teko solve
  b = Teko::buildBlockedMultiVector(b_vec);

  // Build the Teko compatible linear system
  Teko::LinearOp thMat  = Thyra::tpetraLinearOp<double>(range, domain, Mat);
  Teko::LinearOp thZero = Thyra::tpetraLinearOp<double>(range, domain, zeroMat);
  A                     = Thyra::block2x2(thMat, thZero, thZero, thMat);  // build an upper triangular 2x2
}

int extract_int(std::istringstream &iss) {
  std::string s;
  iss >> s;
  if (s != "") {
    return stoi(s);
  } else
    return -1;
}

// Generates proc-by-proc block gid lists
template <class LO, class GO, class NO>
std::vector<std::vector<GO> > read_block_gids(std::string partitionFile, RCP<const Tpetra::Map<LO, GO, NO> > &rowMap) {
  using V   = Tpetra::Vector<LO, LO, GO, NO>;
  using CRS = Tpetra::CrsMatrix<LO, LO, GO, NO>;

  RCP<V> pfile = Tpetra::MatrixMarket::Reader<CRS>::readVectorFile(partitionFile, rowMap->getComm(), rowMap);

  int num_blocks = 1 + pfile->normInf();

  if (!rowMap->getComm()->getRank())
    std::cout << "Reading partition file: Found " << num_blocks << " blocks" << std::endl;

  std::vector<std::vector<GO> > block_gids(num_blocks);
  auto vv = pfile->get1dView();
  for (LO i = 0; i < (LO)vv.size(); i++) {
    LO block_id = vv[i];
    block_gids[block_id].push_back(rowMap->getGlobalElement(i));
  }

  return block_gids;
}

template <class SC, class LO, class GO, class NO>
void ReadSplittingFromDisk(const std::string &partitionFile,
                           const RCP<Tpetra::CrsMatrix<SC, LO, GO, NO> > &crsMat,
                           RCP<Tpetra::Operator<SC, LO, GO, NO> > &A) {
  // The format for the partition file is a number of lines of the form:
  // proc_id block_id gid0, gid1,...
  //
  // Blocks do not need to be uniquely owned by a processor.

  //  Each rank is going to read the file, one at a time (to avoid hammering on the disk)
  auto comm                                        = crsMat->getRowMap()->getComm();
  RCP<const Tpetra::Map<LO, GO, NO> > rowmap       = crsMat->getRowMap();
  std::vector<std::vector<GO> > my_blocks_and_gids = read_block_gids<LO, GO, NO>(partitionFile, rowmap);

  RCP<Teko::TpetraHelpers::BlockedTpetraOperator> rA = rcp(new Teko::TpetraHelpers::BlockedTpetraOperator(my_blocks_and_gids, crsMat));

  A = rA;
}

template <class SC, class LO, class GO, class NO>
int solve_thyra(RCP<Tpetra::CrsMatrix<SC, LO, GO, NO> > &crsMat, const std::string &xmlFile) {
  typedef Tpetra::CrsMatrix<SC> TP_Crs;

  typedef Thyra::MultiVectorBase<SC> MV;
  typedef Thyra::LinearOpBase<SC> OP;
  auto comm = crsMat->getRowMap()->getComm();

  RCP<Stratimikos::DefaultLinearSolverBuilder> linearSolverBuilder = Teuchos::rcp(new Stratimikos::DefaultLinearSolverBuilder);

  /////////////////////////////////////////////////////////
  // Build the Thyra operators
  /////////////////////////////////////////////////////////
  Teko::LinearOp A;
  Teko::MultiVector x, b;
  GenerateDefault_2x2_Splitting<SC, LO, GO, NO>(crsMat, A, x, b);

  /////////////////////////////////////////////////////////
  // Build the preconditioner
  /////////////////////////////////////////////////////////

  // build an InverseLibrary
  RCP<Teko::InverseLibrary> invLib;
  RCP<Teko::InverseFactory> inverse;
  if (xmlFile == "") {
    invLib  = Teko::InverseLibrary::buildFromParameterList(*buildLibPL(), linearSolverBuilder);
    inverse = invLib->getInverseFactory("Gauss-Seidel");
  } else {
    Teuchos::ParameterList xmlList;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFile, Teuchos::Ptr<Teuchos::ParameterList>(&xmlList), *comm);
    invLib  = Teko::InverseLibrary::buildFromParameterList(xmlList, linearSolverBuilder);
    inverse = invLib->getInverseFactory("MyTekoPreconditioner");
  }
  // build the inverse factory needed by the example preconditioner

  // build the preconditioner from the jacobian
  Teko::LinearOp prec = Teko::buildInverse(*inverse, A);

  // Setup the Belos solver
  /////////////////////////////////////////////////////////

  Teuchos::ParameterList belosList;
  belosList.set("Num Blocks", 200);              // Maximum number of blocks in Krylov factorization
  belosList.set("Block Size", 1);                // Blocksize to be used by iterative solver
  belosList.set("Maximum Iterations", 200);      // Maximum number of iterations allowed
  belosList.set("Maximum Restarts", 1);          // Maximum number of restarts allowed
  belosList.set("Convergence Tolerance", 1e-5);  // Relative convergence tolerance requested
  belosList.set("Verbosity", 33);                // Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails );
  belosList.set("Output Frequency", 1);
  belosList.set("Output Style", 1);

  RCP<Belos::LinearProblem<double, MV, OP> > problem = rcp(new Belos::LinearProblem<double, MV, OP>(A, x, b));
  problem->setLeftPrec(prec);
  problem->setProblem();  // should check the return type!!!

  RCP<Belos::SolverManager<double, MV, OP> > solver = rcp(new Belos::BlockGmresSolMgr<double, MV, OP>(problem, rcpFromRef(belosList)));

  //
  // Perform solve
  //
  Belos::ReturnType ret = solver->solve();

  if (ret != Belos::Converged) {
    std::cout << std::endl
              << "ERROR:  Belos did not converge!" << std::endl;
    return -1;
  }

  return 0;
}

template <class SC, class LO, class GO, class NO>
int solve_tpetra(RCP<Tpetra::CrsMatrix<SC, LO, GO, NO> > &crsMat, RCP<Tpetra::Vector<SC, LO, GO, NO> > &b, RCP<Tpetra::MultiVector<SC, LO, GO, NO> > &coords, const std::string &xmlFile, const std::string &partitionFile) {
  // typedef Tpetra::Map<>               TP_Map;
  // typedef Tpetra::Vector<SC>          TP_Vec;
  typedef Tpetra::CrsMatrix<SC> TP_Crs;
  typedef Thyra::PreconditionerFactoryBase<SC> Base;

  typedef Thyra::MultiVectorBase<SC> MV;
  typedef Thyra::LinearOpBase<SC> OP;
  //  typedef Tpetra::Vector<SC,LO,GO,NO>  MV;
  //    typedef Tpetra::Operator<SC,LO,GO,NO>     OP;
  auto comm = crsMat->getRowMap()->getComm();

  // tell Stratimikos => Teko about MueLu
  RCP<Stratimikos::DefaultLinearSolverBuilder> linearSolverBuilder = Teuchos::rcp(new Stratimikos::DefaultLinearSolverBuilder);
  Stratimikos::enableMueLu<SC, LO, GO, NO>(*linearSolverBuilder);

  /////////////////////////////////////////////////////////
  // Build the Thyra operators
  /////////////////////////////////////////////////////////
  RCP<Tpetra::Operator<SC, LO, GO, NO> > A;
  RCP<Tpetra::Vector<SC, LO, GO, NO> > x = rcp(new Tpetra::Vector<SC, LO, GO, NO>(b->getMap()));
  ReadSplittingFromDisk<SC, LO, GO, NO>(partitionFile, crsMat, A);
  x->putScalar(Teuchos::ScalarTraits<SC>::zero());

  RCP<Thyra::MultiVectorBase<SC> > xt = Thyra::createVector(x);
  RCP<Thyra::MultiVectorBase<SC> > bt = Thyra::createVector(b);

  /////////////////////////////////////////////////////////
  // Build the preconditioner
  /////////////////////////////////////////////////////////

  // build an InverseLibrary
  RCP<Teko::InverseLibrary> invLib;
  RCP<Teko::InverseFactory> inverse;
  if (xmlFile == "") {
    invLib  = Teko::InverseLibrary::buildFromParameterList(*buildLibPL(), linearSolverBuilder);
    inverse = invLib->getInverseFactory("Gauss-Seidel");
  } else {
    Teuchos::ParameterList xmlList;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFile, Teuchos::Ptr<Teuchos::ParameterList>(&xmlList), *comm);

    // Add coordinates if we need to
    if (!coords.is_null()) {
      // FIXME:  Please do this more generally
      auto &sublist = xmlList.sublist("MueluScalar").sublist("user data");
      sublist.set("Coordinates", Xpetra::toXpetra(coords));
    }

    invLib  = Teko::InverseLibrary::buildFromParameterList(xmlList, linearSolverBuilder);
    inverse = invLib->getInverseFactory("MyTekoPreconditioner");
  }

  // build the inverse factory needed by the example preconditioner
  Teko::TpetraHelpers::InverseFactoryOperator ifo(inverse);
  ifo.buildInverseOperator(A);

  // build the preconditioner from the jacobian
  Teko::LinearOp At = Thyra::tpetraLinearOp<SC, LO, GO, NO>(Thyra::tpetraVectorSpace<SC, LO, GO, NO>(A->getRangeMap()), Thyra::tpetraVectorSpace<SC, LO, GO, NO>(A->getDomainMap()), A);

  Teko::LinearOp prec = Teko::buildInverse(*inverse, At);

  // Setup the Belos solver
  /////////////////////////////////////////////////////////

  Teuchos::ParameterList belosList;
  belosList.set("Num Blocks", 200);              // Maximum number of blocks in Krylov factorization
  belosList.set("Block Size", 1);                // Blocksize to be used by iterative solver
  belosList.set("Maximum Iterations", 200);      // Maximum number of iterations allowed
  belosList.set("Maximum Restarts", 1);          // Maximum number of restarts allowed
  belosList.set("Convergence Tolerance", 1e-5);  // Relative convergence tolerance requested
  belosList.set("Verbosity", 33);                // Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails );
  belosList.set("Output Frequency", 1);
  belosList.set("Output Style", 1);

  RCP<Belos::LinearProblem<double, MV, OP> > problem = rcp(new Belos::LinearProblem<double, MV, OP>(At, xt, bt));
  problem->setLeftPrec(prec);
  problem->setProblem();  // should check the return type!!!

  RCP<Belos::SolverManager<double, MV, OP> > solver = rcp(new Belos::BlockGmresSolMgr<double, MV, OP>(problem, rcpFromRef(belosList)));

  //
  // Perform solve
  //
  Belos::ReturnType ret = solver->solve();

  if (ret != Belos::Converged) {
    std::cout << std::endl
              << "ERROR:  Belos did not converge!" << std::endl;
    return -1;
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
  std::string xmlFile = "";
  clp.setOption("xml", &xmlFile, "read Tekko parameters from an xml file");
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
    crsMat = Tpetra::MatrixMarket::Reader<TP_Crs>::readSparseFile(matrixFile, Tpetra::getDefaultComm());
  } else {
    RCP<const TP_Map> rowMap = Tpetra::MatrixMarket::Reader<TP_Crs>::readMapFile(rowMapFile, Tpetra::getDefaultComm());
    RCP<const TP_Map> colMap;
    crsMat = Tpetra::MatrixMarket::Reader<TP_Crs>::readSparseFile(matrixFile, rowMap, colMap, rowMap, rowMap);

    if (rhsFile != "") {
      rhs = Tpetra::MatrixMarket::Reader<TP_Crs>::readVectorFile(rhsFile, rowMap->getComm(), rowMap);
    }

    if (coordsFile != "") {
      RCP<const TP_Map> coordsMap = rowMap;
      if (coordsMapFile != "")
        coordsMap = Tpetra::MatrixMarket::Reader<TP_Crs>::readMapFile(coordsMapFile, rowMap->getComm());
      coords = Tpetra::MatrixMarket::Reader<TP_Crs>::readDenseFile(coordsFile, coordsMap->getComm(), coordsMap);
    }
  }

  // Sanity checks
  if (rhs.is_null()) throw std::runtime_error("rhs is null");
  if (crsMat.is_null()) throw std::runtime_error("crsMat is null");

  auto comm = crsMat->getRowMap()->getComm();

  /////////////////////////////////////////////////////////
  // Build the Thyra operators
  /////////////////////////////////////////////////////////
  Teko::LinearOp A;
  Teko::MultiVector x, b;

  int rv;
  if (partitionFile == "")
    rv = solve_thyra(crsMat, xmlFile);
  else
    rv = solve_tpetra<SC, LO, GO, NO>(crsMat, rhs, coords, xmlFile, partitionFile);

  return rv;
}

RCP<Teuchos::ParameterList> buildLibPL() {
  RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList());

  {
    Teuchos::ParameterList &sub_jac = pl->sublist("Jacobi");
    sub_jac.set("Type", "Block Jacobi");
    sub_jac.set("Inverse Type", "Ifpack2");

    Teuchos::ParameterList &sub_gs = pl->sublist("Gauss-Seidel");
    sub_gs.set("Type", "Block Gauss-Seidel");
    sub_gs.set("Use Upper Triangle", true);
    sub_gs.set("Inverse Type", "Ifpack2");
  }
  return pl;
}
