// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include "Teuchos_Assert.hpp"
#include "Teuchos_LocalTestingHelpers.hpp"
#include <Teuchos_StandardCatchMacros.hpp>
#include "Teuchos_ParameterList.hpp"

// Xpetra
#include "Xpetra_Map.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_IO.hpp"
#include "Xpetra_MatrixSplitting.hpp"
#include "Xpetra_RegionAMG_def.hpp"
#include "Xpetra_RegionHandler_def.hpp"
#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraCrsMatrix.hpp"
#endif
#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraCrsMatrix.hpp"
#endif

// Epetra
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"

// =========== //
// main driver //
// =========== //

int main(int argc, char* argv[]) {
  typedef double scalar_type;
  typedef int local_ordinal_type;
  typedef int global_ordinal_type;
  typedef scalar_type Scalar;
  typedef local_ordinal_type LocalOrdinal;
  typedef global_ordinal_type GlobalOrdinal;
  typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType Node;

  typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> multivector_type;
  typedef Xpetra::MatrixSplitting<Scalar, LocalOrdinal, GlobalOrdinal, Node, Xpetra::UseTpetra, false> tpetra_splitting;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm CommEpetra(MPI_COMM_WORLD);
#else
  Epetra_SerialComm CommEpetra;
#endif

  // Here we create the linear problem
  //
  //   Matrix * LHS = RHS
  //
  // with Matrix arising from a 5-point formula discretization.

  TEUCHOS_TEST_FOR_EXCEPT_MSG(argc < 4, "\nInvalid name for input matrix and output file\n");

  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  // process command line arguments
  const char* xmlFileName     = argv[1];
  const char* matrixFileName  = argv[2];
  const char* mappingFileName = argv[3];

  Teuchos::ParameterList xmlParams;
  Teuchos::ParameterList mueluParams;
  Teuchos::updateParametersFromXmlFile(argv[1], Teuchos::inoutArg(mueluParams));

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  if (CommEpetra.MyPID() == 0)
    std::cout << "Number of processors: " << CommEpetra.NumProc() << std::endl;

  // SplittingDriver
  /*Xpetra::SplittingDriver<Scalar, LocalOrdinal, GlobalOrdinal, Node> driver("node.txt", comm);
        Teuchos::Array<GlobalOrdinal> elementlist = driver.GetGlobalRowMap();
        driver.printInactive();
        Xpetra::MatrixSplitting<Scalar,LocalOrdinal,GlobalOrdinal,Node,Xpetra::UseTpetra> xpetraWrapper( argv[1], argv[2], comm );
        std::string output_file="A_write.mtx";
        xpetraWrapper.writeGlobalMatrix();
        xpetraWrapper.writeRegionMatrices();*/

  // Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A;
  // A = Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Read(argv[2], Xpetra::UseTpetra, comm);

  // Create the RegionHandler to deal with mappings of nodes to regions etc.
  Teuchos::RCP<Xpetra::RegionHandler<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regionHandler = Teuchos::rcp(new Xpetra::RegionHandler<Scalar, LocalOrdinal, GlobalOrdinal, Node>(mappingFileName, comm));
  Teuchos::Array<GlobalOrdinal> elementlist                                                     = regionHandler->GetGlobalRowMap();
  std::size_t num_total_elements                                                                = regionHandler->GetNumGlobalElements();
  std::size_t num_total_regions                                                                 = regionHandler->GetNumTotalRegions();

  // Read and split the matrix
  Teuchos::RCP<tpetra_splitting> matrixSplitting = Teuchos::rcp(new tpetra_splitting(matrixFileName, regionHandler, comm));

  // Create region-wise AMG hierarchy
  int max_num_levels    = 4;
  int coarsening_factor = 3;
  Xpetra::RegionAMG<Scalar, LocalOrdinal, GlobalOrdinal, Node> preconditioner(matrixSplitting, regionHandler, comm, mueluParams, max_num_levels, coarsening_factor);

  //  // Setup vectors for test problem
  //  Teuchos::RCP<multivector_type> X = Xpetra::MultiVectorFactory< Scalar, LocalOrdinal, GlobalOrdinal, Node >::Build(preconditioner.getDomainMap(), 1) ;
  //  Teuchos::RCP<multivector_type> Y = Xpetra::MultiVectorFactory< Scalar, LocalOrdinal, GlobalOrdinal, Node >::Build(preconditioner.getRangeMap(), 1) ;
  //  X->randomize();
  //  Y->putScalar((scalar_type) 0.0);
  //
  //  // Apply the preconditioner
  //  preconditioner.apply(*X,*Y);
  //
  //  // Output result to screen
  //  Y->describe(*out, Teuchos::VERB_EXTREME);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return (EXIT_SUCCESS);
}
