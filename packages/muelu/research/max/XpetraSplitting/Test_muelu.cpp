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

// MueLu
#include <MueLu.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <MueLu_Utilities.hpp>

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

  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix;
  typedef Xpetra::MatrixSplitting<Scalar, LocalOrdinal, GlobalOrdinal, Node> MatrixSplitting;
  typedef Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Xpetra::EpetraNode> EpCrsMatrix;

  TEUCHOS_TEST_FOR_EXCEPT_MSG(argc < 2, "\nInvalid name for input matrix\n");

  int numGlobalElements = 1;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // Create Xpetra map
  Teuchos::RCP<const Xpetra::Map<int, GlobalOrdinal, Xpetra::EpetraNode> > xpetraMap;
  xpetraMap = Xpetra::MapFactory<int, GlobalOrdinal, Xpetra::EpetraNode>::Build(Xpetra::UseEpetra, numGlobalElements, 0, comm);

  // Import matrix from an .mtx file into an Xpetra wrapper for an Epetra matrix
  Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > xpetraMatrix = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read(argv[1], Xpetra::UseEpetra, comm);
  // Export matrix from an Xpetra wrapper into an .mtx file
  Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write("A_write.mtx", *xpetraMatrix);

  Teuchos::RCP<MatrixSplitting> xpetraMatrixSplitting;

  Teuchos::ParameterList xmlParams;
  Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Hierarchy = MueLu::CreateXpetraPreconditioner((Teuchos::RCP<Matrix>)xpetraMatrixSplitting, xmlParams);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return (EXIT_SUCCESS);
}
