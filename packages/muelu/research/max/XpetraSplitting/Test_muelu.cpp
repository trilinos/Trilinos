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

int main(int argc, char* argv[])
{


  typedef double                                      scalar_type;
  typedef int                                         local_ordinal_type;
  typedef int                                         global_ordinal_type;
  typedef scalar_type         												Scalar;
  typedef local_ordinal_type  												LocalOrdinal;
  typedef global_ordinal_type 												GlobalOrdinal;
	typedef KokkosClassic::DefaultNode::DefaultNodeType Node;

	typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> 											Matrix;
	typedef Xpetra::MatrixSplitting<Scalar, LocalOrdinal, GlobalOrdinal, Node>							MatrixSplitting;
	typedef Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Xpetra::EpetraNode> EpCrsMatrix;

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
  
	TEUCHOS_TEST_FOR_EXCEPT_MSG(argc<2, "\nInvalid name for input matrix\n");

	int numGlobalElements = 1;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
	if (CommEpetra.MyPID() == 0)
		std::cout<<"Number of processors: "<<CommEpetra.NumProc()<<std::endl;

	//Create Xpetra map
  Teuchos::RCP<const Xpetra::Map<int,GlobalOrdinal,Xpetra::EpetraNode> > xpetraMap;
	xpetraMap = Xpetra::MapFactory<int,GlobalOrdinal,Xpetra::EpetraNode>::Build(Xpetra::UseEpetra, numGlobalElements, 0, comm); 

	//Import matrix from an .mtx file into an Xpetra wrapper for an Epetra matrix
	Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > xpetraMatrix = Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Read(argv[1], Xpetra::UseEpetra, comm);
	//Export matrix from an Xpetra wrapper into an .mtx file
	Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("A_write.mtx", *xpetraMatrix);

	Teuchos::RCP<MatrixSplitting> xpetraMatrixSplitting;
  
	Teuchos::ParameterList xmlParams;
	Teuchos::RCP<MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Hierarchy = MueLu::CreateXpetraPreconditioner( (Teuchos::RCP<Matrix>)xpetraMatrixSplitting, xmlParams );

	//Teuchos::RCP<const Xpetra::Map<int,GlobalOrdinal,Xpetra::EpetraNode> > epmap = Xpetra::MapFactory<int,GlobalOrdinal,Xpetra::EpetraNode>::createUniformContigMap(Xpetra::UseEpetra, 100, comm);	

	const Epetra_Map epetraMap = Xpetra::toEpetra(xpetraMap);


  Teuchos::ParameterList GaleriList;
 
/*
  try
  {
    Map = CreateMap("Cartesian2D", Comm, GaleriList);
    Matrix = CreateCrsMatrix("Laplace2D", Map, GaleriList);
    Epetra_Vector ExactSolution(*Map); ExactSolution.Random();
    Epetra_Vector LHS(*Map); LHS.PutScalar(0.0);
    Epetra_Vector RHS(*Map);

    Matrix->Multiply(false, ExactSolution, RHS);

    Epetra_LinearProblem Problem(Matrix, &LHS, &RHS);

    // at this point any object that understand Epetra_LinearProblem can be
    // used, for example AztecOO, Amesos. IFPACK and ML can be used to define a
    // preconditioner for Matrix. Here we use a simple solver, based on
    // LAPACK, that is meant for simple testing only.
    
    Solve(Problem);

    // and we compute the norm of the true residual. 
    double ResidualNorm = ComputeNorm(Matrix, &LHS, &RHS);

    if (Comm.MyPID() == 0)
      cout << ResidualNorm << endl;

    delete Map;
    delete Matrix;
  }
  catch (Galeri::Exception& rhs)
  {
    if (Comm.MyPID() == 0)
      rhs.Print();
    exit(EXIT_FAILURE);
  }

*/

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}

