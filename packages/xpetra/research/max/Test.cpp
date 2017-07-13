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
  
	TEUCHOS_TEST_FOR_EXCEPT_MSG(argc<3, "\nInvalid name for input matrix and output file\n");

	Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
	if (CommEpetra.MyPID() == 0)
		std::cout<<"Number of processors: "<<CommEpetra.NumProc()<<std::endl;

	//SplittingDriver
	Xpetra::SplittingDriver<Scalar, LocalOrdinal, GlobalOrdinal, Node> driver("node.txt", comm);
	Teuchos::Array<GlobalOrdinal> elementlist = driver.GetGlobalRowMap();
	Teuchos::Array<GlobalOrdinal> elementlist_region1 = driver.GetLocalRowMap(1);
	/*driver.printView();
	driver.printNodesToRegion();*/
	driver.printInactive();
	std::cout<<"PID: "<<comm->getRank()<<" beginning map: "<<elementlist[0]<<std::endl;
	std::cout<<"PID: "<<comm->getRank()<<" end of map: "<<*(elementlist.end()-1)<<std::endl;
	if( elementlist_region1.size()>0 )
		std::cout<<"PID: "<<comm->getRank()<<" num local nodes in region 1: "<<elementlist_region1.size()<<" min_index= "<<*(std::min_element(elementlist_region1.begin(), elementlist_region1.end()))<<" max_index= "<<*(std::max_element(elementlist_region1.begin(), elementlist_region1.end()))<<std::endl;
	else
		std::cout<<"PID: "<<comm->getRank()<<" num local nodes in region 1: "<<elementlist_region1.size()<<std::endl;

	Xpetra::MatrixSplitting<Scalar,LocalOrdinal,GlobalOrdinal,Node,Xpetra::UseEpetra> xpetraWrapper( argv[1], argv[2], comm );
	std::string output_file="A_write.mtx";
	xpetraWrapper.write(output_file.c_str());

	#ifdef HAVE_MPI
 	  MPI_Finalize();
	#endif

	return(EXIT_SUCCESS);
}

