// 
//  Test whether Zoltan's NUM_GLOBAL_PARTITIONS parameter is respected
//  by Isorropia.
//

#include <Isorropia_ConfigDefs.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_EpetraRedistributor.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

#include <ispatest_epetra_utils.hpp>
#include <ispatest_lbeval_utils.hpp>


#ifdef HAVE_EPETRA
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#ifdef HAVE_EPETRAEXT
#include <EpetraExt_CrsMatrixIn.h>
#endif
#endif

#include <iostream>
#include <sstream>
#include <string>

#include <Teuchos_CommandLineProcessor.hpp>

int main(int argc, char** argv) 
{
  int rc=0, fail = 0;  
#ifdef HAVE_EPETRAEXT
  bool verbose = false;
  int numProcs = 1;
  int localProc = 0;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  const Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  const Epetra_SerialComm Comm;
#endif

#ifdef HAVE_ISORROPIA_ZOLTAN
  if (numProcs < 2){
    std::cerr << "This test requires at least 2 processes" << std::endl;
    exit(1);
  }
#else
  std::cerr << "This test requires Zoltan" << std::endl;
  exit(1);
#endif

  // if (getenv("DEBUGME")){
  //   std::cerr << localProc << " gdb test_simple.exe " << getpid() << std::endl;
  //   sleep(15);
  // }

  Teuchos::CommandLineProcessor clp(false,true);

  // --f=fileName provides a different matrix market file for input
  // --v will print out the partitioning (small files only)

  std::string *inputFile = new std::string("simple.mtx");

  clp.setOption( "f", inputFile, 
                "Name of input matrix market file");
  clp.setOption( "v", "q", &verbose, 
                "Display matrix before and after partitioning.");

  Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return =
    clp.parse(argc,argv);

  if( parse_return == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED){
    exit(0);
  }
  if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) {
    exit(0);
  }
  
  const char *fname = inputFile->c_str();

  // Read in the matrix market file and distribute its rows across the
  // processes.

  Epetra_CrsMatrix *matrixPtr;
  rc = EpetraExt::MatrixMarketFileToCrsMatrix(fname, Comm, matrixPtr);
  if (rc){
    if (localProc == 0){
      std::cerr << "error reading input file" << std::endl;
    }
    exit(1);
  }

  //   Test hypergraph partitioning, with num_global_partitions equal
  //     to a value less than the number of partitions.

  Teuchos::RCP<Epetra_CrsMatrix> matrix = Teuchos::rcp(matrixPtr);

  if (verbose){
    ispatest::show_matrix("Before load balancing", matrix->Graph(), Comm);
  }

  Teuchos::ParameterList params;

  // Set the Zoltan parameters for this problem

  Teuchos::ParameterList &sublist = params.sublist("Zoltan");

  int num_partitions = numProcs - 1;
  
  sublist.set("DEBUG_LEVEL", "1");
  sublist.set("LB_METHOD", "HYPERGRAPH");
  sublist.set("LB_APPROACH", "PARTITION");
  sublist.set("PHG_CUT_OBJECTIVE", "CONNECTIVITY");  // "cutl"

  std::ostringstream os;
  os << num_partitions;
  std::string s = os.str();

  sublist.set("NUM_GLOBAL_PARTITIONS", s.c_str());

#if 0
  int num_local = 2;
  os << num_local;
  s = os.str();
  sublist.set("NUM_LOCAL_PARTITIONS", s.c_str());
#endif

  // Perform hyperedge partitioning with Zoltan 

  Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner;

  Teuchos::RCP<const Epetra_RowMatrix> rm = matrix;
  partitioner = Isorropia::Epetra::create_partitioner(rm, params);

  // Create a Redistributor based on the partitioning

  Isorropia::Epetra::Redistributor rd(partitioner);

  // Redistribute the matrix

  Teuchos::RCP<Epetra_CrsMatrix> newMatrix = rd.redistribute(*matrix);

  if (verbose)
    ispatest::show_matrix("After load balancing", newMatrix.get()->Graph(), Comm);

#else
  std::cerr << "test_simple : currently can only test "
         << "with Epetra and EpetraExt enabled." << std::endl;
  fail =  1;
#endif


#if 0
  // Test that new matrix is a valid matrix

  fail = ispatest::test_matrix_vector_multiply(*newMatrix);

  if (fail && !localProc){
    std::cerr << "The rebalanced matrix is not a valid matrix" << std::endl;
  }
#endif

#ifdef HAVE_MPI
  MPI_Finalize();
#endif


  if (localProc == 0){
    if (fail)
      std::cout << "FAIL" << std::endl;
    else
      std::cout << "PASS" << std::endl;
  }

  return fail;
}
