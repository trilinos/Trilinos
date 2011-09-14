#include <mpi.h>
#include <iostream>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Vector.hpp>
#include <MatrixMarket_Tpetra.hpp>

using Teuchos::RCP;
using namespace std;

/////////////////////////////////////////////////////////////////////////////
// Program to demonstrate use of Zoltan2 to partition a TPetra matrix 
// (read from a MatrixMarket file).
// Usage:
//     a.out --inputFile=filename.mtx [--echo] [--verbose]
// Karen Devine, September 2011
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
// Eventually want to use Teuchos unit tests to vary LocalOrdinal and
// GlobalOrdinal.  For now, we set them at compile time.
typedef int LocalOrdinal;
typedef long GlobalOrdinal;
typedef double Scalar;
typedef Kokkos::DefaultNode::DefaultNodeType Node;
typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal> SparseMatrix;
typedef Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal> Vector;

/////////////////////////////////////////////////////////////////////////////
int main(int narg, char** arg)
{
  std::string inputFile;  // Matrix Market file to read
  bool verbose = false;       // Verbosity of output
  bool echo = false;          // Echo the read-in matrix back?

  // Read run-time options.
  Teuchos::CommandLineProcessor cmdp (false, true);
  cmdp.setOption("inputFile", &inputFile,
                 "Name of the Matrix Market sparse matrix file to read.");
  cmdp.setOption("echo", "noecho", &echo,
                 "Whether to echo the read-in matrix back to stdout on Rank 0 "
                 "in Matrix Market format.  Symmetric storage will have been "
                 "expanded, so the result will not be identical to the input "
                 "file, though the matrix represented will be the same.");
  cmdp.setOption("verbose", "quiet", &verbose,
                 "Print messages and results.");
  cmdp.parse(narg, arg);

  ////// Establish MPI session.
  Teuchos::GlobalMPISession mpiSession(&narg, &arg, &cout);
  RCP<const Teuchos::Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  int me = comm->getRank();

  ////// Read Matrix-Market matrix using Tpetra utilities.
  // Need a node for the MatrixMarket reader.
  Teuchos::ParameterList defaultParameters;
  RCP<Node> node = rcp(new Node(defaultParameters));

  RCP<SparseMatrix> origMatrix =
    Tpetra::MatrixMarket::Reader<SparseMatrix>::readSparseFile(
                                       inputFile, comm, node, 
                                       true, false, true);
  if (echo) {
    // Just a sanity check.
    Tpetra::MatrixMarket::Writer<SparseMatrix>::writeSparse(cout, 
                                       origMatrix, verbose);
  }

  ////// Create a vector to use with the matrix.
  RCP<Vector> origVector, origProd;
  origProd   = Tpetra::createVector<Scalar,LocalOrdinal,GlobalOrdinal>(
                                    origMatrix->getRangeMap());
  origVector = Tpetra::createVector<Scalar,LocalOrdinal,GlobalOrdinal>(
                                    origMatrix->getRangeMap());
  origVector->randomize();

  ////// Specify problem parameters

  ////// Create and solve partitioning problem

  ////// Redstribute matrix and vector into new matrix and vector.
  //Zoltan2::PartitioningProblem<LocalOrdinal, GlobalOrdinal> 
//    partitioningProblem(A, params, configuration);

 // partitioningProblem.solve();


  ////// Verify that redistribution is "correct"; perform matvec with 
  ////// original and redistributed matrices/vectors and compare norms.
  int testReturn = 0;

  origMatrix->apply(*origVector, *origProd);

  Scalar origNorm, newNorm;
  origNorm = origProd->norm2();
  if (me == 0)
    cout << "Norm of Original matvec prod:  " << origNorm << endl;

  return testReturn;
}
