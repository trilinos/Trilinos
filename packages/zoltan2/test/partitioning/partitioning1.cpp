#include <mpi.h>
#include <iostream>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Vector.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Zoltan2_PartitioningProblem.hpp>

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
// Eventually want to use Teuchos unit tests to vary z2TestLO and
// GO.  For now, we set them at compile time.
typedef int z2TestLO;
typedef long z2TestGO;
typedef double Scalar;
typedef Kokkos::DefaultNode::DefaultNodeType Node;
typedef Tpetra::CrsMatrix<Scalar, z2TestLO, z2TestGO> SparseMatrix;
typedef Tpetra::Vector<Scalar, z2TestLO, z2TestGO> Vector;

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
  origProd   = Tpetra::createVector<Scalar,z2TestLO,z2TestGO>(
                                    origMatrix->getRangeMap());
  origVector = Tpetra::createVector<Scalar,z2TestLO,z2TestGO>(
                                    origMatrix->getRangeMap());
  origVector->randomize();

  ////// Specify problem parameters
  Teuchos::ParameterList params;
  params.set("APPROACH", "PARTITION");
  params.set("METHOD", "GRAPH");
  params.set("GRAPH_PACKAGE", "PTSCOTCH");

  ////// Create and solve partitioning problem
  Zoltan2::PartitioningProblem<Scalar,z2TestLO,z2TestGO,z2TestLO,z2TestGO,Node>
                               problem(*origMatrix, params);
  problem.solve();
  problem.redistribute();

  ////// Redistribute matrix and vector into new matrix and vector.

  ////// Verify that redistribution is "correct"; perform matvec with 
  ////// original and redistributed matrices/vectors and compare norms.
  int testReturn = 0;

  origMatrix->apply(*origVector, *origProd);

  Scalar origNorm;
  origNorm = origProd->norm2();
  if (me == 0)
    cout << "Norm of Original matvec prod:  " << origNorm << endl;

  return testReturn;
}
