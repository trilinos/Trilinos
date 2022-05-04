#include <iostream>
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosGmresPolySolMgr.hpp"
// I/O for Harwell-Boeing files
#define HIDE_TPETRA_INOUT_IMPLEMENTATIONS
#include <Tpetra_MatrixIO.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_MultiVector_decl.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Kokkos_View.hpp>

using namespace Teuchos;
using Tpetra::Operator;
using Tpetra::MultiVector;
using Tpetra::CrsMatrix;
using Teuchos::tuple;

int main(int argc, char *argv[]) {
  typedef double                           ST;
  typedef Tpetra::MultiVector<ST>          MV;
  typedef Belos::MultiVecTraits<ST,MV>    MVT;
  
  Tpetra::ScopeGuard tpetraScope(&argc,&argv);
  // Adapting from Hello World to print out matrix market file name
  std::string matrix_file("simple.mtx");
  std::string vector_file("");
  // std::cout << "hello world" << std::endl

  // Part C: Teuchos Command Line Option to print out mtx file name
  Teuchos::CommandLineProcessor cmdp(false,true);
  //cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("matrix_name",&matrix_file,"Filename for test matrix.  Acceptable file extensions: *.hb,*.mtx,*.triU,*.triS");
  cmdp.setOption("vector_name",&vector_file,"Filename for test vector.  Acceptable file extensions: *.hb,*.mtx,*.triU,*.triS");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
     return EXIT_FAILURE;
  }
  std::cout << "The matrix file loaded is " << matrix_file << std::endl;
  std::cout << "The vector file loaded is " << vector_file << std::endl;
   
   // Part A: Reads in a Tpetra Multivector from a Matrix Market File
   int numrhs = 1;                           // number of right-hand sides to solve for
   RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
   //RCP<MultiVector<ST> > V; 
   //V = Tpetra::MatrixMarket::Reader<MultiVector<ST> >::readSparseFile(filename,comm); 
   //RCP<const Tpetra::Map<> > map = V->getData();
   RCP<CrsMatrix<ST> > A;
   A = Tpetra::MatrixMarket::Reader<CrsMatrix<ST> >::readSparseFile(matrix_file,comm);
   //Tpetra::Utils::readHBMatrix(filename,comm,A); //This line causes core to dump
   RCP<const Tpetra::Map<> > map = A->getDomainMap();
   RCP<MV> B, X;
   X = rcp( new MV(map,numrhs) );
   B = rcp( new MV(map,numrhs) );
   MVT::MvInit( *X, 0.0 );
   MVT::MvInit( *B, 1.0 );
   
   //RCP<MV> V;
   //V = Tpetra::MatrixMarket::Reader<MV>::readDenseFile(vector_file,comm,map);//TODO MATLAB vector files need to be dense 
   // RCP<const Tpetra::Map<> > map = V->getData();
   // Part B: Prints the multivector to scree
   MVT::MvPrint(*B,std::cout);
   
   return 0;
}
