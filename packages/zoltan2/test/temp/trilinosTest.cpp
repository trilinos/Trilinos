// @HEADER
// ***********************************************************************
//                Copyright message goes here. 
// ***********************************************************************
// @HEADER

#include <Teuchos_DefaultComm.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <MatrixMarket_Tpetra.hpp>

#include <string>
#include <iostream>
#include <exception>

using std::string;
using std::cerr;
using std::endl;
using Teuchos::RCP;

typedef Kokkos::DefaultNode::DefaultNodeType node_t;
typedef float scalar_t;
typedef int lno_t;
typedef int gno_t;

typedef Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> tcrsMatrix_t;


int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // Error: it's a pattern matrix.  (A graph, there are non non-zeros.)
  // packages/tpetra/inout/MatrixMarket_Tpetra.hpp line 1003
  //
  // string fname("commanche_dual.mtx"); 

  // Crash: 
  string fname("USAir97.mtx");

  //string fname("simple.mtx"); // This file is read without error

  RCP<tcrsMatrix_t> M;
  RCP<Kokkos::DefaultNode::DefaultNodeType> dnode
    = Kokkos::DefaultNode::getDefaultNode();

  try{
    M = Tpetra::MatrixMarket::Reader<tcrsMatrix_t>::readSparseFile(
      fname, comm, dnode, true, true);
  }
  catch (std::invalid_argument &e){
    cerr << e.what() << endl;
  }
  catch (...){
    cerr << "error" << endl;
  }
}
