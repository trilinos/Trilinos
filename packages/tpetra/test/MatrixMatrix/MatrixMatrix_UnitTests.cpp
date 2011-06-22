#include "Kokkos_SerialNode.hpp"
#include "Tpetra_MatrixMatrix.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrixMultiplyOp.hpp"
#include "Tpetra_Import.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"
#include <cmath>

namespace {
  static const double defaultEpsilon = 1e-10;
  bool verbose = false;
  std::string matnamesFile;

  using Teuchos::null;
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::tuple;
  using Tpetra::global_size_t;
  using Teuchos::Comm;
  using Tpetra::CrsMatrix;
  using Tpetra::Map;
  using Teuchos::Array;
  using Tpetra::Vector;
  using Tpetra::CrsMatrixMultiplyOp;
  using Tpetra::DefaultPlatform;
  using Tpetra::MatrixMarket::Reader;
  using Teuchos::ArrayView;
  using Teuchos::FancyOStream;
  using Teuchos::ParameterList;
  using Kokkos::SerialNode;
  using Tpetra::RowMatrixTransposer;

  typedef CrsMatrix<double, int, int, SerialNode> Matrix_t;


TEUCHOS_STATIC_SETUP(){
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.setOption("matnames-file", &matnamesFile, 
    "A file containing a list of matricies we'll import", true);
  clp.setOption("v", "not-verbose", &verbose, 
    "Whether or not to use verbose output");
}



typedef struct add_test_results_struct{
  double correctNorm;
  double norm1;
  double epsilon;
} add_test_results;

typedef struct mult_test_results_struct{
  double epsilon;
  double cNorm;
  double compNorm;
} mult_test_results;

template<class CrsMatrix_t> double getNorm(RCP<CrsMatrix_t> matrix){
  double mySum = 0;
  Array<int> inds(matrix->getNodeMaxNumRowEntries());
  Array<double> vals(matrix->getNodeMaxNumRowEntries());
  for(int i =0; ((size_t)i)<matrix->getNodeNumRows(); ++i){
    size_t numRowEnts = matrix->getNumEntriesInLocalRow(i);
    ArrayView<const int> indsView = inds();
    ArrayView<const double> valsView = vals();
    matrix->getLocalRowView(i, indsView, valsView);
    for(size_t j=0; ((size_t)j)<numRowEnts; ++j){
      mySum += valsView[j]*valsView[j];
    }
  }
  double totalSum = 0;
  Teuchos::reduceAll(*(matrix->getComm()), Teuchos::REDUCE_SUM, 1, &mySum, &totalSum);
  return sqrt(totalSum);

}

template<class Ordinal>
add_test_results regular_add_test(
    RCP<Matrix_t > A,
    RCP<Matrix_t > B,
    bool AT,
    bool BT,
    RCP<Matrix_t > C,
    RCP<const Comm<Ordinal> > comm)
{
  add_test_results toReturn;
  toReturn.correctNorm = getNorm(C);

  RCP<const Map<int,int, Kokkos::SerialNode> > rowmap = AT ? A->getDomainMap() : A->getRowMap();
  RCP<Matrix_t> computedC = rcp( new Matrix_t(rowmap, 1));

  Tpetra::MatrixMatrix::Add(*A, AT, 1.0, *B, BT, 1.0, computedC);

  computedC->fillComplete(A->getDomainMap(), A->getRangeMap());
  toReturn.norm1 = getNorm(computedC);
  toReturn.epsilon = fabs(toReturn.correctNorm - toReturn.norm1);
  

  return toReturn;

}

template<class Ordinal>
add_test_results add_into_test(
    RCP<Matrix_t > A,
    RCP<Matrix_t > B,
    bool AT,
    RCP<Matrix_t > C,
    RCP<const Comm<Ordinal> > comm)
{
  add_test_results toReturn;
  toReturn.correctNorm = getNorm(C);

  RCP<const Map<int,int, Kokkos::SerialNode> > rowmap = AT ? A->getDomainMap() : A->getRowMap();
  RCP<Matrix_t> computedC = rcp( new Matrix_t(rowmap, 1));
  Tpetra::MatrixMatrix::Add(*A, AT, 1.0, *B, 1.0);
  B->fillComplete();
  toReturn.norm1 = getNorm(B);
  toReturn.epsilon = fabs(toReturn.correctNorm - toReturn.norm1);

  return toReturn;
}

template<class Ordinal>
mult_test_results multiply_test(
  const std::string& name,
  RCP<Matrix_t> A,
  RCP<Matrix_t> B,
  bool AT,
  bool BT,
  RCP<Matrix_t> C,
  RCP<const Comm<Ordinal> > comm,
  FancyOStream& out)
{

  typedef Map<int, int, SerialNode> Map_t;
  RCP<const Map_t> map = A->getRowMap();

  RCP<Matrix_t> computedC = rcp( new Matrix_t(map, 1));

  Tpetra::MatrixMatrix::Multiply(*A, AT, *B, BT, *computedC, false);
  computedC->globalAssemble();
  /*Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile(
    name+"_calculated.mtx",computedC);
  Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile(
    name+"_real.mtx",C);*/
   
  double cNorm = getNorm(C);
  Tpetra::MatrixMatrix::Add(*C, false, -1.0, *computedC, 1.0);
  if(!AT && !BT){
    computedC->fillComplete(B->getDomainMap(), A->getRangeMap());
  }
  else{
    // We actually should test for other cases like AT*B and BT*A etc,
    // and specify the appropriate domain and range maps but I don't have
    // time for that now.
    // KLN 13/06/2011
    computedC->fillComplete();
  }
  double compNorm = getNorm(computedC);
  mult_test_results results;
  results.epsilon = compNorm/cNorm;
  results.cNorm = cNorm;
  results.compNorm = compNorm;
  return results;
}


TEUCHOS_UNIT_TEST(Tpetra_MatMat, test_find_rows){
  RCP<const Comm<int> > comm = DefaultPlatform::getDefaultPlatform().getComm();
  int numprocs = comm->getSize();
  int localproc = comm->getRank();
  int numlocalrows = 2;
  global_size_t numglobalrows = numprocs*numlocalrows;
  RCP<const Map<int> > rowmap = Tpetra::createUniformContigMap<int,int>(numglobalrows,comm);
  CrsMatrix<double, int> matrix(rowmap, numglobalrows);

  Array<int> cols(numglobalrows);
  Array<double> vals(numglobalrows);

  for(size_t j=0; j<numglobalrows; ++j) {
    cols[j] = j;
    vals[j] = 1.0;
  }

  RCP<const Map<int> > colmap = Tpetra::createNonContigMap<int,int>(cols(), comm);

  for(int i=0; i<numlocalrows; ++i) {
    Array<int> row(1,localproc*numlocalrows+i);
    matrix.insertGlobalValues(
      row[0], row.view(0,1),  vals.view(i,1) );
  }

  matrix.fillComplete();

  typedef Kokkos::DefaultNode::DefaultNodeType DNode;
  typedef Kokkos::DefaultKernels<double, int, DNode>::SparseOps SpMatOps;

  RCP<const Map<int, int, DNode> > map_rows = 
    Tpetra::MMdetails::find_rows_containing_cols<double, int, int, DNode, SpMatOps>(matrix, colmap);

  TEST_EQUALITY(map_rows->getNodeNumElements(), numglobalrows);

}

TEUCHOS_UNIT_TEST(Tpetra_MatMat, operations_test){
  RCP<const Comm<int> > comm = DefaultPlatform::getDefaultPlatform().getComm();
  ParameterList defaultParameters;
  RCP<SerialNode> node = rcp(new SerialNode(defaultParameters));
  Teuchos::RCP<Teuchos::ParameterList> matrixSystems = 
    Teuchos::getParametersFromXmlFile(matnamesFile);
  for(
    Teuchos::ParameterList::ConstIterator it = matrixSystems->begin();
    it != matrixSystems->end();
    ++it)
  {
	  TEST_FOR_EXCEPTION(!it->second.isList(), std::runtime_error,
      "All top level items in the matrix "
	    "file names list must be ParameterLists! In otherwords, you always "
      "need to have matricies "
	    "encapsulated within a matrixsystem" << std::endl <<
      "Bad tag's name: " << it->first << 
      "Type name: " << it->second.getAny().typeName() << 
      std::endl << std::endl);
    Teuchos::ParameterList currentSystem = matrixSystems->sublist(it->first); 
    std::string name = currentSystem.name();
    std::string A_file = currentSystem.get<std::string>("A");
    std::string B_file = currentSystem.get<std::string>("B");
    std::string C_file = currentSystem.get<std::string>("C");
    bool AT = currentSystem.get<bool>("TransA");
    bool BT = currentSystem.get<bool>("TransB");
    double epsilon = currentSystem.get<double>("epsilon", defaultEpsilon);
    std::string op = currentSystem.get<std::string>("op");

    RCP<Matrix_t > A = Reader<Matrix_t >::readSparseFile(A_file, comm, node);
    RCP<Matrix_t > B = Reader<Matrix_t >::readSparseFile(B_file, comm, node);
    RCP<Matrix_t > C = Reader<Matrix_t >::readSparseFile(C_file, comm, node);
 

    TEST_FOR_EXCEPTION(op != "multiply" && op != "add", std::runtime_error,
      "Unrecognized Matrix Operation: " << op << "!" << std::endl);
  
    if(op == "multiply"){
      if(verbose){
        out << "Running multiply test for " << currentSystem.name() << std::endl;
      }
      mult_test_results results = multiply_test(name, A,B,AT,BT,C,comm, out);
      if(verbose){
        out << "Results:" <<std::endl;
        out << "\tEpsilon: " << results.epsilon << std::endl;
        out << "\tcNorm: " << results.cNorm << std::endl;
        out << "\tcompNorm: " << results.compNorm << std::endl;
      }
      TEST_COMPARE(results.epsilon, <, epsilon)
    }
    else if(op == "add"){
      if(verbose){
        out << "Running add test for " << currentSystem.name() << std::endl;
      }
      add_test_results results = regular_add_test(A,B,AT,BT,C,comm);
      TEST_COMPARE(results.epsilon, <, epsilon)
      out << "Regular Add Test Results: " << std::endl;
      out << "\tCorrect Norm: " << results.correctNorm << std::endl;
      out << "\tNorm 1: " << results.norm1 << std::endl;
      out << "\tEpsilon: " << results.epsilon << std::endl;
      B = Reader<Matrix_t >::readSparseFile(B_file, comm, node, false);
      if(!BT){
        results = add_into_test(A,B,AT,C,comm);
        TEST_COMPARE(results.epsilon, <, epsilon)
        out << "Add Into Test Results: " << std::endl;
        out << "\tCorrect Norm: " << results.correctNorm << std::endl;
        out << "\tNorm 1: " << results.norm1 << std::endl;
        out << "\tEpsilon: " << results.epsilon << std::endl;
      }
    }
  }   
}

TEUCHOS_UNIT_TEST(Tpetra_MatMat, sparse_dot_test){
  Array<double> uVal = tuple<double>(4,8,1,6);
  Array<double> vVal = tuple<double>(3,2,4,50);
  Array<int> uInd = tuple<int>(0,5,7,9);
  Array<int> vInd = tuple<int>(0,9,10,11);
  TEST_EQUALITY_CONST(
    Tpetra::MMdetails::sparsedot(uVal(), uInd(), vVal(), vInd()), 24);
}

/*
 * This test was created at the request of Chris Seifert to verify
 * that some inexplicable behaviour in MueLue was not due to a faulty
 * assumption in the Matrix Matrix Multiply Kernel.
 * KLN 15/06/2011
 */
TEUCHOS_UNIT_TEST(Tpetra_MatMat, range_row_test){
  RCP<const Comm<int> > comm = DefaultPlatform::getDefaultPlatform().getComm();
  ParameterList defaultParameters;
  RCP<SerialNode> node = rcp(new SerialNode(defaultParameters));
  int numProcs = comm->getSize();
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //THIS NUMBER MUST BE EVEN SO THAT WHEN I CALCULATE THE NUMBER
  //OF ROWS IN THE DOMAIN MAP I DON'T ENCOUNTER ANY 
  //WEIRD RESULTS DUE TO INTEGER DIVISION
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  int numRowsPerProc = 4;
  int rank = comm->getRank();
  global_size_t globalNumRows = numRowsPerProc*numProcs;



//Create identity matrix
  RCP<const Map<int,int,SerialNode> > identityRowMap = 
    Tpetra::createUniformContigMapWithNode<int,int,SerialNode>(
      globalNumRows, comm, node);
  RCP<CrsMatrix<double,int,int,SerialNode> > identityMatrix = 
    Tpetra::createCrsMatrix<double,int,int,SerialNode>(identityRowMap, 1);
  for(
    ArrayView<const int>::iterator it = 
      identityRowMap->getNodeElementList().begin();
    it != identityRowMap->getNodeElementList().end();
    ++it)
  {
    Array<int> col(1,*it);
    Array<double> val(1,1.0);
    identityMatrix->insertGlobalValues(*it, col(), val());
  }
  identityMatrix->fillComplete();


//Create "B"
  Array<int> myRows = tuple<int>(
    rank*numRowsPerProc,
    rank*numRowsPerProc+1,
    rank*numRowsPerProc+2,
    rank*numRowsPerProc+3);
  Array<int> rangeElements;
  if(rank == 0){
    rangeElements = tuple<int>(
      (numProcs-1)*numRowsPerProc+1,
      (numProcs-1)*numRowsPerProc+2,
      (numProcs-1)*numRowsPerProc,
      (numProcs-1)*numRowsPerProc+3);
  }
  else{
    rangeElements = tuple<int>(
      (rank-1)*numRowsPerProc+1,
      (rank-1)*numRowsPerProc+2,
      (rank-1)*numRowsPerProc,
      (rank-1)*numRowsPerProc+3);
  }

  RCP<const Map<int,int,SerialNode> > bRowMap =
    Tpetra::createNonContigMapWithNode<int,int,SerialNode>(myRows, comm, node);
  RCP<const Map<int,int,SerialNode> > bRangeMap =
    Tpetra::createNonContigMapWithNode<int,int,SerialNode>(rangeElements, comm, node);
  //We divide by 2 to make the matrix tall and "skinny"
  RCP<const Map<int,int,SerialNode> > bDomainMap = 
    Tpetra::createUniformContigMapWithNode<int,int,SerialNode>(
      globalNumRows/2, comm, node);
  RCP<CrsMatrix<double,int,int,SerialNode> > bMatrix = 
    Tpetra::createCrsMatrix<double,int,int,SerialNode>(bRowMap, 1);
  for(
    ArrayView<const int>::iterator it = 
      bRowMap->getNodeElementList().begin();
    it != bRowMap->getNodeElementList().end();
    ++it)
  {
    Array<int> col(1,(*it)/2);
    Array<double> val(1,3.0);
    bMatrix->insertGlobalValues(*it, col(), val());
  }
  bMatrix->fillComplete(bDomainMap, bRangeMap); 

  out << "Regular I*P" << std::endl;
  mult_test_results results = multiply_test(
    "Different Range and Row Maps", 
    identityMatrix,
    bMatrix,
    false,
    false,
    bMatrix,
    comm,
    out);
  if(verbose){
    out << "Results:" <<std::endl;
    out << "\tEpsilon: " << results.epsilon << std::endl;
    out << "\tcNorm: " << results.cNorm << std::endl;
    out << "\tcompNorm: " << results.compNorm << std::endl;
  }
  TEST_COMPARE(results.epsilon, <, defaultEpsilon)



  RCP<const Map<int,int,SerialNode> > identity2RowMap = 
    Tpetra::createUniformContigMapWithNode<int, int, SerialNode>(
      globalNumRows/2, comm,node);
  RCP<CrsMatrix<double,int,int,SerialNode> > identity2 = 
    Tpetra::createCrsMatrix<double,int,int,SerialNode>(identity2RowMap,1);
  for(
    ArrayView<const int>::iterator it = 
      identity2RowMap->getNodeElementList().begin();
    it != identity2RowMap->getNodeElementList().end();
    ++it)
  {
    Array<int> col(1,*it);
    Array<double> val(1,1.0);
    identity2->insertGlobalValues(*it, col(), val());
  }
  identity2->fillComplete();


  RCP<const Map<int,int,SerialNode> > bTransRowMap =
    Tpetra::createUniformContigMapWithNode<int,int,SerialNode>(globalNumRows/2,comm,node);

  RCP<CrsMatrix<double,int,int,SerialNode> > bTrans = 
    Tpetra::createCrsMatrix<double,int,int,SerialNode>(bTransRowMap, 1);
  Array<int> bTransRangeElements;
  if(rank == 0){
    bTransRangeElements = tuple<int>(
      (numProcs-1)*(numRowsPerProc/2)+1,
      (numProcs-1)*(numRowsPerProc/2));
  }
  else{
    bTransRangeElements = tuple<int>(
      (rank-1)*(numRowsPerProc/2)+1,
      (rank-1)*(numRowsPerProc/2));
  }
  out << bTransRangeElements << std::endl;
  RCP<const Map<int,int,SerialNode> > bTransRangeMap = 
    Tpetra::createNonContigMapWithNode<int,int,SerialNode>(bTransRangeElements, comm, node);
  RCP<const Map<int,int,SerialNode> > bTransDomainMap =
    Tpetra::createUniformContigMapWithNode<int,int,SerialNode>(globalNumRows,comm,node);
  Tpetra::MatrixMatrix::Multiply(*identity2,false,*bMatrix, true, *bTrans, false);
  bTrans->fillComplete(bTransDomainMap, bTransRangeMap);

  RCP<CrsMatrix<double,int,int,SerialNode> > bTransTest = 
    Tpetra::createCrsMatrix<double,int,int,SerialNode>(bTransRowMap, 1);

  for(
    ArrayView<const int>::iterator it = 
      bRowMap->getNodeElementList().begin();
    it != bRowMap->getNodeElementList().end();
    ++it)
  {
    Array<int> col(1,*it);
    Array<double> val(1,3.0);
    bTransTest->insertGlobalValues((*it)/2, col(), val());
  }
  bTransTest->fillComplete(bTransDomainMap, bTransRangeMap);


  out << "Regular I*P^T" << std::endl;
  RCP<CrsMatrix<double,int,int,SerialNode> > bTransDiff = 
    Tpetra::createCrsMatrix<double,int,int,SerialNode>(bTransRowMap, 1);
  Tpetra::MatrixMatrix::Add<double,int,int,SerialNode>(*bTransTest, false, -1.0, *bTrans, false, 1.0,bTransDiff);
  bTransDiff->fillComplete(bTransDomainMap, bDomainMap);
  double diffNorm = getNorm(bTransDiff);
  double realNorm = getNorm(bTransTest);
  double calcEpsilon = diffNorm/realNorm;
    
  out << "B" << std::endl;


  if(verbose){
    out << "Results:" <<std::endl;
    out << "\tEpsilon: " << calcEpsilon<< std::endl;
    out << "\treal norm: " << realNorm<< std::endl;
    out << "\tcompNorm: " << diffNorm<< std::endl;
  }
  TEST_COMPARE(calcEpsilon, <, defaultEpsilon)

}


} //namespace Tpetra

