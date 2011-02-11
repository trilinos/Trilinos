#include "Tpetra_MatrixMatrix.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"

namespace Tpetra{

static const double defaultEpsilon = 1e-10;
bool verbose = false;
std::string matnamesFile;
bool write_result_hb = false;

TEUCHOS_STATIC_SETUP(){
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.setOption("matnames-file", &matnamesFile, 
    "A file containing a list of matricies we'll import", true);
  clp.setOption("writeresults", "no-write-results", &write_result_hb, 
    "Whether or not to write the resutling matricies to hb files");
  clp.setOption("v", "not-verbose", &verbose, 
    "Whether or not to use verbose output");
}

template<class Ordinal>
int add_test(
    RCP<CrsMatrix<double,int> > A,
    RCP<CrsMatrix<double,int> > B,
    RCP<CrsMatrix<double,int> > C,
    bool AT,
    bool BT,
    double epsilon,
    RCP<const Comm<Ordinal> > comm,
    bool verbose)
{
  typedef Kokkos::DefaultNode::DefaultNodeType DNode;

  //int localProc = comm->getRank();

  RCP<CrsMatrix<double,int> > computedC = null;
  RCP<const Map<int> > rowmap = AT ? A->getDomainMap() : A->getRowMap();

  computedC = rcp( new CrsMatrix<double,int>(rowmap, 1));

  MatrixMatrix::Add(*A, false, 1.0, *B, false, 1.0, computedC);

  computedC->fillComplete(C->getDomainMap(), C->getRangeMap());

  MatrixMatrix::Add(*C, false, -1.0, *computedC, 1.0);

  double calculated_euc_norm = computedC->getEuclideanNorm();
  double c_euc_norm = C->getEuclideanNorm();
  double resultVal1 = calculated_euc_norm/c_euc_norm;

  MatrixMatrix::Add(*A, false, 1.0, *B, 1.0);

  MatrixMatrix::Add(*C, false, -1.0, *B, 1.0);

  calculated_euc_norm = B->getEuclideanNorm();
  double resultVal2 = calculated_euc_norm/c_euc_norm;



  if (resultVal1 < epsilon && resultVal2 < epsilon) {
/*    if (localProc == 0 && verbose) {
      std::cout << "Test Passed" << std::endl;
      std::cout << "||(A+B)-C||/||C|| " << resultVal1 <<std::endl ;
      std::cout << "||(A+B)-C||/||C|| " << resultVal2 <<std::endl << std::endl;
    }*/
    return 0;
  }
  else {
    /*if (localProc == 0) {
      std::cout << "Test Failed: " << std::endl;
      std::cout << "||(A+B)-C||/||C|| " << resultVal1 <<std::endl ;
      std::cout << "||(A+B)-C||/||C|| " << resultVal2 <<std::endl << std::endl;
    }*/
    return -1;
  }

}

template<class Ordinal>
int multiply_test(
  RCP<CrsMatrix<double, int> > A,
  RCP<CrsMatrix<double, int> > B,
  RCP<CrsMatrix<double, int> > C_check,
  bool AT,
  bool BT,
  double epsilon,
  RCP<const Comm<Ordinal> > comm,
  bool verbose)
{

  //int localProc = comm->getRank();
  RCP<CrsMatrix<double,int> > computedC = null;
  //RCP<const Map<int> > rowmap = AT ? A->getDomainMap() : A->getRowMap();

  computedC = rcp( new CrsMatrix<double,int>(C_check->getRowMap(), 1));

  RCP<const CrsMatrix<double,int> > constA = A;
  RCP<const CrsMatrix<double,int> > constB = B;
  typedef Kokkos::DefaultNode::DefaultNodeType DNode;


  MatrixMatrix::Multiply(*constA, AT, *constB, BT, *computedC);

//  std::cout << "A: " << *A << std::endl << "B: "<<*B<<std::endl<<"C: "<<*C<<std::endl;
  //if (result_mtx_to_file) {
   // EpetraExt::RowMatrixToMatrixMarketFile("result.mtx", *C);
  //}

  

  MatrixMatrix::Add(*C_check, false, -1.0, *computedC, 1.0);

  double c_check_euc_norm = C_check->getEuclideanNorm();
  double c_euc_norm = computedC->getEuclideanNorm();

  double diff_result = c_euc_norm/c_check_euc_norm;

  int return_code =0;
  if (diff_result < epsilon) {
/*    if (localProc == 0 && verbose) {
      std::cout << "Test Passed" << std::endl;
      std::cout << "||A*B-C||/||C|| " << diff_result << std::endl << std::endl;
    }*/
  }
  else {
    return_code = -1;
    /*if (localProc == 0) {
      std::cout << "Test Failed: ||A*B-C||/||C|| " << diff_result << std::endl << std::endl;
    }*/
  }

  return(return_code);

}


template<class Ordinal>
int run_test(RCP<const Comm<Ordinal> > comm,
             Teuchos::ParameterList matrixSystem,
             bool result_mtx_to_file,
             bool verbose)
{
  std::string A_file = matrixSystem.get<std::string>("A");
  std::string B_file = matrixSystem.get<std::string>("B");
  std::string C_file = matrixSystem.get<std::string>("C");
  bool AT = matrixSystem.get<bool>("TransA");
  bool BT = matrixSystem.get<bool>("TransB");
  double epsilon = matrixSystem.get<double>("epsilon", defaultEpsilon);
  std::string op = matrixSystem.get<std::string>("op");


  int localProc = comm->getRank();


  RCP<CrsMatrix<double,int> > A = null;
  RCP<CrsMatrix<double,int> > B = null;
  RCP<CrsMatrix<double,int> > C_check = null;

  Utils::readHBMatrix(A_file, comm, Kokkos::DefaultNode::getDefaultNode(), A);
  Utils::readHBMatrix(B_file, comm, Kokkos::DefaultNode::getDefaultNode(), B);
  Utils::readHBMatrix(C_file, comm, Kokkos::DefaultNode::getDefaultNode(), C_check);

  if(op == "multiply"){
    if(localProc == 0 && verbose){
      std::cout << "Running multiply test for " << matrixSystem.name() << 
        std::endl;
    }
    return multiply_test(A,B,C_check,AT,BT,epsilon,comm,verbose);
  }
  else if(op == "add"){
    if(localProc == 0 && verbose){
      std::cout << "Running add test for " << matrixSystem.name() << 
        std::endl;
    }
    return add_test(A,B,C_check,AT,BT,epsilon,comm,verbose);
  }
  else{
    if(localProc == 0 && verbose){
      std::cout<< "Unrecognize matrix operation: " << op << ".";
    }
    return -1;
  }


}


TEUCHOS_UNIT_TEST(Tpetra_MatMat, TwoProcTest){

  typedef Map<int>                                       Map;
  typedef CrsMatrix<double,int>                          CrsMatrix;
  typedef CrsMatrix::mat_vec_type                       MatVec;
  typedef CrsMatrix::node_type                          DNode;
   
  RCP<const Comm<int> > comm = DefaultPlatform::getDefaultPlatform().getComm();

  int thisproc = comm->getRank();
  int numprocs = comm->getSize();

  if(numprocs <2){ return; }


  //RCP< Teuchos::FancyOStream > out = Teuchos::fancyOStream(rcp(&std::cout,false),"",0,false,10,false,true);
  MMdebug::debug_stream = rcpFromRef(out);
  MMdebug::debug_level  = Teuchos::VERB_NONE;

  // set up a row-std::map with 2 global elements,
  // 1 on each proc.
  const int numGlobalRows = 2;
  ArrayRCP<int> myrow(1);
  if (thisproc == 1) myrow[0] = 7;
  else               myrow[0] = 3;
  RCP<const Map> rowmap = rcp(new Map(numGlobalRows, myrow(), 0, comm));

  //set up a domain-std::map with columns 0 - 4 on proc 0,
  //and columns 5 - 9 on proc 1.
  const int numGlobalCols = 10;
  const int numMyCols = 5;
  ArrayRCP<int> mycols(numGlobalCols);
  for(int i=0; i<numGlobalCols; ++i) {
    mycols[i] = i;
  }
  RCP<const Map> domainmap = 
    rcp(new Map(numGlobalCols, mycols(thisproc*numMyCols,numMyCols), 0, comm));

  // now create matrices A, B and C with rowmap; the second argument is just the suggested allocation size
  RCP<CrsMatrix> A = rcp(new CrsMatrix(rowmap, 1));
  A->setObjectLabel("Factor Matrix A");
  RCP<CrsMatrix> C = rcp(new CrsMatrix(rowmap, 1));
  C->setObjectLabel("Product matrix C");

  ArrayRCP<double> coefs(numGlobalCols);
  for(int i=0; i<numGlobalCols; ++i) {
    coefs[i] = 1.0*i;
  }

  A->insertGlobalValues(
    myrow[0], 
    mycols(thisproc*numMyCols, numMyCols), 
    coefs(thisproc*numMyCols, numMyCols));

  A->fillComplete(domainmap, rowmap);
  // A->describe(*out, Teuchos::VERB_EXTREME);

  MatrixMatrix::Multiply(*A, false, *A, true, *C);
  C->describe(out, Teuchos::VERB_EXTREME);


  TEST_EQUALITY(C->getGlobalNumEntries(), 4);


}

TEUCHOS_UNIT_TEST(Tpetra_MatMat, test_find_rows){
  RCP<const Comm<int> > comm = DefaultPlatform::getDefaultPlatform().getComm();
  int numprocs = comm->getSize();
  int localproc = comm->getRank();
  int numlocalrows = 2;
  global_size_t numglobalrows = numprocs*numlocalrows;
  RCP<Map<int> > rowmap = 
    rcp(new Map<int>(numglobalrows, 0, comm));
  CrsMatrix<double, int> matrix(rowmap, numglobalrows);

  Array<int> cols(numglobalrows);
  Array<double> vals(numglobalrows);

  for(size_t j=0; j<numglobalrows; ++j) {
    cols[j] = j;
    vals[j] = 1.0;
  }

  RCP<Map<int> > colmap = 
    rcp(new Map<int>(-1, cols(), 0, comm));

  for(int i=0; i<numlocalrows; ++i) {
    Array<int> row(1,localproc*numlocalrows+i);
    matrix.insertGlobalValues(
      row[0], row.view(0,1),  vals.view(i,1) );
  }

  matrix.fillComplete();

  typedef Kokkos::DefaultNode::DefaultNodeType DNode;
  typedef Kokkos::DefaultKernels<double, int, DNode>::SparseOps SpMatOps;

  RCP<const Map<int> > map_rows = 
    MMdetails::find_rows_containing_cols<double, int, int, DNode, SpMatOps>(matrix, colmap);

  /*if (map_rows->getNodeNumElements() != numglobalrows) {
    if(localproc ==0){
      std::cout << "Error in test_find_rows" << std::endl <<
      "Num elements found: " << map_rows->getNodeNumElements() << 
      std::endl <<
      "Num global rows: " << numglobalrows << std::endl;
    }
    return(-1);
  }*/
  TEST_EQUALITY(map_rows->getNodeNumElements(), numglobalrows);

}

TEUCHOS_UNIT_TEST(Tpetra_MatMat, operations_test){
  RCP<const Comm<int> > comm = DefaultPlatform::getDefaultPlatform().getComm();
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
      
    TEST_EQUALITY(
    Tpetra::run_test<int>(
      comm, 
      matrixSystems->sublist(it->first), 
      write_result_hb, 
      verbose), 0);
  }
}


} //namespace Tpetra

