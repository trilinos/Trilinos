#include "Tpetra_MatrixMatrix.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace Tpetra{
/*extern
template<
  class Scalar,
  class LocalOrdinal,
  class GlobalOrdinal,
  class Node,
  class SpMatVec,
  class SpMatSlv>
RCP<Map<LocalOrdinal, GlobalOrdinal, Node> > find_row_containing_cols(
  RCP<const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv> > M,
  RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > colmap);
}*/

static const double defaultEpsilon = 1e-10;

template<class Ordinal>
int test_find_rows(RCP<const Comm<Ordinal> > Comm)
{
  int numprocs = Comm->getSize();
  int localproc = Comm->getRank();
  int numlocalrows = 2;
  global_size_t numglobalrows = numprocs*numlocalrows;
  RCP<Map<int> > rowmap = 
    rcp(new Map<int>(numglobalrows, 0, Comm));
  CrsMatrix<double, int> matrix(rowmap, numglobalrows);

  Array<int> cols(numglobalrows);
  Array<double> vals(numglobalrows);

  for(size_t j=0; j<numglobalrows; ++j) {
    cols[j] = j;
    vals[j] = 1.0;
  }

  RCP<Map<int> > colmap = 
    rcp(new Map<int>(-1, cols(), 0, Comm));

  for(int i=0; i<numlocalrows; ++i) {
    Array<int> row(1,localproc*numlocalrows+i);
    matrix.insertGlobalValues(
      row[0], row.view(0,1),  vals.view(i,1) );
  }

  matrix.fillComplete();

  typedef Kokkos::DefaultNode::DefaultNodeType DNode;
  typedef typename Kokkos::DefaultKernels<double, int, DNode>::SparseOps SpMatOps;

  RCP<const Map<int> > map_rows = 
    MMdetails::find_rows_containing_cols<double, int, int, DNode, SpMatOps>(matrix, colmap);

  if (map_rows->getNodeNumElements() != numglobalrows) {
    if(localproc ==0){
      std::cout << "Error in test_find_rows" << std::endl <<
      "Num elements found: " << map_rows->getNodeNumElements() << 
      std::endl <<
      "Num global rows: " << numglobalrows << std::endl;
    }
    return(-1);
  }


  return(0);
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

  int localProc = comm->getRank();

  RCP<CrsMatrix<double,int> > computedC = null;
  RCP<const Map<int> > rowmap = AT ? A->getDomainMap() : A->getRowMap();

  computedC = rcp( new CrsMatrix<double,int>(rowmap, 1));

  MatrixMatrix::Add(*A, false, 1.0, *B, false, 1.0, *computedC);

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
    if (localProc == 0 && verbose) {
      std::cout << "Test Passed" << std::endl;
      std::cout << "||(A+B)-C||/||C|| " << resultVal1 <<std::endl ;
      std::cout << "||(A+B)-C||/||C|| " << resultVal2 <<std::endl << std::endl;
    }
    return 0;
  }
  else {
    if (localProc == 0) {
      std::cout << "Test Failed: " << std::endl;
      std::cout << "||(A+B)-C||/||C|| " << resultVal1 <<std::endl ;
      std::cout << "||(A+B)-C||/||C|| " << resultVal2 <<std::endl << std::endl;
    }
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

  int localProc = comm->getRank();
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
    if (localProc == 0 && verbose) {
      std::cout << "Test Passed" << std::endl;
      std::cout << "||A*B-C||/||C|| " << diff_result << std::endl << std::endl;
    }
  }
  else {
    return_code = -1;
    if (localProc == 0) {
      std::cout << "Test Failed: ||A*B-C||/||C|| " << diff_result << std::endl << std::endl;
    }
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

template<class Ordinal>
int two_proc_test(RCP<const Comm<Ordinal> > Comm,
  bool verbose)
{
  typedef Map<int>                                       Map;
  typedef CrsMatrix<double,int>                          CrsMatrix;
  typedef typename CrsMatrix::mat_vec_type                       MatVec;
  typedef typename CrsMatrix::node_type                          DNode;

  (void)verbose;
  int thisproc = Comm->getRank();
  int numprocs = Comm->getSize();
  int err =0;

  //only run this test on 2 procs
  if (numprocs != 2) return(0);

  RCP< Teuchos::FancyOStream > out = Teuchos::fancyOStream(rcp(&std::cout,false),"",0,false,10,false,true);
  MMdebug::debug_stream = out;
  MMdebug::debug_level  = Teuchos::VERB_NONE;

  // set up a row-std::map with 2 global elements,
  // 1 on each proc.
  const int numGlobalRows = 2;
  ArrayRCP<int> myrow(1);
  if (thisproc == 1) myrow[0] = 7;
  else               myrow[0] = 3;
  RCP<const Map> rowmap = rcp(new Map(numGlobalRows, myrow(), 0, Comm));

  //set up a domain-std::map with columns 0 - 4 on proc 0,
  //and columns 5 - 9 on proc 1.
  const int numGlobalCols = 10;
  const int numMyCols = 5;
  ArrayRCP<int> mycols(numGlobalCols);
  for(int i=0; i<numGlobalCols; ++i) {
    mycols[i] = i;
  }
  RCP<const Map> domainmap = rcp(new Map(numGlobalCols, mycols(thisproc*numMyCols,numMyCols), 0, Comm));

  // now create matrices A, B and C with rowmap; the second argument is just the suggested allocation size
  RCP<CrsMatrix> A = rcp(new CrsMatrix(rowmap, 1));
  A->setObjectLabel("Factor Matrix A");
  RCP<CrsMatrix> C = rcp(new CrsMatrix(rowmap, 1));
  C->setObjectLabel("Product matrix C");

  ArrayRCP<double> coefs(numGlobalCols);
  for(int i=0; i<numGlobalCols; ++i) {
    coefs[i] = 1.0*i;
  }

  A->insertGlobalValues(myrow[0], mycols(thisproc*numMyCols, numMyCols), coefs(thisproc*numMyCols, numMyCols));

  A->fillComplete(domainmap, rowmap);
  // A->describe(*out, Teuchos::VERB_EXTREME);

  MatrixMatrix::Multiply(*A, false, *A, true, *C);
  C->describe(*out, Teuchos::VERB_EXTREME);

  if (C->getGlobalNumEntries() != 4) {
    err += 1;
  }

  return(err);
}
/*
int time_matrix_matrix_multiply(Epetra_Comm& Comm, bool verbose)
{

  const int magic_num = 3000;
  // 2009/02/23: rabartl: If you are going to do a timing test you need to
  // make this number adjustable form the command-line and you need to put in
  // a real test that compares against hard numbers for pass/fail.

  int localn = magic_num/Comm.NumProc();

  Epetra_CrsMatrix* A = create_crsmatrix(Comm.NumProc(),
                                                Comm.MyPID(),
                                                localn);

  Epetra_CrsMatrix* C = new Epetra_CrsMatrix(Copy, A->RowMap(), 0);

  Epetra_Time timer(Comm);
  double start_time = timer.WallTime();

  int err = EpetraExt::MatrixMatrix::Multiply(*A, false, *A, false, *C);

  int globaln = localn*Comm.NumProc();
  if (verbose) {
    std::cout << "size: " << globaln << "x"<<globaln<<", C = A*A, time: "
       << timer.WallTime()-start_time << std::endl;
  }

  C->FillComplete();

  start_time = timer.WallTime();

  err = EpetraExt::MatrixMatrix::Multiply(*A, false, *A, false, *C);

  if (verbose) {
    std::cout << "size: " << globaln << "x"<<globaln<<", C = A*A, time: "
       << timer.WallTime()-start_time << " (C already Filled)\n" <<std::endl;
  }

  delete C;

  C = new Epetra_CrsMatrix(Copy, A->RowMap(), 0);

  start_time = timer.WallTime();

  err = EpetraExt::MatrixMatrix::Multiply(*A, false, *A, true, *C);

  if (verbose) {
    std::cout << "size: " << globaln << "x"<<globaln<<", C = A*A^T, time: "
       << timer.WallTime()-start_time << std::endl;
  }

  C->FillComplete();

  start_time = timer.WallTime();

  err = EpetraExt::MatrixMatrix::Multiply(*A, false, *A, true, *C);

  if (verbose) {
    std::cout << "size: " << globaln << "x"<<globaln<<", C = A*A^T, time: "
       << timer.WallTime()-start_time << " (C already Filled)\n" <<std::endl;
  }

  delete C;

  C = new Epetra_CrsMatrix(Copy, A->RowMap(), 0);

  start_time = timer.WallTime();

  err = EpetraExt::MatrixMatrix::Multiply(*A, true, *A, false, *C);

  if (verbose) {
    std::cout << "size: " << globaln << "x"<<globaln<<", C = A^T*A, time: "
       << timer.WallTime()-start_time << std::endl;
  }

  C->FillComplete();

  start_time = timer.WallTime();

  err = EpetraExt::MatrixMatrix::Multiply(*A, true, *A, false, *C);

  if (verbose) {
    std::cout << "size: " << globaln << "x"<<globaln<<", C = A^T*A, time: "
       << timer.WallTime()-start_time << " (C already Filled)\n"<<std::endl;
  }

  delete C;

  C = new Epetra_CrsMatrix(Copy, A->RowMap(), 0);

  start_time = timer.WallTime();

  err = EpetraExt::MatrixMatrix::Multiply(*A, true, *A, true, *C);

  if (verbose) {
    std::cout << "size: " << globaln << "x"<<globaln<<", C = A^T*A^T, time: "
       << timer.WallTime()-start_time << std::endl;
  }

  C->FillComplete();

  start_time = timer.WallTime();

  err = EpetraExt::MatrixMatrix::Multiply(*A, true, *A, true, *C);

  if (verbose) {
    std::cout << "size: " << globaln << "x"<<globaln<<", C = A^T*A^T, time: "
       << timer.WallTime()-start_time << " (C already Filled)\n" <<std::endl;
  }

  delete C;

  delete A;

  return(err);
}
*/
template<
  class Scalar,
  class LocalOrdinal,
  class GlobalOrdinal,
  class Node,
  class SpMatVec,
  class SpMatSlv,
  class CommOrdinal>
RCP<CrsMatrix<LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv> > 
create_crsmatrix(
  RCP<const Comm<CommOrdinal> > comm,
  size_t local_n,
  bool callFillComplete,
  bool symmetric)
{
  int numProcs = comm->getSize();
  global_size_t global_num_rows = numProcs*local_n;

  RCP<Map<LocalOrdinal, GlobalOrdinal, Node> > rowmap = 
    rcp(global_num_rows, local_n, 0, comm);

  size_t nnz_per_row = 9;
  RCP<CrsMatrix<LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv> > matrix =
    rcp(new CrsMatrix<LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv>(rowmap, nnz_per_row));

  // Add  rows one-at-a-time
  ArrayRCP<Scalar> negOne(1, -ScalarTraits<Scalar>::one());
  ArrayRCP<Scalar> posTwo(1, 2*ScalarTraits<Scalar>::one());
  ArrayRCP<Scalar> val_L(1, symmetric ? negOne : 0.5*negOne);

  GlobalOrdinal GO1 = OrdinalTraits<GlobalOrdinal>::one();
  GlobalOrdinal GO0 = OrdinalTraits<GlobalOrdinal>::zero();

  for (int i=0; i<local_n; i++) {
    ArrayRCP<GlobalOrdinal> GlobalRow(1, matrix->getRowMap()->getGlobalElement(i));
    ArrayRCP<GlobalOrdinal> RowLess1(1,GlobalRow - GO1);
    ArrayRCP<GlobalOrdinal> RowPlus1(1,GlobalRow + GO1);
    ArrayRCP<GlobalOrdinal> RowLess5(1,GlobalRow - (5*GO1));
    ArrayRCP<GlobalOrdinal> RowPlus5(1,GlobalRow + (5*GO1));
    ArrayRCP<GlobalOrdinal> RowLess9(1,GlobalRow - (9*GO1));
    ArrayRCP<GlobalOrdinal> RowPlus9(1,GlobalRow + (9*GO1));
    ArrayRCP<GlobalOrdinal> RowLess24(1,GlobalRow - (24*GO1));
    ArrayRCP<GlobalOrdinal> RowPlus24(1,GlobalRow + (24*GO1));
    ArrayRCP<GlobalOrdinal> RowLess48(1,GlobalRow - (48*GO1));
    ArrayRCP<GlobalOrdinal> RowPlus48(1,GlobalRow + (48*GO1));

//    if (!symmetric) RowLess5 -= 2;

    if (RowLess48>=GO0) {
      matrix->insertGlobalValues(GlobalRow, RowLess48(), val_L());
    }
    if (RowLess24>=GO0) {
      matrix->insertGlobalValues(GlobalRow, RowLess24(), val_L());
    }
    if (RowLess9>=GO0) {
      matrix->insertGlobalValues(GlobalRow, RowLess9(), val_L());
    }
    if (RowLess5>=GO0) {
      matrix->insertGlobalValues(GlobalRow, RowLess5(), val_L());
    }
    if (RowLess1>=GO0) {
      matrix->insertGlobalValues(GlobalRow, RowLess1(), val_L());
    }
    if (RowPlus1<global_num_rows) {
      matrix->insertGlobalValues(GlobalRow, RowPlus1(), negOne());
    }
    if (RowPlus5<global_num_rows) {
      matrix->insertGlobalValues(GlobalRow, RowPlus5(), negOne());
    }
    if (RowPlus9<global_num_rows) {
      matrix->insertGlobalValues(GlobalRow, RowPlus9(), negOne());
    }
    if (RowPlus24<global_num_rows) {
      matrix->insertGlobalValues(GlobalRow, RowPlus24(), negOne());
    }
    if (RowPlus48<global_num_rows) {
      matrix->insertGlobalValues(GlobalRow, RowPlus48(), negOne());
    }

    matrix->insertGlobalValues(GlobalRow, GlobalRow(), posTwo());
  }

  if (callFillComplete) {
    matrix->fillComplete();
  }

  return(matrix);
}

} //namespace

/*
template<class Ordinal>
int time_matrix_matrix_multiply(RCP<Comm<Ordinal> > Comm,
  bool verbose);
*/
/////////////////////////////////////
//Global variable!!!!
//char* path;
////////////////////////////////////

int main(int argc, char* argv[]) {
  Teuchos::CommandLineProcessor clp;
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  bool verbose = false;
  std::string matnamesFile;
  bool write_result_hb = false;
  clp.setOption("matnames-file", &matnamesFile, "A file containing a list of matricies we'll import", true);
  clp.setOption("writeresults", "no-write-results", &write_result_hb, "Whether or not to write the resutling matricies to hb files");
  clp.setOption("v", "not-verbose", &verbose, "Whether or not to use verbose output");
  clp.parse(argc, argv);


  int localProc = comm->getRank();

  int err = Tpetra::two_proc_test(comm, verbose);
  if (err != 0) {
    if (localProc == 0 && verbose) {
      std::cerr << "two_proc_test returned err=="<<err<<std::endl;
    }
    return(err);
  }

  err = Tpetra::test_find_rows(comm);
  if (err != 0 && comm->getRank() ==0) {
    if (localProc == 0 && verbose) {
      std::cerr << "test_find_rows returned err=="<<err<<std::endl;
    }
    return(err);
  } 
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
      
    if( Tpetra::run_test<int>(
      comm, 
      matrixSystems->sublist(it->first), 
      write_result_hb, 
      verbose) != 0)
    {
      err = -1;
    }

  }


  //THIS SHOULD BECOME A TEMPLATED UNIT TEST
  /*RCP<CrsMatrix<double, int> > D = 
     create_crsmatrix<double, int, int,
     Kokkos::DefaultNode::DefaultNodeType,
     Kokkos::DefaultSparseMultiply<double,int,Kokkos::DefaultNode::DefaultNodeType>,
     Kokkos::DefaultSparseSolve<double,int,Kokkos::DefaultNode::DefaultNodeType> >  (comm, (size_t)10, true, false);

  std::cout << "D: \n"  << *D << std::endl;*/

  //EpetraExt::MatrixMatrix::Add(*D, true, 0.5, *D, 0.5);

//  std::cout << "symm D: \n"  << *D << std::endl;

  //delete D;

  /*
  if (err == 0) {
    err = time_matrix_matrix_multiply(comm, verbose);
  }*/

  int global_err = err;

  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &err, &global_err);

  return(global_err);

}
