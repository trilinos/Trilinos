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

static const double defaultEpsilon = 1e-7;

template<class Ordinal>
int run_test(RCP<const Comm<Ordinal> > Comm,
  Teuchos::ParameterList matrixSystem,
  bool result_mtx_to_file=false,
  bool verbose=false);

template<class Ordinal>
int simple_add_test(RCP<const Comm<Ordinal> > Comm, bool verbose);

template<class CommOrdinal>
int two_proc_test(RCP<const Comm<CommOrdinal> > Comm,
  bool verbose=false);

template<class Ordinal>
int test_find_rows(RCP<const Comm<Ordinal> > Comm);

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
  bool symmetric);

template<class Ordinal>
int simple_add_test(RCP<const Comm<Ordinal> > Comm, bool verbose){
  int localProc = Comm->getRank();
  RCP<CrsMatrix<double,int> > A = null;
  RCP<CrsMatrix<double,int> > B = null;
  RCP<CrsMatrix<double,int> > C = null;

  if (localProc == 0 && verbose) {
    std::cout << "Before reading" << std::endl;
  }
  Utils::readHBMatrix("matrices/addA3.hb", Comm, Kokkos::DefaultNode::getDefaultNode(), A);
  Utils::readHBMatrix("matrices/addB3.hb", Comm, Kokkos::DefaultNode::getDefaultNode(), B);
  Utils::readHBMatrix("matrices/addC3.hb", Comm, Kokkos::DefaultNode::getDefaultNode(), C);
  if (localProc == 0 && verbose) {
    std::cout << "after reading" << std::endl;
  }

  typedef Kokkos::DefaultNode::DefaultNodeType DNode;

  

  MatrixMatrix<
    double, 
    int,
    int,
    DNode,
    typename Kokkos::DefaultKernels<double,int,DNode>::SparseOps>::
  Add(A, false, 1.0, B, 1.0);

  MatrixMatrix<
    double, 
    int,
    int,
    DNode,
    typename Kokkos::DefaultKernels<double,int,DNode>::SparseOps>::
  Add(C, false, -1.0, B, 1.0);

  double calculated_euc_norm = B->getEuclideanNorm();
  double c_euc_norm = C->getEuclideanNorm();
  double resultVal = calculated_euc_norm/c_euc_norm;

  int return_code =0;
  if (resultVal < defaultEpsilon) {
    if (localProc == 0 && verbose) {
      std::cout << "Test Passed" << std::endl;
      std::cout << "||(A+B)-C||/||C|| " << resultVal << std::endl;
    }
  }
  else {
    return_code = -1;
    if (localProc == 0) {
      std::cout << "Test Failed: ||(A+B)-C||/||C|| " << resultVal << std::endl;
    }
  }

  return(return_code);


}

template<class Ordinal>
int test_find_rows(RCP<const Comm<Ordinal> > Comm)
{
  int numprocs = Comm->getSize();
  int localproc = Comm->getRank();
  int numlocalrows = 2;
  global_size_t numglobalrows = numprocs*numlocalrows;
  RCP<Map<int> > rowmap = 
    rcp(new Map<int>(numglobalrows, 0, Comm));
  RCP<CrsMatrix<double, int> > matrix = 
    rcp(new CrsMatrix<double, int>(rowmap, numglobalrows));

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
    matrix->insertGlobalValues(
      row[0], row.view(0,1),  vals.view(i,1) );
  }

  matrix->fillComplete();

  typedef Kokkos::DefaultNode::DefaultNodeType DNode;

  RCP<const Map<int> > map_rows = 
    MatrixMatrix<
      double, 
      int,
      int,
      DNode,
      typename Kokkos::DefaultKernels<double,int,DNode>::SparseOps
      >::find_rows_containing_cols(
      matrix.getConst(), colmap);

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

/*int expand_name_list(const char* newname,
                     const char**& names,
                     int& alloc_len,
                     int& num_names)
{
  int offset = num_names;
  if (offset >= alloc_len) {
    int alloc_increment = 8;
    const char** newlist = new const char*[alloc_len+alloc_increment];
    for(int i=0; i<offset; ++i) {
      newlist[i] = names[i];
    }
    delete [] names;
    names = newlist;
    alloc_len += alloc_increment;
    for(int i=offset; i<alloc_len; ++i) {
      names[i] = NULL;
    }
  }

  names[offset] = newname;
  ++num_names;
  return(0);
}

template<class Ordinal>
int broadcast_name(RCP<const Comm<Ordinal> > Comm, const char*& name)
{
  if (Comm->getSize() < 2) return(0);

  Ordinal len;
  int localProc = Comm->getRank();
  if (localProc == 0) {
    len = (Ordinal)strlen(name)+1;
    
    //MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
	Teuchos::broadcast(*Comm, 0, 1, &len);
    //MPI_Bcast((void*)name, len, MPI_CHAR, 0, MPI_COMM_WORLD);
	Comm->broadcast(0, len, (char*)name);

  }
  else {
    //MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
	Teuchos::broadcast(*Comm, 0, 1, &len);
    name = new char[len];
    //MPI_Bcast((void*)name, len, MPI_CHAR, 0, MPI_COMM_WORLD);
	Comm->broadcast(0, len, (char*)name);
  }

  return(0);
}

template<class Odrinal>
int read_input_file(RCP<const Comm<Odrinal> > Comm,
                    const char* input_file_name,
                    const char**& filenames,
                    int& numfiles,
                    int& numfilenames_allocated)
{
  int local_err = 0, global_err = 0;
  std::ifstream* infile = NULL;
  int pathlen = path != 0 ? (int)strlen(path): 0;

  if (Comm->getRank() == 0) {
    char* full_name = NULL;
    int filenamelen = input_file_name != 0 ? (int)strlen(input_file_name) : 0;

    full_name = new char[pathlen+filenamelen+2];
    if (path != 0) {
      sprintf(full_name, "%s/%s",path,input_file_name);
    }
    else {
      sprintf(full_name, "%s", input_file_name);
    }

    infile = new std::ifstream(full_name);
    if (!(*infile)) {
      local_err = -1;
      delete infile;
    }
    delete [] full_name;
  }

  Teuchos::reduceAll(*Comm, Teuchos::REDUCE_SUM, 1, &local_err, &global_err);
  if (global_err != 0) {
    return(global_err);
  }


  if (Comm->getRank() == 0) {
    int linelen = 512;
    char* line = NULL;

    std::ifstream& ifile = *infile;
    while(!ifile.eof()) {
      line = new char[pathlen+1+linelen];
      if (pathlen>0) {
        sprintf(line,"%s/",path);
        ifile.getline(&(line[pathlen+1]), linelen);
      }
      else {
        ifile.getline(line, linelen);
      }

      if (ifile.fail()) {
	delete [] line;
        break;
      }
      if (strchr(line, '#') == NULL) {
        expand_name_list(line, filenames, numfilenames_allocated, numfiles);
      }
      else {
        delete [] line;
      }
    }

    Teuchos::broadcast(*Comm, 0, 1, &numfiles);
    for(int i=0; i<numfiles; ++i) {
      broadcast_name(Comm, filenames[i]);
    }

    delete infile;
  }
  else {
    Teuchos::broadcast(*Comm, 1, &numfiles);
    filenames = new const char*[numfiles];
    numfilenames_allocated = numfiles;
    for(int i=0; i<numfiles; ++i) {
      broadcast_name(Comm, filenames[i]);
    }
  }
  
  return(0);
}
*/

template<class Ordinal>
int run_test(RCP<const Comm<Ordinal> > Comm,
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


  int localProc = Comm->getRank();

  if (localProc == 0 && verbose) {
    std::cout << "Testing C=A"<<AT<<"*B"<<BT<< "; A:" << A_file
              << ", B:" << B_file << ", C:" << C_file << std::endl;
  }

  RCP<CrsMatrix<double,int> > A = null;
  RCP<CrsMatrix<double,int> > B = null;
  RCP<CrsMatrix<double,int> > C = null;
  RCP<CrsMatrix<double,int> > C_check = null;

/*
  Utils::readMatrixMarketMatrix(A_file, Comm, Kokkos::DefaultNode::getDefaultNode(), A);
  Utils::readMatrixMarketMatrix(B_file, Comm, Kokkos::DefaultNode::getDefaultNode(), B);
  Utils::readMatrixMarketMatrix(C_file, Comm, Kokkos::DefaultNode::getDefaultNode(), C_check);*/
  Utils::readHBMatrix(A_file, Comm, Kokkos::DefaultNode::getDefaultNode(), A);
  Utils::readHBMatrix(B_file, Comm, Kokkos::DefaultNode::getDefaultNode(), B);
  Utils::readHBMatrix(C_file, Comm, Kokkos::DefaultNode::getDefaultNode(), C_check);


  /*RCP<Teuchos::basic_FancyOStream<char> > outstream = 
    Teuchos::getFancyOStream(rcpFromRef(std::cout));
  B->describe(*outstream, Teuchos::VERB_EXTREME);*/
  RCP<const Map<int> > rowmap = AT ? A->getDomainMap() : A->getRowMap();

  C = rcp( new CrsMatrix<double,int>(rowmap, 1));

  RCP<const CrsMatrix<double,int> > constA = A;
  RCP<const CrsMatrix<double,int> > constB = B;
  typedef Kokkos::DefaultNode::DefaultNodeType DNode;


  MatrixMatrix<
    double, 
    int,
    int,
    DNode,
    typename Kokkos::DefaultKernels<double,int,DNode>::SparseOps>::
  Multiply(constA, AT, constB, BT, C);

//  std::cout << "A: " << *A << std::endl << "B: "<<*B<<std::endl<<"C: "<<*C<<std::endl;
  //if (result_mtx_to_file) {
   // EpetraExt::RowMatrixToMatrixMarketFile("result.mtx", *C);
  //}

  
  double c_euc_norm = C_check->getEuclideanNorm();

  MatrixMatrix<
    double, 
    int,
    int,
    DNode,
    typename Kokkos::DefaultKernels<double,int,DNode>::SparseOps>::
  Add(C, false, -1.0, C_check, 1.0);

  double c_check_euc_norm = C_check->getEuclideanNorm();

  double diff_result = c_check_euc_norm/c_euc_norm;

  int return_code =0;
  if (diff_result < epsilon) {
    if (localProc == 0 && verbose) {
      std::cout << "Test Passed" << std::endl;
      std::cout << "||A*B-C||/||C|| " << diff_result << std::endl;
    }
  }
  else {
    return_code = -1;
    if (localProc == 0) {
      std::cout << "Test Failed: ||A*B-C||/||C|| " << diff_result << std::endl;
    }
  }

  return(return_code);
}

template<class Ordinal>
int two_proc_test(RCP<const Comm<Ordinal> > Comm,
  bool verbose)
{
  typedef Map<int>                                       Map;
  typedef CrsMatrix<double,int>                          CrsMatrix;
  typedef typename CrsMatrix::mat_vec_type                       MatVec;
  typedef typename CrsMatrix::node_type                          DNode;
  typedef MatrixMatrix< double, int, int, DNode, MatVec> MatMat;

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

  MatMat::Multiply(A, false, A, true, C);
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

  int err = Tpetra::two_proc_test(comm, verbose);
  if (err != 0) {
    std::cerr << "two_proc_test returned err=="<<err<<std::endl;
    return(err);
  }

  err = Tpetra::simple_add_test(comm, verbose);
  if(err != 0){
    std::cerr << "two_proc_test returned err=="<<err<<std::endl;
    return(err);
  }


  err = Tpetra::test_find_rows(comm);
  if (err != 0 && comm->getRank() ==0) {
    std::cerr << "test_find_rows returned err=="<<err<<std::endl;
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

/**
 * Makes sure the multiply doesn't work if given bad matrices
 */
  Teuchos::ParameterList wrongList;
  wrongList.set("A", "matrices/wrong_m.hb");
  wrongList.set("B", "matrices/wrong_tce.hb");
  wrongList.set("C", "matrices/wrong_d.hb");
  wrongList.set("TransA", false);
  wrongList.set("TransB", true);
  if(
    Tpetra::run_test<int>(
      comm, 
      wrongList, 
      write_result_hb, 
      verbose) == 0 )
  {
    err = -1;
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
