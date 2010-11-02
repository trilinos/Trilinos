#include "Tpetra_MatrixMatrix.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

/*namespace Tpetra{
extern
template<
  class Scalar,
  class LocalOrdinal,
  class GlobalOrdinal,
  class Node,
  class SpMatVec,
  class SpMatSlv>
Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> > find_row_containing_cols(
  Teuchos::RCP<const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv> > M,
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > colmap);
}*/

/*template<class Ordinal>
int read_input_file(Teuchos::RCP<const Teuchos::Comm<Ordinal> > Comm,
  const char* input_file_name,
  const char**& filenames,
  int& numfiles,
  int& numfilenames_allocated);

template<class Ordinal>
int read_matrix_file_names(Teuchos::RCP<Teuchos::Comm<Ordinal> > Comm,
  const char* input_file_name,
  char*& A_file,
  bool& transA,
  char*& B_file,
  bool& transB,
  char*& C_file);

template<class Ordinal>
int broadcast_name(Teuchos::RCP<const Teuchos::Comm<Ordinal> > Comm, const char*& name);

template<class Ordinal, class LocalOrdinal, class GlobalOrdinal, class Node>
int create_maps(
  Teuchos::RCP<Teuchos::Comm<Ordinal> > Comm,
  const char* input_file_name,
  Teuchos::RCP<Tpetra::Map<LocalOrdinal, Global, Node> > row_map,
  Teuchos::RCP<Tpetra::Map<LocalOrdinal, Global, Node> >col_map,
  Teuchos::RCP<Tpetra::Map<LocalOrdinal, Global, Node> > range_map,
  Teuchos::RCP<Tpetra::Map<LocalOrdinal, Global, Node> > domain_map);

template<
  class Ordinal,
  class Scalar,
  class LocalOrdinal,
  class GlobalOrdinal,
  class Node,
  class SpMatVec,
  class SpMatSlv>
int read_matrix(
  const char* filename,
  Teuchos::RCP<Teuchos::Comm<Ordinal> > Comm,
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal, Global, Node> > rowmap,
  Teuchos::RCP<Tpetra::Map<LocalOrdinal, Global, Node> > colmap,
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal, Global, Node> > rangemap,
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal, Global, Node> > domainmap,
  Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv> > mat);
*/
template<class Ordinal>
int run_test(Teuchos::RCP<const Teuchos::Comm<Ordinal> > Comm,
  Teuchos::ParameterList matrixSystem,
  bool result_mtx_to_file=false,
  bool verbose=false);

template<class CommOrdinal>
int two_proc_test(Teuchos::RCP<const Teuchos::Comm<CommOrdinal> > Comm,
  bool verbose=false);

template<class Ordinal>
int test_find_rows(Teuchos::RCP<const Teuchos::Comm<Ordinal> > Comm);

template<
  class Scalar,
  class LocalOrdinal,
  class GlobalOrdinal,
  class Node,
  class SpMatVec,
  class SpMatSlv,
  class CommOrdinal>
Teuchos::RCP<Tpetra::CrsMatrix<LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv> > 
create_crsmatrix(
  Teuchos::RCP<const Teuchos::Comm<CommOrdinal> > comm,
  size_t local_n,
  bool callFillComplete,
  bool symmetric);
/*
template<class Ordinal>
int time_matrix_matrix_multiply(Teuchos::RCP<Teuchos::Comm<Ordinal> > Comm,
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
  /*int write = 0;
  bool path_specified = false;
  char* input_file = NULL;
  bool input_file_specified = false;

  if (comm->getRank()==0) {
    for(int ii=0; ii<argc; ++ii) {
      if (!strcmp("-write_result", argv[ii])) write_result_mtx = true;
      if (!strcmp("-v", argv[ii])) verbose = true;
      if (!strcmp("-i", argv[ii])) {
        input_file = argv[ii+1];
        input_file_specified = true;
      }
      if (!strcmp("-d",argv[ii])) {
        path = argv[ii+1];
        path_specified = true;
      }
    }
    write = write_result_mtx ? 1 : 0;
  }

  Teuchos::broadcast(*comm, 0, 1, &write);
  if (write) write_result_mtx = true;

  if (!path_specified) {
    path = new char[32];
    sprintf(path, "%s", ".");
  }*/

  int err = two_proc_test(comm, verbose);
  if (err != 0) {
    std::cerr << "two_proc_test returned err=="<<err<<std::endl;
    return(err);
  }

/*  if (!input_file_specified) {
    input_file = new char[64];
    sprintf(input_file, "%s", "infiles");
  }

  const char** filenames = NULL;
  int numfiles = 0;
  int numfilenames_allocated = 0;

  err = read_input_file(comm, input_file,
                        filenames, numfiles, numfilenames_allocated);
  if (err != 0) {
    if (path_specified) path_specified = false;
    sprintf(path, "%s", "./MatrixMatrix");
    read_input_file(comm, input_file,
                    filenames, numfiles, numfilenames_allocated);
  }
*/
  err = test_find_rows(comm);
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
      
    if( run_test<int>(
      comm, 
      matrixSystems->sublist(it->first), 
      write_result_hb, 
      verbose) != 0)
    {
      err = -1;
    }

  }


  Teuchos::ParameterList wrongList;
  wrongList.set("A", "wrong_m.hb");
  wrongList.set("B", "wrong_tce.hb");
  wrongList.set("C", "wrong_d.hb");
  wrongList.set("TransA", false);
  wrongList.set("TransB", true);
  if(
    run_test<int>(
      comm, 
      wrongList, 
      write_result_hb, 
      verbose) == 0 )
  {
    err = -1;
  }

/*
  for(int i=0; i<numfiles; ++i) {
    err = run_test(comm, matnamesFile, write_result_mtx, verbose);
  //  delete [] filenames[i];
    if (err != 0) break;
  }*/
/*
  for(int j=numfiles; j<numfilenames_allocated; ++j) {
    delete [] filenames[j];
  }

  delete [] filenames;

  if (!input_file_specified) delete [] input_file;
  if (!path_specified) delete [] path;*/

  //THIS SHOULD BECOME A TEMPLATED UNIT TEST
  /*Teuchos::RCP<Tpetra::CrsMatrix<double, int> > D = 
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

template<class Ordinal>
int test_find_rows(Teuchos::RCP<const Teuchos::Comm<Ordinal> > Comm)
{
  int numprocs = Comm->getSize();
  int localproc = Comm->getRank();
  int numlocalrows = 2;
  Tpetra::global_size_t numglobalrows = numprocs*numlocalrows;
  Teuchos::RCP<Tpetra::Map<int> > rowmap = 
    Teuchos::rcp(new Tpetra::Map<int>(numglobalrows, 0, Comm));
  Teuchos::RCP<Tpetra::CrsMatrix<double, int> > matrix = 
    Teuchos::rcp(new Tpetra::CrsMatrix<double, int>(rowmap, numglobalrows));

  Teuchos::Array<int> cols(numglobalrows);
  Teuchos::Array<double> vals(numglobalrows);

  for(size_t j=0; j<numglobalrows; ++j) {
    cols[j] = j;
    vals[j] = 1.0;
  }

  Teuchos::RCP<Tpetra::Map<int> > colmap = 
    Teuchos::rcp(new Tpetra::Map<int>(-1, cols(), 0, Comm));

  for(int i=0; i<numlocalrows; ++i) {
    Teuchos::Array<int> row(1,localproc*numlocalrows+i);
    matrix->insertGlobalValues(
      row[0], row.view(0,1),  vals.view(i,1) );
  }

  matrix->fillComplete();

  typedef Kokkos::DefaultNode::DefaultNodeType DNode;

  Teuchos::RCP<const Tpetra::Map<int> > map_rows = 
    Tpetra::MatrixMatrix<
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
int broadcast_name(Teuchos::RCP<const Teuchos::Comm<Ordinal> > Comm, const char*& name)
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
int read_input_file(Teuchos::RCP<const Teuchos::Comm<Odrinal> > Comm,
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
int run_test(Teuchos::RCP<const Teuchos::Comm<Ordinal> > Comm,
             Teuchos::ParameterList matrixSystem,
             bool result_mtx_to_file,
             bool verbose)
{
  std::string A_file = matrixSystem.get<std::string>("A");
  std::string B_file = matrixSystem.get<std::string>("B");
  std::string C_file = matrixSystem.get<std::string>("C");
  bool AT = matrixSystem.get<bool>("TransA");
  bool BT = matrixSystem.get<bool>("TransB");


  int localProc = Comm->getRank();

  if (localProc == 0 && verbose) {
    std::cout << "Testing C=A"<<AT<<"*B"<<BT<< "; A:" << A_file
              << ", B:" << B_file << ", C:" << C_file << std::endl;
  }

  Teuchos::RCP<Tpetra::CrsMatrix<double,int> > A = Teuchos::null;
  Teuchos::RCP<Tpetra::CrsMatrix<double,int> > B = Teuchos::null;
  Teuchos::RCP<Tpetra::CrsMatrix<double,int> > C = Teuchos::null;
  Teuchos::RCP<Tpetra::CrsMatrix<double,int> > C_check = Teuchos::null;

/*
  Tpetra::Utils::readMatrixMarketMatrix(A_file, Comm, Kokkos::DefaultNode::getDefaultNode(), A);
  Tpetra::Utils::readMatrixMarketMatrix(B_file, Comm, Kokkos::DefaultNode::getDefaultNode(), B);
  Tpetra::Utils::readMatrixMarketMatrix(C_file, Comm, Kokkos::DefaultNode::getDefaultNode(), C_check);*/
  Tpetra::Utils::readHBMatrix(A_file, Comm, Kokkos::DefaultNode::getDefaultNode(), A);
  Tpetra::Utils::readHBMatrix(B_file, Comm, Kokkos::DefaultNode::getDefaultNode(), B);
  Tpetra::Utils::readHBMatrix(C_file, Comm, Kokkos::DefaultNode::getDefaultNode(), C_check);


  /*Teuchos::RCP<Teuchos::basic_FancyOStream<char> > outstream = 
    Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
  B->describe(*outstream, Teuchos::VERB_EXTREME);*/
  Teuchos::RCP<const Tpetra::Map<int> > rowmap = AT ? A->getDomainMap() : A->getRowMap();

  C = Teuchos::rcp( new Tpetra::CrsMatrix<double,int>(rowmap, 1));

  Teuchos::RCP<const Tpetra::CrsMatrix<double,int> > constA = A;
  Teuchos::RCP<const Tpetra::CrsMatrix<double,int> > constB = B;
  typedef Kokkos::DefaultNode::DefaultNodeType DNode;


  Tpetra::MatrixMatrix<
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

  

  Tpetra::MatrixMatrix<
    double, 
    int,
    int,
    DNode,
    typename Kokkos::DefaultKernels<double,int,DNode>::SparseOps>::
  Add(C, false, -1.0, C_check, 1.0);

  double euc_norm = C_check->getEuclideanNorm();

  int return_code = 0;


  if (euc_norm < 1.e-13) {
    if (localProc == 0 && verbose) {
      std::cout << "Test Passed" << std::endl;
      std::cout << "euc_norm " << euc_norm << std::endl;
    }
  }
  else {
    return_code = -1;
    if (localProc == 0) {
      std::cout << "Test Failed, euc_norm = " << euc_norm << std::endl;
    }
  }

  return(return_code);
}
/*
template<class Ordinal>
int read_matrix_file_names(Teuchos::RCP<Teuchos::Comm<Ordinal> > Comm,
                           const char* input_file_name,
                           char*& A_file,
                           bool& transA,
                           char*& B_file,
                           bool& transB,
                           char*& C_file)
{
  int pathlen = path!=0 ? (int)strlen(path) : 0;

  if (Comm->getRank()==0) {
    std::ifstream infile(input_file_name);
    if (!infile) {
      std::cout << "error opening input file " << input_file_name << std::endl;
      return(-1);
    }

    int linelen = 512;
    char line[512];

    infile.getline(line, linelen);
    if (!infile.eof()) {
      if (strchr(line, '#') == NULL) {
        A_file = new char[pathlen+strlen(line)+2];
        sprintf(A_file, "%s/%s",path,line);
      }
    }

    infile.getline(line, linelen);
    if (!infile.eof()) {
      if (!strcmp(line, "TRANSPOSE")) {
        transA = true;
      }
      else transA = false;
    }

    infile.getline(line, linelen);
    if (!infile.eof()) {
      if (strchr(line, '#') == NULL) {
        B_file = new char[pathlen+strlen(line)+2];
        sprintf(B_file, "%s/%s",path,line);
      }
    }

    infile.getline(line, linelen);
    if (!infile.eof()) {
      if (!strcmp(line, "TRANSPOSE")) {
        transB = true;
      }
      else transB = false;
    }

    infile.getline(line, linelen);
    if (!infile.eof()) {
      if (strchr(line, '#') == NULL) {
        C_file = new char[pathlen+strlen(line)+2];
        sprintf(C_file, "%s/%s", path, line);
      }
    }

    broadcast_name(Comm, (const char*&)A_file);
    broadcast_name(Comm, (const char*&)B_file);
    broadcast_name(Comm, (const char*&)C_file);
    int len = transA ? 1 : 0;
	Teuchos::broadcast(Comm, 0, 1, &len);
    len = transB ? 1 : 0;
	Teuchos::broadcast(Comm, 0, 1, &len);
	
  }
  else {
    broadcast_name(Comm, (const char*&)A_file);
    broadcast_name(Comm, (const char*&)B_file);
    broadcast_name(Comm, (const char*&)C_file);
    int len = 0;
	Teuchos::broadcast(Comm, 0, 1, &len);
    transA = len==1 ? true : false;
	Teuchos::broadcast(Comm, 0, 1, &len);
    transB = len==1 ? true : false;
  }

  return(0);
}

int create_maps(Epetra_Comm& Comm,
                const char* input_file_name,
                Epetra_Map*& row_map,
                Epetra_Map*& col_map,
                Epetra_Map*& range_map,
                Epetra_Map*& domain_map)
{
  return( EpetraExt::MatrixMarketFileToBlockMaps(input_file_name,
                                         Comm,
                                         (Epetra_BlockMap*&)row_map,
                                         (Epetra_BlockMap*&)col_map,
                                         (Epetra_BlockMap*&)range_map,
                                         (Epetra_BlockMap*&)domain_map)
  );
}

int read_matrix(const char* filename,
                Epetra_Comm& Comm,
                const Epetra_Map* rowmap,
                Epetra_Map* colmap,
                const Epetra_Map* rangemap,
                const Epetra_Map* domainmap,
                Epetra_CrsMatrix*& mat)
{
  (void)Comm;
  int err = EpetraExt::MatrixMarketFileToCrsMatrix(filename, *rowmap, *colmap,
                                                   *rangemap, *domainmap, mat);

  return(err);
}
*/
template<class Ordinal>
int two_proc_test(Teuchos::RCP<const Teuchos::Comm<Ordinal> > Comm,
  bool verbose)
{
  (void)verbose;
  int thisproc = Comm->getRank();
  int numprocs = Comm->getSize();
  int err =0;

  //only run this test on 2 procs
  if (numprocs != 2) return(0);

  //set up a row-std::map with 2 global elements,
  //1 on each proc.
  int numGlobalRows = 2;
  Teuchos::ArrayRCP<int> myrow(1,3);
  if (thisproc == 1) myrow[0] = 7;
  Teuchos::RCP<const Tpetra::Map<int> > rowmap = Teuchos::rcp(new Tpetra::Map<int>(numGlobalRows, myrow(), 0, Comm));

  //set up a domain-std::map with columns 0 - 4 on proc 0,
  //and columns 5 - 9 on proc 1.
  int numGlobalCols = 10;
  int numMyCols = 5;
  Teuchos::ArrayRCP<int> mycols(numGlobalCols);
  int i;
  for(i=0; i<numGlobalCols; ++i) {
    mycols[i] = i;
  }

  Teuchos::RCP<const Tpetra::Map<int> > domainmap = Teuchos::rcp(new Tpetra::Map<int>(numGlobalCols, mycols(thisproc*numMyCols,numMyCols), 0, Comm));

  //now create matrices A, B and C with rowmap.
  Teuchos::RCP<Tpetra::CrsMatrix<double,int> > A = Teuchos::rcp(new Tpetra::CrsMatrix<double, int>(rowmap, numGlobalCols));
  Teuchos::RCP<Tpetra::CrsMatrix<double,int> > B = Teuchos::rcp(new Tpetra::CrsMatrix<double, int>(rowmap, numGlobalCols));
  Teuchos::RCP<Tpetra::CrsMatrix<double,int> > C = Teuchos::rcp(new Tpetra::CrsMatrix<double, int>(rowmap, numGlobalCols));

  Teuchos::ArrayRCP<double> coefs(numGlobalCols);
  for(i=0; i<numGlobalCols; ++i) {
    coefs[i] = 1.0*i;
  }

  A->insertGlobalValues(myrow[0], mycols(thisproc*numMyCols, numMyCols), coefs(thisproc*numMyCols, numMyCols));

  B->insertGlobalValues(myrow[0], mycols(thisproc*numMyCols, numMyCols), coefs(thisproc*numMyCols, numMyCols));

  A->fillComplete(domainmap, rowmap);
  B->fillComplete(domainmap, rowmap);

  Teuchos::RCP<const Tpetra::CrsMatrix<double,int> > constA = A;
  Teuchos::RCP<const Tpetra::CrsMatrix<double,int> > constB = B;
  typedef Kokkos::DefaultNode::DefaultNodeType DNode;

  Tpetra::MatrixMatrix<
    double, 
    int,
    int,
    DNode,
    typename Kokkos::DefaultKernels<double,int,DNode>::SparseOps>::
  Multiply(constA, false, constB, true, C);

  //std::cout << "two_proc_test, A: "<<std::endl;
  //std::cout << A << std::endl;

  //std::cout << "two_proc_test, B: "<<std::endl;
  //std::cout << B << std::endl;

  //std::cout << "two_proc_test, C: "<<std::endl;
  //std::cout << C << std::endl;

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
Teuchos::RCP<Tpetra::CrsMatrix<LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv> > 
create_crsmatrix(
  Teuchos::RCP<const Teuchos::Comm<CommOrdinal> > comm,
  size_t local_n,
  bool callFillComplete,
  bool symmetric)
{
  int numProcs = comm->getSize();
  Tpetra::global_size_t global_num_rows = numProcs*local_n;

  Teuchos::RCP<Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > rowmap = 
    Teuchos::rcp(global_num_rows, local_n, 0, comm);

  size_t nnz_per_row = 9;
  Teuchos::RCP<Tpetra::CrsMatrix<LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv> > matrix =
    Teuchos::rcp(new Tpetra::CrsMatrix<LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv>(rowmap, nnz_per_row));

  // Add  rows one-at-a-time
  Teuchos::ArrayRCP<Scalar> negOne(1, -Teuchos::ScalarTraits<Scalar>::one());
  Teuchos::ArrayRCP<Scalar> posTwo(1, 2*Teuchos::ScalarTraits<Scalar>::one());
  Teuchos::ArrayRCP<Scalar> val_L(1, symmetric ? negOne : 0.5*negOne);

  GlobalOrdinal GO1 = Teuchos::OrdinalTraits<GlobalOrdinal>::one();
  GlobalOrdinal GO0 = Teuchos::OrdinalTraits<GlobalOrdinal>::zero();

  for (int i=0; i<local_n; i++) {
    Teuchos::ArrayRCP<GlobalOrdinal> GlobalRow(1, matrix->getRowMap()->getGlobalElement(i));
    Teuchos::ArrayRCP<GlobalOrdinal> RowLess1(1,GlobalRow - GO1);
    Teuchos::ArrayRCP<GlobalOrdinal> RowPlus1(1,GlobalRow + GO1);
    Teuchos::ArrayRCP<GlobalOrdinal> RowLess5(1,GlobalRow - (5*GO1));
    Teuchos::ArrayRCP<GlobalOrdinal> RowPlus5(1,GlobalRow + (5*GO1));
    Teuchos::ArrayRCP<GlobalOrdinal> RowLess9(1,GlobalRow - (9*GO1));
    Teuchos::ArrayRCP<GlobalOrdinal> RowPlus9(1,GlobalRow + (9*GO1));
    Teuchos::ArrayRCP<GlobalOrdinal> RowLess24(1,GlobalRow - (24*GO1));
    Teuchos::ArrayRCP<GlobalOrdinal> RowPlus24(1,GlobalRow + (24*GO1));
    Teuchos::ArrayRCP<GlobalOrdinal> RowLess48(1,GlobalRow - (48*GO1));
    Teuchos::ArrayRCP<GlobalOrdinal> RowPlus48(1,GlobalRow + (48*GO1));

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


