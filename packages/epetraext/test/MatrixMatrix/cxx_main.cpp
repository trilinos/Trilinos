
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include <Epetra_ConfigDefs.h>

#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#endif

#include <Epetra_SerialComm.h>
#include <Epetra_Time.h>
#include <Epetra_Import.h>
#include <Epetra_Map.h>
#include <Epetra_LocalMap.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <EpetraExt_MatrixMatrix.h>

#include <EpetraExt_BlockMapIn.h>
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_RowMatrixOut.h>

namespace EpetraExt {
extern
Epetra_Map* find_rows_containing_cols(const Epetra_CrsMatrix& M,
                                      const Epetra_Map* colmap);
}

int read_input_file(Epetra_Comm& Comm,
                    const char* input_file_name,
                    const char**& filenames,
                    int& numfiles,
                    int& numfilenames_allocated);

int read_matrix_file_names(Epetra_Comm& Comm,
                           const char* input_file_name,
                           char*& A_file,
                           bool& transA,
                           char*& B_file,
                           bool& transB,
                           char*& C_file);

int broadcast_name(Epetra_Comm& Comm, const char*& name);

int create_maps(Epetra_Comm& Comm,
                const char* input_file_name,
                Epetra_Map*& row_map,
                Epetra_Map*& col_map,
                Epetra_Map*& range_map,
                Epetra_Map*& domain_map);

int read_matrix(const char* filename,
                Epetra_Comm& Comm,
                const Epetra_Map* rowmap,
                Epetra_Map* colmap,
                const Epetra_Map* rangemap,
                const Epetra_Map* domainmap,
                Epetra_CrsMatrix*& mat);

int run_test(Epetra_Comm& Comm, const char* filename,
             bool result_mtx_to_file=false,
             bool verbose=false);

int test_find_rows(Epetra_Comm& Comm);

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << "-i file_name" << std::endl;
    std::cout << "  (where file_name contains a list of input-files)"<<std::endl;
    return(-1);
  }

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  bool write_result_mtx = false;
  bool verbose = false;
  int write = 0;
  char* input_file = NULL;
  bool input_file_specified = false;

  if (Comm.MyPID()==0) {
    for(int ii=0; ii<argc; ++ii) {
      if (!strcmp("-write_result", argv[ii])) write_result_mtx = true;
      if (!strcmp("-v", argv[ii])) verbose = true;
      if (!strcmp("-i", argv[ii])) {
        input_file = argv[ii+1];
        input_file_specified = true;
      }
    }
    write = write_result_mtx ? 1 : 0;
  }
#ifdef EPETRA_MPI
  MPI_Bcast(&write, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (write) write_result_mtx = true;
#endif

  if (!input_file_specified) {
    input_file = new char[16];
    sprintf(input_file, "./infiles");
  }

  const char** filenames = NULL;
  int numfiles = 0;
  int numfilenames_allocated = 0;

  int err = read_input_file(Comm, input_file,
                            filenames, numfiles, numfilenames_allocated);
  if (err != 0) {
    std::cout << "read_input_file returned err=="<<err<<std::endl;
    return(err);
  }

  err = test_find_rows(Comm);
  if (err != 0) {
    std::cout << "test_find_rows returned err=="<<err<<endl;
    return(err);
  }

  for(int i=0; i<numfiles; ++i) {
    err = run_test(Comm, filenames[i], write_result_mtx, verbose);
    delete [] filenames[i];
    if (err != 0) break;
  }

  for(int j=numfiles; j<numfilenames_allocated; ++j) {
    delete [] filenames[j];
  }

  delete [] filenames;

  if (!input_file_specified) delete [] input_file;

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return(err);
}

int test_find_rows(Epetra_Comm& Comm)
{
  int numprocs = Comm.NumProc();
  int localproc = Comm.MyPID();
  int numlocalrows = 2;
  int numglobalrows = numprocs*numlocalrows;
  Epetra_Map rowmap(numlocalrows*numprocs, 0, Comm);
  Epetra_CrsMatrix matrix(Copy, rowmap, numglobalrows);

  int err = 0;
  int* cols = new int[numglobalrows];
  double*vals = new double[numglobalrows];

  for(int j=0; j<numglobalrows; ++j) {
    cols[j] = j;
    vals[j] = 1.0;
  }

  for(int i=0; i<numlocalrows; ++i) {
    int row = localproc*numlocalrows+i;
    err = matrix.InsertGlobalValues(row, numglobalrows, vals, cols);
    if (err != 0) {
      return(err);
    }
  }

  err = matrix.FillComplete();
  if (err != 0) {
    return(err);
  }

  Epetra_Map* map_rows = EpetraExt::find_rows_containing_cols(matrix,
							      &(matrix.ColMap()));

  if (map_rows->NumMyElements() != numglobalrows) {
    return(-1);
  }

  delete map_rows;
  delete [] cols;
  delete [] vals;

  return(0);
}

int expand_name_list(const char* newname,
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
  }

  names[offset] = newname;
  ++num_names;
  return(0);
}

int broadcast_name(Epetra_Comm& Comm, const char*& name)
{
  if (Comm.NumProc() < 2) return(0);

#ifdef EPETRA_MPI
  int len;
  int localProc = Comm.MyPID();
  if (localProc == 0) {
    len = strlen(name)+1;
    
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast((void*)name, len, MPI_CHAR, 0, MPI_COMM_WORLD);
  }
  else {
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    name = new char[len];
    MPI_Bcast((void*)name, len, MPI_CHAR, 0, MPI_COMM_WORLD);
  }

#endif
  return(0);
}

int read_input_file(Epetra_Comm& Comm,
                    const char* input_file_name,
                    const char**& filenames,
                    int& numfiles,
                    int& numfilenames_allocated)
{
  if (Comm.MyPID() == 0) {
    ifstream infile(input_file_name);
    if (!infile) {
      std::cout << "ERROR opening file "<<input_file_name << std::endl;
      return(-1);
    }

    int linelen = 256;
    char* line = NULL;

    while(!infile.eof()) {
      line = new char[linelen];
      infile.getline(line, linelen);
      if (infile.fail()) {
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

#ifdef EPETRA_MPI
    MPI_Bcast(&numfiles, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    for(int i=0; i<numfiles; ++i) {
      broadcast_name(Comm, filenames[i]);
    }
  }
  else {
#ifdef EPETRA_MPI
    MPI_Bcast(&numfiles, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    filenames = new const char*[numfiles];
    for(int i=0; i<numfiles; ++i) {
      broadcast_name(Comm, filenames[i]);
    }
  }
  
  return(0);
}

int run_test(Epetra_Comm& Comm, const char* filename,
             bool result_mtx_to_file,
             bool verbose)
{
  char* A_file = NULL;
  char AT[3]; AT[0] = '^'; AT[1] = 'T'; AT[2] = '\0';
  char* B_file = NULL;
  char BT[3]; BT[0] = '^'; BT[1] = 'T'; BT[2] = '\0';
  char* C_file = NULL;
  bool transA, transB;

  int err = read_matrix_file_names(Comm, filename, A_file, transA,
                                   B_file, transB, C_file);
  if (err != 0) {
    std::cout << "Error, read_matrix_file_names returned " << err << std::endl;
    return(err);
  }

  if (!transA) AT[0] = '\0';
  if (!transB) BT[0] = '\0';

  int localProc = Comm.MyPID();

  if (localProc == 0 && verbose) {
    std::cout << "Testing C=A"<<AT<<"*B"<<BT<< "; A:" << A_file
              << ", B:" << B_file << ", C:" << C_file << std::endl;
  }

  Epetra_CrsMatrix* A = NULL;
  Epetra_CrsMatrix* B = NULL;
  Epetra_CrsMatrix* C = NULL;
  Epetra_CrsMatrix* C_check = NULL;

  Epetra_Map* A_row_map = NULL;
  Epetra_Map* A_col_map = NULL;
  Epetra_Map* A_range_map = NULL;
  Epetra_Map* A_domain_map = NULL;
  err = create_maps(Comm, A_file, A_row_map, A_col_map, A_range_map, A_domain_map);
  if (err != 0) {
    std::cout << "create_maps A returned err=="<<err<<std::endl;
    return(err);
  }

  Epetra_Map* B_row_map = NULL;
  Epetra_Map* B_col_map = NULL;
  Epetra_Map* B_range_map = NULL;
  Epetra_Map* B_domain_map = NULL;
  err = create_maps(Comm, B_file, B_row_map, B_col_map, B_range_map, B_domain_map);
  if (err != 0) {
    std::cout << "create_maps A returned err=="<<err<<std::endl;
    return(err);
  }

  err = read_matrix(A_file, Comm, A_row_map, A_col_map,
                    A_range_map, A_domain_map, A);
  delete [] A_file;
  if (err != 0) {
    std::cout << "read_matrix A returned err=="<<err<<std::endl;
    return(err);
  }

  err = read_matrix(B_file, Comm, B_row_map, B_col_map,
                    B_range_map, B_domain_map, B);
  delete [] B_file;
  if (err != 0) {
    std::cout << "read_matrix B returned err=="<<err<<std::endl;
    return(-1);
  }

  const Epetra_Map* rowmap = transA ? &(A->DomainMap()) : &(A->RowMap());
  const Epetra_Map* Cdomainmap = transB ? &(B->RangeMap()) : B_domain_map;

  C = new Epetra_CrsMatrix(Copy, *rowmap, 1);

  err = EpetraExt::MatrixMatrix::Multiply(*A, transA, *B, transB, *C);
  if (err != 0) {
    std::cout << "err "<<err<<" from MatrixMatrix::Multiply"<<std::endl;
    return(err);
  }

  if (result_mtx_to_file) {
    EpetraExt::RowMatrixToMatrixMarketFile("result.mtx", *C);
  }

  Epetra_Map* Cck_row_map = NULL;
  Epetra_Map* Cck_col_map = NULL;
  Epetra_Map* Cck_range_map = NULL;
  Epetra_Map* Cck_domain_map = NULL;
  err = create_maps(Comm, C_file, Cck_row_map, Cck_col_map,
                    Cck_range_map, Cck_domain_map);
  if (err != 0) {
    std::cout << "create_maps C returned err=="<<err<<std::endl;
    return(err);
  }

  err = read_matrix(C_file, Comm, Cck_row_map, Cck_col_map,
                     Cck_range_map, Cck_domain_map, C_check);
  delete [] C_file;
  if (err != 0) {
    std::cout << "read_matrix C returned err=="<<err<<std::endl;
    return(-1);
  }

  EpetraExt::MatrixMatrix::Add(*C, false, -1.0, *C_check, 1.0);

  double inf_norm = C_check->NormInf();

  int return_code = 0;

  if (inf_norm < 1.e-13) {
    if (localProc == 0 && verbose) {
      std::cout << "Test Passed" << std::endl;
    }
  }
  else {
    return_code = -1;
    if (localProc == 0) {
      std::cout << "Test Failed, inf_norm = " << inf_norm << std::endl;
    }
  }

  delete A;
  delete B;
  delete C;
  delete C_check;

  delete A_row_map;
  delete A_col_map;
  delete A_range_map;
  delete A_domain_map;

  delete B_row_map;
  delete B_col_map;
  delete B_range_map;
  delete B_domain_map;

  delete Cck_row_map;
  delete Cck_col_map;
  delete Cck_range_map;
  delete Cck_domain_map;

  return(return_code);
}

int read_matrix_file_names(Epetra_Comm& Comm,
                           const char* input_file_name,
                           char*& A_file,
                           bool& transA,
                           char*& B_file,
                           bool& transB,
                           char*& C_file)
{
  if (Comm.MyPID()==0) {
    ifstream infile(input_file_name);
    if (!infile) {
      std::cout << "error opening input file " << input_file_name << std::endl;
      return(-1);
    }

    char line[256];

    infile.getline(line, 256);
    if (!infile.eof()) {
      if (strchr(line, '#') == NULL) {
        A_file = new char[strlen(line)+1];
        sprintf(A_file, line);
      }
    }

    infile.getline(line, 256);
    if (!infile.eof()) {
      if (!strcmp(line, "TRANSPOSE")) {
        transA = true;
      }
      else transA = false;
    }

    infile.getline(line, 256);
    if (!infile.eof()) {
      if (strchr(line, '#') == NULL) {
        B_file = new char[strlen(line)+1];
        sprintf(B_file, line);
      }
    }

    infile.getline(line, 256);
    if (!infile.eof()) {
      if (!strcmp(line, "TRANSPOSE")) {
        transB = true;
      }
      else transB = false;
    }

    infile.getline(line, 256);
    if (!infile.eof()) {
      if (strchr(line, '#') == NULL) {
        C_file = new char[strlen(line)+1];
        sprintf(C_file, line);
      }
    }

    broadcast_name(Comm, (const char*&)A_file);
    broadcast_name(Comm, (const char*&)B_file);
    broadcast_name(Comm, (const char*&)C_file);
    int len = transA ? 1 : 0;
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    len = transB ? 1 : 0;
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  else {
    broadcast_name(Comm, (const char*&)A_file);
    broadcast_name(Comm, (const char*&)B_file);
    broadcast_name(Comm, (const char*&)C_file);
    int len = 0;
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    transA = len==1 ? true : false;
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
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
  int err = EpetraExt::MatrixMarketFileToCrsMatrix(filename, *rowmap, *colmap,
                                                   *rangemap, *domainmap, mat);

  return(err);
}

