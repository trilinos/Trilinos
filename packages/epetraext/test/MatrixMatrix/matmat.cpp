
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#include <mpi.h>
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
#include <EpetraExt_BlockMapOut.h>
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_RowMatrixOut.h>

int read_input_file(Epetra_Comm& Comm,
                    const char* input_file_name,
                    char*& A_file,
                    bool& transA,
                    char*& B_file,
                    bool& transB,
                    char*& C_file);

int create_maps(Epetra_Comm& Comm,
                const char* input_file_name,
                Epetra_Map*& row_map,
                Epetra_Map*& col_map,
                Epetra_Map*& domain_map);

int read_matrix(const char* filename,
                Epetra_Comm& Comm,
                const Epetra_Map* rowmap,
                Epetra_Map* colmap,
                const Epetra_Map* domainmap,
                Epetra_CrsMatrix*& mat);

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " input_file_name" << std::endl;
    return(-1);
  }

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

#else
  Epetra_SerialComm Comm;
#endif

  int localProc = Comm.MyPID();

  char* A_file = NULL;
  char* B_file = NULL;
  char* C_file = NULL;
  bool transA, transB;

  int err = read_input_file(Comm, argv[1], A_file, transA, B_file, transB, C_file);
  if (err != 0) {
    std::cout << "read_input_file returned err=="<<err<<std::endl;
    return(err);
  }

  if (localProc == 0) {
    std::cout << "A_file: " << A_file << ", transA: " << transA << std::endl;
    std::cout << "B_file: " << B_file << ", transB: " << transB << std::endl;
    std::cout << "C_file: " << C_file << std::endl;
  }

  Epetra_CrsMatrix* A = NULL;
  Epetra_CrsMatrix* B = NULL;
  Epetra_CrsMatrix* C = NULL;
  Epetra_CrsMatrix* C_check = NULL;

  Epetra_Map* A_row_map = NULL;
  Epetra_Map* A_col_map = NULL;
  Epetra_Map* A_domain_map = NULL;
  err = create_maps(Comm, A_file, A_row_map, A_col_map, A_domain_map);
  if (err != 0) {
    std::cout << "create_maps A returned err=="<<err<<std::endl;
    return(err);
  }

  Epetra_Map* B_row_map = NULL;
  Epetra_Map* B_col_map = NULL;
  Epetra_Map* B_domain_map = NULL;
  err = create_maps(Comm, B_file, B_row_map, B_col_map, B_domain_map);
  if (err != 0) {
    std::cout << "create_maps A returned err=="<<err<<std::endl;
    return(err);
  }

  err = read_matrix(A_file, Comm, A_row_map, A_col_map, A_domain_map, A);
  if (err != 0) {
    std::cout << "read_matrix A returned err=="<<err<<std::endl;
    return(err);
  }

  //std::cout << "******** A *********" << std::endl << *A << std::endl;

  err = read_matrix(B_file, Comm, B_row_map, B_col_map, B_domain_map, B);
  if (err != 0) {
    std::cout << "read_matrix B returned err=="<<err<<std::endl;
    return(-1);
  }

  //std::cout << "******** B *********" << std::endl << *B << std::endl;

  const Epetra_Map* rowmap = transA ? &(A->DomainMap()) : &(A->RowMap());
  const Epetra_Map* Cdomainmap = transB ? &(B->RangeMap()) : B_domain_map;

  C = new Epetra_CrsMatrix(Copy, *rowmap, 1);

  err = EpetraExt::MatrixMatrix::Multiply(*A, transA, *B, transB, *C);
  if (err != 0) {
    std::cout << "err "<<err<<" from MatrixMatrix::Multiply"<<std::endl;
    return(err);
  }

  std::cout << std::endl << "********* C *********" <<std::endl<< *C << std::endl;

  err = read_matrix(C_file, Comm, rowmap, NULL, Cdomainmap, C_check);
  if (err != 0) {
    std::cout << "read_matrix C returned err=="<<err<<std::endl;
    return(-1);
  }

  std::cout << std::endl << "********* C_check *********" <<std::endl<< *C_check << std::endl;

  EpetraExt::MatrixMatrix::Add(*C, false, -1.0, *C_check, 1.0);

  double inf_norm = C_check->NormInf();

  if (inf_norm == 0.0) {
    if (localProc == 0) {
      std::cout << "Test Passed" << std::endl;
    }
  }
  else {
    if (localProc == 0) {
      std::cout << "Test Failed" << std::endl;
    }
  }

  delete [] A_file;
  delete [] B_file;

  delete A;
  delete B;
  delete C;
  delete C_check;

  delete A_row_map;
  delete A_col_map;
  delete A_domain_map;
  delete B_row_map;
  delete B_col_map;
  delete B_domain_map;

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return(0);
}

int read_input_file(Epetra_Comm& Comm,
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
  }

  if (Comm.NumProc() > 1) {
    if (Comm.MyPID() == 0) {
      int len = strlen(A_file)+1;
#ifdef EPETRA_MPI
      MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(A_file, len, MPI_CHAR, 0, MPI_COMM_WORLD);
      len = strlen(B_file)+1;
      MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(B_file, len, MPI_CHAR, 0, MPI_COMM_WORLD);
      len = strlen(C_file)+1;
      MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(C_file, len, MPI_CHAR, 0, MPI_COMM_WORLD);
      len = transA ? 1 : 0;
      MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
      len = transB ? 1 : 0;
      MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    }
    else {
      int len = 0;
#ifdef EPETRA_MPI
      MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
      A_file = new char[len];
      MPI_Bcast(A_file, len, MPI_CHAR, 0, MPI_COMM_WORLD);
      MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
      B_file = new char[len];
      MPI_Bcast(B_file, len, MPI_CHAR, 0, MPI_COMM_WORLD);
      MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
      C_file = new char[len];
      MPI_Bcast(C_file, len, MPI_CHAR, 0, MPI_COMM_WORLD);
      MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
      transA = len==1 ? true : false;
      MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
      transB = len==1 ? true : false;
#endif
    }
  }

  return(0);
}

int create_maps(Epetra_Comm& Comm,
                const char* input_file_name,
                Epetra_Map*& row_map,
                Epetra_Map*& col_map,
                Epetra_Map*& domain_map)
{
  row_map = NULL;
  col_map = NULL;
  domain_map = NULL;

  int numRows, numCols;
  ifstream infile(input_file_name);
  if (!infile) {
    std::cout << "error opening file " << input_file_name << std::endl;
    return(-1);
  }

  infile >> numRows;
  infile >> numCols;

  row_map = new Epetra_Map(numRows, 0, Comm);
  int num_map_cols = 0, insertPoint, foundOffset;
  int allocLen = numCols;
  int* map_cols = new int[allocLen];

  while(!infile.eof()) {
    int row, col;
    double val;
    infile >> row;
    infile >> col;
    infile >> val;
    if (!infile.eof()) {
      if (row_map->MyGID(row)) {
        foundOffset = Epetra_Util_binary_search(col, map_cols, num_map_cols,
                                                insertPoint);
        if (foundOffset < 0) {
          Epetra_Util_insert(col, insertPoint, map_cols,
                             num_map_cols, allocLen);
        }
      }
    }
  }
 
  col_map = new Epetra_Map(-1, num_map_cols, map_cols, 0, Comm);
  domain_map = new Epetra_Map(numCols, 0, Comm);

  delete [] map_cols;
  return(0);
}

int read_matrix(const char* filename,
                Epetra_Comm& Comm,
                const Epetra_Map* rowmap,
                Epetra_Map* colmap,
                const Epetra_Map* domainmap,
                Epetra_CrsMatrix*& mat)
{
  int numRows, numCols;
  ifstream infile(filename);
  if (!infile) {
    std::cout << "error opening file " << filename << std::endl;
    return(-1);
  }

  infile >> numRows;
  infile >> numCols;

  if (colmap == NULL) {
    mat = new Epetra_CrsMatrix(Copy, *rowmap, 1);
  }
  else {
    mat = new Epetra_CrsMatrix(Copy, *rowmap, *colmap, 1);
  }

  while(!infile.eof()) {
    int row, col;
    double val;
    infile >> row;
    infile >> col;
    infile >> val;
    if (!infile.eof()) {
      if (rowmap->MyGID(row)) {
        mat->InsertGlobalValues(row, 1, &val, &col);
      }
    }
  }

  mat->FillComplete(*domainmap, *rowmap);

  return(0);
}

