/*
//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/


#include <string>
#include <vector>
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
#include <Epetra_Export.h>
#include <Epetra_Map.h>
#include <Epetra_LocalMap.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <EpetraExt_MMHelpers.h>
#include <EpetraExt_MatrixMatrix.h>

#include <EpetraExt_BlockMapIn.h>
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_RowMatrixOut.h>

namespace EpetraExt {
extern
Epetra_Map* find_rows_containing_cols(const Epetra_CrsMatrix& M,
                                      const Epetra_Map& column_map);
}

int read_input_file(Epetra_Comm& Comm,
                    const std::string& input_file_name,
                    std::vector<std::string>& filenames);

int read_matrix_file_names(Epetra_Comm& Comm,
                           const std::string& input_file_name,
                           std::string& A_file,
                           bool& transA,
                           std::string& B_file,
                           bool& transB,
                           std::string& C_file);

int broadcast_name(Epetra_Comm& Comm, std::string& name);

int create_maps(Epetra_Comm& Comm,
                const std::string& input_file_name,
                Epetra_Map*& row_map,
                Epetra_Map*& col_map,
                Epetra_Map*& range_map,
                Epetra_Map*& domain_map);

int read_matrix(const std::string& filename,
                Epetra_Comm& Comm,
                const Epetra_Map* rowmap,
                Epetra_Map* colmap,
                const Epetra_Map* rangemap,
                const Epetra_Map* domainmap,
                Epetra_CrsMatrix*& mat);

int run_test(Epetra_Comm& Comm,
             const std::string& filename,
	     bool do_FillComplete,
             bool result_mtx_to_file=false,
             bool verbose=false);

int two_proc_test(Epetra_Comm& Comm,
                  bool verbose=false);

int test_find_rows(Epetra_Comm& Comm);

Epetra_CrsMatrix* create_epetra_crsmatrix(int numProcs,
                                          int localProc,
                                          int local_n,
                                          bool callFillComplete = true,
                                          bool symmetric = true);

int time_matrix_matrix_multiply(Epetra_Comm& Comm,
                                bool verbose);

int test_drumm1(Epetra_Comm& Comm);

/////////////////////////////////////
//Global variable!!!!
std::string path;
////////////////////////////////////

int main(int argc, char** argv) {

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  bool write_result_mtx = false;
  bool verbose = false;
  int write = 0, inp_specified = 0;
  bool path_specified = false;
  std::string input_file;
  bool input_file_specified = false;

  if (Comm.MyPID()==0) {
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
    inp_specified = input_file_specified ? 1 : 0;
  }
#ifdef EPETRA_MPI
  MPI_Bcast(&write, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (write) write_result_mtx = true;
  MPI_Bcast(&inp_specified, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (inp_specified) input_file_specified = true;
#endif

  if (!path_specified) {
    path = ".";
  }

  int err = two_proc_test(Comm, verbose);
  if (err != 0) {
    std::cout << "two_proc_test returned err=="<<err<<std::endl;
    return(err);
  }

  std::vector<std::string> filenames;

  if (input_file_specified) {
    filenames.push_back(std::string(input_file));
  }
  else {
    input_file = "infiles";
    std::cout << "specific input-file not specified, getting list of input-files from file 'infiles'." << std::endl;

    err = read_input_file(Comm, input_file, filenames);
    if (err != 0) {
      if (path_specified) path_specified = false;
      path = "./MatrixMatrix";
      read_input_file(Comm, input_file, filenames);
    }
  }

  err = test_find_rows(Comm);
  if (err != 0) {
    std::cout << "test_find_rows returned err=="<<err<<std::endl;
    return(err);
  }

  err = test_drumm1(Comm);
  if (err != 0) {
    std::cout << "test_drumm1 returned err=="<<err<<std::endl;
    return(err);
  }

  for(size_t i=0; i<filenames.size(); ++i) {
    err = run_test(Comm, filenames[i], true, write_result_mtx, verbose);
    if (err != 0) break;
    err = run_test(Comm, filenames[i], false, write_result_mtx, verbose);
    if (err != 0) break;
  }

  Epetra_CrsMatrix* D = create_epetra_crsmatrix(Comm.NumProc(),
                                                Comm.MyPID(), 10,
                                                true, false);

//  std::cout << "D: \n"  << *D << std::endl;

  EpetraExt::MatrixMatrix::Add(*D, true, 0.5, *D, 0.5);

//  std::cout << "symm D: \n"  << *D << std::endl;

  delete D;

  if (err == 0) {
    err = time_matrix_matrix_multiply(Comm, verbose);
  }

  int global_err = err;

#ifdef EPETRA_MPI
  MPI_Allreduce(&err, &global_err, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Finalize();
#endif

  return(global_err);
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

  Epetra_Map colmap(-1, numglobalrows, cols, 0, Comm);

  for(int i=0; i<numlocalrows; ++i) {
    int row = localproc*numlocalrows+i;
    err = matrix.InsertGlobalValues(row, 1, &(vals[i]), &row);
    if (err != 0) {
      return(err);
    }
  }

  err = matrix.FillComplete();
  if (err != 0) {
    return(err);
  }

  Epetra_Map* map_rows = EpetraExt::find_rows_containing_cols(matrix, colmap);

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
    for(int i=offset; i<alloc_len; ++i) {
      names[i] = NULL;
    }
  }

  names[offset] = newname;
  ++num_names;
  return(0);
}

int broadcast_name(Epetra_Comm& Comm, std::string& name)
{
  if (Comm.NumProc() < 2) return(0);

#ifdef EPETRA_MPI
  int len = name.size();
  int localProc = Comm.MyPID();
  if (localProc == 0) {
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast((void*)name.c_str(), len+1, MPI_CHAR, 0, MPI_COMM_WORLD);
  }
  else {
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    name.resize(len+1, ' ');
    MPI_Bcast((void*)name.c_str(), len+1, MPI_CHAR, 0, MPI_COMM_WORLD);
  }

#endif
  return(0);
}

int read_input_file(Epetra_Comm& Comm,
                    const std::string& input_file_name,
                    std::vector<std::string>& filenames)
{
  int local_err = 0, global_err = 0;
  std::ifstream* infile = NULL;
  int pathlen = path.size();

  if (Comm.MyPID() == 0) {
    std::string full_name = path.empty() ? input_file_name : path+"/"+input_file_name;

    infile = new std::ifstream(full_name.c_str());
    if (!(*infile)) {
      local_err = -1;
      delete infile;
    }
  }

  Comm.SumAll(&local_err, &global_err, 1);

  if (global_err != 0) {
    return(global_err);
  }


  if (Comm.MyPID() == 0) {
    const int linelen = 512;
    char line[linelen];

    std::ifstream& ifile = *infile;
    while(!ifile.eof()) {
      if (pathlen>0) {
        sprintf(line,"%s/",path.c_str());
        ifile.getline(&(line[pathlen+1]), linelen);
      }
      else {
        ifile.getline(line, linelen);
      }

      if (ifile.fail()) {
        break;
      }
      if (strchr(line, '#') == NULL) {
        filenames.push_back(std::string(line));
      }
    }

    int numfiles = filenames.size();
#ifdef EPETRA_MPI
    MPI_Bcast(&numfiles, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    for(int i=0; i<numfiles; ++i) {
      broadcast_name(Comm, filenames[i]);
    }

    delete infile;
  }
  else {
    int numfiles = 0;
#ifdef EPETRA_MPI
    MPI_Bcast(&numfiles, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    filenames.resize(numfiles);
    for(int i=0; i<numfiles; ++i) {
      broadcast_name(Comm, filenames[i]);
    }
  }
  
  return(0);
}

int run_test(Epetra_Comm& Comm,
             const std::string& filename,
	     bool do_FillComplete,
             bool result_mtx_to_file,
             bool verbose)
{
  std::string A_file;
  char AT[3]; AT[0] = '^'; AT[1] = 'T'; AT[2] = '\0';
  std::string B_file;
  char BT[3]; BT[0] = '^'; BT[1] = 'T'; BT[2] = '\0';
  std::string C_file;
  bool transA, transB;

  if(!Comm.MyPID()) cout<<"Testing: "<<filename<<endl;

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
  if (err != 0) {
    std::cout << "read_matrix A returned err=="<<err<<std::endl;
    return(err);
  }

  err = read_matrix(B_file, Comm, B_row_map, B_col_map,
                    B_range_map, B_domain_map, B);
  if (err != 0) {
    std::cout << "read_matrix B returned err=="<<err<<std::endl;
    return(-1);
  }

  const Epetra_Map* rowmap    = transA ? &(A->DomainMap()) : &(A->RowMap());
  const Epetra_Map* domainMap = transB ? &(B->RangeMap()) : &(B->DomainMap());
  const Epetra_Map* rangeMap  = transA ? &(A->DomainMap()) : &(A->RangeMap());


  C = new Epetra_CrsMatrix(Copy, *rowmap, 1);

  if(C->Comm().MyPID()) printf("transA = %d transB = %d\n",(int)transA,(int)transB);

  err = EpetraExt::MatrixMatrix::Multiply(*A, transA, *B, transB, *C,do_FillComplete);
  if (err != 0) {
    std::cout << "err "<<err<<" from MatrixMatrix::Multiply"<<std::endl;
    return(err);
  }

  if(!do_FillComplete) C->FillComplete(*domainMap,*rangeMap);


//  std::cout << "A: " << *A << std::endl << "B: "<<*B<<std::endl<<"C: "<<*C<<std::endl;
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
  if (err != 0) {
    std::cout << "read_matrix C returned err=="<<err<<std::endl;
    return(-1);
  }

  //now subtract our check-matrix from our result matrix (C) and
  //the infinity-norm of that should be zero.
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
      std::cout << "Test Failed ("<<filename<<"), inf_norm = " << inf_norm << std::endl;
    }
Comm.Barrier();
std::cout << "C"<<std::endl;
std::cout << *C<<std::endl;
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
                           const std::string& input_file_name,
                           std::string& A_file,
                           bool& transA,
                           std::string& B_file,
                           bool& transB,
                           std::string& C_file)
{
  if (Comm.MyPID()==0) {
    std::ifstream infile(input_file_name.c_str());
    if (!infile) {
      std::cout << "error opening input file " << input_file_name << std::endl;
      return(-1);
    }

    const int linelen = 512;
    char line[linelen];

    infile.getline(line, linelen);
    if (!infile.eof()) {
      if (strchr(line, '#') == NULL) {
        A_file = path+"/"+line;
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
        B_file = path+"/"+line;
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
        C_file = path+"/"+line;
      }
    }

    broadcast_name(Comm, A_file);
    broadcast_name(Comm, B_file);
    broadcast_name(Comm, C_file);
#ifdef EPETRA_MPI
    int len = transA ? 1 : 0;
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    len = transB ? 1 : 0;
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  }
  else {
    broadcast_name(Comm, A_file);
    broadcast_name(Comm, B_file);
    broadcast_name(Comm, C_file);
#ifdef EPETRA_MPI
    int len = 0;
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    transA = len==1 ? true : false;
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    transB = len==1 ? true : false;
#endif
  }

  return(0);
}

int create_maps(Epetra_Comm& Comm,
                const std::string& input_file_name,
                Epetra_Map*& row_map,
                Epetra_Map*& col_map,
                Epetra_Map*& range_map,
                Epetra_Map*& domain_map)
{
  return( EpetraExt::MatrixMarketFileToBlockMaps(input_file_name.c_str(),
                                         Comm,
                                         (Epetra_BlockMap*&)row_map,
                                         (Epetra_BlockMap*&)col_map,
                                         (Epetra_BlockMap*&)range_map,
                                         (Epetra_BlockMap*&)domain_map)
  );
}

int read_matrix(const std::string& filename,
                Epetra_Comm& Comm,
                const Epetra_Map* rowmap,
                Epetra_Map* colmap,
                const Epetra_Map* rangemap,
                const Epetra_Map* domainmap,
                Epetra_CrsMatrix*& mat)
{
  (void)Comm;
  int err = EpetraExt::MatrixMarketFileToCrsMatrix(filename.c_str(), *rowmap,
                                                   *rangemap, *domainmap, mat);

  return(err);
}

int two_proc_test(Epetra_Comm& Comm,
                  bool verbose)
{
  (void)verbose;
  int thisproc = Comm.MyPID();
  int numprocs = Comm.NumProc();

  //only run this test on 2 procs
  if (numprocs != 2) return(0);

  //set up a row-std::map with 2 global elements,
  //1 on each proc.
  int numGlobalRows = 2;
  int numMyRows = 1;
  int myrow = 3;
  if (thisproc == 1) myrow = 7;
  Epetra_Map rowmap(numGlobalRows, numMyRows, &myrow, 0, Comm);

  //set up a domain-map with columns 0 - 4 on proc 0,
  //and columns 5 - 9 on proc 1.
  int numGlobalCols = 10;
  int numMyCols = 5;
  int* mycols = new int[numGlobalCols];
  int i;
  for(i=0; i<numGlobalCols; ++i) {
    mycols[i] = i;
  }

  Epetra_Map domainmap(numGlobalCols, numMyCols, &(mycols[thisproc*numMyCols]),
                       0, Comm);

  //now create matrices A, B and C with rowmap.
  Epetra_CrsMatrix A(Copy, rowmap, 10);
  Epetra_CrsMatrix B(Copy, rowmap, 10);
  Epetra_CrsMatrix C(Copy, rowmap, 10);

  double* coefs = new double[numGlobalCols];
  for(i=0; i<numGlobalCols; ++i) {
    coefs[i] = 1.0*i;
  }

  int numCols = numGlobalCols - 2;
  int offset = 0;
  if (thisproc == 1) offset = 2;

  int err = A.InsertGlobalValues(myrow, numCols, &coefs[offset], &mycols[offset]);

  err += B.InsertGlobalValues(myrow, numMyCols, &(coefs[thisproc*numMyCols]),
                       &(mycols[thisproc*numMyCols]));

  err += A.FillComplete(domainmap, rowmap);
  err += B.FillComplete(domainmap, rowmap);

  err += EpetraExt::MatrixMatrix::Multiply(A, false, B, true, C);

  //std::cout << "two_proc_test, A: "<<std::endl;
  //std::cout << A << std::endl;

  //std::cout << "two_proc_test, B: "<<std::endl;
  //std::cout << B << std::endl;

  //std::cout << "two_proc_test, C: "<<std::endl;
  //std::cout << C << std::endl;

  if (C.NumGlobalNonzeros() != 4) {
    err += 1;
  }

  delete [] coefs;
  delete [] mycols;

  return(err);
}

int time_matrix_matrix_multiply(Epetra_Comm& Comm, bool verbose)
{

  const int magic_num = 3000;
  // 2009/02/23: rabartl: If you are going to do a timing test you need to
  // make this number adjustable form the command-line and you need to put in
  // a real test that compares against hard numbers for pass/fail.

  int localn = magic_num/Comm.NumProc();

  Epetra_CrsMatrix* A = create_epetra_crsmatrix(Comm.NumProc(),
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

Epetra_CrsMatrix* create_epetra_crsmatrix(int numProcs,
                                          int localProc,
                                          int local_n,
                                          bool callFillComplete,
                                          bool symmetric)
{
  (void)localProc;
#ifdef EPETRA_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  int global_num_rows = numProcs*local_n;

  Epetra_Map rowmap(global_num_rows, local_n, 0, comm);

  int nnz_per_row = 9;
  Epetra_CrsMatrix* matrix =
    new Epetra_CrsMatrix(Copy, rowmap, nnz_per_row);

  // Add  rows one-at-a-time
  double negOne = -1.0;
  double posTwo = 2.0;
  double val_L = symmetric ? negOne : -0.5;

  for (int i=0; i<local_n; i++) {
    int GlobalRow = matrix->GRID(i);
    int RowLess1 = GlobalRow - 1;
    int RowPlus1 = GlobalRow + 1;
    int RowLess5 = GlobalRow - 5;
    int RowPlus5 = GlobalRow + 5;
    int RowLess9 = GlobalRow - 9;
    int RowPlus9 = GlobalRow + 9;
    int RowLess24 = GlobalRow - 24;
    int RowPlus24 = GlobalRow + 24;
    int RowLess48 = GlobalRow - 48;
    int RowPlus48 = GlobalRow + 48;

//    if (!symmetric) RowLess5 -= 2;

    if (RowLess48>=0) {
      matrix->InsertGlobalValues(GlobalRow, 1, &val_L, &RowLess48);
    }
    if (RowLess24>=0) {
      matrix->InsertGlobalValues(GlobalRow, 1, &val_L, &RowLess24);
    }
    if (RowLess9>=0) {
      matrix->InsertGlobalValues(GlobalRow, 1, &val_L, &RowLess9);
    }
    if (RowLess5>=0) {
      matrix->InsertGlobalValues(GlobalRow, 1, &val_L, &RowLess5);
    }
    if (RowLess1>=0) {
      matrix->InsertGlobalValues(GlobalRow, 1, &val_L, &RowLess1);
    }
    if (RowPlus1<global_num_rows) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowPlus1);
    }
    if (RowPlus5<global_num_rows) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowPlus5);
    }
    if (RowPlus9<global_num_rows) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowPlus9);
    }
    if (RowPlus24<global_num_rows) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowPlus24);
    }
    if (RowPlus48<global_num_rows) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowPlus48);
    }

    matrix->InsertGlobalValues(GlobalRow, 1, &posTwo, &GlobalRow);
  }

  if (callFillComplete) {
    int err = matrix->FillComplete();
    if (err != 0) {
      std::cout << "create_epetra_matrix: error in matrix.FillComplete()"
         <<std::endl;
    }
  }

  return(matrix);
}

int test_drumm1(Epetra_Comm& Comm)
{
  int size = Comm.NumProc();
  if (size != 2) return 0;

  int rank = Comm.MyPID();

  int indexBase = 0;
  int numGlobalElements = 2;
  
  Epetra_Map emap(numGlobalElements, indexBase, Comm);

  Epetra_CrsMatrix A(Copy, emap, 0);

  // 2x2 matrix:
  //   3 4
  //   1 2
  std::vector<std::vector<double> > vals(numGlobalElements);
  vals[0].push_back(3); vals[0].push_back(4);
  vals[1].push_back(1); vals[1].push_back(2);

  std::vector<int> indices;
  indices.push_back(0); indices.push_back(1);

  for (int row=0; row<numGlobalElements; ++row) {
    if ( A.MyGRID(row) )
      A.InsertGlobalValues(row, numGlobalElements, &(vals[row][0]), &indices[0]);
  }

  A.FillComplete();

  Epetra_CrsMatrix B(Copy, emap, 0);
  EpetraExt::MatrixMatrix::Multiply(A, true, A, false, B);

  // B = Transpose(A) x A should be
  //  10 14
  //  14 20
  int idx[2];
  int tmp;
  double val[2];

  //for this little test, global_row == rank:
  B.ExtractGlobalRowCopy(rank, 2, tmp, val, idx);

  int test_result = 0;

  if (rank == 0) {
    if (idx[0] == 0 && val[0] != 10.0) test_result = 1;
    if (idx[1] == 0 && val[1] != 10.0) test_result = 1;
    if (idx[0] == 1 && val[0] != 14.0) test_result = 1;
    if (idx[1] == 1 && val[1] != 14.0) test_result = 1;
  }
  else {
    if (idx[0] == 0 && val[0] != 14.0) test_result = 1;
    if (idx[1] == 0 && val[1] != 14.0) test_result = 1;
    if (idx[0] == 1 && val[0] != 20.0) test_result = 1;
    if (idx[1] == 1 && val[1] != 20.0) test_result = 1;
  }

  int global_test_result = 0;
  Comm.SumAll(&test_result, &global_test_result, 1);

  return global_test_result;
}

