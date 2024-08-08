//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
#include <cstdio>

#include <ctime>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <limits.h>
#include <cmath>

#include <Kokkos_Core.hpp>
#include <Kokkos_MV.hpp>
#include <KokkosSparse_spmv.hpp>
#include <Kokkos_CrsMatrix.hpp>
#include <Kokkos_Timer.hpp>

typedef Kokkos::DefaultExecutionSpace execution_space;

#ifdef INT64
typedef long long int LocalOrdinalType;
#else
typedef int LocalOrdinalType;
#endif
template <typename ScalarType, typename OrdinalType>
int SparseMatrix_MatrixMarket_read(const char* filename, OrdinalType& nrows, OrdinalType& ncols, OrdinalType& nnz,
                                   ScalarType*& values, OrdinalType*& rowPtr, OrdinalType*& colInd) {
  FILE* file = fopen(filename, "r");
  char line[512];
  line[0]         = '%';
  int count       = -1;
  char* symmetric = NULL;
  int nlines;

  /*  //determine the size of the file to allocate memory
    fseek(file,0,SEEK_END);
    long long int fileSize = ftell(file);
    fseek(file,0,SEEK_SET);

    //now create array to take the input data from file.
    char* buffer = (char*) malloc(sizeof(char)*fileSize);

    //read content of file to memory
    newFileSize = fread (buffer,1,fileSize,file);*/

  while (line[0] == '%') {
    fgets(line, 511, file);
    count++;
    if (count == 0) symmetric = strstr(line, "symmetric");
  }
  rewind(file);
  for (int i = 0; i < count; i++) fgets(line, 511, file);
  fscanf(file, "%i", &nrows);
  fscanf(file, "%i", &ncols);
  fscanf(file, "%i", &nlines);
  printf("Matrix dimension: %i %i %i %s\n", nrows, ncols, nlines, symmetric ? "Symmetric" : "General");
  if (symmetric)
    nnz = nlines * 2;
  else
    nnz = nlines;

  OrdinalType* colIndtmp            = new OrdinalType[nnz];
  OrdinalType* rowIndtmp            = new OrdinalType[nnz];
  double* valuestmp                 = new double[nnz];
  OrdinalType* priorEntrySameRowInd = new OrdinalType[nnz];
  OrdinalType* lastEntryWithRowInd  = new OrdinalType[nrows];
  for (int i = 0; i < nrows; i++) lastEntryWithRowInd[i] = -1;
  nnz = 0;
  for (int ii = 0; ii < nlines; ii++) {
    fscanf(file, "%i %i %le", &rowIndtmp[nnz], &colIndtmp[nnz], &valuestmp[nnz]);
    priorEntrySameRowInd[nnz]               = lastEntryWithRowInd[rowIndtmp[nnz] - 1];
    lastEntryWithRowInd[rowIndtmp[nnz] - 1] = nnz;
    if ((symmetric) && (rowIndtmp[nnz] != colIndtmp[nnz])) {
      nnz++;
      rowIndtmp[nnz]                          = colIndtmp[nnz - 1];
      colIndtmp[nnz]                          = rowIndtmp[nnz - 1];
      valuestmp[nnz]                          = valuestmp[nnz - 1];
      priorEntrySameRowInd[nnz]               = lastEntryWithRowInd[rowIndtmp[nnz] - 1];
      lastEntryWithRowInd[rowIndtmp[nnz] - 1] = nnz;
    }
    nnz++;
  }

  values = new ScalarType[nnz];
  colInd = new OrdinalType[nnz];
  rowPtr = new OrdinalType[nrows + 1];

  int pos = 0;
  for (int row = 0; row < nrows; row++) {
    int j       = lastEntryWithRowInd[row];
    rowPtr[row] = pos;
    while (j > -1) {
      values[pos] = valuestmp[j];
      colInd[pos] = colIndtmp[j] - 1;
      j           = priorEntrySameRowInd[j];
      pos++;
    }
  }
  rowPtr[nrows] = pos;
  printf("Number of Non-Zeros: %i\n", pos);
  delete[] valuestmp;
  delete[] colIndtmp;
  delete[] rowIndtmp;
  delete[] priorEntrySameRowInd;
  delete[] lastEntryWithRowInd;
  return nnz;
}

template <typename ScalarType, typename OrdinalType>
int SparseMatrix_ExtractBinaryGraph(const char* filename, OrdinalType& nrows, OrdinalType& ncols, OrdinalType& nnz,
                                    ScalarType*& values, OrdinalType*& rowPtr, OrdinalType*& colInd) {
  nnz = SparseMatrix_MatrixMarket_read<ScalarType, OrdinalType>(filename, nrows, ncols, nnz, values, rowPtr, colInd);
  char* filename_row = new char[strlen(filename) + 5];
  char* filename_col = new char[strlen(filename) + 5];
  strcpy(filename_row, filename);
  strcpy(filename_col, filename);
  strcat(filename_row, "_row");
  strcat(filename_row, "_col");
  FILE* RowFile = fopen(filename_row, "w");
  FILE* ColFile = fopen(filename_col, "w");

  fwrite(rowPtr, sizeof(OrdinalType), nrows + 1, RowFile);
  fwrite(colInd, sizeof(OrdinalType), nnz, ColFile);
  fclose(RowFile);
  fclose(ColFile);
  return nnz;
}

template <typename ScalarType, typename OrdinalType>
int SparseMatrix_ReadBinaryGraph(const char* filename, OrdinalType& nrows, OrdinalType& ncols, OrdinalType& nnz,
                                 ScalarType*& values, OrdinalType*& rowPtr, OrdinalType*& colInd) {
  char* filename_descr = new char[strlen(filename) + 7];
  strcpy(filename_descr, filename);
  strcat(filename_descr, "_descr");
  FILE* file = fopen(filename_descr, "r");
  char line[512];
  line[0]         = '%';
  int count       = -1;
  char* symmetric = NULL;
  int nlines;

  while (line[0] == '%') {
    fgets(line, 511, file);
    count++;
    if (count == 0) symmetric = strstr(line, "symmetric");
  }
  rewind(file);
  for (int i = 0; i < count; i++) fgets(line, 511, file);
  fscanf(file, "%i", &nrows);
  fscanf(file, "%i", &ncols);
  fscanf(file, "%i", &nlines);
  printf("Matrix dimension: %i %i %i %s\n", nrows, ncols, nlines, symmetric ? "Symmetric" : "General");
  if (symmetric)
    nnz = nlines * 2;
  else
    nnz = nlines;
  fclose(file);

  char* filename_row = new char[strlen(filename) + 5];
  char* filename_col = new char[strlen(filename) + 5];
  strcpy(filename_row, filename);
  strcpy(filename_col, filename);
  strcat(filename_row, "_row");
  strcat(filename_col, "_col");
  FILE* RowFile = fopen(filename_row, "r");
  FILE* ColFile = fopen(filename_col, "r");

  values = new ScalarType[nnz];
  rowPtr = new OrdinalType[nrows + 1];
  colInd = new OrdinalType[nnz];

  fread(rowPtr, sizeof(OrdinalType), nrows + 1, RowFile);
  fread(colInd, sizeof(OrdinalType), nnz, ColFile);
  fclose(RowFile);
  fclose(ColFile);
  return nnz;
}

template <typename ScalarType, typename OrdinalType>
int SparseMatrix_generate(OrdinalType nrows, OrdinalType ncols, OrdinalType& nnz, OrdinalType varianz_nel_row,
                          OrdinalType width_row, ScalarType*& values, OrdinalType*& rowPtr, OrdinalType*& colInd) {
  rowPtr = new OrdinalType[nrows + 1];

  OrdinalType elements_per_row = nnz / nrows;
  srand(13721);
  rowPtr[0] = 0;
  for (int row = 0; row < nrows; row++) {
    int varianz     = (1.0 * rand() / INT_MAX - 0.5) * varianz_nel_row;
    rowPtr[row + 1] = rowPtr[row] + elements_per_row + varianz;
  }
  nnz    = rowPtr[nrows];
  values = new ScalarType[nnz];
  colInd = new OrdinalType[nnz];
  for (int row = 0; row < nrows; row++) {
    for (int k = rowPtr[row]; k < rowPtr[row + 1]; k++) {
      int pos = (1.0 * rand() / INT_MAX - 0.5) * width_row + row;
      if (pos < 0) pos += ncols;
      if (pos >= ncols) pos -= ncols;
      colInd[k] = pos;
      values[k] = 100.0 * rand() / INT_MAX - 50.0;
    }
  }
  return nnz;
}

template <typename Scalar>
int test_crs_matrix_test(LocalOrdinalType numRows, LocalOrdinalType numCols, LocalOrdinalType nnz,
                         LocalOrdinalType numVecs, LocalOrdinalType test, const char* filename, const bool binaryfile) {
  typedef Kokkos::CrsMatrix<Scalar, LocalOrdinalType, execution_space, void, int> matrix_type;
  typedef typename Kokkos::MultiVectorDynamic<Scalar, execution_space>::type mv_type;
  typedef typename Kokkos::MultiVectorDynamic<Scalar, execution_space>::random_read_type mv_random_read_type;
  typedef typename mv_type::HostMirror h_mv_type;

  Scalar* val           = NULL;
  LocalOrdinalType* row = NULL;
  LocalOrdinalType* col = NULL;

  srand(17312837);
  if (filename == NULL)
    nnz = SparseMatrix_generate<Scalar, LocalOrdinalType>(numRows, numCols, nnz, nnz / numRows * 0.2, numRows * 0.01,
                                                          val, row, col);
  else if (!binaryfile)
    nnz = SparseMatrix_MatrixMarket_read<Scalar, LocalOrdinalType>(filename, numRows, numCols, nnz, val, row, col);
  else
    nnz = SparseMatrix_ReadBinaryGraph<Scalar, LocalOrdinalType>(filename, numRows, numCols, nnz, val, row, col);

  matrix_type A("CRS::A", numRows, numCols, nnz, val, row, col, false);

  mv_type x("X", numCols, numVecs);
  mv_random_read_type t_x(x);
  mv_type y("Y", numRows, numVecs);
  h_mv_type h_x                                                = Kokkos::create_mirror_view(x);
  h_mv_type h_y                                                = Kokkos::create_mirror_view(y);
  h_mv_type h_y_compare                                        = Kokkos::create_mirror(y);
  typename matrix_type::StaticCrsGraphType::HostMirror h_graph = Kokkos::create_mirror(A.graph);
  typename matrix_type::values_type::HostMirror h_values       = Kokkos::create_mirror_view(A.values);

  // Kokkos::deep_copy(h_graph.row_map,A.graph.row_map);
  for (LocalOrdinalType k = 0; k < numVecs; k++) {
    // h_a(k) = (Scalar) (1.0*(rand()%40)-20.);
    for (LocalOrdinalType i = 0; i < numCols; i++) {
      h_x(i, k) = (Scalar)(1.0 * (rand() % 40) - 20.);
      h_y(i, k) = (Scalar)(1.0 * (rand() % 40) - 20.);
    }
  }
  for (LocalOrdinalType i = 0; i < numRows; i++) {
    LocalOrdinalType start = h_graph.row_map(i);
    LocalOrdinalType end   = h_graph.row_map(i + 1);
    for (LocalOrdinalType j = start; j < end; j++) {
      h_values(j) = h_graph.entries(j) + i;
    }
    for (LocalOrdinalType k = 0; k < numVecs; k++) h_y_compare(i, k) = 0;
    for (LocalOrdinalType j = start; j < end; j++) {
      Scalar val           = h_graph.entries(j) + i;
      LocalOrdinalType idx = h_graph.entries(j);
      for (LocalOrdinalType k = 0; k < numVecs; k++) h_y_compare(i, k) += val * h_x(idx, k);
    }
  }

  Kokkos::deep_copy(x, h_x);
  Kokkos::deep_copy(y, h_y);
  Kokkos::deep_copy(A.graph.entries, h_graph.entries);
  Kokkos::deep_copy(A.values, h_values);

#ifdef NEWKERNEL
  KokkosSparse::spmv("N", 1.0, A, x, 0.0, y);
#else
  Kokkos::MV_Multiply(y, A, x);
#endif
  execution_space().fence();
  Kokkos::deep_copy(h_y, y);
  Scalar error[numVecs];
  Scalar sum[numVecs];
  for (LocalOrdinalType k = 0; k < numVecs; k++) {
    error[k] = 0;
    sum[k]   = 0;
  }
  for (LocalOrdinalType i = 0; i < numRows; i++)
    for (LocalOrdinalType k = 0; k < numVecs; k++) {
      error[k] += (h_y_compare(i, k) - h_y(i, k)) * (h_y_compare(i, k) - h_y(i, k));
      sum[k] += h_y_compare(i, k) * h_y_compare(i, k);
      // prLocalOrdinalTypef("%i %i %lf %lf
      // %lf\n",i,k,h_y_compare(i,k),h_y(i,k),h_x(i,k));
    }

  // for(LocalOrdinalType i=0;i<A.nnz;i++) prLocalOrdinalTypef("%i
  // %lf\n",h_graph.entries(i),h_values(i));
  LocalOrdinalType num_errors = 0;
  double total_error          = 0;
  double total_sum            = 0;
  for (LocalOrdinalType k = 0; k < numVecs; k++) {
    num_errors += (error[k] / (sum[k] == 0 ? 1 : sum[k])) > 1e-5 ? 1 : 0;
    total_error += error[k];
    total_sum += sum[k];
  }

  LocalOrdinalType loop = 10;
  Kokkos::Timer timer;
  for (LocalOrdinalType i = 0; i < loop; i++)
#ifdef NEWKERNEL
    KokkosSparse::spmv("N", 1.0, A, x, 0.0, y);
#else
    Kokkos::MV_Multiply(y, A, x);
#endif
  execution_space().fence();
  double time = timer.seconds();
  double matrix_size =
      1.0 * ((nnz * (sizeof(Scalar) + sizeof(LocalOrdinalType)) + numRows * sizeof(LocalOrdinalType))) / 1024 / 1024;
  double vector_size      = 2.0 * numRows * numVecs * sizeof(Scalar) / 1024 / 1024;
  double vector_readwrite = (nnz + numCols) * numVecs * sizeof(Scalar) / 1024 / 1024;

  double problem_size = matrix_size + vector_size;
  printf("%i %i %i %i %6.2lf MB %6.2lf GB/s %6.2lf GFlop/s %6.3lf ms %i\n", nnz, numRows, numCols, numVecs,
         problem_size, (matrix_size + vector_readwrite) / time * loop / 1024, 2.0 * nnz * numVecs * loop / time / 1e9,
         time / loop * 1000, num_errors);
  return (int)total_error;
}

template <typename Scalar>
int test_crs_matrix_test_singlevec(int numRows, int numCols, int nnz, int test, const char* filename,
                                   const bool binaryfile) {
  typedef Kokkos::CrsMatrix<Scalar, int, execution_space, void, int> matrix_type;
  typedef typename Kokkos::View<Scalar*, Kokkos::LayoutLeft, execution_space> mv_type;
  typedef typename Kokkos::View<Scalar*, Kokkos::LayoutLeft, execution_space, Kokkos::MemoryRandomAccess>
      mv_random_read_type;
  typedef typename mv_type::HostMirror h_mv_type;

  Scalar* val = NULL;
  int* row    = NULL;
  int* col    = NULL;

  srand(17312837);
  if (filename == NULL)
    nnz = SparseMatrix_generate<Scalar, int>(numRows, numCols, nnz, nnz / numRows * 0.2, numRows * 0.01, val, row, col);
  else if (!binaryfile)
    nnz = SparseMatrix_MatrixMarket_read<Scalar, int>(filename, numRows, numCols, nnz, val, row, col);
  else
    nnz = SparseMatrix_ReadBinaryGraph<Scalar, int>(filename, numRows, numCols, nnz, val, row, col);

  matrix_type A("CRS::A", numRows, numCols, nnz, val, row, col, false);

  mv_type x("X", numCols);
  mv_random_read_type t_x(x);
  mv_type y("Y", numRows);
  h_mv_type h_x                                                = Kokkos::create_mirror_view(x);
  h_mv_type h_y                                                = Kokkos::create_mirror_view(y);
  h_mv_type h_y_compare                                        = Kokkos::create_mirror(y);
  typename matrix_type::StaticCrsGraphType::HostMirror h_graph = Kokkos::create_mirror(A.graph);
  typename matrix_type::values_type::HostMirror h_values       = Kokkos::create_mirror_view(A.values);

  // Kokkos::deep_copy(h_graph.row_map,A.graph.row_map);
  // h_a(k) = (Scalar) (1.0*(rand()%40)-20.);
  for (int i = 0; i < numCols; i++) {
    h_x(i) = (Scalar)(1.0 * (rand() % 40) - 20.);
    h_y(i) = (Scalar)(1.0 * (rand() % 40) - 20.);
  }
  for (int i = 0; i < numRows; i++) {
    int start = h_graph.row_map(i);
    int end   = h_graph.row_map(i + 1);
    for (int j = start; j < end; j++) {
      h_values(j) = h_graph.entries(j) + i;
    }
    h_y_compare(i) = 0;
    for (int j = start; j < end; j++) {
      Scalar val = h_graph.entries(j) + i;
      int idx    = h_graph.entries(j);
      h_y_compare(i) += val * h_x(idx);
    }
  }

  Kokkos::deep_copy(x, h_x);
  Kokkos::deep_copy(y, h_y);
  Kokkos::deep_copy(A.graph.entries, h_graph.entries);
  Kokkos::deep_copy(A.values, h_values);
  /*for(int i=0;i<numRows;i++)
          for(int k = 0; k<numVecs; k++) {
    //error[k]+=(h_y_compare(i,k)-h_y(i,k))*(h_y_compare(i,k)-h_y(i,k));
    printf("%i %i %lf %lf %lf\n",i,k,h_y_compare(i,k),h_y(i,k),h_x(i,k));
          }*/
  typename Kokkos::CrsMatrix<Scalar, int, execution_space, void, int>::values_type x1("X1", numCols);
  typename Kokkos::CrsMatrix<Scalar, int, execution_space, void, int>::values_type y1("Y1", numRows);
#ifdef NEWKERNEL
  KokkosSparse::spmv("N", 1.0, A, x1, 0.0, y1);
#else
  Kokkos::MV_Multiply(y1, A, x1);
#endif

#ifdef NEWKERNEL
  KokkosSparse::spmv("N", 1.0, A, x, 0.0, y);
#else
  Kokkos::MV_Multiply(y, A, x);
#endif
  execution_space().fence();
  Kokkos::deep_copy(h_y, y);
  Scalar error = 0;
  Scalar sum   = 0;
  for (int i = 0; i < numRows; i++) {
    error += (h_y_compare(i) - h_y(i)) * (h_y_compare(i) - h_y(i));
    sum += h_y_compare(i) * h_y_compare(i);
    // printf("%i %i %lf %lf %lf\n",i,k,h_y_compare(i,k),h_y(i,k),h_x(i,k));
  }

  // for(int i=0;i<A.nnz;i++) printf("%i %lf\n",h_graph.entries(i),h_values(i));
  int num_errors     = 0;
  double total_error = 0;
  double total_sum   = 0;
  num_errors += (error / (sum == 0 ? 1 : sum)) > 1e-5 ? 1 : 0;
  total_error += error;
  total_sum += sum;

  int loop = 100;
  Kokkos::Timer timer;

  for (int i = 0; i < loop; i++)
#ifdef NEWKERNEL
    KokkosSparse::spmv("N", 1.0, A, x, 0.0, y);
#else
    Kokkos::MV_Multiply(y, A, x);
#endif
  execution_space().fence();
  double time             = timer.seconds();
  double matrix_size      = 1.0 * ((nnz * (sizeof(Scalar) + sizeof(int)) + numRows * sizeof(int))) / 1024 / 1024;
  double vector_size      = 2.0 * numRows * sizeof(Scalar) / 1024 / 1024;
  double vector_readwrite = (nnz + numCols) * sizeof(Scalar) / 1024 / 1024;

  double problem_size = matrix_size + vector_size;
  printf("%i %i %i %i %6.2lf MB %6.2lf GB/s %6.2lf GFlop/s %6.3lf ms %i\n", nnz, numRows, numCols, 1, problem_size,
         (matrix_size + vector_readwrite) / time * loop / 1024, 2.0 * nnz * loop / time / 1e9, time / loop * 1000,
         num_errors);
  return (int)total_error;
}

int test_crs_matrix_type(int numrows, int numcols, int nnz, int numVecs, int type, int test, const char* filename,
                         const bool binaryfile) {
  if (numVecs == 1)
    return test_crs_matrix_test_singlevec<double>(numrows, numcols, nnz, test, filename, binaryfile);
  else
    return test_crs_matrix_test<double>(numrows, numcols, nnz, numVecs, test, filename, binaryfile);
}

int main(int argc, char** argv) {
  long long int size = 110503;  // a prime number
  int numVecs        = 4;
  int test           = -1;
  int type           = -1;
  char* filename     = NULL;
  bool binaryfile    = false;

  for (int i = 0; i < argc; i++) {
    if ((strcmp(argv[i], "-s") == 0)) {
      size = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-v") == 0)) {
      numVecs = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "--test") == 0)) {
      test = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "--type") == 0)) {
      type = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-f") == 0)) {
      filename = argv[++i];
      continue;
    }
    if ((strcmp(argv[i], "-fb") == 0)) {
      filename   = argv[++i];
      binaryfile = true;
      continue;
    }
  }

  Kokkos::initialize(argc, argv);

  int numVecsList[10] = {1, 2, 3, 4, 5, 8, 11, 15, 16, 17};
  int maxNumVecs      = numVecs == -1 ? 17 : numVecs;
  int numVecIdx       = 0;
  if (numVecs == -1) numVecs = numVecsList[numVecIdx++];

  int total_errors = 0;
  while (numVecs <= maxNumVecs) {
    total_errors += test_crs_matrix_type(size, size, size * 10, numVecs, type, test, filename, binaryfile);
    if (numVecs < maxNumVecs)
      numVecs = numVecsList[numVecIdx++];
    else
      numVecs++;
  }

  if (total_errors == 0)
    printf("Kokkos::MultiVector Test: Passed\n");
  else
    printf("Kokkos::MultiVector Test: Failed\n");

  Kokkos::finalize();
}
