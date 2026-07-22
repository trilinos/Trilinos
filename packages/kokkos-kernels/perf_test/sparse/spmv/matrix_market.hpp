// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef MATRIX_MARKET_HPP_
#define MATRIX_MARKET_HPP_

#include <cstdio>
#include <cstddef>
#include <array>
#include <memory>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <limits.h>

namespace Impl {
template <typename OrdinalType>
void SparseGraph_SortRows(OrdinalType nrows, OrdinalType* rowPtr, OrdinalType* colInd) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int row = 0; row < nrows; row++) {
    OrdinalType row_start = rowPtr[row];
    OrdinalType row_end   = rowPtr[row + 1];
    for (OrdinalType i = row_start; i < row_end - 1; i++) {
      for (OrdinalType j = row_end - 1; j > i; j--) {
        if (colInd[j] < colInd[j - 1]) {
          int idx       = colInd[j];
          colInd[j]     = colInd[j - 1];
          colInd[j - 1] = idx;
        }
      }
    }
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int row = 0; row < nrows; row++) {
    OrdinalType row_start = rowPtr[row];
    OrdinalType row_end   = rowPtr[row + 1];
    for (OrdinalType i = row_start; i < row_end - 1; i++) {
      if (colInd[i + 1] < colInd[i]) printf("Error Not sorted %i %i | %i %i\n", row, i, colInd[i], colInd[i + 1]);
    }
  }
}
}  // namespace Impl

template <typename ScalarType, typename OrdinalType>
int SparseMatrix_MatrixMarket_read(const char* filename, OrdinalType& nrows, OrdinalType& ncols, OrdinalType& nnz,
                                   ScalarType*& values, OrdinalType*& rowPtr, OrdinalType*& colInd, bool sort,
                                   OrdinalType idx_offset = 0) {
  auto file = std::unique_ptr<FILE, decltype(&fclose)>(fopen(filename, "r"), &fclose);
  if (!file) {
    fprintf(stderr, "SparseMatrix_MatrixMarket_read: could not open \"%s\"\n", filename);
    return -1;
  }
  std::array<char, 512> line{};
  line[0]         = '%';
  int count       = -1;
  char* symmetric = NULL;
  char* pattern   = NULL;
  int nlines;

  while (line[0] == '%') {
    if (!fgets(line.data(), static_cast<int>(line.size()), file.get())) {
      fprintf(stderr, "SparseMatrix_MatrixMarket_read: unexpected EOF in header of \"%s\"\n", filename);
      return -1;
    }
    count++;
    if (count == 0) {
      symmetric = strstr(line.data(), "symmetric");
      pattern   = strstr(line.data(), "pattern");
    }
  }
  rewind(file.get());
  for (int i = 0; i < count; i++) {
    if (!fgets(line.data(), static_cast<int>(line.size()), file.get())) {
      fprintf(stderr, "SparseMatrix_MatrixMarket_read: unexpected EOF rewinding \"%s\"\n", filename);
      return -1;
    }
  }
  if (fscanf(file.get(), "%i", &nrows) != 1 || fscanf(file.get(), "%i", &ncols) != 1 ||
      fscanf(file.get(), "%i", &nlines) != 1) {
    fprintf(stderr, "SparseMatrix_MatrixMarket_read: invalid dimension line in \"%s\"\n", filename);
    return -1;
  }
  if (nrows <= 0 || ncols <= 0 || nlines < 0) {
    fprintf(stderr, "SparseMatrix_MatrixMarket_read: non-positive dimensions in \"%s\"\n", filename);
    return -1;
  }
  printf("Matrix dimension: %i %i %i %s %s\n", nrows, ncols, nlines, symmetric ? "Symmetric" : "General",
         pattern ? "Pattern" : "Real");

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
    if (pattern) {
      fscanf(file.get(), "%i %i", &rowIndtmp[nnz], &colIndtmp[nnz]);
      valuestmp[nnz] = (1.0 * (ii % nrows)) / ncols;
    } else
      fscanf(file.get(), "%i %i %le", &rowIndtmp[nnz], &colIndtmp[nnz], &valuestmp[nnz]);
    if (ii < 10 || ii > nlines - 10)
      printf("Read: %i %i %i %le\n", nnz, rowIndtmp[nnz], colIndtmp[nnz], valuestmp[nnz]);
    rowIndtmp[nnz] -= idx_offset;
    colIndtmp[nnz] -= idx_offset;
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

  size_t min_span = nrows + 1;
  size_t max_span = 0;
  size_t ave_span = 0;
  for (int row = 0; row < nrows; row++) {
    int min = nrows + 1;
    int max = 0;
    for (int i = rowPtr[row]; i < rowPtr[row + 1]; i++) {
      if (colInd[i] < min) min = colInd[i];
      if (colInd[i] > max) max = colInd[i];
    }
    if (rowPtr[row + 1] > rowPtr[row]) {
      size_t span = max - min;
      if (span < min_span) min_span = span;
      if (span > max_span) max_span = span;
      ave_span += span;
    } else
      min_span = 0;
  }

  printf("%lu Spans: %lu %lu %lu\n", (size_t)nnz, min_span, max_span, ave_span / nrows);
  if (sort) Impl::SparseGraph_SortRows<OrdinalType>(nrows, rowPtr, colInd);
  return nnz;
}

template <typename ScalarType, typename OrdinalType>
int SparseMatrix_WriteBinaryFormat(const char* filename, OrdinalType& nrows, OrdinalType& ncols, OrdinalType& nnz,
                                   ScalarType*& values, OrdinalType*& rowPtr, OrdinalType*& colInd, bool sort,
                                   OrdinalType idx_offset = 0) {
  nnz = SparseMatrix_MatrixMarket_read<ScalarType, OrdinalType>(filename, nrows, ncols, nnz, values, rowPtr, colInd,
                                                                sort, idx_offset);
  if (nnz < 0) return -1;

  char* filename_row   = new char[strlen(filename) + 5];
  char* filename_col   = new char[strlen(filename) + 5];
  char* filename_vals  = new char[strlen(filename) + 6];
  char* filename_descr = new char[strlen(filename) + 7];
  strcpy(filename_row, filename);
  strcpy(filename_col, filename);
  strcpy(filename_vals, filename);
  strcpy(filename_descr, filename);
  strcat(filename_row, "_row");
  strcat(filename_col, "_col");
  strcat(filename_vals, "_vals");
  strcat(filename_descr, "_descr");
  FILE* RowFile   = fopen(filename_row, "w");
  FILE* ColFile   = fopen(filename_col, "w");
  FILE* ValsFile  = fopen(filename_vals, "w");
  FILE* DescrFile = fopen(filename_descr, "w");
  if (!RowFile || !ColFile || !ValsFile || !DescrFile) {
    fprintf(stderr, "SparseMatrix_WriteBinaryFormat: could not open output file(s) for \"%s\"\n", filename);
    if (RowFile) fclose(RowFile);
    if (ColFile) fclose(ColFile);
    if (ValsFile) fclose(ValsFile);
    if (DescrFile) fclose(DescrFile);
    return -1;
  }

  auto file = std::unique_ptr<FILE, decltype(&fclose)>(fopen(filename, "r"), &fclose);
  if (!file) {
    fprintf(stderr, "SparseMatrix_WriteBinaryFormat: could not reopen \"%s\" for description copy\n", filename);
    fclose(RowFile);
    fclose(ColFile);
    fclose(ValsFile);
    fclose(DescrFile);
    return -1;
  }
  std::array<char, 512> line{};
  line[0]   = '%';
  int count = -1;
  // char* symmetric = NULL;
  // int nlines;

  while (line[0] == '%') {
    if (!fgets(line.data(), static_cast<int>(line.size()), file.get())) {
      fprintf(stderr, "SparseMatrix_WriteBinaryFormat: unexpected EOF in \"%s\"\n", filename);
      fclose(RowFile);
      fclose(ColFile);
      fclose(ValsFile);
      fclose(DescrFile);
      return -1;
    }
    count++;
    // if(count==0) symmetric=strstr(line,"symmetric");

    if (line[0] == '%')
      fprintf(DescrFile, "%s", line.data());
    else
      fprintf(DescrFile, "%i %i %i\n", nrows, ncols, nnz);
  }
  fprintf(DescrFile, "\n");

  fwrite(rowPtr, sizeof(OrdinalType), nrows + 1, RowFile);
  fwrite(colInd, sizeof(OrdinalType), nnz, ColFile);
  fwrite(values, sizeof(ScalarType), nnz, ValsFile);

  fclose(RowFile);
  fclose(ColFile);
  fclose(ValsFile);
  fclose(DescrFile);

  size_t min_span = nrows + 1;
  size_t max_span = 0;
  size_t ave_span = 0;
  for (int row = 0; row < nrows; row++) {
    int min = nrows + 1;
    int max = 0;
    for (int i = rowPtr[row]; i < rowPtr[row + 1]; i++) {
      if (colInd[i] < min) min = colInd[i];
      if (colInd[i] > max) max = colInd[i];
    }
    if (rowPtr[row + 1] > rowPtr[row]) {
      size_t span = max - min;
      if (span < min_span) min_span = span;
      if (span > max_span) max_span = span;
      ave_span += span;
    } else
      min_span = 0;
  }
  printf("%lu Spans: %lu %lu %lu\n", (size_t)nnz, min_span, max_span, ave_span / nrows);

  return nnz;
}

template <typename ScalarType, typename OrdinalType>
int SparseMatrix_ReadBinaryFormat(const char* filename, OrdinalType& nrows, OrdinalType& ncols, OrdinalType& nnz,
                                  ScalarType*& values, OrdinalType*& rowPtr, OrdinalType*& colInd) {
  char* filename_descr = new char[strlen(filename) + 7];
  strcpy(filename_descr, filename);
  strcat(filename_descr, "_descr");
  auto file = std::unique_ptr<FILE, decltype(&fclose)>(fopen(filename_descr, "r"), &fclose);
  if (!file) {
    fprintf(stderr, "SparseMatrix_ReadBinaryFormat: could not open \"%s\"\n", filename_descr);
    delete[] filename_descr;
    return -1;
  }
  std::array<char, 512> line{};
  line[0]         = '%';
  int count       = -1;
  char* symmetric = NULL;
  // int nlines;

  while (line[0] == '%') {
    if (!fgets(line.data(), static_cast<int>(line.size()), file.get())) {
      fprintf(stderr, "SparseMatrix_ReadBinaryFormat: unexpected EOF in \"%s\"\n", filename_descr);
      delete[] filename_descr;
      return -1;
    }
    count++;
    if (count == 0) symmetric = strstr(line.data(), "symmetric");
  }
  rewind(file.get());
  for (int i = 0; i < count; i++) {
    if (!fgets(line.data(), static_cast<int>(line.size()), file.get())) {
      fprintf(stderr, "SparseMatrix_ReadBinaryFormat: unexpected EOF rewinding \"%s\"\n", filename_descr);
      delete[] filename_descr;
      return -1;
    }
  }
  if (fscanf(file.get(), "%i", &nrows) != 1 || fscanf(file.get(), "%i", &ncols) != 1 ||
      fscanf(file.get(), "%i", &nnz) != 1) {
    fprintf(stderr, "SparseMatrix_ReadBinaryFormat: invalid dimension line in \"%s\"\n", filename_descr);
    delete[] filename_descr;
    return -1;
  }
  if (nrows <= 0 || ncols <= 0 || nnz < 0) {
    fprintf(stderr, "SparseMatrix_ReadBinaryFormat: non-positive dimensions in \"%s\"\n", filename_descr);
    delete[] filename_descr;
    return -1;
  }
  printf("Matrix dimension: %i %i %i %s\n", nrows, ncols, nnz, symmetric ? "Symmetric" : "General");

  delete[] filename_descr;

  char* filename_row  = new char[strlen(filename) + 5];
  char* filename_col  = new char[strlen(filename) + 5];
  char* filename_vals = new char[strlen(filename) + 6];
  strcpy(filename_row, filename);
  strcpy(filename_col, filename);
  strcpy(filename_vals, filename);
  strcat(filename_row, "_row");
  strcat(filename_col, "_col");
  strcat(filename_vals, "_vals");
  FILE* RowFile  = fopen(filename_row, "r");
  FILE* ColFile  = fopen(filename_col, "r");
  FILE* ValsFile = fopen(filename_vals, "r");

  const bool read_values = (ValsFile != NULL);

  if (!RowFile || !ColFile) {
    fprintf(stderr, "SparseMatrix_ReadBinaryFormat: could not open row/col binary for \"%s\"\n", filename);
    if (RowFile) fclose(RowFile);
    if (ColFile) fclose(ColFile);
    if (ValsFile) fclose(ValsFile);
    delete[] filename_row;
    delete[] filename_col;
    delete[] filename_vals;
    return -1;
  }

  values = new ScalarType[nnz];
  rowPtr = new OrdinalType[nrows + 1];
  colInd = new OrdinalType[nnz];

  if (fread(rowPtr, sizeof(OrdinalType), nrows + 1, RowFile) != static_cast<size_t>(nrows + 1) ||
      fread(colInd, sizeof(OrdinalType), nnz, ColFile) != static_cast<size_t>(nnz)) {
    fprintf(stderr, "SparseMatrix_ReadBinaryFormat: short read row/col data for \"%s\"\n", filename);
    delete[] values;
    delete[] rowPtr;
    delete[] colInd;
    fclose(RowFile);
    fclose(ColFile);
    if (ValsFile) fclose(ValsFile);
    delete[] filename_row;
    delete[] filename_col;
    delete[] filename_vals;
    return -1;
  }

  fclose(RowFile);
  fclose(ColFile);
  if (read_values) {
    if (fread(values, sizeof(ScalarType), nnz, ValsFile) != static_cast<size_t>(nnz)) {
      fprintf(stderr, "SparseMatrix_ReadBinaryFormat: short read values for \"%s\"\n", filename);
      delete[] values;
      delete[] rowPtr;
      delete[] colInd;
      fclose(ValsFile);
      delete[] filename_row;
      delete[] filename_col;
      delete[] filename_vals;
      return -1;
    }
    fclose(ValsFile);
  } else {
    for (int i = 0; i < nnz; i++) values[i] = 0.001 * (rand() % 1000);
  }

  delete[] filename_row;
  delete[] filename_col;
  delete[] filename_vals;

  size_t min_span = nrows + 1;
  size_t max_span = 0;
  size_t ave_span = 0;
  for (int row = 0; row < nrows; row++) {
    int min = nrows + 1;
    int max = 0;
    for (int i = rowPtr[row]; i < rowPtr[row + 1]; i++) {
      if (colInd[i] < min) min = colInd[i];
      if (colInd[i] > max) max = colInd[i];
    }
    if (rowPtr[row + 1] > rowPtr[row]) {
      size_t span = max - min;
      if (span < min_span) min_span = span;
      if (span > max_span) max_span = span;
      ave_span += span;
    } else
      min_span = 0;
  }
  printf("%lu Spans: %lu %lu %lu\n", (size_t)nnz, min_span, max_span, ave_span / nrows);

  return nnz;
}

#endif /* MATRIX_MARKET_HPP_ */
