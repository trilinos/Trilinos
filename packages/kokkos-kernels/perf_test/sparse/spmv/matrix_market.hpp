/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef MATRIX_MARKET_HPP_
#define MATRIX_MARKET_HPP_

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <limits.h>

namespace Impl {
template< typename OrdinalType>
void SparseGraph_SortRows(OrdinalType nrows, OrdinalType* rowPtr, OrdinalType* colInd) {
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int row = 0; row < nrows; row++) {
    OrdinalType row_start = rowPtr[row];
    OrdinalType row_end = rowPtr[row+1];
    for(OrdinalType i = row_start; i < row_end-1; i++) {
      for(OrdinalType j = row_end-1; j > i; j--) {
        if(colInd[j] < colInd[j-1]) {
          int idx = colInd[j];
          colInd[j] = colInd[j-1];
          colInd[j-1] = idx;
        }
      }
    }
  }

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int row = 0; row < nrows; row++) {
    OrdinalType row_start = rowPtr[row];
    OrdinalType row_end = rowPtr[row+1];
    for(OrdinalType i = row_start; i < row_end-1; i++) {
      if(colInd[i+1] < colInd[i]) printf("Error Not sorted %i %i | %i %i\n",row,i,colInd[i],colInd[i+1]);
    }
  }
}
}

template< typename ScalarType , typename OrdinalType>
int SparseMatrix_MatrixMarket_read(const char* filename, OrdinalType &nrows, OrdinalType &ncols, OrdinalType &nnz,
                                   ScalarType* &values, OrdinalType* &rowPtr, OrdinalType* &colInd, bool sort, OrdinalType idx_offset  = 0)
{
  FILE* file = fopen(filename,"r");
  char line[512];
  line[0]='%';
  int count=-1;
  char* symmetric = NULL;
  char* pattern = NULL;
  int nlines;

  while(line[0]=='%')
  {
          fgets(line,511,file);
          count++;
          if(count==0) {
            symmetric=strstr(line,"symmetric");
            pattern=strstr(line,"pattern");
          }
  }
  rewind(file);
  for(int i=0;i<count;i++)
          fgets(line,511,file);
  fscanf(file,"%i",&nrows);
  fscanf(file,"%i",&ncols);
  fscanf(file,"%i",&nlines);
  printf("Matrix dimension: %i %i %i %s %s\n",nrows,ncols,nlines,symmetric?"Symmetric":"General",pattern?"Pattern":"Real");

  if(symmetric) nnz=nlines*2;
  else nnz=nlines;

  OrdinalType* colIndtmp = new OrdinalType[nnz];
  OrdinalType* rowIndtmp = new OrdinalType[nnz];
  double* valuestmp = new double[nnz];
  OrdinalType* priorEntrySameRowInd = new OrdinalType[nnz];
  OrdinalType* lastEntryWithRowInd = new OrdinalType[nrows];
  for(int i=0;i<nrows;i++) lastEntryWithRowInd[i]=-1;
  nnz=0;
  for(int ii=0;ii<nlines;ii++)
  {

          if(pattern) {
            fscanf(file,"%i %i",&rowIndtmp[nnz],&colIndtmp[nnz]);
            valuestmp[nnz] = (1.0*(ii%nrows))/ncols;
          }
          else
            fscanf(file,"%i %i %le",&rowIndtmp[nnz],&colIndtmp[nnz],&valuestmp[nnz]);
          if(ii<10||ii>nlines-10)
            printf("Read: %i %i %i %le\n",nnz,rowIndtmp[nnz],colIndtmp[nnz],valuestmp[nnz]);
          rowIndtmp[nnz]-= idx_offset;
          colIndtmp[nnz]-= idx_offset;
          priorEntrySameRowInd[nnz] = lastEntryWithRowInd[rowIndtmp[nnz]-1];
          lastEntryWithRowInd[rowIndtmp[nnz]-1]=nnz;
          if((symmetric) && (rowIndtmp[nnz]!=colIndtmp[nnz]))
          {
            nnz++;
            rowIndtmp[nnz]=colIndtmp[nnz-1];
            colIndtmp[nnz]=rowIndtmp[nnz-1];
            valuestmp[nnz]=valuestmp[nnz-1];
            priorEntrySameRowInd[nnz] = lastEntryWithRowInd[rowIndtmp[nnz]-1];
            lastEntryWithRowInd[rowIndtmp[nnz]-1]=nnz;
          }

          nnz++;
  }

  values = new ScalarType[nnz];
  colInd = new OrdinalType[nnz];
  rowPtr = new OrdinalType[nrows+1];

  int pos = 0;
  for(int row=0;row<nrows;row++)
  {
        int j = lastEntryWithRowInd[row];
        rowPtr[row]=pos;
    while(j>-1)
    {
        values[pos] = valuestmp[j];
        colInd[pos] = colIndtmp[j]-1;
        j = priorEntrySameRowInd[j];
        pos++;
    }
  }
  rowPtr[nrows]=pos;

  printf("Number of Non-Zeros: %i\n",pos);
  delete [] valuestmp;
  delete [] colIndtmp;
  delete [] rowIndtmp;
  delete [] priorEntrySameRowInd;
  delete [] lastEntryWithRowInd;

  size_t min_span = nrows+1;
  size_t max_span = 0;
  size_t ave_span = 0;
  for(int row=0; row<nrows;row++) {
    int min = nrows+1; int max = 0;
    for(int i=rowPtr[row]; i<rowPtr[row+1]; i++) {
      if(colInd[i]<min) min = colInd[i];
      if(colInd[i]>max) max = colInd[i];
    }
    if(rowPtr[row+1]>rowPtr[row]) {
      size_t span = max-min;
      if(span<min_span) min_span = span;
      if(span>max_span) max_span = span;
      ave_span += span;
    } else min_span = 0;
  }

  printf("%lu Spans: %lu %lu %lu\n",(size_t) nnz,min_span,max_span,ave_span/nrows);
  if(sort)
    Impl::SparseGraph_SortRows<OrdinalType>(nrows,rowPtr,colInd);
  return nnz;
}

template< typename ScalarType , typename OrdinalType>
int SparseMatrix_WriteBinaryFormat(const char* filename, OrdinalType &nrows, OrdinalType &ncols, OrdinalType &nnz,
                                    ScalarType* &values, OrdinalType* &rowPtr, OrdinalType* &colInd, bool sort, OrdinalType idx_offset = 0)
{
  nnz = SparseMatrix_MatrixMarket_read<ScalarType,OrdinalType>(filename,nrows,ncols,nnz,values,rowPtr,colInd,sort, idx_offset);


  char * filename_row = new char[strlen(filename)+5];
  char * filename_col = new char[strlen(filename)+5];
  char * filename_vals = new char[strlen(filename)+6];
  char * filename_descr = new char[strlen(filename)+7];
  strcpy(filename_row,filename);
  strcpy(filename_col,filename);
  strcpy(filename_vals,filename);
  strcpy(filename_descr,filename);
  strcat(filename_row,"_row");
  strcat(filename_col,"_col");
  strcat(filename_vals,"_vals");
  strcat(filename_descr,"_descr");
  FILE* RowFile = fopen(filename_row,"w");
  FILE* ColFile = fopen(filename_col,"w");
  FILE* ValsFile = fopen(filename_vals,"w");
  FILE* DescrFile = fopen(filename_descr,"w");


  FILE* file = fopen(filename,"r");
  char line[512];
  line[0]='%';
  int count=-1;
  //char* symmetric = NULL;
  //int nlines;

  while(line[0]=='%')
  {
          fgets(line,511,file);
          line[511] = 0;
          count++;
          //if(count==0) symmetric=strstr(line,"symmetric");

          if(line[0]=='%')
            fprintf ( DescrFile , "%s",line);
          else
            fprintf ( DescrFile , "%i %i %i\n",nrows,ncols,nnz);
  }
  fprintf ( DescrFile , "\n");

  fwrite ( rowPtr, sizeof(OrdinalType), nrows+1, RowFile);
  fwrite ( colInd, sizeof(OrdinalType), nnz, ColFile);
  fwrite ( values, sizeof(ScalarType), nnz, ValsFile);

  fclose(RowFile);
  fclose(ColFile);
  fclose(ValsFile);
  fclose(DescrFile);

  size_t min_span = nrows+1;
  size_t max_span = 0;
  size_t ave_span = 0;
  for(int row=0; row<nrows;row++) {
    int min = nrows+1; int max = 0;
    for(int i=rowPtr[row]; i<rowPtr[row+1]; i++) {
      if(colInd[i]<min) min = colInd[i];
      if(colInd[i]>max) max = colInd[i];
    }
    if(rowPtr[row+1]>rowPtr[row]) {
      size_t span = max-min;
      if(span<min_span) min_span = span;
      if(span>max_span) max_span = span;
      ave_span += span;
    } else min_span = 0;
  }
  printf("%lu Spans: %lu %lu %lu\n",(size_t) nnz,min_span,max_span,ave_span/nrows);

  return nnz;
}

template< typename ScalarType , typename OrdinalType>
int SparseMatrix_ReadBinaryFormat(const char* filename, OrdinalType &nrows, OrdinalType &ncols, OrdinalType &nnz, ScalarType* &values, OrdinalType* &rowPtr, OrdinalType* &colInd)
{
  char * filename_descr = new char[strlen(filename)+7];
  strcpy(filename_descr,filename);
  strcat(filename_descr,"_descr");
  FILE* file = fopen(filename_descr,"r");
  char line[512];
  line[0]='%';
  int count=-1;
  char* symmetric = NULL;
  //int nlines;

  while(line[0]=='%')
  {
          fgets(line,511,file);
          count++;
          if(count==0) symmetric=strstr(line,"symmetric");
  }
  rewind(file);
  for(int i=0;i<count;i++)
          fgets(line,511,file);
  fscanf(file,"%i",&nrows);
  fscanf(file,"%i",&ncols);
  fscanf(file,"%i",&nnz);
  printf("Matrix dimension: %i %i %i %s\n",nrows,ncols,nnz,symmetric?"Symmetric":"General");

  fclose(file);

  char * filename_row = new char[strlen(filename)+5];
  char * filename_col = new char[strlen(filename)+5];
  char * filename_vals = new char[strlen(filename)+6];
  strcpy(filename_row,filename);
  strcpy(filename_col,filename);
  strcpy(filename_vals,filename);
  strcat(filename_row,"_row");
  strcat(filename_col,"_col");
  strcat(filename_vals,"_vals");
  FILE* RowFile = fopen(filename_row,"r");
  FILE* ColFile = fopen(filename_col,"r");
  FILE* ValsFile = fopen(filename_vals,"r");

  bool read_values = false; 
  if(ValsFile == NULL) 
    read_values = false;

  values = new ScalarType[nnz];
  rowPtr = new OrdinalType[nrows+1];
  colInd = new OrdinalType[nnz];

  fread ( rowPtr, sizeof(OrdinalType), nrows+1, RowFile);
  fread ( colInd, sizeof(OrdinalType), nnz, ColFile);
  
  if(read_values)

  fclose(RowFile);
  fclose(ColFile);
  if(read_values) {
    fread ( values, sizeof(ScalarType), nnz, ValsFile);
    fclose(ValsFile);
  } else {
    for(int i = 0; i<nnz; i++) 
      values[i] = 0.001*(rand()%1000);
  }

  size_t min_span = nrows+1;
  size_t max_span = 0;
  size_t ave_span = 0;
  for(int row=0; row<nrows;row++) {
    int min = nrows+1; int max = 0;
    for(int i=rowPtr[row]; i<rowPtr[row+1]; i++) {
      if(colInd[i]<min) min = colInd[i];
      if(colInd[i]>max) max = colInd[i];
    }
    if(rowPtr[row+1]>rowPtr[row]) {
      size_t span = max-min;
      if(span<min_span) min_span = span;
      if(span>max_span) max_span = span;
      ave_span += span;
    } else min_span = 0;
  }
  printf("%lu Spans: %lu %lu %lu\n",(size_t) nnz,min_span,max_span,ave_span/nrows);


  return nnz;
}

#endif /* MATRIX_MARKET_HPP_ */
