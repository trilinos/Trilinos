// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//
/////////////////////////////////////////////////////////////////////////////
// Written by Seher Acer, 2019
// This code reads a mtx file and writes a binary file in CRS format:
// [nrows nnzs rowPtr[0 nrows] colInd[0 nnzs]
// 1) It symmetrizes the input, regardless of it is already symmetric or not
// 2) It removes the diagonal entries
// 3) The column indices are sorted in increasing order for each row
// 4) Index base for the output file is 0
/////////////////////////////////////////////////////////////////////////////
// The corresponding file reader exists in "readMatrixFromBinaryFile.hpp"
// Note: This is research code. We do not guarantee it works in all cases. 
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <vector>
#include <algorithm>
#include <queue>

typedef long long ord_type;

void findLargestComponent(ord_type seed,
			  ord_type inn, ord_type *inRowPtr, ord_type *inColInd,
			  ord_type &outn, ord_type *outRowPtr, ord_type *outColInd)
{
  
  std::vector<int> mark(inn, -1);
  std::queue<ord_type> active;
  ord_type i, v, lastSearch, root, popCnt = 0, ncc = 0;

  root = seed;
  lastSearch = 0;

  while(popCnt < inn) {
    active.push(root);
    mark[root] = ncc;
  
    while(active.empty() == false) {

      i = active.front();
      active.pop();
      popCnt ++;

      for(ord_type j = inRowPtr[i]; j < inRowPtr[i+1]; ++j) {
        v = inColInd[j];
        if(mark[v] == -1) {
          mark[v] = ncc;
          active.push(v);
        }
      }
    }

    //search for an unmarked vertex to be the next root
    for(ord_type i_r = lastSearch; i_r < inn; ++i_r) {
      if(mark[i_r] == -1) {
        root = i_r;
        lastSearch = i_r;
        break;
      }
    }

    //increase component id
    ncc++;
  }

  // find component sizes
  std::vector<ord_type> compSize(ncc, 0);
  for(ord_type i_c = 0; i_c < inn; ++i_c) {
    if(mark[i_c] == -1) {
      std::cout << "Vertex " << i_c << " is untouched,, Exiting!\n";
      exit(1);
    }
    else
      compSize[mark[i_c]]++;
  }

  // find the largest component
  ord_type maxSize = 0;
  ord_type maxId = -1; 
  for(ord_type i_l = 0; i_l < ncc; ++i_l) {
    if(compSize[i_l] > maxSize) {
      maxSize = compSize[i_l];
      maxId = i_l;
    }
  }

  // renumber the vertices
  outn = 0;
  for(ord_type i_v = 0; i_v < inn; ++i_v) {
    if(mark[i_v] == maxId)
      mark[i_v] = outn++;
    else
      mark[i_v] = -1;
  }

  if(outn != maxSize) {
    std::cout << "The number of vertices in this component: " << outn << " does not match the maximum component size: " << maxSize << "\n";
    exit(1);
  }

  // write the largest component to the output arrays
  ord_type ptr = 0;
  outn = 0;
  outRowPtr[outn] = 0;
  for(ord_type ii = 0; ii < inn; ii++) {
    if(mark[ii] > -1){
      for(ord_type j = inRowPtr[ii]; j < inRowPtr[ii+1]; ++j) {
        v = inColInd[j];
        if(mark[v] == -1) {
          std::cout << "Neighbor " << v << " of " << ii  << " is found to be absent in the component\n";
          exit(1);
        }
        outColInd[ptr++] = mark[v];
      }
      outn++;
      outRowPtr[outn] = ptr;
    }
  }

  if(outn != maxSize) {
    std::cout << "The number of vertices written: " << outn << " does not match the maximum component size: " << maxSize << "\n";
    exit(1);
  }

}

int main(int argc, char* argv[])
{
  
  std::string line;
  ord_type nrows, nnzs, r, c = 0;

  std::ifstream in(argv[1]);

  do{
    getline(in, line);
  }
  while(line[0] == '%');

  std::stringstream ss1(line);
  ss1 >> nrows >> nrows >> nnzs;
  std::cout << "Number of Rows " << nrows << " Number of Columns " << nrows << std::endl;
  std::vector<bool> diag(nrows, false);

  ord_type *rowPtr = new ord_type[nrows+2];
  for(ord_type i = 0; i < nrows+2; i++)
    rowPtr[i] = 0;

  while(getline(in, line)) {

    std::stringstream ss2(line);
    ss2 >> r >> c;
    r--;
    c--;

    if(r != c) {
      rowPtr[r+2]++;
      rowPtr[c+2]++;
    }
  }
  in.close();

  for(ord_type j0 = 0; j0 < nrows; j0++)
    rowPtr[j0+2]++;


  for(ord_type k = 2; k < nrows+2; k++)
    rowPtr[k] += rowPtr[k-1];

  ord_type *colInd = new ord_type[rowPtr[nrows+1]];

  // re-read from the beginning, skip the intro
  in.open(argv[1]);
  do {
    getline(in, line);
  }
  while(line[0] == '%');


  while(getline(in, line)) {

    std::stringstream ss2(line);
    ss2 >> r >> c;
    r--;
    c--;

    if(r != c) {
      colInd[rowPtr[r+1]++] = c;
      colInd[rowPtr[c+1]++] = r;
    }

  }
  in.close();


  for(ord_type i0 = 0; i0 < nrows; i0++)
    colInd[rowPtr[i0+1]++] = i0;


  ord_type *rowPtrNew = new ord_type[nrows+1];
  ord_type *colIndNew = new ord_type[rowPtr[nrows+1]];

  rowPtrNew[0] = 0;
  ord_type ptr = 0;
  ord_type prev = -1;
  for(ord_type i1 = 0; i1 < nrows; i1++) {

    ord_type deg = rowPtr[i1+1] - rowPtr[i1];
    if(deg > 0) {

      std::sort(&colInd[rowPtr[i1]], &colInd[rowPtr[i1+1]]);
      colIndNew[ptr++] = colInd[rowPtr[i1]];
      prev = colInd[rowPtr[i1]];
      for(ord_type j1 = rowPtr[i1]+1; j1 < rowPtr[i1+1]; j1++) {
        if(colInd[j1] != prev) {
          colIndNew[ptr++] = colInd[j1];
          prev = colInd[j1];
        }
      }
    }

    rowPtrNew[i1+1] = ptr;
  }

  for(ord_type i2 = 0; i2 < nrows; i2++) {
    for(ord_type j2 = rowPtrNew[i2]; j2 < rowPtrNew[i2+1]; ++j2)
      if(colIndNew[j2] == i2)
        diag[i2] = true;

    if(diag[i2] == false)
      std::cout << "ROW " << i2 << " misses diagonal\n";
  }

  std::cout << argv[1] << " " << nrows << " " << ptr << " "; 

  ord_type newnrows = -1;
  findLargestComponent(0, //seed
      nrows, rowPtrNew, colIndNew,
      newnrows, rowPtr, colInd );
  ptr = rowPtr[newnrows]; //new number of nonzeros 

  ord_type deg, max = 0;
  for(ord_type i3 = 0; i3 < nrows; i3++) {
    deg = rowPtrNew[i3+1] - rowPtrNew[i3];
    if(deg > max)
      max = deg;
  }

  // write into the output file
  std::ofstream out(argv[2], std::ios::out | std::ios::binary);
  out.write((char*)&newnrows, sizeof(ord_type));
  out.write((char*)&ptr, sizeof(ord_type));

  out.write((char*)rowPtr, sizeof(ord_type)*(newnrows+1));
  out.write((char*)colInd, sizeof(ord_type)*(ptr));

  out.close();

  std::cout << newnrows << " " << ptr << " " << max << "\n";

  delete [] rowPtr;
  delete [] colInd;

  delete [] rowPtrNew;
  delete [] colIndNew;

  return 0;
}
