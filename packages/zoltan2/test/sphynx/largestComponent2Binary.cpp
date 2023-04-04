// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact Seher Acer        (sacer@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
// ***********************************************************************
//
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
    for(ord_type i = lastSearch; i < inn; ++i) {
      if(mark[i] == -1) {
        root = i;
        lastSearch = i;
        break;
      }
    }

    //increase component id
    ncc++;
  }

  // find component sizes
  std::vector<ord_type> compSize(ncc, 0);
  for(ord_type i = 0; i < inn; ++i) {
    if(mark[i] == -1) {
      std::cout << "Vertex " << i << " is untouched,, Exiting!\n";
      exit(1);
    }
    else
      compSize[mark[i]]++;
  }

  // find the largest component
  ord_type maxSize = 0;
  ord_type maxId = -1; 
  for(ord_type i = 0; i < ncc; ++i) {
    if(compSize[i] > maxSize) {
      maxSize = compSize[i];
      maxId = i;
    }
  }

  // renumber the vertices
  outn = 0;
  for(ord_type i = 0; i < inn; ++i) {
    if(mark[i] == maxId)
      mark[i] = outn++;
    else
      mark[i] = -1;
  }

  if(outn != maxSize) {
    std::cout << "The number of vertices in this component: " << outn << " does not match the maximum component size: " << maxSize << "\n";
    exit(1);
  }

  // write the largest component to the output arrays
  ord_type ptr = 0;
  outn = 0;
  outRowPtr[outn] = 0;
  for(ord_type i = 0; i < inn; i++) {
    if(mark[i] > -1){
      for(ord_type j = inRowPtr[i]; j < inRowPtr[i+1]; ++j) {
        v = inColInd[j];
        if(mark[v] == -1) {
          std::cout << "Neighbor " << v << " of " << i  << " is found to be absent in the component\n";
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
  ord_type nrows, nnzs, r, c, ndd = 0;

  std::ifstream in(argv[1]);

  do{
    getline(in, line);
  }
  while(line[0] == '%');

  std::stringstream ss1(line);
  ss1 >> nrows >> nrows >> nnzs;

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

  for(ord_type i = 0; i < nrows; i++)
    rowPtr[i+2]++;


  for(ord_type i = 2; i < nrows+2; i++)
    rowPtr[i] += rowPtr[i-1];

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


  for(ord_type i = 0; i < nrows; i++)
    colInd[rowPtr[i+1]++] = i;


  ord_type *rowPtrNew = new ord_type[nrows+1];
  ord_type *colIndNew = new ord_type[rowPtr[nrows+1]];

  rowPtrNew[0] = 0;
  ord_type ptr = 0;
  ord_type prev = -1;
  for(ord_type i = 0; i < nrows; i++) {

    ord_type deg = rowPtr[i+1] - rowPtr[i];
    if(deg > 0) {

      std::sort(&colInd[rowPtr[i]], &colInd[rowPtr[i+1]]);
      colIndNew[ptr++] = colInd[rowPtr[i]];
      prev = colInd[rowPtr[i]];
      for(ord_type j = rowPtr[i]+1; j < rowPtr[i+1]; j++) {
        if(colInd[j] != prev) {
          colIndNew[ptr++] = colInd[j];
          prev = colInd[j];
        }
      }
    }

    rowPtrNew[i+1] = ptr;
  }

  ord_type diagcnt = 0;
  for(ord_type i = 0; i < nrows; i++) {
    for(ord_type j = rowPtrNew[i]; j < rowPtrNew[i+1]; ++j)
      if(colIndNew[j] == i)
        diag[i] = true;

    if(diag[i] == false)
      std::cout << "ROW " << i << " misses diagonal\n";
  }

  std::cout << argv[1] << " " << nrows << " " << ptr << " "; 

  ord_type newnrows = -1;
  findLargestComponent(0, //seed
      nrows, rowPtrNew, colIndNew,
      newnrows, rowPtr, colInd );
  ptr = rowPtr[newnrows]; //new number of nonzeros 

  ord_type deg, max = 0;
  for(ord_type i = 0; i < nrows; i++) {
    deg = rowPtrNew[i+1] - rowPtrNew[i];
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
