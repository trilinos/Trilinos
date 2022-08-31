/*
//@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

#include <fstream>

#include "Teuchos_CommHelpers.hpp"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Trilinos_Util.h"

Epetra_CrsMatrix*
read_matrix_mm(const std::string& mm_file,
               const Epetra_Comm& comm)
{
  int my_proc = comm.MyPID();

  int num_global_rows = 0;
  int nnz_per_row = 0;

  std::ifstream* infile = NULL;
  infile = new std::ifstream(mm_file.c_str());
  if (infile == NULL || !*infile) {
    throw std::runtime_error("Failed to open file "+mm_file);
  }

  //first skip over the file header, which has
  //lines beginning with '%'.
  std::string line;
  do {
    getline(*infile, line);
  } while(line[0] == '%');

  //now get the matrix dimensions.

  int numrows, numcols, nnz;
  std::istringstream isstr(line);
  isstr >> numrows >> numcols >> nnz;

  //make sure we successfully read the three ints from that line.
  if (isstr.fail()) {
    throw std::runtime_error("Failed to parse matrix-market header.");
  }

  if (my_proc == 0) {
    num_global_rows = numrows;
    nnz_per_row = nnz/numrows;
  }

  comm.Broadcast(&num_global_rows, 1, 0);
  comm.Broadcast(&nnz_per_row, 1, 0);

  const int indexBase = 0;
  Epetra_Map rowmap(num_global_rows, indexBase, comm);

  Epetra_CrsMatrix* A = new Epetra_CrsMatrix(Copy, rowmap, nnz_per_row);

  Teuchos::Array<int> col;
  Teuchos::Array<double> coef;

  int irow=0, icol=0;
  int g_row=-1, last_row=-1;
  double val=0;

  while(!infile->eof()) {
    getline(*infile, line);
    std::istringstream isstr(line);
    isstr >> irow >> icol >> val;
  
    if (isstr.fail()) continue;
    if (!rowmap.MyGID(irow-1)) continue;

    g_row = irow-1;
    if (g_row != last_row) {
      if (col.size() > 0) {
        A->InsertGlobalValues(last_row, col.size(), &coef[0], &col[0] );
        col.clear();
        coef.clear();
      }
      last_row = g_row;
    }
    col.push_back(icol-1);
    coef.push_back(val);
  }

  if (col.size() > 0) {
    A->InsertGlobalValues(g_row, col.size(), &coef[0], &col[0]);
  }

  A->FillComplete();
  delete infile;
  return A;
}

Epetra_Vector*
read_vector_mm(const std::string& mm_file,
               const Epetra_Comm& comm)
{
  int my_proc = comm.MyPID();

  int num_global_rows = 0;

  std::ifstream* infile = NULL;
  if (my_proc == 0) {
    infile = new std::ifstream(mm_file.c_str());
    if (infile == NULL || !*infile) {
      throw std::runtime_error("Failed to open file "+mm_file);
    }

    //first skip over the file header, which has
    //lines beginning with '%'.
    std::string line;
    do {
      getline(*infile, line);
    } while(line[0] == '%');

    //now get the matrix dimensions.

    int numrows, numcols;
    std::istringstream isstr(line);
    isstr >> numrows >> numcols;

    //make sure we successfully read the ints from that line.
    if (isstr.fail()) {
      throw std::runtime_error("Failed to parse matrix-market header.");
    }

    num_global_rows = numrows;
  }

  comm.Broadcast(&num_global_rows, 1, 0);

  const int indexBase = 0;
  Epetra_Map rowmap(num_global_rows, indexBase, comm);

  Epetra_Vector* b = new Epetra_Vector(rowmap, 1);

  if (my_proc == 0) {
    int irow=0, icol=0;
    double val=0;

    std::string line;
    std::ifstream& in = *infile;
    while(!infile->eof()) {
      getline(in, line);
      std::istringstream isstr(line);
      isstr >> val;
    
      if (isstr.fail()) continue;

      b->ReplaceGlobalValue(irow++, icol, val);
    }
  }
  delete infile;
  return b;
}

void read_matrix_hb(const std::string& hb_file,
                    const Epetra_Comm& Comm,
                    Epetra_CrsMatrix*& A,
                    Epetra_Vector*& b)
{
  Epetra_Map* Map = NULL;
  Epetra_Vector* x = NULL;
  Epetra_Vector* xexact = NULL;
  Trilinos_Util_ReadHb2Epetra(const_cast<char*>(hb_file.c_str()), Comm, Map,
                             A, x, b, xexact);
  delete Map;
  delete x;
  delete xexact;
}

