//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER

#include <EpetraExt_ConfigDefs.h>
#include <EpetraExt_MMHelpers.h>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>

namespace EpetraExt {

CrsMatrixStruct::CrsMatrixStruct()
 : numRows(0), numEntriesPerRow(NULL), indices(NULL), values(NULL),
      remote(NULL), numRemote(0), rowMap(NULL), colMap(NULL),
      domainMap(NULL), importColMap(NULL), importMatrix(NULL)
{
}

CrsMatrixStruct::~CrsMatrixStruct()
{
  deleteContents();
}

void CrsMatrixStruct::deleteContents()
{
  numRows = 0;
  delete [] numEntriesPerRow; numEntriesPerRow = NULL;
  delete [] indices; indices = NULL;
  delete [] values; values = NULL;
  delete [] remote; remote = NULL;
  numRemote = 0;
  delete importMatrix;
}

int dumpCrsMatrixStruct(const CrsMatrixStruct& M)
{
  cout << "proc " << M.rowMap->Comm().MyPID()<<endl;
  cout << "numRows: " << M.numRows<<endl;
  for(int i=0; i<M.numRows; ++i) {
    for(int j=0; j<M.numEntriesPerRow[i]; ++j) {
      if (M.remote[i]) {
        cout << "  *"<<M.rowMap->GID(i)<<"   "
             <<M.importColMap->GID(M.indices[i][j])<<"   "<<M.values[i][j]<<endl;
      }
      else {
        cout << "   "<<M.rowMap->GID(i)<<"   "
             <<M.colMap->GID(M.indices[i][j])<<"   "<<M.values[i][j]<<endl;
      }
    }
  }
  return(0);
}

CrsWrapper_Epetra_CrsMatrix::CrsWrapper_Epetra_CrsMatrix(Epetra_CrsMatrix& epetracrsmatrix)
 : ecrsmat_(epetracrsmatrix)
{
}

CrsWrapper_Epetra_CrsMatrix::~CrsWrapper_Epetra_CrsMatrix()
{
}

const Epetra_Map&
CrsWrapper_Epetra_CrsMatrix::RowMap() const
{
  return ecrsmat_.RowMap();
}

bool CrsWrapper_Epetra_CrsMatrix::Filled()
{
  return ecrsmat_.Filled();
}

int
CrsWrapper_Epetra_CrsMatrix::InsertGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices)
{
  return ecrsmat_.InsertGlobalValues(GlobalRow, NumEntries, Values, Indices);
}

int
CrsWrapper_Epetra_CrsMatrix::SumIntoGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices)
{
  return ecrsmat_.SumIntoGlobalValues(GlobalRow, NumEntries, Values, Indices);
}


//------------------------------------

CrsWrapper_GraphBuilder::CrsWrapper_GraphBuilder(const Epetra_Map& emap)
 : graph_(),
   rowmap_(emap),
   max_row_length_(0)
{
  int num_rows = emap.NumMyElements();
  int* rows = emap.MyGlobalElements();

  for(int i=0; i<num_rows; ++i) {
    graph_[rows[i]] = new std::set<int>;
  }
}

CrsWrapper_GraphBuilder::~CrsWrapper_GraphBuilder()
{
  std::map<int,std::set<int>*>::iterator
    iter = graph_.begin(), iter_end = graph_.end();
  for(; iter!=iter_end; ++iter) {
    delete iter->second;
  }

  graph_.clear();
}

bool CrsWrapper_GraphBuilder::Filled()
{
  return false;
}

int
CrsWrapper_GraphBuilder::InsertGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices)
{
  std::map<int,std::set<int>*>::iterator
    iter = graph_.find(GlobalRow);

  if (iter == graph_.end()) return(-1);

  std::set<int>& cols = *(iter->second);

  for(int i=0; i<NumEntries; ++i) {
    cols.insert(Indices[i]);
  }

  int row_length = cols.size();
  if (row_length > max_row_length_) max_row_length_ = row_length;

  return(0);
}

int
CrsWrapper_GraphBuilder::SumIntoGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices)
{
  return InsertGlobalValues(GlobalRow, NumEntries, Values, Indices);
}

std::map<int,std::set<int>*>&
CrsWrapper_GraphBuilder::get_graph()
{
  return graph_;
}

void insert_matrix_locations(CrsWrapper_GraphBuilder& graphbuilder,
                              Epetra_CrsMatrix& C)
{
  int max_row_length = graphbuilder.get_max_row_length();
  if (max_row_length < 1) return;

  std::vector<int> indices(max_row_length);
  int* indices_ptr = &indices[0];
  std::vector<double> zeros(max_row_length, 0.0);
  double* zeros_ptr = &zeros[0];

  std::map<int,std::set<int>*>& graph = graphbuilder.get_graph();

  std::map<int,std::set<int>*>::iterator
    iter = graph.begin(), iter_end = graph.end();

  for(; iter!=iter_end; ++iter) {
    int row = iter->first;
    std::set<int>& cols = *(iter->second);
    int num_entries = cols.size();

    std::set<int>::iterator
      col_iter = cols.begin(), col_end = cols.end();
    for(int j=0; col_iter!=col_end; ++col_iter, ++j) {
      indices_ptr[j] = *col_iter;
    }

    C.InsertGlobalValues(row, num_entries, zeros_ptr, indices_ptr);
  }
}

}//namespace EpetraExt

