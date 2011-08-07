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

void pack_outgoing_rows(const Epetra_CrsMatrix& mtx,
                        const std::vector<int>& proc_col_ranges,
                        std::vector<int>& send_rows,
                        std::vector<int>& rows_per_send_proc)
{
  const Epetra_Map& rowmap = mtx.RowMap();
  int numrows = mtx.NumMyRows();
  const Epetra_CrsGraph& graph = mtx.Graph();
  int rowlen = 0;
  int* col_indices = NULL;
  int num_col_ranges = proc_col_ranges.size()/2;
  rows_per_send_proc.resize(num_col_ranges);
  send_rows.clear();
  for(int nc=0; nc<num_col_ranges; ++nc) {
    int first_col = proc_col_ranges[nc*2];
    int last_col = proc_col_ranges[nc*2+1];
    int num_send_rows = 0;
    for(int i=0; i<numrows; ++i) {
      int grow = rowmap.GID(i);
      if (mtx.Filled()) {
        const Epetra_Map& colmap = mtx.ColMap();
        graph.ExtractMyRowView(i, rowlen, col_indices);
        if (rowlen > 0) {
          int begin = colmap.GID(col_indices[0]);
          int end = colmap.GID(col_indices[rowlen-1]);
          if (first_col <= end && last_col >= begin) {
            ++num_send_rows;
            send_rows.push_back(grow);
          }
        }
      }
      else {
        graph.ExtractGlobalRowView(grow, rowlen, col_indices);
        for(int j=0; j<rowlen; ++j) {
          if (first_col <= col_indices[j] && last_col >= col_indices[j]) {
            ++num_send_rows;
            send_rows.push_back(grow);
            break;
          }
        }
      }
    }
    rows_per_send_proc[nc] = num_send_rows;
  }
}

std::pair<int,int> get_col_range(const Epetra_CrsMatrix& mtx)
{
  std::pair<int,int> col_range;
  if (mtx.Filled()) {
    col_range = get_col_range(mtx.ColMap());
  }
  else {
    const Epetra_Map& row_map = mtx.RowMap();
    col_range.first = row_map.MaxMyGID();
    col_range.second = row_map.MinMyGID();
    int rowlen = 0;
    int* col_indices = NULL;
    const Epetra_CrsGraph& graph = mtx.Graph();
    for(int i=0; i<row_map.NumMyElements(); ++i) {
      graph.ExtractGlobalRowView(row_map.GID(i), rowlen, col_indices);
      for(int j=0; j<rowlen; ++j) {
        if (col_indices[j] < col_range.first) col_range.first = col_indices[j];
        if (col_indices[j] > col_range.second) col_range.second = col_indices[j];
      }
    }
  }

  return col_range;
}

}//namespace EpetraExt

