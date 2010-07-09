//@HEADER
// ************************************************************************
// 
//               Tpetra: Templated Linear Algebra Services Package 
//                 Copyright (2008) Sandia Corporation
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
// ************************************************************************
//@HEADER

#include "Teuchos_VerboseObject.hpp"
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_CrsMatrix.hpp"

#ifndef TPETRA_MATRIXMATRIX_DEF_HPP
#define TPETRA_MATRIXMATRIX_DEF_HPP

#ifdef DOXYGEN_USE_ONLY
//#include "Tpetra_MMHelpers_decl.hpp"
#endif

/*! \file Tpetra_MMHelpers_def.hpp 

    The implementations for the MatrixMatrix helper classes. 
 */
namespace Tpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatVec, class SpMatSlv>
CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv>::CrsMatrixStruct()
 : numRows(0), numEntriesPerRow(NULL), indices(NULL), values(NULL),
      remote(NULL), numRemote(0), rowMap(NULL), colMap(NULL),
      domainMap(NULL), importColMap(NULL), importMatrix(NULL)
{
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatVec, class SpMatSlv>
CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv>::~CrsMatrixStruct()
{
  deleteContents();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatVec, class SpMatSlv>
void CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv>::deleteContents()
{
  numRows = 0;
/*  delete [] numEntriesPerRow; numEntriesPerRow = NULL;
  delete [] indices; indices = NULL;
  delete [] values; values = NULL;
  delete [] remote; remote = NULL;*/
  numEntriesPerRow.clear();
  indices.clear();
  values.clear();
  remote.clear();
  numRemote = 0;
  //delete importMatrix;
  importMatrix.reset();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatVec, class SpMatSlv>
int dumpCrsMatrixStruct(const CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv>& M)
{
  std::cout << "proc " << M.rowMap->Comm().MyPID()<<std::endl;
  std::cout << "numRows: " << M.numRows<<std::endl;
  for(int i=0; i<M.numRows; ++i) {
    for(int j=0; j<M.numEntriesPerRow[i]; ++j) {
      if (M.remote[i]) {
        std::cout << "  *"<<M.rowMap->GID(i)<<"   "
             <<M.importColMap->GID(M.indices[i][j])<<"   "<<M.values[i][j]<<std::endl;
      }
      else {
        std::cout << "   "<<M.rowMap->GID(i)<<"   "
             <<M.colMap->GID(M.indices[i][j])<<"   "<<M.values[i][j]<<std::endl;
      }
    }
  }
  return(0);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatVec, class SpMatSlv>
CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv>::CrsWrapper_CrsMatrix(CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv>& crsmatrix)
 : crsmat_(crsmatrix)
{
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatVec, class SpMatSlv>
CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv>::~CrsWrapper_CrsMatrix()
{
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatVec, class SpMatSlv>
const Map<LocalOrdinal, GlobalOrdinal, Node>&
CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv>::getRowMap() const
{
  return crsmat_.getRowMap();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatVec, class SpMatSlv>
bool CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv>::isFillComplete()
{
  return crsmat_.isFillComplete();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatVec, class SpMatSlv>
int
CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv>::insertGlobalValues(GlobalOrdinal globalRow, const Teuchos::ArrayView<const GlobalOrdinal> &indices, const Teuchos::ArrayView<const Scalar> &values)
{
  return crsmat_.insertGlobalValues(globalRow, indices, values);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatVec, class SpMatSlv>
int
CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv>::sumIntoGlobalValues(GlobalOrdinal globalRow, const Teuchos::ArrayView<const GlobalOrdinal> &indices, const Teuchos::ArrayView<const Scalar> &values)
{
  return crsmat_.sumIntoGlobalValues(globalRow, indices, values);
}


//------------------------------------

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
CrsWrapper_GraphBuilder<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CrsWrapper_GraphBuilder(const Map<LocalOrdinal, GlobalOrdinal, Node>& map)
 : graph_(),
   rowmap_(map),
   max_row_length_(0)
{
  //size_t num_rows = map.getNodeNumElements();
  //int* rows = map.MyGlobalElements();
  Teuchos::ArrayView<const GlobalOrdinal> rows= map.getNodeElementList();

  //for(int i=0; i<num_rows; ++i) {
  for(int i=0; i<rows.size(); ++i) {
    graph_[rows[i]] = new std::set<int>;
  }
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
CrsWrapper_GraphBuilder<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~CrsWrapper_GraphBuilder()
{
  typename std::map<GlobalOrdinal,std::set<GlobalOrdinal>*>::iterator
    iter = graph_.begin(), iter_end = graph_.end();
  for(; iter!=iter_end; ++iter) {
    delete iter->second;
  }

  graph_.clear();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool CrsWrapper_GraphBuilder<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isFillComplete()
{
  return false;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int
CrsWrapper_GraphBuilder<Scalar, LocalOrdinal, GlobalOrdinal, Node>::insertGlobalValues(GlobalOrdinal globalRow, const Teuchos::ArrayView<const GlobalOrdinal> &indices, const Teuchos::ArrayView<const Scalar> &values)
{
  typename std::map<GlobalOrdinal,std::set<GlobalOrdinal>*>::iterator
    iter = graph_.find(globalRow);

  if (iter == graph_.end()) return(-1);

  std::set<GlobalOrdinal>& cols = *(iter->second);

  //for(int i=0; i<NumEntries; ++i) {
  for(int i=0; i<indices.size(); ++i) {
    cols.insert(indices[i]);
  }

  global_size_t row_length = cols.size();
  if (row_length > max_row_length_) max_row_length_ = row_length;

  return(0);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int
CrsWrapper_GraphBuilder<Scalar, LocalOrdinal, GlobalOrdinal, Node>::sumIntoGlobalValues(GlobalOrdinal globalRow, const Teuchos::ArrayView<const GlobalOrdinal> &indices, const Teuchos::ArrayView<const Scalar> &values)
{
  return insertGlobalValues(globalRow, indices, values);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::map<GlobalOrdinal,std::set<GlobalOrdinal>*>&
CrsWrapper_GraphBuilder<Scalar, LocalOrdinal, GlobalOrdinal, Node>::get_graph()
{
  return graph_;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatVec, class SpMatSlv>
void insert_matrix_locations(CrsWrapper_GraphBuilder<Scalar, LocalOrdinal, GlobalOrdinal, Node>& graphbuilder,
                              CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv>& C)
{
  global_size_t max_row_length = graphbuilder.get_max_row_length();
  if (max_row_length < 1) return;

  std::vector<GlobalOrdinal> indices(max_row_length);
  //ArrayView<GlobalOrdinal> indices_ptr = &indices[0];
  std::vector<Scalar> zeros(max_row_length, 0.0);
  //ArrayView<Scalar> zeros_ptr = &zeros[0];

  std::map<GlobalOrdinal,std::set<GlobalOrdinal>*>& graph = graphbuilder.get_graph();

  typename std::map<GlobalOrdinal,std::set<GlobalOrdinal>*>::iterator
    iter = graph.begin(), iter_end = graph.end();

  for(; iter!=iter_end; ++iter) {
    GlobalOrdinal row = iter->first;
    std::set<GlobalOrdinal>& cols = *(iter->second);
    global_size_t num_entries = cols.size();

    /*std::set<GlobalOrdinal>::iterator
      col_iter = cols.begin(), col_end = cols.end();
    for(global_size_t j=0; col_iter!=col_end; ++col_iter, ++j) {
      indices_ptr[j] = *col_iter;
    }*/

    C.insertGlobalValues(row, Teuchos::ArrayView<GlobalOrdinal>(indices), Teuchos::ArrayView<Scalar>(zeros));
  }
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_CRSMATRIXSTRUCT_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template class CrsMatrixStruct< SCALAR , LO , GO , NODE >;

#define TPETRA_CRSWRAPPER_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template class CrsWrapper< SCALAR , LO , GO , NODE >;

#define TPETRA_CRSWRAPPER_CRSMATRIX_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template class CrsWrapper_CrsMatrix< SCALAR , LO , GO , NODE >;

#define TPETRA_CRSWRAPPER_GRAPHBUILDER_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template class CrsWrapper_GraphBuilder< SCALAR , LO , GO , NODE >;

}
#endif // TPETRA_MATRIXMATRIX_DEF_HPP
