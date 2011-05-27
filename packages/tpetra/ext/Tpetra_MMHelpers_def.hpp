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

#ifndef TPETRA_MMHELPERS_DEF_HPP
#define TPETRA_MMHELPERS_DEF_HPP

#include "Teuchos_VerboseObject.hpp"
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MMHelpers_decl.hpp"

/*! \file Tpetra_MMHelpers_def.hpp 

    The implementations for the MatrixMatrix helper classes. 
 */
namespace Tpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatOps>
CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::CrsMatrixStruct()
{
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatOps>
CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::~CrsMatrixStruct()
{
  deleteContents();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatOps>
void CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::deleteContents()
{
  numRows = 0;
  numEntriesPerRow.clear();
  indices.clear();
  values.clear();
  remote.clear();
  numRemote = 0;
  importMatrix.reset();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatOps>
int dumpCrsMatrixStruct(const CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& M)
{
  std::cout << "proc " << M.rowMap->Comm().MyPID()<<std::endl;
  std::cout << "numRows: " << M.numRows<<std::endl;
  for(LocalOrdinal i=0; i<M.numRows; ++i) {
    for(LocalOrdinal j=0; j<M.numEntriesPerRow[i]; ++j) {
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

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatOps>
CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::CrsWrapper_CrsMatrix(CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& crsmatrix)
 : crsmat_(crsmatrix)
{
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatOps>
CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::~CrsWrapper_CrsMatrix()
{
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatOps>
Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::getRowMap() const
{
  return crsmat_.getRowMap();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatOps>
bool CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::isFillComplete()
{
  return crsmat_.isFillComplete();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatOps>
void
CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::insertGlobalValues(GlobalOrdinal globalRow, const Teuchos::ArrayView<const GlobalOrdinal> &indices, const Teuchos::ArrayView<const Scalar> &values)
{
  crsmat_.insertGlobalValues(globalRow, indices, values);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatOps>
void 
CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::sumIntoGlobalValues(GlobalOrdinal globalRow, const Teuchos::ArrayView<const GlobalOrdinal> &indices, const Teuchos::ArrayView<const Scalar> &values)
{
  crsmat_.sumIntoGlobalValues(globalRow, indices, values);
}


//------------------------------------

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
CrsWrapper_GraphBuilder<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CrsWrapper_GraphBuilder(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& map)
 : graph_(),
   rowmap_(map),
   max_row_length_(0)
{
  Teuchos::ArrayView<const GlobalOrdinal> rows= map->getNodeElementList();

  for(int i=0; i<rows.size(); ++i) {
    graph_[rows[i]] = new std::set<GlobalOrdinal>;
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
void
CrsWrapper_GraphBuilder<Scalar, LocalOrdinal, GlobalOrdinal, Node>::insertGlobalValues(GlobalOrdinal globalRow, const Teuchos::ArrayView<const GlobalOrdinal> &indices, const Teuchos::ArrayView<const Scalar> &values)
{
  /*std::cout << "inserting" << std::endl;
  std::cout << "  indices: " << indices << std::endl;
  std::cout << "  values: " << values << std::endl;*/
  typename std::map<GlobalOrdinal,std::set<GlobalOrdinal>*>::iterator
    iter = graph_.find(globalRow);

  TEST_FOR_EXCEPTION(iter == graph_.end(), std::runtime_error,
  Teuchos::typeName(*this) << "::insertGlobalValues could not find row " << globalRow << " in the graph. Super bummer man. Hope you figure it out.\n");

  std::set<GlobalOrdinal>& cols = *(iter->second);

  for(int i=0; i<indices.size(); ++i) {
    cols.insert(indices[i]);
  }

  global_size_t row_length = cols.size();
  if (row_length > max_row_length_) max_row_length_ = row_length;

  
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
CrsWrapper_GraphBuilder<Scalar, LocalOrdinal, GlobalOrdinal, Node>::sumIntoGlobalValues(GlobalOrdinal globalRow, const Teuchos::ArrayView<const GlobalOrdinal> &indices, const Teuchos::ArrayView<const Scalar> &values)
{
  insertGlobalValues(globalRow, indices, values);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::map<GlobalOrdinal,std::set<GlobalOrdinal>*>&
CrsWrapper_GraphBuilder<Scalar, LocalOrdinal, GlobalOrdinal, Node>::get_graph()
{
  return graph_;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatOps>
void insert_matrix_locations(CrsWrapper_GraphBuilder<Scalar, LocalOrdinal, GlobalOrdinal, Node> & graphbuilder,
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& C)
{
  global_size_t max_row_length = graphbuilder.get_max_row_length();
  if (max_row_length < 1) return;

  Array<GlobalOrdinal> indices(max_row_length);
  Array<Scalar> zeros(max_row_length, 0.0);

  typedef std::map<GlobalOrdinal,std::set<GlobalOrdinal>*> Graph;
  typedef typename Graph::iterator    GraphIter;
  typedef std::set<GlobalOrdinal>     Set;
  typedef typename Set::iterator      SetIter;
  Graph& graph = graphbuilder.get_graph();

  const GraphIter iter_end = graph.end();
  for(GraphIter iter=graph.begin(); iter!=iter_end; ++iter) {
    const GlobalOrdinal row = iter->first;
    const std::set<GlobalOrdinal> &cols = *(iter->second);
    // "copy" entries out of set into contiguous array storage
    const size_t num_entries = std::copy(cols.begin(), cols.end(), indices.begin()) - indices.begin();
    // insert zeros into the result matrix at the appropriate locations
    C.insertGlobalValues(row, indices(0,num_entries), zeros(0,num_entries));
  }
}

} //namespace Tpetra

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

#endif // TPETRA_MMHELPERS_DEF_HPP
