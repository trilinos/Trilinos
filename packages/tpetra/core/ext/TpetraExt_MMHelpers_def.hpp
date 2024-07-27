// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_MMHELPERS_DEF_HPP
#define TPETRA_MMHELPERS_DEF_HPP

#include "TpetraExt_MMHelpers_decl.hpp"
#include "Teuchos_VerboseObject.hpp"

/*! \file TpetraExt_MMHelpers_def.hpp

    The implementations for the MatrixMatrix helper classes.
 */
namespace Tpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CrsMatrixStruct()
{
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~CrsMatrixStruct()
{
  deleteContents();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
deleteContents ()
{
  importMatrix.reset();
  origMatrix = Teuchos::null;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockCrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BlockCrsMatrixStruct(const LocalOrdinal blocksize_)
  : blocksize(blocksize_)
{
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockCrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~BlockCrsMatrixStruct()
{
  deleteContents();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockCrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
deleteContents ()
{
  importMatrix.reset();
  origMatrix = Teuchos::null;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int dumpCrsMatrixStruct (const CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& M)
{
  std::cout << "proc " << M.rowMap->Comm().MyPID()<<std::endl;
  std::cout << "numRows: " << M.numRows<<std::endl;
  for(LocalOrdinal i=0; i<M.numRows; ++i) {
    for(LocalOrdinal j=0; j<M.numEntriesPerRow[i]; ++j) {
      std::cout << "   "<<M.rowMap->GID(i)<<"   "
                <<M.colMap->GID(M.indices[i][j])<<"   "<<M.values[i][j]<<std::endl;
    }
  }

  return 0;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
CrsWrapper_CrsMatrix (CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& crsmatrix)
 : crsmat_ (crsmatrix)
{
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~CrsWrapper_CrsMatrix()
{
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getRowMap() const
{
  return crsmat_.getRowMap();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
isFillComplete ()
{
  return crsmat_.isFillComplete ();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
insertGlobalValues (GlobalOrdinal globalRow,
                    const Teuchos::ArrayView<const GlobalOrdinal> &indices,
                    const Teuchos::ArrayView<const Scalar> &values)
{
  crsmat_.insertGlobalValues (globalRow, indices, values);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
sumIntoGlobalValues (GlobalOrdinal globalRow,
                     const Teuchos::ArrayView<const GlobalOrdinal> &indices,
                     const Teuchos::ArrayView<const Scalar> &values)
{
  crsmat_.sumIntoGlobalValues (globalRow, indices, values);
}



template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
CrsWrapper_GraphBuilder<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
CrsWrapper_GraphBuilder (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& map)
 : graph_(),
   rowmap_(map),
   max_row_length_(0)
{
  Teuchos::ArrayView<const GlobalOrdinal> rows = map->getLocalElementList ();
  const LocalOrdinal numRows = static_cast<LocalOrdinal> (rows.size ());
  for (LocalOrdinal i = 0; i < numRows; ++i) {
    graph_[rows[i]] = new std::set<GlobalOrdinal>;
  }
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
CrsWrapper_GraphBuilder<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
~CrsWrapper_GraphBuilder ()
{
  typename std::map<GlobalOrdinal,std::set<GlobalOrdinal>*>::iterator
    iter = graph_.begin(), iter_end = graph_.end();
  for (; iter != iter_end; ++iter) {
    delete iter->second;
  }
  graph_.clear ();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool CrsWrapper_GraphBuilder<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isFillComplete()
{
  return false;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
CrsWrapper_GraphBuilder<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
insertGlobalValues (GlobalOrdinal globalRow,
                    const Teuchos::ArrayView<const GlobalOrdinal> &indices,
                    const Teuchos::ArrayView<const Scalar> &/* values */)
{
  typename std::map<GlobalOrdinal,std::set<GlobalOrdinal>*>::iterator
    iter = graph_.find (globalRow);

  TEUCHOS_TEST_FOR_EXCEPTION(
    iter == graph_.end(), std::runtime_error,
    "Tpetra::CrsWrapper_GraphBuilder::insertGlobalValues could not find row "
    << globalRow << " in the graph. Super bummer man. Hope you figure it out.");

  std::set<GlobalOrdinal>& cols = * (iter->second);

  for (typename Teuchos::ArrayView<const GlobalOrdinal>::size_type i = 0;
       i < indices.size (); ++i) {
    cols.insert (indices[i]);
  }

  const global_size_t row_length = static_cast<global_size_t> (cols.size ());
  if (row_length > max_row_length_) {
    max_row_length_ = row_length;
  }
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
CrsWrapper_GraphBuilder<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
sumIntoGlobalValues (GlobalOrdinal globalRow,
                     const Teuchos::ArrayView<const GlobalOrdinal> &indices,
                     const Teuchos::ArrayView<const Scalar> &values)
{
  insertGlobalValues (globalRow, indices, values);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::map<GlobalOrdinal,std::set<GlobalOrdinal>*>&
CrsWrapper_GraphBuilder<Scalar, LocalOrdinal, GlobalOrdinal, Node>::get_graph ()
{
  return graph_;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
insert_matrix_locations (CrsWrapper_GraphBuilder<Scalar, LocalOrdinal, GlobalOrdinal, Node>& graphbuilder,
                         CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C)
{
  global_size_t max_row_length = graphbuilder.get_max_row_length();
  if (max_row_length < 1) return;

  Teuchos::Array<GlobalOrdinal> indices(max_row_length);
  Teuchos::Array<Scalar> zeros(max_row_length, Teuchos::ScalarTraits<Scalar>::zero());

  typedef std::map<GlobalOrdinal,std::set<GlobalOrdinal>*> Graph;
  typedef typename Graph::iterator GraphIter;
  Graph& graph = graphbuilder.get_graph ();

  const GraphIter iter_end = graph.end ();
  for (GraphIter iter = graph.begin (); iter != iter_end; ++iter) {
    const GlobalOrdinal row = iter->first;
    const std::set<GlobalOrdinal>& cols = * (iter->second);
    // "copy" entries out of set into contiguous array storage
    const size_t num_entries = std::copy (cols.begin (), cols.end (), indices.begin ()) - indices.begin ();
    // insert zeros into the result matrix at the appropriate locations
    C.insertGlobalValues (row, indices (0, num_entries), zeros (0, num_entries));
  }
}

} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_CRSMATRIXSTRUCT_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template class CrsMatrixStruct< SCALAR , LO , GO , NODE >;

#define TPETRA_BLOCKCRSMATRIXSTRUCT_INSTANT(SCALAR,LO,GO,NODE) \
  \
 template class BlockCrsMatrixStruct< SCALAR , LO , GO , NODE >;

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
