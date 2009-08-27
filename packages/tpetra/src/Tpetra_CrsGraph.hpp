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

#ifndef TPETRA_CRSGRAPH_HPP
#define TPETRA_CRSGRAPH_HPP

// TODO: filter column indices first in insertLocalIndices()
// TODO: filter column indices first in insertGlobalIndices()

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_CrsGraph.hpp>

#include <Teuchos_Describable.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_getRawPtr.hpp>
#include <Teuchos_NullIteratorTraits.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_RowGraph.hpp"
#include "Tpetra_Util.hpp"

namespace Tpetra 
{

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration
  template <class S, class LO, class GO, class N>
  class CrsMatrix;
#endif

  //! \brief A class for constructing and using sparse compressed index graphs with row access.
  /*! This class is templated on \c LocalOrdinal and \c GlobalOrdinal. If the \c GlobalOrdinal is not specified, then 
   *  it takes the same type as the \c LocalOrdinal.
   */
  template<class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class CrsGraph : public RowGraph<LocalOrdinal,GlobalOrdinal,Node> {
    template <class S, class LO, class GO, class N, class SpMV>
    friend class CrsMatrix;

    public: 

      //! @name Constructor/Destructor Methods
      //@{ 

      //! Constructor with fixed number of indices per row.
      CrsGraph(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype = DynamicProfile);

      //! Constructor with variable number of indices per row.
      CrsGraph(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const Teuchos::ArrayView<size_t> &NumEntriesPerRowToAlloc, ProfileType pftype = DynamicProfile);

      //! Constructor with fixed number of indices per row and specified column map.
      CrsGraph(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, size_t maxNumEntriesPerRow, ProfileType pftype = DynamicProfile);

      //! Constructor with variable number of indices per row and specified column map.
      CrsGraph(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, const Teuchos::ArrayView<size_t> &NumEntriesPerRowToAlloc, ProfileType pftype = DynamicProfile);

      // !Destructor.
      virtual ~CrsGraph();

      //@}

      //! @name Insertion/Removal Methods
      //@{ 

      //! allows insertion of global indices, either belonging to this node or intended for another.
      void insertGlobalIndices(GlobalOrdinal row, const Teuchos::ArrayView<const GlobalOrdinal> &indices);

      //! allows insertion of local indices intended for this node.
      void insertLocalIndices(LocalOrdinal row, const Teuchos::ArrayView<const LocalOrdinal> &indices);

      //@}

      //! @name Transformational Methods
      //@{ 

      /*! \brief Signal that data entry is complete, specifying domain and range maps. 
          Off-node entries are distributed, repeated entries are summed, and global indices are transformed to local indices.
          If \c OptimizeStorage is true, then optimizeStorage() is called as well.
       */ 
      void fillComplete(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap, OptimizeOption os = DoOptimizeStorage);

      /*! \brief Signal that data entry is complete. 
          Off-node entries are distributed, repeated entries are summed, and global indices are transformed to local indices.
          If \c OptimizeStorage is true, then optimizeStorage() is called as well.
          \note This method calls fillComplete( getRowMap(), getRowMap(), OptimizeStorage ).
       */
      void fillComplete(OptimizeOption os = DoOptimizeStorage);

      //! \brief Communicate non-local contributions to other nodes.
      void globalAssemble();

      //! \brief Re-allocate the data into contiguous storage.
      void optimizeStorage();

      //! Returns \c true if fillComplete() has been called.
      bool isFillComplete() const;

      //! Returns \c true if optimizeStorage() has been called.
      bool isStorageOptimized() const;

      //! Returns \c true if the graph data was allocated in static data structures.
      ProfileType getProfileType() const;

      //@}

      //! @name Graph Query Methods
      //@{ 

      //! Returns the communicator.
      Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

      //! Returns the Map that describes the row distribution in this graph.
      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRowMap() const;

      //! \brief Returns the Map that describes the column distribution in this graph.
      /*! Cannot be called before fillComplete(), unless the matrix was constructed with a column map. */
      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getColMap() const;

      //! Returns the Map associated with the domain of this graph.
      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getDomainMap() const;

      //! Returns the Map associated with the domain of this graph.
      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRangeMap() const;

      //! Returns the importer associated with this graph.
      Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > getImporter() const;

      //! Returns the exporter associated with this graph.
      Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > getExporter() const;

      //! Returns the number of global rows in the graph.
      GlobalOrdinal getGlobalNumRows() const;

      //! \brief Returns the number of global columns in the graph.
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      GlobalOrdinal getGlobalNumCols() const;

      //! Returns the number of rows owned on the calling node.
      size_t getNodeNumRows() const;

      //! Returns the number of columns connected to the locally owned rows of this graph.
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      size_t getNodeNumCols() const;

      //! Returns the index base for global indices for this graph. 
      GlobalOrdinal getIndexBase() const;

      //! \brief Returns the number of entries in the global matrix. 
      /*! Returns the number of global entries in the associated graph. */
      global_size_t getGlobalNumEntries() const;

      //! \brief Returns the number of entries in the calling image's portion of the matrix. 
      /*! Before fillComplete() is called, this could include duplicated entries. */
      size_t getNodeNumEntries() const;

      //! \brief Returns the current number of entries on this node in the specified global row .
      /*! Throws exception std::runtime_error if the specified global row does not belong to this node. */
      size_t getNumEntriesForGlobalRow(GlobalOrdinal globalRow) const;

      //! Returns the current number of entries on this node in the specified local row.
      /*! Throws exception std::runtime_error if the specified local row is not valid for this node. */
      size_t getNumEntriesForLocalRow(LocalOrdinal localRow) const;

      //! \brief Returns the total number of indices allocated for the graph, across all row owned by this processor.
      /*! Returns zero if <tt>getProfileType() == DynamicProfile</tt>. */
      size_t totalAllocation() const;

      //! \brief Returns the current number of allocated entries for this node in the specified global row .
      /*! Throws exception std::runtime_error if the specified global row does not belong to this node. */
      size_t numAllocatedEntriesForGlobalRow(GlobalOrdinal globalRow) const;

      //! Returns the current number of allocated entries on this node in the specified local row.
      /*! Throws exception std::runtime_error if the specified local row is not valid for this node. */
      size_t numAllocatedEntriesForLocalRow(LocalOrdinal localRow) const;

      //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      global_size_t getNumGlobalDiags() const;

      //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      size_t getNumLocalDiags() const;

      //! \brief Returns the maximum number of entries across all rows/columns on all nodes. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      size_t getGlobalMaxNumRowEntries() const;

      //! \brief Returns the maximum number of entries across all rows/columns on this node. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      size_t getLocalMaxNumRowEntries() const;

      //! \brief Indicates whether the matrix has a well-defined column map. 
      /*! The column map does not exist until after fillComplete(), unless the matrix was constructed with one. */
      bool hasColMap() const; 

      //! \brief Indicates whether the graph is lower triangular.
      bool isLowerTriangular() const;

      //! \brief Indicates whether the graph is upper triangular.
      bool isUpperTriangular() const;

      //! \brief If graph indices are in the local range, this function returns true. Otherwise, this function returns false. */
      bool isLocallyIndexed() const;

      //! \brief If graph indices are in the global range, this function returns true. Otherwise, this function returns false. */
      bool isGloballyIndexed() const;

      //@}

      //! @name Extraction Methods
      //@{ 
          
      //! Extract a list of elements in a specified global row of the graph. Put into pre-allocated storage.
      /*!
        \param LocalRow - (In) Global row number for which indices are desired.
        \param Indices - (Out) Global column indices corresponding to values.
        \param NumIndices - (Out) Number of indices.

         Note: A std::runtime_error exception is thrown indices is not large enough to hold the column indices associated
         with row \c GlobalRow. If \c GlobalRow does not belong to this node, then \c indices is unchanged and \c NumIndices is 
         returned as Teuchos::OrdinalTraits<size_t>::invalid().
       */
      void getGlobalRowCopy(GlobalOrdinal GlobalRow, const Teuchos::ArrayView<GlobalOrdinal> &indices, size_t &NumIndices) const;

      //! Extract a list of elements in a specified local row of the graph. Put into storage allocated by calling routine.
      /*!
        \param LocalRow - (In) Local row number for which indices are desired.
        \param Indices - (Out) Local column indices corresponding to values.
        \param NumIndices - (Out) Number of indices.

         Note: A std::runtime_error exception is thrown indices is not large enough to hold the column indices associated
         with row \c LocalRow. If \c LocalRow is not valid for this node, then \c indices is unchanged and \c NumIndices is 
         returned as Teuchos::OrdinalTraits<size_t>::invalid().

        \pre isLocallyIndexed()==true
       */
      void getLocalRowCopy(LocalOrdinal LocalRow, const Teuchos::ArrayView<LocalOrdinal> &indices, size_t &NumIndices) const;

      //! Get a persisting non-const view of the elements in a specified global row of the graph.
      /*!
        \param GlobalRow - (In) Global row number to get indices.

         Note: If \c GlobalRow does not belong to this node, then returns <tt>Teuchos::null</tt>.

        \pre isGloballyIndexed()==true
       */
      Teuchos::ArrayRCP<GlobalOrdinal> getGlobalRowViewNonConst(GlobalOrdinal GlobalRow);

      //! Get a persisting non-const view of the elements in a specified local row of the graph.
      /*!
        \param LocalRow - (In) Local row number to get indices.

         Note: If \c LocalRow is not valid for this node, then returns <tt>Teuchos::null</tt>.

        \pre isLocallyIndexed()==true
       */
      Teuchos::ArrayRCP<LocalOrdinal> getLocalRowViewNonConst(LocalOrdinal LocalRow);

      //! Get a persisting const view of the elements in a specified global row of the graph.
      /*!
        \param GlobalRow - (In) Global row number to get indices.

         Note: If \c GlobalRow does not belong to this node, then returns <tt>Teuchos::null</tt>.

        \pre isGloballyIndexed()==true
       */
      Teuchos::ArrayRCP<const GlobalOrdinal> getGlobalRowView(GlobalOrdinal GlobalRow) const;

      //! Get a persisting const view of the elements in a specified local row of the graph.
      /*!
        \param LocalRow - (In) Local row number to get indices.

         Note: If \c LocalRow is not valid for this node, then returns <tt>Teuchos::null</tt>.

        \pre isLocallyIndexed()==true
       */
      Teuchos::ArrayRCP<const LocalOrdinal> getLocalRowView(LocalOrdinal LocalRow) const;

      //@}

      //! @name I/O Methods
      //@{ 
      
      //! Prints the graph on the specified stream. This is very verbose.
      void print(std::ostream& os) const;

      // @}

    protected:
      void initNumAlloc(const Teuchos::ArrayView<size_t> &NumEntriesPerRowToAlloc);
      void allocateIndices(bool local);
      void makeIndicesLocal(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap);
      void makeColMap(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap);
      void computeIndexState();
      void sortIndices();
      void removeRedundantIndices();
      void makeImportExport();
      bool noRedundancies() const;
      bool indicesAreSorted() const;
      void indicesAreSorted(bool sorted);
      bool indicesAreAllocated() const;
      void staticAssertions();
      size_t findMyIndex(size_t row, LocalOrdinal ind) const;
      size_t findGlobalIndex(size_t row, GlobalOrdinal ind) const;
      inline size_t RNNZ(size_t row) const;
      inline size_t RNumAlloc(size_t row) const;
      inline typename Teuchos::ArrayRCP<const LocalOrdinal>::iterator getLptr(size_t row) const;
      inline typename Teuchos::ArrayRCP<const GlobalOrdinal>::iterator getGptr(size_t row) const;
      inline typename Teuchos::ArrayRCP<LocalOrdinal>::iterator getLptr(size_t row);
      inline typename Teuchos::ArrayRCP<GlobalOrdinal>::iterator getGptr(size_t row);

      // Tpetra support objects
      Teuchos::RCP<const Teuchos::Comm<int> > comm_;
      Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > rowMap_, colMap_, rangeMap_, domainMap_;
      Teuchos::RCP<Import<LocalOrdinal,GlobalOrdinal,Nod> > importer_;
      Teuchos::RCP<Export<LocalOrdinal,GlobalOrdinal,Nod> > exporter_;

      // Local and Global Counts
      global_size_t globalNumEntries_, globalNumDiags_, globalMaxNumRowEntries_;
      size_t          nodeNumEntries_,   nodeNumDiags_,   nodeMaxNumRowEntries_, nodeNumAllocated_;
    
      // 
      ProfileType pftype_;

      // local data. two graphs: one of GIDs and one of LIDs. 
      // except during makeIndicesLocal(), at most one of these is active, as indicated by isLocallyIndexed() and isGloballyIndexed()
      // we provide the buffers to the graphs
      Kokkos::CrsGraph<GlobalOrdinal,Node> gidGraph_;
      Kokkos::CrsGraph< LocalOrdinal,Node> lidGraph_;

      // non-local data
      std::map<GlobalOrdinal, std::deque<GlobalOrdinal> > nonlocals_;
  };

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::CrsGraph(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype)
  : rowMap_(rowMap), lclGraph_(rowMap.getNode()), pftype_(pftype) {
    staticAssertions();
    TEST_FOR_EXCEPTION(maxNumEntriesPerRow < 1 && maxNumEntriesPerRow != 0, std::runtime_error,
        Teuchos::typeName(*this) << "::CrsGraph(rowMap,maxNumEntriesPerRow): maxNumEntriesPerRow must be non-negative.");
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::CrsGraph(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, 
                                                      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, 
                                                      size_t maxNumEntriesPerRow, ProfileType pftype)
  : rowMap_(rowMap), colMap_(colMap), lclGraph_(rowMap.getNode()), pftype_(pftype) {
    staticAssertions();
    TEST_FOR_EXCEPTION(maxNumEntriesPerRow < 1 && maxNumEntriesPerRow != 0, std::runtime_error,
        Teuchos::typeName(*this) << "::CrsGraph(rowMap,colMap,maxNumEntriesPerRow): maxNumEntriesPerRow must be non-negative.");
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::CrsGraph(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, 
                                                      const Teuchos::ArrayView<size_t> &NumEntriesPerRowToAlloc, ProfileType pftype)
  : rowMap_(rowMap), lclGraph_(rowMap.getNode()), pftype_(pftype) {
    staticAssertions();
    TEST_FOR_EXCEPTION(NumEntriesPerRowToAlloc.size() != getNodeNumRows(), std::runtime_error,
        Teuchos::typeName(*this) << "::CrsGraph(rowMap,NumEntriesPerRowToAlloc): NumEntriesPerRowToAlloc must have as many entries as specified by rowMap for this node.");
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::CrsGraph(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, const Teuchos::ArrayView<size_t> &NumEntriesPerRowToAlloc, ProfileType pftype)
  : rowMap_(rowMap), colMap_(colMap), lclGraph_(rowMap.getNode()), pftype_(pftype) {
    staticAssertions();
    TEST_FOR_EXCEPTION(NumEntriesPerRowToAlloc.size() != getNodeNumRows(), std::runtime_error,
        Teuchos::typeName(*this) << "::CrsGraph(rowMap,colMap,NumEntriesPerRowToAlloc): NumEntriesPerRowToAlloc must have as many entries as specified by rowMap for this node.");
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::allocateIndices(bool local) {
    // this is a protected function, only callable by us. if it was called incorrectly, it is our fault.
    TEST_FOR_EXCEPTION(isLocallyIndexed() && local==false, std::logic_error,
        Teuchos::typeName(*this) << "::allocateIndices(): Internal logic error. Please contact Tpetra team.");
    TEST_FOR_EXCEPTION(isGloballyIndexed() && local==true, std::logic_error,
        Teuchos::typeName(*this) << "::allocateIndices(): Internal logic error. Please contact Tpetra team.");
    if (indicesAreAllocated()) {
      return;
    }
    indicesAreLocal_  =  local;
    indicesAreGlobal_ = !local;
    if (getNodeNumRows() > 0) {
      if (getProfileType() == StaticProfile()) {
        typename Teuchos::ArrayRCP<LocalOrdinal>::iterator    
            ntaptr = rowNumToAlloc_.begin(),
            ntaend = rowNumToAlloc_.end();
        size_t sofar = 0, r = 0;
        if (local) {
          colLIndsPtrs_ = Teuchos::arcp<typename Teuchos::ArrayRCP<LocalOrdinal>::iterator>(getNodeNumRows()+1);
          if (totalAllocation()) {
            contigColLInds_ = Teuchos::arcp<LocalOrdinal>(totalAllocation());
            while (ntaptr != ntaend) {
              colLIndsPtrs_[r++] = contigColLInds_.persistingView(sofar,*ntaptr).begin();
              sofar += (*ntaptr++);
            }
            colLIndsPtrs_[r] = contigColLInds_.end();
          }
          else {
            std::fill(colLIndsPtrs_.begin(),colLIndsPtrs_.end(),
                      Teuchos::NullIteratorTraits<typename Teuchos::ArrayRCP<LocalOrdinal>::iterator>::getNull());
          }
        }
        else {
          colGIndsPtrs_ = Teuchos::arcp<typename Teuchos::ArrayRCP<GlobalOrdinal>::iterator>(getNodeNumRows()+1);
          if (totalAllocation()) {
            contigColGInds_ = Teuchos::arcp<GlobalOrdinal>(totalAllocation());
            while (ntaptr != ntaend) {
              colGIndsPtrs_[r++] = contigColGInds_.persistingView(sofar,*ntaptr).begin();
              sofar += (*ntaptr++);
            }
            colGIndsPtrs_[r] = contigColGInds_.end();
          }
          else {
            std::fill(colGIndsPtrs_.begin(),colGIndsPtrs_.end(),
                      Teuchos::NullIteratorTraits<typename Teuchos::ArrayRCP<GlobalOrdinal>::iterator>::getNull());
          }
        }
#ifdef HAVE_TPETRA_DEBUG
        TEST_FOR_EXCEPTION(sofar != totalAllocation(), std::logic_error, 
            Teuchos::typeName(*this) << "::allocateIndices(): Internal logic error. Please contact Tpetra team.");
#endif
      }
      else {
        if (local)  // allocate local indices in colLInds_
        {
          typename Teuchos::ArrayRCP<LocalOrdinal>::iterator 
              ntaptr = rowNumToAlloc_.begin(),
              ntaend = rowNumToAlloc_.end();
          colLInds_ = Teuchos::arcp<Teuchos::ArrayRCP<LocalOrdinal> >(getNodeNumRows());
          typename Teuchos::ArrayRCP<Teuchos::ArrayRCP<LocalOrdinal> >::iterator rinds  = colLInds_.begin();
          for (; ntaptr != ntaend; ++ntaptr, ++rinds) {
            if (*ntaptr > 0) {
              (*rinds) = Teuchos::arcp<LocalOrdinal>(*ntaptr);
            }
          }
        }
        else        // allocate global indices in colGInds_
        {
          typename Teuchos::ArrayRCP<LocalOrdinal>::iterator 
              ntaptr = rowNumToAlloc_.begin(),
              ntaend = rowNumToAlloc_.end();
          colGInds_ = Teuchos::arcp<Teuchos::ArrayRCP<GlobalOrdinal> >(getNodeNumRows());
          typename Teuchos::ArrayRCP<Teuchos::ArrayRCP<GlobalOrdinal> >::iterator rinds  = colGInds_.begin();
          for (; ntaptr != ntaend; ++ntaptr, ++rinds) {
            if (*ntaptr > 0) {
              (*rinds) = Teuchos::arcp<GlobalOrdinal>(*ntaptr);
            }
          }
        }
      }
    }
    rowNumToAlloc_ = Teuchos::null;
    indicesAreAllocated_ = true;    
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::initNumAlloc(const Teuchos::ArrayView<size_t> &NumEntriesPerRowToAlloc)
  {
    typename Teuchos::ArrayView<size_t>::iterator nnz    = NumEntriesPerRowToAlloc.begin(),
                                                  nnzend = NumEntriesPerRowToAlloc.end();
    typename Teuchos::ArrayRCP<LocalOrdinal>::iterator ralloc = rowNumToAlloc_.begin();
    totalNumAllocated_ = 0;
    for (; nnz != nnzend; ++nnz, ++ralloc) {
      if (*nnz > 0) {
        (*ralloc) = *nnz;
        totalNumAllocated_ += *nnz;
      }
      else if (*nnz != 0) {
        TEST_FOR_EXCEPTION(true, std::runtime_error,
            Teuchos::typeName(*this) << "::initNumAlloc(NumEntriesPerRowToAlloc): all entries of NumEntriesPerRowToAlloc must be non-negative.");
      }
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::~CrsGraph()
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  GlobalOrdinal CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getGlobalNumRows() const
  { return rowMap_.getGlobalNumEntries(); }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  GlobalOrdinal CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getGlobalNumCols() const
  {
    TEST_FOR_EXCEPTION(isFillComplete() != true, std::runtime_error,
      Teuchos::typeName(*this) << ": cannot call getGlobalNumCols() until fillComplete() is called.");
    return getDomainMap().getGlobalNumEntries();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNodeNumRows() const
  { return rowMap_.getNumLocalElements(); }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNodeNumCols() const
  {
    TEST_FOR_EXCEPTION(hasColMap() != true, std::runtime_error,
      Teuchos::typeName(*this) << ": cannot call getNodeNumCols() without a column map.");
    return colMap_.getNumLocalElements();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNumLocalDiags() const
  { 
    TEST_FOR_EXCEPTION(!isFillComplete(), std::runtime_error,
      Teuchos::typeName(*this) << ": cannot call getNumLocalDiags() until fillComplete() has been called.");
    return numLocalDiags_; 
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  global_size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNumGlobalDiags() const
  {
    TEST_FOR_EXCEPTION(!isFillComplete(), std::runtime_error,
      Teuchos::typeName(*this) << ": cannot call getNumGlobalDiags() until fillComplete() has been called.");
    return numGlobalDiags_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<LocalOrdinal> CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getLocalRowView(LocalOrdinal LocalRow) const {
    TEST_FOR_EXCEPTION(isGloballyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowView(): local indices do not exist.");
    TEST_FOR_EXCEPTION(getRowMap().isNodeLocalElement(LocalRow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowView(LocalRow): LocalRow (== " << LocalRow << ") is not valid on this node.");
    const size_t rnnz = RNNZ(LocalRow);
    if (indicesAreAllocated() == false || rnnz == 0) {
      indices = Teuchos::ArrayView<const LocalOrdinal>(Teuchos::null);
      return;
    }
    // RNNZ > 0, so there must be something allocated
    if (isStorageOptimized() == true || (getProfileType() == StaticProfile)) {
      indices = Teuchos::arrayView<const LocalOrdinal>(Teuchos::getRawPtr(colLIndsPtrs_[LocalRow]),rnnz);
    }
    else {
      indices = colLInds_[LocalRow](0,rnnz);
    }
    return;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<const GlobalOrdinal> CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getGlobalRowView(GlobalOrdinal GlobalRow) const {
    TEST_FOR_EXCEPTION(isLocallyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalRowView(): global indices do not exist.");
    const size_t lrow = getRowMap().getLocalIndex(GlobalRow);
    TEST_FOR_EXCEPTION(lrow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalRowView(GlobalRow): GlobalRow (== " << GlobalRow << ") does not belong to this node.");
    const size_t rnnz = RNNZ(lrow);
    if (indicesAreAllocated() == false || rnnz == 0) {
      indices = Teuchos::ArrayView<const GlobalOrdinal>(Teuchos::null);
      return;
    }
    // RNNZ > 0, so there must be something allocated
    if (isStorageOptimized() == true || (getProfileType() == StaticProfile)) {
      indices = Teuchos::arrayView<const GlobalOrdinal>(Teuchos::getRawPtr(colGIndsPtrs_[lrow]),rnnz);
    }
    else {
      indices = colGInds_[lrow](0,rnnz);
    }
    return;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<LocalOrdinal> CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getLocalRowViewNonConst(LocalOrdinal LocalRow) { 
    TEST_FOR_EXCEPTION(isGloballyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowViewNonConst(): local indices do not exist.");
    TEST_FOR_EXCEPTION(getRowMap().isNodeLocalElement(LocalRow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowViewNonConst(LocalRow): LocalRow (== " << LocalRow << ") is not valid on this node.");
    const size_t rnnz = RNNZ(LocalRow);
    if (indicesAreAllocated() == false || rnnz == 0) {
      indices = Teuchos::ArrayView<LocalOrdinal>(Teuchos::null);
      return;
    }
    // RNNZ > 0, so there must be something allocated
    if (isStorageOptimized() == true || (getProfileType() == StaticProfile)) {
      indices = Teuchos::arrayView<LocalOrdinal>(Teuchos::getRawPtr(colLIndsPtrs_[LocalRow]),rnnz);
    }
    else {
      indices = colLInds_[LocalRow](0,rnnz);
    }
    return;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<GlobalOrdinal> CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getGlobalRowViewNonConst(GlobalOrdinal GlobalRow) {
    TEST_FOR_EXCEPTION(isLocallyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalRowViewNonConst(): global indices do not exist.");
    const size_t lrow = getRowMap().getLocalIndex(GlobalRow);
    TEST_FOR_EXCEPTION(lrow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalRowViewNonConst(GlobalRow,...): GlobalRow (== " << GlobalRow << ") does not belong to this node.");
    const size_t rnnz = RNNZ(lrow);
    if (indicesAreAllocated() == false || rnnz == 0) {
      indices = Teuchos::ArrayView<GlobalOrdinal>(Teuchos::null);
      return;
    }
    // RNNZ > 0, so there must be something allocated
    if (isStorageOptimized() == true || (getProfileType() == StaticProfile)) {
      indices = Teuchos::arrayView<GlobalOrdinal>(Teuchos::getRawPtr(colGIndsPtrs_[lrow]),rnnz);
    }
    else {
      indices = colGInds_[lrow](0,rnnz);
    }
    return;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getLocalRowCopy(LocalOrdinal LocalRow, const Teuchos::ArrayView<LocalOrdinal> &indices, size_t &NumIndices) const
  {
    // can only do this if we have local indices
    TEST_FOR_EXCEPTION(isGloballyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowCopy(): local indices do not exist.");
    TEST_FOR_EXCEPTION(getRowMap().isNodeLocalElement(LocalRow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowCopy(LocalRow,...): LocalRow (== " << LocalRow << ") is not valid on this node.");
    size_t rnnz = RNNZ(LocalRow);
    TEST_FOR_EXCEPTION(indices.size() < rnnz, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowCopy(): specified storage (size==" << indices.size() 
        << ") is not large enough to hold all entries for this row (rnnz=" << rnnz << ").");
    // use one of the view routines to get the proper view, then copy it over
    NumIndices = rnnz;
    Teuchos::ArrayView<const LocalOrdinal> lview;
    getLocalRowView(LocalRow,lview);
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(lview.size() != rnnz, std::logic_error,
        typeName(*this) << "::getLocalRowCopy(): Internal logic error. Please contact Tpetra team.");
#endif
    std::copy(lview.begin(),lview.end(),indices.begin());
    return;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getGlobalRowCopy(GlobalOrdinal GlobalRow, const Teuchos::ArrayView<GlobalOrdinal> &indices, size_t &NumIndices) const
  {
    using Teuchos::ArrayView;
    // because the data is copied into user storage, we can always do this
    size_t lrow = getRowMap().getLocalIndex(GlobalRow);
    TEST_FOR_EXCEPTION(lrow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalRowCopy(GlobalRow,...): GlobalRow (== " << GlobalRow << ") does not belong to this node.");
    size_t rnnz = RNNZ(lrow);
    TEST_FOR_EXCEPTION(indices.size() < rnnz, std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalRowCopy(): specified storage (size==" << indices.size() 
        << ") is not large enough to hold all entries for this row (rnnz=" << rnnz << ").");
    // use one of the view routines to get the proper view, then copy it over
    NumIndices = rnnz;
    if (isLocallyIndexed()) {
      ArrayView<const LocalOrdinal> lview;
      getLocalRowView(lrow,lview);
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION(lview.size() != rnnz, std::logic_error,
          typeName(*this) << "::getGlobalRowCopy(): Internal logic error. Please contact Tpetra team.");
#endif
      // copy and convert
      typename ArrayView<      GlobalOrdinal>::iterator i=indices.begin();
      typename ArrayView<const LocalOrdinal >::iterator l=lview.begin();
      while (l!=lview.end())
      {
        (*i++) = getColMap().getGlobalIndex(*l++);
      }
    }
    else if (isGloballyIndexed()) {
      ArrayView<const GlobalOrdinal> gview;
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION(gview.size() != rnnz, std::logic_error,
          typeName(*this) << "::getGlobalRowCopy(): Internal logic error. Please contact Tpetra team.");
#endif
      getGlobalRowView(GlobalRow,gview);
      std::copy(gview.begin(), gview.end(), indices.begin());
    }
    return;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  GlobalOrdinal CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getGlobalMaxNumRowEntries() const
  {
    TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalMaxRowEntries(): global quantities not computed until fillComplete() has been called.");
    return globalMaxNumRowEntries_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getRowMap() const
  { return rowMap_; }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getColMap() const
  { 
    TEST_FOR_EXCEPTION(hasColMap() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getColMap(): graph does not have a column map yet; call fillComplete() first.");
    return colMap_; 
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getDomainMap() const
  {
    TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getDomainMap(): graph does not have a domain map yet; call fillComplete() first.");
    return domainMap_; 
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getRangeMap() const
  {
    TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getRangeMap(): graph does not have a range map yet; call fillComplete() first.");
    return rangeMap_; 
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> >
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getImporter() const
  { 
    TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getImporter(): cannot get importer until fillComplete() has been called.");
    return importer_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> >
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getExporter() const
  {
    TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getExporter(): cannot get exporter until fillComplete() has been called.");
    return exporter_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::hasColMap() const
  { return hasColMap_; }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::isStorageOptimized() const
  { return storageOptimized_; }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ProfileType CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getProfileType() const
  { return profile_; }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  GlobalOrdinal CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getGlobalNumEntries() const
  { return numGlobalEntries_; }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::isFillComplete() const
  { return fillComplete_; }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::isUpperTriangular() const
  { return upperTriangular_; }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::isLowerTriangular() const
  { return lowerTriangular_; }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::isLocallyIndexed() const
  { return indicesAreLocal_; }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::isGloballyIndexed() const
  { return indicesAreGlobal_; }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNumEntriesForGlobalRow(GlobalOrdinal globalRow) const
  {
    using Teuchos::OrdinalTraits;
    LocalOrdinal rlid = getRowMap().getLocalIndex(globalRow);
    if (rlid == OrdinalTraits<LocalOrdinal>::invalid()) return OrdinalTraits<size_t>::invalid();
    return RNNZ(rlid);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNumEntriesForLocalRow(LocalOrdinal localRow) const
  {
    using Teuchos::OrdinalTraits;
    if (!getRowMap().isNodeLocalElement(localRow)) return OrdinalTraits<size_t>::invalid();
    return RNNZ(localRow);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::totalAllocation() const
  {
    return totalNumAllocated_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::numAllocatedEntriesForGlobalRow(GlobalOrdinal globalRow) const
  {
    using Teuchos::OrdinalTraits;
    LocalOrdinal rlid = getRowMap().getLocalIndex(globalRow);
    if (rlid == OrdinalTraits<LocalOrdinal>::invalid()) return OrdinalTraits<size_t>::invalid();
    return RNumAlloc(rlid);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::numAllocatedEntriesForLocalRow(LocalOrdinal localRow) const
  {
    using Teuchos::OrdinalTraits;
    if (!getRowMap().isNodeLocalElement(localRow)) return OrdinalTraits<size_t>::invalid();
    return RNumAlloc(localRow);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNodeNumEntries() const
  { return numLocalEntries_; }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getLocalMaxNumRowEntries() const
  { return localMaxNumEntries_; }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::insertLocalIndices(LocalOrdinal lrow, const Teuchos::ArrayView<const LocalOrdinal> &indices)
  {
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(isStorageOptimized() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertLocalIndices(): cannot insert new indices after optimizeStorage() has been called.");
    TEST_FOR_EXCEPTION(isGloballyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertLocalIndices(): graph indices are global; use insertGlobalIndices().");
    TEST_FOR_EXCEPTION(hasColMap() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::insertLocalIndices(): cannot insert local indices without a column map.");
    TEST_FOR_EXCEPTION(getRowMap().isNodeLocalElement(lrow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::insertLocalIndices(): row does not belong to this node.");
    if (indicesAreAllocated() == false) {
      bool local = true;
      allocateIndices(local);
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION(indicesAreAllocated() == false, std::logic_error, 
          Teuchos::typeName(*this) << "::insertLocalIndices(): Internal logic error. Please contact Tpetra team.");
#endif
    }

    indicesAreSorted_ = false;
    noRedundancies_ = false;
    haveGlobalConstants_ = false;

    // add to allocated space
    size_t rowNE = rowNumEntries_[lrow],
           rowNA = RNumAlloc(lrow),
           toAdd = indices.size();
    if (rowNE+toAdd > rowNA) {
      TEST_FOR_EXCEPTION(getProfileType() == StaticProfile, std::runtime_error,
          Teuchos::typeName(*this) << "::insertLocalIndices(): new indices exceed statically allocated graph structure.");
      TPETRA_EFFICIENCY_WARNING(true,std::runtime_error,
          "::insertLocalIndices(): Pre-allocated space has been exceeded, requiring new allocation. To improve efficiency, suggest larger allocation.");
      // increase allocation to necessary amount, copy previous data, reassign ArrayRCP
      ArrayView<const LocalOrdinal> curInds;
      getLocalRowView(lrow,curInds);
      rowNA = rowNE+toAdd;             
      ArrayRCP<LocalOrdinal> newInds = Teuchos::arcp<LocalOrdinal>(rowNA);
      std::copy(curInds.begin(), curInds.end(), newInds.begin());
      colLInds_[lrow] = newInds;
    }
    // check the local indices against the column map; only add ones that are defined
    typename ArrayView<const LocalOrdinal>::iterator srcind = indices.begin();
    typename ArrayRCP<LocalOrdinal>::iterator dstbeg, dstind;
    dstbeg = getLptr(lrow);
    dstind = dstbeg + rowNE;
    while (srcind != indices.end()) 
    {
      if (getColMap().isNodeLocalElement(*srcind)) {
        (*dstind++) = (*srcind);
      }
      ++srcind;
    }
    rowNumEntries_[lrow] = dstind - dstbeg;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::insertGlobalIndices(GlobalOrdinal grow, const Teuchos::ArrayView<const GlobalOrdinal> &indices)
  {
    using Teuchos::OrdinalTraits;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    TEST_FOR_EXCEPTION(isLocallyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertGlobalIndices(): graph indices are local; use insertLocalIndices().");
    if (indicesAreAllocated() == false) {
      bool local = false;
      allocateIndices(local);
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION(indicesAreAllocated() == false, std::logic_error, 
          Teuchos::typeName(*this) << "::insertGlobalIndices(): Internal logic error. Please contact Tpetra team.");
#endif
    }

    indicesAreSorted_ = false;
    noRedundancies_ = false;
    haveGlobalConstants_ = false;

    LocalOrdinal lrow = getRowMap().getLocalIndex(grow);
    if (lrow != OrdinalTraits<LocalOrdinal>::invalid()) {
      // add to allocated space
      size_t rowNE = rowNumEntries_[lrow],
             rowNA = RNumAlloc(lrow),
             toAdd = indices.size();
      if (rowNE+toAdd > rowNA) {
        TEST_FOR_EXCEPTION(getProfileType() == StaticProfile, std::runtime_error,
            Teuchos::typeName(*this) << "::insertGlobalIndices(): new indices exceed statically allocated graph structure.");
        TPETRA_EFFICIENCY_WARNING(true,std::runtime_error,
            "::insertGlobalIndices(): Pre-allocated space has been exceeded, requiring new allocation. To improve efficiency, suggest larger allocation.");
        // increase allocation to necessary amount, copy previous data, reassign ArrayRCP
        ArrayView<const GlobalOrdinal> curInds;
        getGlobalRowView(grow,curInds);
        rowNA = rowNE+toAdd;
        ArrayRCP<GlobalOrdinal> newInds = Teuchos::arcp<GlobalOrdinal>(rowNA);
        std::copy(curInds.begin(), curInds.end(), newInds.begin());
        colGInds_[lrow] = newInds;
      }
      typename ArrayRCP<GlobalOrdinal>::iterator dstbeg, dstind;
      dstbeg = getGptr(lrow);
      dstind = dstbeg + rowNE;
      if (hasColMap()) {
        // check the global indices against the column map; only add ones that are defined
        typename ArrayView<const GlobalOrdinal>::iterator srcind = indices.begin();
        while (srcind != indices.end()) 
        {
          if (getColMap().isNodeGlobalElement(*srcind)) {
            (*dstind++) = (*srcind);
          }
          ++srcind;
        }
        rowNumEntries_[lrow] = dstind - dstbeg;
      }
      else {
        std::copy( indices.begin(), indices.end(), dstind );
        rowNumEntries_[lrow] += toAdd;
      }
    }
    else {
      // a nonlocal
      for (typename Teuchos::ArrayView<const GlobalOrdinal>::iterator i=indices.begin(); i != indices.end(); ++i) 
      {
        nonlocals_[grow].push_back(*i);
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Teuchos::Comm<int> > 
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getComm() const
  { return comm_; }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  GlobalOrdinal CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getIndexBase() const
  { return rowMap_.getIndexBase(); }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::indicesAreSorted() const
  { return indicesAreSorted_; }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::noRedundancies() const
  { return noRedundancies_; }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::indicesAreSorted(bool sorted) 
  { indicesAreSorted_ = sorted; }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::fillComplete(OptimizeOption os) {
    fillComplete(getRowMap(),getRowMap(),os);
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::optimizeStorage()
  {
    // optimizeStorage will perform two functions:
    // 1) create a single allocation of memory
    // 2) pack data in that allocation
    // if getProfileType() == StaticProfile, then 1) has already been done
    if (isStorageOptimized() == true) return;
    TEST_FOR_EXCEPTION(isFillComplete() == false || indicesAreSorted() == false || noRedundancies() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::optimizeStorage(): fillComplete() must be called before optimizeStorage().");

    // 1) allocate single memory block
    const size_t nlrs = getNodeNumRows();
    if (nlrs > 0) {
      if (getProfileType() == DynamicProfile) {
        // need to allocate storage, create pointers, copy data, and delete old data
        const size_t nE = getNodeNumEntries();
        colLIndsPtrs_ = Teuchos::arcp<typename Teuchos::ArrayRCP<LocalOrdinal>::iterator>(nlrs+1);
        if (nE > 0) {
          contigColLInds_ = Teuchos::arcp<LocalOrdinal>(nE);
          size_t sofar = 0;
          for (size_t r=0; r<nlrs; ++r) {
            size_t rne = rowNumEntries_[r];
            colLIndsPtrs_[r] = contigColLInds_.persistingView(sofar,rne).begin();
            std::copy(colLInds_[r].begin(),colLInds_[r].begin()+rne,colLIndsPtrs_[r]);
            colLInds_[r] = Teuchos::null;
            sofar += rne;
          }
          colLIndsPtrs_[nlrs] = contigColLInds_.end();
        }
        else {
          std::fill(colLIndsPtrs_.begin(),colLIndsPtrs_.end(),
                    Teuchos::NullIteratorTraits<typename Teuchos::ArrayRCP<LocalOrdinal>::iterator>::getNull());
        }
        colLInds_ = Teuchos::null;
      }
      else {
        // storage is already allocated and pointers are set; just need to pack
        // need to change pointers and move data
        // remember, these aren't just pointers, but also bounds-checked 
        const size_t nE = getNodeNumEntries();
        if (nE > 0) {
          size_t sofar = 0;
          typename Teuchos::ArrayRCP<LocalOrdinal>::iterator newptr, oldptr;
          for (size_t r=0; r<nlrs; ++r) {
            size_t rne = rowNumEntries_[r];
            newptr = contigColLInds_.persistingView(sofar,rne).begin();
            oldptr = colLIndsPtrs_[r];
            if (newptr != oldptr) {
              for (size_t j=0; j<rne; ++j) {
                newptr[j] = oldptr[j];
              }
              colLIndsPtrs_[r] = newptr;
            }
            sofar += rne;
          }
          colLIndsPtrs_[nlrs] = contigColLInds_.persistingView(sofar,0).begin();
        }
      }
    }
    rowNumEntries_ = Teuchos::null;
    totalNumAllocated_ = getNodeNumEntries();
    storageOptimized_ = true;
#ifdef HAVE_TPETRA_DEBUG 
    // check that we only have the data structures that we need
    TEST_FOR_EXCEPTION( rowNumEntries_ != Teuchos::null, std::logic_error,
        Teuchos::typeName(*this) << "::optimizeStorage(): Internal logic error. Please contact Tpetra team.");
    TEST_FOR_EXCEPTION( colLInds_ != Teuchos::null, std::logic_error,
        Teuchos::typeName(*this) << "::optimizeStorage(): Internal logic error. Please contact Tpetra team.");
#endif
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::globalAssemble()
  {
    using Teuchos::arcp;
    using Teuchos::Array;
    using Teuchos::SerialDenseMatrix;
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::deque;
    using std::pair;
    using std::make_pair;
    typedef typename std::map<GlobalOrdinal,std::deque<GlobalOrdinal> >::const_iterator NLITER;
    int numImages = Teuchos::size(*getComm());
    int myImageID = Teuchos::rank(*getComm());
    // Determine if any nodes have global entries to share
    size_t MyNonlocals = nonlocals_.size(), 
           MaxGlobalNonlocals;
    Teuchos::reduceAll<size_t>(*getComm(),Teuchos::REDUCE_MAX,MyNonlocals,&MaxGlobalNonlocals);
    if (MaxGlobalNonlocals == 0) return;  // no entries to share

    // compute a list of NLRs from nonlocals_ and use it to compute:
    //      IdsAndRows: a vector of (id,row) pairs
    //          NLR2Id: a map from NLR to the Id that owns it
    // globalNeighbors: a global graph of connectivity between images: globalNeighbors(i,j) indicates that j sends to i
    //         sendIDs: a list of all images I send to
    //         recvIDs: a list of all images I receive from (constructed later)
    Array<pair<int,GlobalOrdinal> > IdsAndRows;
    std::map<GlobalOrdinal,int> NLR2Id;
    SerialDenseMatrix<int,char> globalNeighbors;
    Array<int> sendIDs, recvIDs;
    {
      // nonlocals_ contains the entries we are holding for all non-local rows
      // we want a list of the rows for which we have data
      Array<GlobalOrdinal> NLRs;
      std::set<GlobalOrdinal> setOfRows;
      for (NLITER iter = nonlocals_.begin(); iter != nonlocals_.end(); ++iter)
      {
        setOfRows.insert(iter->first);
      }
      // copy the elements in the set into an Array
      NLRs.resize(setOfRows.size());
      std::copy(setOfRows.begin(), setOfRows.end(), NLRs.begin());

      // get a list of ImageIDs for the non-local rows (NLRs)
      Array<int> NLRIds(NLRs.size());
      {
        bool invalidGIDs = getRowMap().getRemoteIndexList(NLRs(),NLRIds());
        char lclerror = ( invalidGIDs ? 1 : 0 );
        char gblerror;
        Teuchos::reduceAll(*getComm(),Teuchos::REDUCE_MAX,lclerror,&gblerror);
        TEST_FOR_EXCEPTION(gblerror, std::runtime_error,
            Teuchos::typeName(*this) << "::globalAssemble(): non-local entries correspond to invalid rows.");
      }

      // build up a list of neighbors, as well as a map between NLRs and Ids
      // localNeighbors[i] != 0 iff I have data to send to image i
      // put NLRs,Ids into an array of pairs
      IdsAndRows.reserve(NLRs.size());
      Array<char> localNeighbors(numImages,0);
      typename Array<GlobalOrdinal>::const_iterator nlr;
      typename Array<int>::const_iterator id;
      for (nlr = NLRs.begin(), id = NLRIds.begin();
           nlr != NLRs.end(); ++nlr, ++id) 
      {
        NLR2Id[*nlr] = *id;
        localNeighbors[*id] = 1;
        IdsAndRows.push_back(make_pair<int,GlobalOrdinal>(*id,*nlr));
      }
      for (int j=0; j<numImages; ++j)
      {
        if (localNeighbors[j]) {
          sendIDs.push_back(j);
        }
      }
      // sort IdsAndRows, by Ids first, then rows
      std::sort(IdsAndRows.begin(),IdsAndRows.end());
      // gather from other nodes to form the full graph
      globalNeighbors.shapeUninitialized(numImages,numImages);
      Teuchos::gatherAll(*getComm(),numImages,localNeighbors.getRawPtr(),numImages*numImages,globalNeighbors.values());
      // globalNeighbors at this point contains (on all images) the
      // connectivity between the images. 
      // globalNeighbors(i,j) != 0 means that j sends to i/that i receives from j
    }

    ////////////////////////////////////////////////////////////////////////////////////// 
    // FIGURE OUT WHO IS SENDING TO WHOM AND HOW MUCH
    // DO THIS IN THE PROCESS OF PACKING ALL OUTGOING DATA ACCORDING TO DESTINATION ID
    ////////////////////////////////////////////////////////////////////////////////////// 

    // loop over all columns to know from which images I can expect to receive something
    for (int j=0; j<numImages; ++j)
    {
      if (globalNeighbors(myImageID,j)) {
        recvIDs.push_back(j);
      }
    }
    size_t numRecvs = recvIDs.size();

    // we know how many we're sending to already
    // form a contiguous list of all data to be sent
    // track the number of entries for each ID
    Array<pair<GlobalOrdinal,GlobalOrdinal> > IJSendBuffer;
    Array<size_t> sendSizes(sendIDs.size(), 0);
    size_t numSends = 0;
    for (typename Array<pair<int,GlobalOrdinal> >::const_iterator IdAndRow = IdsAndRows.begin();
         IdAndRow != IdsAndRows.end(); ++IdAndRow) 
    {
      int            id = IdAndRow->first;
      GlobalOrdinal row = IdAndRow->second;
      // have we advanced to a new send?
      if (sendIDs[numSends] != id) {
        numSends++;
        TEST_FOR_EXCEPTION(sendIDs[numSends] != id, std::logic_error, Teuchos::typeName(*this) << "::globalAssemble(): internal logic error. Contact Tpetra team.");
      }
      // copy data for row into contiguous storage
      for (typename deque<GlobalOrdinal>::const_iterator j = nonlocals_[row].begin(); j != nonlocals_[row].end(); ++j)
      {
        IJSendBuffer.push_back( pair<GlobalOrdinal,GlobalOrdinal>(row,*j) );
        sendSizes[numSends]++;
      }
    }
    if (IdsAndRows.size() > 0) {
      numSends++; // one last increment, to make it a count instead of an index
    }
    TEST_FOR_EXCEPTION(Teuchos::as<typename Array<int>::size_type>(numSends) != sendIDs.size(), std::logic_error, Teuchos::typeName(*this) << "::globalAssemble(): internal logic error. Contact Tpetra team.");

    // don't need this data anymore
    nonlocals_.clear();

    ////////////////////////////////////////////////////////////////////////////////////// 
    // TRANSMIT SIZE INFO BETWEEN SENDERS AND RECEIVERS
    ////////////////////////////////////////////////////////////////////////////////////// 
    // perform non-blocking sends: send sizes to our recipients
    Array<RCP<Teuchos::CommRequest> > sendRequests;
    for (int s=0; s < numSends ; ++s)
    {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      sendRequests.push_back( Teuchos::isend<int,int>(*getComm(),rcp(&sendSizes[s],false),sendIDs[s]) );
    }
    // perform non-blocking receives: receive sizes from our senders
    Array<RCP<Teuchos::CommRequest> > recvRequests;
    Array<size_t> recvSizes(numRecvs);
    for (int r=0; r < numRecvs; ++r) {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      recvRequests.push_back( Teuchos::ireceive(*getComm(),rcp(&recvSizes[r],false),recvIDs[r]) );
    }
    // wait on all 
    if (!sendRequests.empty()) {
      Teuchos::waitAll(*getComm(),sendRequests());
    }
    if (!recvRequests.empty()) {
      Teuchos::waitAll(*getComm(),recvRequests());
    }
    Teuchos::barrier(*getComm());
    sendRequests.clear();
    recvRequests.clear();

    ////////////////////////////////////////////////////////////////////////////////////
    // NOW SEND/RECEIVE ALL ROW DATA
    ////////////////////////////////////////////////////////////////////////////////////
    // from the size info, build the ArrayViews into IJSendBuffer
    Array<ArrayView<pair<GlobalOrdinal,GlobalOrdinal> > > sendBuffers(numSends,Teuchos::null);
    {
      size_t cur = 0;
      for (size_t s=0; s<numSends; ++s) {
        sendBuffers[s] = IJSendBuffer(cur,sendSizes[s]);
        cur += sendSizes[s];
      }
    }
    // perform non-blocking sends
    for (size_t s=0; s < numSends ; ++s)
    {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      ArrayRCP<pair<GlobalOrdinal,GlobalOrdinal> > tmparcp = arcp(sendBuffers[s].getRawPtr(),0,sendBuffers[s].size(),false);
      sendRequests.push_back( Teuchos::isend<int,pair<GlobalOrdinal,GlobalOrdinal> >(*getComm(),tmparcp,sendIDs[s]) );
    }
    // calculate amount of storage needed for receives
    // setup pointers for the receives as well
    size_t totalRecvSize = std::accumulate(recvSizes.begin(),recvSizes.end(),0);
    Array<pair<GlobalOrdinal,GlobalOrdinal> > IJRecvBuffer(totalRecvSize);
    // from the size info, build the ArrayViews into IJRecvBuffer
    Array<ArrayView<pair<GlobalOrdinal,GlobalOrdinal> > > recvBuffers(numRecvs,Teuchos::null);
    {
      size_t cur = 0;
      for (size_t r=0; r<numRecvs; ++r) {
        recvBuffers[r] = IJRecvBuffer(cur,recvSizes[r]);
        cur += recvSizes[r];
      }
    }
    // perform non-blocking recvs
    for (size_t r=0; r < numRecvs ; ++r)
    {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      ArrayRCP<pair<GlobalOrdinal,GlobalOrdinal> > tmparcp = arcp(recvBuffers[r].getRawPtr(),0,recvBuffers[r].size(),false);
      recvRequests.push_back( Teuchos::ireceive(*getComm(),tmparcp,recvIDs[r]) );
    }
    // perform waits
    if (!sendRequests.empty()) {
      Teuchos::waitAll(*getComm(),sendRequests());
    }
    if (!recvRequests.empty()) {
      Teuchos::waitAll(*getComm(),recvRequests());
    }
    Teuchos::barrier(*getComm());
    sendRequests.clear();
    recvRequests.clear();


    ////////////////////////////////////////////////////////////////////////////////////
    // NOW PROCESS THE RECEIVED ROW DATA
    ////////////////////////////////////////////////////////////////////////////////////
    // TODO: instead of adding one entry at a time, add one row at a time.
    //       this requires resorting; they arrived sorted by sending node, so that entries could be non-contiguous if we received
    //       multiple entries for a particular row from different processors.
    //       it also requires restoring the data, which may make it not worth the trouble.
    for (typename Array<pair<GlobalOrdinal,GlobalOrdinal> >::const_iterator ij = IJRecvBuffer.begin(); ij != IJRecvBuffer.end(); ++ij)
    {
      insertGlobalIndices(ij->first, Teuchos::arrayView(&ij->second,1));
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::fillComplete(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, 
                                                               const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap, 
                                                               OptimizeOption os) {
    using Teuchos::Array;
    using Teuchos::ArrayRCP;

    Teuchos::ArrayView<const GlobalOrdinal> myGlobalEntries = rowMap_.getMyGlobalEntries();
    domainMap_ = domainMap;
    rangeMap_  = rangeMap;

    if (Teuchos::size(*getComm()) > 1) {
      globalAssemble();
    }
    else {
      TEST_FOR_EXCEPTION(nonlocals_.size() > 0, std::runtime_error,
          Teuchos::typeName(*this) << "::fillComplete(): cannot have non-local entries on a serial run. Invalid entries were submitted to the CrsMatrix.");
    }

    makeIndicesLocal(domainMap, rangeMap);  // create column map, allocate/reapportion space for local indices, transform global indices to local indices
    sortIndices();                          // sort local indices
    removeRedundantIndices();               // remove redundant indices, count non-zero entries, count diagonals, check triangularity
    makeImportExport();                     // create import and export object

    // compute global constants using computed local constants
    if (haveGlobalConstants_ == false) {
      GlobalOrdinal lcl[2], gbl[2];
      lcl[0] = numLocalEntries_;
      lcl[1] = numLocalDiags_;
      Teuchos::reduceAll<int,GlobalOrdinal>(*getComm(),Teuchos::REDUCE_SUM,2,lcl,gbl);
      numGlobalEntries_ = gbl[0]; 
      numGlobalDiags_   = gbl[1];
      Teuchos::reduceAll<int,GlobalOrdinal>(*getComm(),Teuchos::REDUCE_MAX,localMaxNumEntries_,&globalMaxNumRowEntries_);
      haveGlobalConstants_ = true;
    }

    // mark transformation as successfully completed
    fillComplete_ = true;

#ifdef HAVE_TPETRA_DEBUG 
    // check that we only have the data structures that we need
    TEST_FOR_EXCEPTION( contigColGInds_ != Teuchos::null, std::logic_error,
        Teuchos::typeName(*this) << "::fillComplete(): Internal logic error. Please contact Tpetra team.");
    TEST_FOR_EXCEPTION( colGInds_ != Teuchos::null, std::logic_error,
        Teuchos::typeName(*this) << "::fillComplete(): Internal logic error. Please contact Tpetra team.");
    TEST_FOR_EXCEPTION( indicesAreAllocated() && rowNumToAlloc_ != Teuchos::null, std::logic_error,
        Teuchos::typeName(*this) << "::fillComplete(): Internal logic error. Please contact Tpetra team.");
    TEST_FOR_EXCEPTION( colGIndsPtrs_ != Teuchos::null, std::logic_error,
        Teuchos::typeName(*this) << "::fillComplete(): Internal logic error. Please contact Tpetra team.");
#endif
    if (os == OptimizeStorage) optimizeStorage();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::makeIndicesLocal(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap)
  {
    using Teuchos::ArrayRCP;
    using Teuchos::NullIteratorTraits;
    computeIndexState(); // Update index state by checking isLocallyIndexed/Global on all nodes
    TEST_FOR_EXCEPTION(isLocallyIndexed() && isGloballyIndexed(), std::logic_error,
        Teuchos::typeName(*this) << "::makeIndicesLocal(): indices are marked as both global and local.");

    makeColMap(domainMap, rangeMap); // If user has not prescribed column map, create one from indices

    // Transform indices to local index space
    const size_t nlrs = getNodeNumRows();

    if (isGloballyIndexed() && indicesAreAllocated() && nlrs > 0) {
      // allocate data for local indices
      if (getProfileType() == StaticProfile) {
        colLIndsPtrs_ = Teuchos::arcp<typename Teuchos::ArrayRCP<LocalOrdinal>::iterator>(getNodeNumRows()+1);
        if (contigColGInds_ != Teuchos::null) {
          contigColLInds_ = Teuchos::arcp_reinterpret_cast<LocalOrdinal>(contigColGInds_);
          size_t sofar = 0, r = 0;
          while (r < getNodeNumRows()) {
            // RNumAlloc(r) uses colGIndsPtrs_[r+1] and colGIndsPtrs_[r], can't erase till done
            size_t numalloc = RNumAlloc(r),
                   numentry = RNNZ(r);
            colLIndsPtrs_[r] = contigColLInds_.persistingView(sofar,numalloc).begin();
            for (size_t j=0; j<numentry; ++j) {
              colLIndsPtrs_[r][j] = getColMap().getLocalIndex( colGIndsPtrs_[r][j] );
#ifdef HAVE_TPETRA_DEBUG
              TEST_FOR_EXCEPTION(colLIndsPtrs_[r][j] == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::logic_error,
                  Teuchos::typeName(*this) << ": Internal error in fillComplete(). Please contact Tpetra team.");
#endif
            }
            colGIndsPtrs_[r] = NullIteratorTraits<typename ArrayRCP<GlobalOrdinal>::iterator>::getNull();
            ++r;
            sofar += numalloc;
          }
          colLIndsPtrs_[r] = contigColLInds_.end();
#ifdef HAVE_TPETRA_DEBUG
          TEST_FOR_EXCEPTION(sofar != contigColGInds_.size(), std::logic_error,
              Teuchos::typeName(*this) << "::makeIndicesLocal(): Internal logic error. Please contact Tpetra team.");
#endif
        }
        else {
          std::fill(colLIndsPtrs_.begin(),colLIndsPtrs_.end(),
                    Teuchos::NullIteratorTraits<typename Teuchos::ArrayRCP<LocalOrdinal>::iterator>::getNull());
        }
        colGIndsPtrs_ = Teuchos::null;
        contigColGInds_ = Teuchos::null;
      }
      else {
        colLInds_ = Teuchos::arcp<Teuchos::ArrayRCP<LocalOrdinal> >(getNodeNumRows());
        for (size_t r=0; r < nlrs; ++r)
        {
          colLInds_[r] = Teuchos::arcp_reinterpret_cast<LocalOrdinal>(colGInds_[r]);
          typename ArrayRCP<GlobalOrdinal>::const_iterator 
            cindG = colGInds_[r].begin(),
            cendG = colGInds_[r].begin() + rowNumEntries_[r];
          typename ArrayRCP<LocalOrdinal>::iterator cindL = colLInds_[r].begin();
          for (; cindG != cendG; ++cindG, ++cindL)
          {
            (*cindL) = colMap_.getLocalIndex(*cindG);
#ifdef HAVE_TPETRA_DEBUG
            TEST_FOR_EXCEPTION((*cindL) == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::logic_error,
                Teuchos::typeName(*this) << ": Internal error in fillComplete(). Please contact Tpetra team.");
#endif
          }
          colGInds_[r] = Teuchos::null;   // data is owned by colLInds_[r] now
        }
        colGInds_ = Teuchos::null;
      }
    }
    indicesAreLocal_  = true;
    indicesAreGlobal_ = false;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::sortIndices()
  {
    using Teuchos::ArrayView;
    TEST_FOR_EXCEPT(isGloballyIndexed()==true);   // this should be called only after makeIndicesLocal()
    if (indicesAreSorted()) return;
    if (indicesAreAllocated()) {                 // are there any indices?
      const size_t nlrs = getNodeNumRows();
      for (size_t r=0; r < nlrs; ++r)
      {
        size_t rnnz = RNNZ(r);
        typename Teuchos::ArrayRCP<LocalOrdinal>::iterator rowbeg = getLptr(r);
        std::sort(rowbeg, rowbeg + rnnz);
      }
    }
    indicesAreSorted_ = true;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::makeColMap(
                                      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, 
                                      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap) {
    using Teuchos::ArrayRCP;
    using Teuchos::Array;
    using Teuchos::ArrayView;
    typedef Teuchos::OrdinalTraits<GlobalOrdinal> GOT;
    const size_t nlrs = getNodeNumRows();
    if (hasColMap()) return;
    computeIndexState();
    TEST_FOR_EXCEPTION(isLocallyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::makeColMap(): indices must still be global when making the column map.");

    // ultimate: list of indices for column map
    Array<GlobalOrdinal> myColumns;
    if (isGloballyIndexed() == true) { // otherwise, we have no indices to bother with
      // Construct two lists of columns for which we have a non-zero
      // Local GIDs: Column GIDs which are locally present on the Domain map
      // Remote GIDs: Column GIDs which are not locally present present on the Domain map
      //
      // instead of list of local GIDs, we will use a vector<bool> LocalGID, where LocalGID[LID]
      // indicates that, for LID=DomainMap.LID(GID), GID is local. on modern compilers, vector<bool>
      // should be a partial specialization of vector<> which represents each entry via a single
      // bit, instead of a full byte.
      const LocalOrdinal LINV = Teuchos::OrdinalTraits<LocalOrdinal>::invalid();
      size_t numLocalColGIDs = 0, numRemoteColGIDs = 0;
      //
      // intitial: partitioning into local and remote
      std::vector<bool>       GIDisLocal(domainMap.getNumLocalElements());
      std::set<GlobalOrdinal> RemoteGIDSet;
      for (size_t r=0; r < nlrs; ++r)
      {
        typename ArrayRCP<GlobalOrdinal>::const_iterator cind, cend;
        cind = getGptr(r);
        cend = cind + rowNumEntries_[r];
        for (; cind != cend; ++cind) {
          GlobalOrdinal gid = (*cind);
          LocalOrdinal lid = domainMap.getLocalIndex(gid);
          if (lid != LINV) {
            bool alreadyFound = GIDisLocal[lid];
            if (!alreadyFound) {
              GIDisLocal[lid] = true;
              ++numLocalColGIDs;
            }
          }
          else {
            std::pair<typename std::set<GlobalOrdinal>::iterator, bool> ip;
            ip = RemoteGIDSet.insert(gid);
            if (ip.second == true) { // gid did not exist in the set and was actually inserted
              ++numRemoteColGIDs;
            }
          }
        }
      }

      // Possible short-circuit for serial scenario
      // If the all domain GIDs are present as column indices, then set ColMap=DomainMap
      // By construction, LocalGIDs \subset DomainGIDs
      // If we have
      //   * Number of remote GIDs is 0, so that ColGIDs == LocalGIDs,
      // and
      //   * Number of local GIDs is number of domain GIDs
      // then 
      //   * LocalGIDs \subset DomainGIDs && size(LocalGIDs) == size(DomainGIDs) => DomainGIDs == LocalGIDs == ColGIDs
      // on this node. 
      // We will concern ourselves only with the special case of a serial DomainMap, obliviating the need for a communication.
      // If 
      //   * DomainMap has a serial comm
      // then we can set Column map as Domain map and return. Benefit: this graph won't need an Import object later.
      // 
      // Note, for a serial domain map, there can be no RemoteGIDs, because there are no remote nodes. 
      // Likely explanations for this are:
      //  * user submitted erroneous column indices
      //  * user submitted erroneous domain map
      if (Teuchos::size(*domainMap.getComm()) == 1) {
        TEST_FOR_EXCEPTION(numRemoteColGIDs != 0, std::runtime_error,
            Teuchos::typeName(*this) << "::makeColMap(): Some column IDs are not in the domain map." << std::endl 
            << "Either these column IDs are invalid or the domain map is invalid." << std::endl
            << "Remember, for a rectangular matrix, the domain map must be passed to fillComplete().");
        if (numLocalColGIDs == domainMap.getNumLocalElements()) {
          colMap_ = domainMap;
          hasColMap_ = true;
          return;
        }
      }

      // Now, populate myColumns() with a list of all column GIDs. 
      // Put local GIDs at the front: they correspond to "same" and "permuted" entries between the column map and the domain map
      // Put remote GIDs at the back
      myColumns.resize(numLocalColGIDs + numRemoteColGIDs);
      // get pointers into myColumns for each part
      ArrayView<GlobalOrdinal> LocalColGIDs, RemoteColGIDs;
      LocalColGIDs = myColumns(0,numLocalColGIDs);
      RemoteColGIDs = myColumns(numLocalColGIDs,numRemoteColGIDs);

      // Copy the Remote GIDs into myColumns
      std::copy(RemoteGIDSet.begin(), RemoteGIDSet.end(), RemoteColGIDs.begin());
      // We will make a list of correspodning node IDs
      Array<int> RemoteImageIDs(numRemoteColGIDs);
      // Lookup the Remote Nodes IDs in the Domain map
      bool missingID_lcl = domainMap.getRemoteIndexList(RemoteColGIDs, RemoteImageIDs());
      // This error check is crucial: this tells us that one of the remote indices was not present in the domain map. 
      // This means that the Import object cannot be constructed, because of incongruity between the column map and domain map.
      //   * The user has made a mistake in the column indices
      //   * The user has made a mistake w.r.t. the domain map 
      // Same error message as above for serial case.
      char missingID_gbl;
      Teuchos::reduceAll<int,char>(*getComm(),Teuchos::REDUCE_MAX,(missingID_lcl ? 1 : 0),&missingID_gbl);
      TEST_FOR_EXCEPTION(missingID_gbl == 1, std::runtime_error,
          Teuchos::typeName(*this) << "::makeColMap(): Some column IDs are not in the domain map." << std::endl 
          << "Either these column IDs are invalid or the domain map is invalid." << std::endl
          << "Remember, for a rectangular matrix, the domain map must be passed to fillComplete().");
      // Sort External column indices so that all columns coming from a given remote processor are contiguous
      // This obliviates the need for the Distributor associated with the Import from having to reorder data.
      sort2(RemoteImageIDs.begin(), RemoteImageIDs.end(), RemoteColGIDs.begin());

      // Copy the Local GIDs into myColumns. Two cases:
      // (1) If the number of Local column GIDs is the same as the number of Local domain GIDs, we
      //     can simply read the domain GIDs into the front part of ColIndices (see logic above from the serial short circuit case)
      // (2) We step through the GIDs of the DomainMap, checking to see if each domain GID is a column GID.
      //     we want to do this to maintain a consistent ordering of GIDs between the columns and the domain.
      ArrayView<const GlobalOrdinal> mge = domainMap.getMyGlobalEntries();
      if (numLocalColGIDs == domainMap.getNumLocalElements()) {
        std::copy(mge.begin(), mge.end(), LocalColGIDs.begin());
      }
      else {
        size_t numlocalagain = 0;
        for (size_t i=0; i<domainMap.getNumLocalElements(); ++i) {
          if (GIDisLocal[i]) {
            LocalColGIDs[numlocalagain++] = mge[i];
          }
        }
        TEST_FOR_EXCEPTION(numlocalagain != numLocalColGIDs, std::logic_error,
            Teuchos::typeName(*this) << "::makeColMap(): Internal logic error. Please contact Tpetra team.");
      }
    }
    colMap_ = Teuchos::rcp(new Map<LocalOrdinal,GlobalOrdinal,Node>(GOT::invalid(), myColumns, domainMap.getIndexBase(), domainMap.getComm()) );
    hasColMap_ = true;
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::removeRedundantIndices() 
  {
    TEST_FOR_EXCEPT(isGloballyIndexed()==true);    // this should be called only after makeIndicesLocal()
    TEST_FOR_EXCEPT(indicesAreSorted()==false);   // this should be called only after sortIndices()
    TEST_FOR_EXCEPT(isStorageOptimized()==true);  // assumptions about available data structures
    if (noRedundancies()) return;
    const size_t nlrs = getNodeNumRows();
    Teuchos::ArrayView<const GlobalOrdinal> myGlobalEntries = getRowMap().getMyGlobalEntries();
    // reset all local quantities
    upperTriangular_ = true;
    lowerTriangular_ = true;
    localMaxNumEntries_ = 0;
    numLocalEntries_    = 0;
    numLocalDiags_      = 0;
    // indices are already sorted in each row
    if (indicesAreAllocated()) {
      for (size_t r=0; r < nlrs; ++r)
      {
        GlobalOrdinal rgid = myGlobalEntries[r];
        LocalOrdinal rlcid = getColMap().getLocalIndex(rgid);   // determine the local column index for this row, used for delimiting the diagonal
        typename Teuchos::ArrayRCP<LocalOrdinal>::iterator beg, end;
        beg = getLptr(r);
        end = beg + rowNumEntries_[r];
        typename Teuchos::ArrayRCP<LocalOrdinal>::iterator newend, cur;
        newend = beg;
        if (beg != end) {
          if (rlcid == (*newend)) ++numLocalDiags_;
          for (cur = beg + 1; cur != end; ++cur) {
            // it is unique; add it and test it for the diagonal
            // otherwise, it is a dup and should be ignored
            if (*cur != *newend) {
              ++newend;
              (*newend) = (*cur);
              // is this the diagonal?
              if (rlcid == (*newend)) ++numLocalDiags_;
            }
          }
          // because of sorting, smallest column index is (*beg); it indicates upper triangularity
          if ((*beg) < r) upperTriangular_ = false;
          // because of sorting, largest column index is (*newend); it indicates lower triangularity
          if (r < (*newend)) lowerTriangular_ = false;
          // increment newend so that our range is [beg,newend) instead of [beg,newend]
          ++newend;
        }
        // compute num entries for this row, accumulate into numLocalEntries_, update localMaxNumEntries_
        rowNumEntries_[r] = newend - beg;
        localMaxNumEntries_ = std::max( localMaxNumEntries_, rowNumEntries_[r] );
        numLocalEntries_ += rowNumEntries_[r];
      }
    }
    noRedundancies_ = true;
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::makeImportExport()
  {
    TEST_FOR_EXCEPT(hasColMap()==false); // must have column map
    // create import, export
    if (!domainMap_.isSameAs(colMap_)) {
      importer_ = Teuchos::rcp( new Import<LocalOrdinal,GlobalOrdinal,Node>(domainMap_,colMap_) );
    }
    else {
      importer_ = Teuchos::null;
    }
    if (!rangeMap_.isSameAs(rowMap_)) {
      exporter_ = Teuchos::rcp( new Export<LocalOrdinal,GlobalOrdinal,Node>(rowMap_,rangeMap_) );
    }
    else {
      exporter_ = Teuchos::null;
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::computeIndexState()
  {
    char myIndicesAreLocal = 0;
    char myIndicesAreGlobal = 0;
    if(indicesAreLocal_) 
      myIndicesAreLocal = 1;
    if(indicesAreGlobal_) 
      myIndicesAreGlobal = 1;
    char allIndicesAreLocal;
    char allIndicesAreGlobal;
    Teuchos::reduceAll(*getComm(),Teuchos::REDUCE_MAX,myIndicesAreLocal ,&allIndicesAreLocal );
    Teuchos::reduceAll(*getComm(),Teuchos::REDUCE_MAX,myIndicesAreGlobal,&allIndicesAreGlobal );
    indicesAreLocal_  = (allIndicesAreLocal ==1);  // If indices are local on one PE, should be local on all
    indicesAreGlobal_ = (allIndicesAreGlobal==1);  // If indices are global on one PE should be local on all
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::indicesAreAllocated() const
  { return indicesAreAllocated_; }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::RNNZ(size_t row) const
  {
    // if storage is optimized, then numalloc == numnz, and indices are local
    if (isStorageOptimized()) {
      return (colLIndsPtrs_[row+1] - colLIndsPtrs_[row]);
    }
    // otherwise, we are still numnz
    else {
      return rowNumEntries_[row];
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::RNumAlloc(size_t row) const
  {
    if (indicesAreAllocated() == false) {
      return rowNumToAlloc_[row];
    }
    // if static graph or optimize storage, pointers tell us the allocation size
    else if (isStorageOptimized() || (getProfileType() == StaticProfile)) {
      if (isLocallyIndexed()) {
        return colLIndsPtrs_[row+1] - colLIndsPtrs_[row];
      }
      else {
        return colGIndsPtrs_[row+1] - colGIndsPtrs_[row];
      }
    }
    // otherwise, the ArrayRCP knows
    else {
      if (isLocallyIndexed()) {
        return colLInds_[row].size();
      }
      else {
        return colGInds_[row].size();
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::staticAssertions()
  {
    using Teuchos::OrdinalTraits;
    // Assumption: sizeof(GlobalOrdinal) >= sizeof(LocalOrdinal)
    //    This is so that we can store LocalOrdinals in the memory formerly occupied by GlobalOrdinals
    // Assumption: max(GlobalOrdinal) >= max(LocalOrdinal)  and  max(size_t) >= max(LocalOrdinal)
    //    This is so that we can represent any LocalOrdinal as a size_t, and any LocalOrdinal as a GlobalOrdinal
    Teuchos::CompileTimeAssert<sizeof(GlobalOrdinal) < sizeof(LocalOrdinal)> cta_size;
    (void)cta_size;
    // can't call max() with CompileTimeAssert, because it isn't a constant expression; will need to make this a runtime check
    TEST_FOR_EXCEPTION( OrdinalTraits<GlobalOrdinal>::max() < OrdinalTraits<LocalOrdinal>::max(), std::runtime_error, 
        Teuchos::typeName(*this) << ": Object cannot be allocated with stated template arguments: size assumptions are not valid.");
    TEST_FOR_EXCEPTION( OrdinalTraits<size_t>::max() < OrdinalTraits<LocalOrdinal>::max(), std::runtime_error,
        Teuchos::typeName(*this) << ": Object cannot be allocated with stated template arguments: size assumptions are not valid.");
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Teuchos::ArrayRCP<const LocalOrdinal>::iterator 
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getLptr(size_t row) const
  {
    if ((getProfileType() == StaticProfile) || isStorageOptimized()) {
      return colLIndsPtrs_[row];
    }
    else {
      return colLInds_[row].begin();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Teuchos::ArrayRCP<LocalOrdinal>::iterator 
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getLptr(size_t row)
  {
    if ((getProfileType() == StaticProfile) || isStorageOptimized()) {
      return colLIndsPtrs_[row];
    }
    else {
      return colLInds_[row].begin();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Teuchos::ArrayRCP<const GlobalOrdinal>::iterator 
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getGptr(size_t row) const
  {
    if ((getProfileType() == StaticProfile) || isStorageOptimized()) {
      return colGIndsPtrs_[row];
    }
    else {
      return colGInds_[row].begin();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Teuchos::ArrayRCP<GlobalOrdinal>::iterator 
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getGptr(size_t row) 
  {
    if ((getProfileType() == StaticProfile) || isStorageOptimized()) {
      return colGIndsPtrs_[row];
    }
    else {
      return colGInds_[row].begin();
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t 
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::findMyIndex(size_t row, LocalOrdinal ind) const
  {
    typedef typename Teuchos::ArrayRCP<const LocalOrdinal>::iterator IT;
    bool found = true;
    IT rptr, locptr;
    rptr = getLptr(row);
    size_t nE = RNNZ(row);
    if (indicesAreSorted()) {
      // binary search
      std::pair<IT,IT> p = std::equal_range(rptr,rptr+nE,ind);
      if (p.first == p.second) found = false;
      else locptr = p.first;
    }
    else {
      // direct search
      locptr = std::find(rptr,rptr+nE,ind);
      if (locptr == rptr+nE) found = false;
    }
    if (!found) {
      return Teuchos::OrdinalTraits<size_t>::invalid();
    }
    return (locptr - rptr);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t 
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::findGlobalIndex(size_t row, GlobalOrdinal ind) const
  {
    typedef typename Teuchos::ArrayRCP<const GlobalOrdinal>::iterator IT;
    bool found = true;
    IT rptr, locptr;
    rptr = getGptr(row);
    size_t nE = RNNZ(row);
    if (indicesAreSorted()) {
      // binary search
      std::pair<IT,IT> p = std::equal_range(rptr,rptr+nE,ind);
      if (p.first == p.second) found = false;
      else locptr = p.first;
    }
    else {
      // direct search
      locptr = std::find(rptr,rptr+nE,ind);
      if (locptr == rptr+nE) found = false;
    }
    if (!found) {
      return Teuchos::OrdinalTraits<size_t>::invalid();
    }
    return (locptr - rptr);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::print(std::ostream &os) const
  {
    // support this operation whether fillComplete() or not
    using std::endl;
    int myImageID = Teuchos::rank(*getComm());
    if (myImageID == 0)
    {
      os << "Tpetra::CrsGraph, label = " << this->getObjectLabel() << endl;
      os << "Number of global rows   = " << getGlobalNumRows() << endl;
      if (isFillComplete())
      {
        os << "Number of global columns    = " << getGlobalNumCols() << endl;
        os << "Status = fill complete" << endl;
        os << "Number of global nonzeros   = " << getGlobalNumEntries() << endl;
        os << "Global max nonzeros per row = " << getGlobalMaxNumRowEntries() << endl;
      }
      else
      {
        os << "Status = fill not complete" << endl;
      }
    }
    if (isFillComplete())
    {
      for (int pid=0; pid < Teuchos::size(*getComm()); ++pid)
      {
        if (pid == myImageID)
        {
          Teuchos::Array<GlobalOrdinal> indices(getLocalMaxNumRowEntries());
          Teuchos::ArrayView<const GlobalOrdinal> mgis = getRowMap().getMyGlobalEntries();
          os << "% Number of rows on image " << myImageID << " = " << getNodeNumRows() << endl;
          for (size_t r=0; r < getNodeNumRows(); ++r)
          {                                                        
            GlobalOrdinal globalRow = mgis[r];
            size_t rowSize;
            getGlobalRowCopy(globalRow, indices(), rowSize);
            if (rowSize > 0) {
              os << "Row " << globalRow << ":";
              for (size_t j=0; j < rowSize; ++j) {
                os << " " << indices[j];
              }
              os << std::endl;
            }
          }
        }
        Teuchos::barrier(*getComm());
        Teuchos::barrier(*getComm());
        Teuchos::barrier(*getComm());
      }
    }
  }


} // namespace Tpetra

#endif
