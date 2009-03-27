//@HEADER
// ************************************************************************
// 
//               Tpetra: Templated Linear Algebra Services Package 
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
// ************************************************************************
//@HEADER

#ifndef TPETRA_CRSGRAPH_HPP
#define TPETRA_CRSGRAPH_HPP

#include <Teuchos_Describable.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_getRawPtr.hpp>

#include "Tpetra_RowGraph.hpp"
#include "Tpetra_CrsGraphData.hpp"
#include "Tpetra_Util.hpp"

namespace Tpetra 
{

  // forward declaration
  template <class S, class LO, class GO>
  class CrsMatrix;

  //! CrsGraph
  template<class LocalOrdinal, class GlobalOrdinal=LocalOrdinal>
  class CrsGraph : public RowGraph<LocalOrdinal,GlobalOrdinal>
  {
    template <class S, class LO, class GO>
    friend class CrsMatrix;
    public: 

      //! @name Constructor/Destructor Methods
      //@{ 

      //! Constructor with fixed number of indices per row.
      CrsGraph(const Map<LocalOrdinal,GlobalOrdinal> &rowMap, Teuchos_Ordinal maxNumEntriesPerRow, bool staticProfile = false);

      //! Constructor with variable number of indices per row.
      CrsGraph(const Map<LocalOrdinal,GlobalOrdinal> &rowMap, const Teuchos::ArrayView<Teuchos_Ordinal> &NumEntriesPerRowToAlloc, bool staticProfile = false);

      //! Constructor with fixed number of indices per row and specified column map.
      CrsGraph(const Map<LocalOrdinal,GlobalOrdinal> &rowMap, const Map<LocalOrdinal,GlobalOrdinal> &colMap, Teuchos_Ordinal maxNumEntriesPerRow, bool staticProfile = false);

      //! Constructor with variable number of indices per row and specified column map.
      CrsGraph(const Map<LocalOrdinal,GlobalOrdinal> &rowMap, const Map<LocalOrdinal,GlobalOrdinal> &colMap, const Teuchos::ArrayView<Teuchos_Ordinal> &NumEntriesPerRowToAlloc, bool staticProfile = false);

      // !Destructor.
      virtual ~CrsGraph();

      //@}

      //! @name Insertion/Removal Methods
      //@{ 

      //! allows insertion of global indices, either belonging to this node or intended for another.
      void insertGlobalIndices(GlobalOrdinal row, const Teuchos::ArrayView<const GlobalOrdinal> &indices);

      //! allows insertion of local indices intended for this node.
      void insertMyIndices(LocalOrdinal row, const Teuchos::ArrayView<const LocalOrdinal> &indices);

      //@}

      //! @name Transformational Methods
      //@{ 

      /*! \brief Signal that data entry is complete, specifying domain and range maps. 
          Off-node entries are distributed, repeated entries are summed, and global indices are transformed to local indices.
          If \c OptimizeStorage is true, then optimizeStorage() is called as well.
       */ 
      void fillComplete(const Map<LocalOrdinal,GlobalOrdinal> &domainMap, const Map<LocalOrdinal,GlobalOrdinal> &rangeMap, bool OptimizeStorage=false);

      /*! \brief Signal that data entry is complete. 
          Off-node entries are distributed, repeated entries are summed, and global indices are transformed to local indices.
          If \c OptimizeStorage is true, then optimizeStorage() is called as well.
          \note This method calls fillComplete( getRowMap(), getRowMap(), OptimizeStorage ).
       */
      void fillComplete(bool OptimizeStorage=false);

      //! \brief Communicate non-local contributions to other nodes.
      void globalAssemble();

      //! \brief Re-allocate the data into contiguous storage.
      void optimizeStorage();

      //! Returns \c true if fillComplete() has been called.
      bool isFillComplete() const;

      //! Returns \c true if optimizeStorage() has been called.
      bool isStorageOptimized() const;

      //! Returns \c true if the graph data was allocated in static data structures.
      bool isStaticProfile() const;

      //@}

      //! @name Graph Query Methods
      //@{ 

      //! Returns the communicator.
      Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

      //! Returns the Map that describes the row distribution in this graph.
      const Map<LocalOrdinal,GlobalOrdinal> & getRowMap() const;

      //! \brief Returns the Map that describes the column distribution in this graph.
      /*! Cannot be called before fillComplete(), unless the matrix was constructed with a column map. */
      const Map<LocalOrdinal,GlobalOrdinal> & getColMap() const;

      //! Returns the Map associated with the domain of this graph.
      const Map<LocalOrdinal,GlobalOrdinal> & getDomainMap() const;

      //! Returns the Map associated with the domain of this graph.
      const Map<LocalOrdinal,GlobalOrdinal> & getRangeMap() const;

      //! Returns the importer associated with this graph.
      Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal> > getImporter() const;

      //! Returns the exporter associated with this graph.
      Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal> > getExporter() const;

      //! Returns the number of global rows in the graph.
      GlobalOrdinal numGlobalRows() const;

      //! \brief Returns the number of global columns in the graph.
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      GlobalOrdinal numGlobalCols() const;

      //! Returns the number of rows owned on the calling node.
      Teuchos_Ordinal numLocalRows() const;

      //! Returns the number of columns connected to the locally owned rows of this graph.
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      Teuchos_Ordinal numLocalCols() const;

      //! Returns the index base for global indices for this graph. 
      Teuchos_Ordinal getIndexBase() const;

      //! \brief Returns the number of entries in the global matrix. 
      /*! Returns the number of global entries in the associated graph. */
      GlobalOrdinal numGlobalEntries() const;

      //! \brief Returns the number of entries in the calling image's portion of the matrix. 
      /*! Before fillComplete() is called, this could include duplicated entries. */
      Teuchos_Ordinal numMyEntries() const;

      //! \brief Returns the current number of entries on this node in the specified global row .
      /*! Throws exception std::runtime_error if the specified global row does not belong to this node. */
      Teuchos_Ordinal numEntriesForGlobalRow(GlobalOrdinal globalRow) const;

      //! Returns the current number of entries on this node in the specified local row.
      /*! Throws exception std::runtime_error if the specified local row is not valid for this node. */
      Teuchos_Ordinal numEntriesForMyRow(LocalOrdinal localRow) const;

      //! \brief Returns the current number of allocated entries for this node in the specified global row .
      /*! Throws exception std::runtime_error if the specified global row does not belong to this node. */
      Teuchos_Ordinal numAllocatedEntriesForGlobalRow(GlobalOrdinal globalRow) const;

      //! Returns the current number of allocated entries on this node in the specified local row.
      /*! Throws exception std::runtime_error if the specified local row is not valid for this node. */
      Teuchos_Ordinal numAllocatedEntriesForMyRow(LocalOrdinal localRow) const;

      //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      GlobalOrdinal numGlobalDiagonals() const;

      //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      Teuchos_Ordinal numMyDiagonals() const;

      //! \brief Returns the maximum number of entries across all rows/columns on all nodes. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      GlobalOrdinal globalMaxNumRowEntries() const;

      //! \brief Returns the maximum number of entries across all rows/columns on this node. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      Teuchos_Ordinal myMaxNumRowEntries() const;

      //! \brief Indicates whether the matrix has a well-defined column map. 
      /*! The column map does not exist until after fillComplete(), unless the matrix was constructed with one. */
      bool hasColMap() const; 

      //! \brief Indicates whether the graph is lower triangular.
      bool lowerTriangular() const;

      //! \brief Indicates whether the graph is upper triangular.
      bool upperTriangular() const;

      //! \brief If graph indices are in the local range, this function returns true. Otherwise, this function returns false. */
      bool indicesAreLocal() const;

      //! \brief If graph indices are in the global range, this function returns true. Otherwise, this function returns false. */
      bool indicesAreGlobal() const;

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
         returned as Teuchos::OrdinalTraits<Teuchos_Ordinal>::invalid().
       */
      void extractGlobalRowCopy(GlobalOrdinal GlobalRow, Teuchos::ArrayView<GlobalOrdinal> indices, Teuchos_Ordinal &NumIndices) const;

      //! Extract a list of elements in a specified local row of the graph. Put into storage allocated by calling routine.
      /*!
        \param LocalRow - (In) Local row number for which indices are desired.
        \param Indices - (Out) Local column indices corresponding to values.
        \param NumIndices - (Out) Number of indices.

         Note: A std::runtime_error exception is thrown indices is not large enough to hold the column indices associated
         with row \c LocalRow. If \c LocalRow is not valid for this node, then \c indices is unchanged and \c NumIndices is 
         returned as Teuchos::OrdinalTraits<Teuchos_Ordinal>::invalid().

        \pre indicesAreLocal()==true
       */
      void extractMyRowCopy(Teuchos_Ordinal LocalRow, Teuchos::ArrayView<LocalOrdinal> indices, Teuchos_Ordinal& NumIndices) const;

      //! Get a non-persisting view of the elements in a specified global row of the graph.
      /*!
        \param GlobalRow - (In) Global row number to get indices.
        \param Indices - (Out) Global column indices corresponding to values.

         Note: If \c GlobalRow does not belong to this node, then \c indices is unchanged and \c NumIndices is 
         returned as Teuchos::OrdinalTraits<Teuchos_Ordinal>::invalid().

        \pre indicesAreLocal()==false
       */
      void extractGlobalRowView(GlobalOrdinal GlobalRow, Teuchos::ArrayView<GlobalOrdinal> &indices);

      //! Get a view of the elements in a specified local row of the graph.
      /*!
        \param LocalRow - (In) Local row number to get indices.
        \param Indices - (Out) Local column indices corresponding to values.

         Note: If \c LocalRow is not valid for this node, then \c indices is unchanged and \c NumIndices is 
         returned as Teuchos::OrdinalTraits<Teuchos_Ordinal>::invalid().

        \pre indicesAreLocal()==true
       */
      void extractMyRowView(LocalOrdinal LocalRow, Teuchos::ArrayView<LocalOrdinal> &indices);

      //! Get a non-persisting view of the elements in a specified global row of the graph.
      /*!
        \param GlobalRow - (In) Global row number to get indices.
        \param Indices - (Out) Global column indices corresponding to values.

         Note: If \c GlobalRow does not belong to this node, then \c indices is unchanged and \c NumIndices is 
         returned as Teuchos::OrdinalTraits<Teuchos_Ordinal>::invalid().

        \pre indicesAreLocal()==false
       */
      void extractGlobalRowConstView(GlobalOrdinal GlobalRow, Teuchos::ArrayView<const GlobalOrdinal> &indices) const;

      //! Get a view of the elements in a specified local row of the graph.
      /*!
        \param LocalRow - (In) Local row number to get indices.
        \param Indices - (Out) Local column indices corresponding to values.

         Note: If \c LocalRow is not valid for this node, then \c indices is unchanged and \c NumIndices is 
         returned as Teuchos::OrdinalTraits<Teuchos_Ordinal>::invalid().

        \pre indicesAreLocal()==true
       */
      void extractMyRowConstView(LocalOrdinal LocalRow, Teuchos::ArrayView<const LocalOrdinal> &indices) const;

      //@}

      //! @name I/O Methods
      //@{ 
      
      //! Prints the graph on the specified stream. This is very verbose.
      void print(std::ostream& os) const;

      // @}

    protected:
      void initNumAlloc(const Teuchos::ArrayView<Teuchos_Ordinal> &NumEntriesPerRowToAlloc);
      void allocateIndices(bool local);
      void makeIndicesLocal(const Map<LocalOrdinal,GlobalOrdinal> &domainMap, const Map<LocalOrdinal,GlobalOrdinal> &rangeMap);
      void makeColMap(const Map<LocalOrdinal,GlobalOrdinal> &domainMap, const Map<LocalOrdinal,GlobalOrdinal> &rangeMap);
      void computeIndexState();
      void sortIndices();
      void removeRedundantIndices();
      void makeImportExport();
      bool indicesAreSorted() const;
      void indicesAreSorted(bool sorted);
      bool indicesAreAllocated() const;
      void staticAssertions();
      inline Teuchos_Ordinal RNNZ(Teuchos_Ordinal row) const;
      inline Teuchos_Ordinal RNumAlloc(Teuchos_Ordinal row) const;
      inline typename Teuchos::ArrayRCP<const LocalOrdinal>::iterator getLptr(Teuchos_Ordinal row) const;
      inline typename Teuchos::ArrayRCP<const GlobalOrdinal>::iterator getGptr(Teuchos_Ordinal row) const;
      inline typename Teuchos::ArrayRCP<LocalOrdinal>::iterator getLptr(Teuchos_Ordinal row);
      inline typename Teuchos::ArrayRCP<GlobalOrdinal>::iterator getGptr(Teuchos_Ordinal row);

      Teuchos::RCP<CrsGraphData<LocalOrdinal,GlobalOrdinal> > graphData_;
      std::map<GlobalOrdinal, std::deque<GlobalOrdinal> > nonlocals_;
  };

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  template <class LocalOrdinal, class GlobalOrdinal>
  CrsGraph<LocalOrdinal,GlobalOrdinal>::CrsGraph(const Map<LocalOrdinal,GlobalOrdinal> &rowMap, Teuchos_Ordinal maxNumEntriesPerRow, bool staticProfile)
  : graphData_(Teuchos::rcp(new CrsGraphData<LocalOrdinal,GlobalOrdinal>(rowMap)))
  {
    staticAssertions();
    graphData_->staticProfile_ = staticProfile;
    TEST_FOR_EXCEPTION(maxNumEntriesPerRow < 1 && maxNumEntriesPerRow != 0, std::runtime_error,
        Teuchos::typeName(*this) << "::CrsGraph(rowMap,maxNumEntriesPerRow): maxNumEntriesPerRow must be non-negative.");
    std::fill(graphData_->rowNumToAlloc_.begin(), graphData_->rowNumToAlloc_.end(), maxNumEntriesPerRow);
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  CrsGraph<LocalOrdinal,GlobalOrdinal>::CrsGraph(const Map<LocalOrdinal,GlobalOrdinal> &rowMap, const Map<LocalOrdinal,GlobalOrdinal> &colMap, Teuchos_Ordinal maxNumEntriesPerRow, bool staticProfile)
  : graphData_(Teuchos::rcp(new CrsGraphData<LocalOrdinal,GlobalOrdinal>(rowMap,colMap)))
  {
    staticAssertions();
    graphData_->staticProfile_ = staticProfile;
    TEST_FOR_EXCEPTION(maxNumEntriesPerRow < 1 && maxNumEntriesPerRow != 0, std::runtime_error,
        Teuchos::typeName(*this) << "::CrsGraph(rowMap,colMap,maxNumEntriesPerRow): maxNumEntriesPerRow must be non-negative.");
    std::fill(graphData_->rowNumToAlloc_.begin(), graphData_->rowNumToAlloc_.end(), maxNumEntriesPerRow);
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  CrsGraph<LocalOrdinal,GlobalOrdinal>::CrsGraph(const Map<LocalOrdinal,GlobalOrdinal> &rowMap, const Teuchos::ArrayView<Teuchos_Ordinal> &NumEntriesPerRowToAlloc, bool staticProfile)
  : graphData_(Teuchos::rcp(new CrsGraphData<LocalOrdinal,GlobalOrdinal>(rowMap)))
  {
    staticAssertions();
    graphData_->staticProfile_ = staticProfile;
    TEST_FOR_EXCEPTION(NumEntriesPerRowToAlloc.size() != numLocalRows(), std::runtime_error,
        Teuchos::typeName(*this) << "::CrsGraph(rowMap,NumEntriesPerRowToAlloc): NumEntriesPerRowToAlloc must have as many entries as specified by rowMap for this node.");
    initNumAlloc(NumEntriesPerRowToAlloc);
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  CrsGraph<LocalOrdinal,GlobalOrdinal>::CrsGraph(const Map<LocalOrdinal,GlobalOrdinal> &rowMap, const Map<LocalOrdinal,GlobalOrdinal> &colMap, const Teuchos::ArrayView<Teuchos_Ordinal> &NumEntriesPerRowToAlloc, bool staticProfile)
  : graphData_(Teuchos::rcp(new CrsGraphData<LocalOrdinal,GlobalOrdinal>(rowMap,colMap)))
  {
    staticAssertions();
    graphData_->staticProfile_ = staticProfile;
    TEST_FOR_EXCEPTION(NumEntriesPerRowToAlloc.size() != numLocalRows(), std::runtime_error,
        Teuchos::typeName(*this) << "::CrsGraph(rowMap,colMap,NumEntriesPerRowToAlloc): NumEntriesPerRowToAlloc must have as many entries as specified by rowMap for this node.");
    initNumAlloc(NumEntriesPerRowToAlloc);
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  void CrsGraph<LocalOrdinal,GlobalOrdinal>::allocateIndices(bool local)
  {
    // this is a protected function, only callable by us. if it was called incorrectly, it is our fault.
    TEST_FOR_EXCEPTION(indicesAreLocal() && local==false, std::logic_error,
        Teuchos::typeName(*this) << "::allocateIndices(): Internal logic error. Please contact Tpetra team.");
    TEST_FOR_EXCEPTION(indicesAreGlobal() && local==true, std::logic_error,
        Teuchos::typeName(*this) << "::allocateIndices(): Internal logic error. Please contact Tpetra team.");
    if (indicesAreAllocated()) {
      return;
    }
    graphData_->indicesAreLocal_  =  local;
    graphData_->indicesAreGlobal_ = !local;
    if (isStaticProfile()) {
      typename Teuchos::Array<LocalOrdinal>::iterator    
        ntaptr = graphData_->rowNumToAlloc_.begin(),
        ntaend = graphData_->rowNumToAlloc_.end();
      Teuchos_Ordinal ntasum = std::accumulate(ntaptr, ntaend, 0);
      Teuchos_Ordinal sofar = 0, r = 0;
      if (local) {
        graphData_->colLIndsPtrs_.resize(numLocalRows()+1);
        if (ntasum) {
          graphData_->contigColLInds_ = Teuchos::arcp<LocalOrdinal>(ntasum);
          while (ntaptr != ntaend) {
            graphData_->colLIndsPtrs_[r++] = graphData_->contigColLInds_.persistingView(sofar,*ntaptr).begin();
            sofar += (*ntaptr++);
          }
          graphData_->colLIndsPtrs_[r] = graphData_->contigColLInds_.end();
        }
      }
      else {
        graphData_->colGIndsPtrs_.resize(numLocalRows()+1);
        if (ntasum) {
          graphData_->contigColGInds_ = Teuchos::arcp<GlobalOrdinal>(ntasum);
          while (ntaptr != ntaend) {
            graphData_->colGIndsPtrs_[r++] = graphData_->contigColGInds_.persistingView(sofar,*ntaptr).begin();
            sofar += (*ntaptr++);
          }
          graphData_->colGIndsPtrs_[r] = graphData_->contigColGInds_.end();
        }
      }
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION(sofar != ntasum, std::logic_error, 
          Teuchos::typeName(*this) << "::allocateIndices(): Internal logic error. Please contact Tpetra team.");
#endif
    }
    else {
      if (local)  // allocate local indices in graphData_->colLInds_
      {
        typename Teuchos::Array<LocalOrdinal>::iterator 
          ntaptr = graphData_->rowNumToAlloc_.begin(),
          ntaend = graphData_->rowNumToAlloc_.end();
        graphData_->colLInds_.resize(numLocalRows());
        typename Teuchos::Array<Teuchos::ArrayRCP<LocalOrdinal> >::iterator rinds  = graphData_->colLInds_.begin();
        for (; ntaptr != ntaend; ++ntaptr, ++rinds) {
          if (*ntaptr > 0) {
            (*rinds) = Teuchos::arcp<LocalOrdinal>(*ntaptr);
          }
        }
      }
      else        // allocate global indices in graphData_->colGInds_
      {
        typename Teuchos::Array<LocalOrdinal>::iterator 
          ntaptr = graphData_->rowNumToAlloc_.begin(),
          ntaend = graphData_->rowNumToAlloc_.end();
        graphData_->colGInds_.resize(numLocalRows());
        typename Teuchos::Array<Teuchos::ArrayRCP<GlobalOrdinal> >::iterator rinds  = graphData_->colGInds_.begin();
        for (; ntaptr != ntaend; ++ntaptr, ++rinds) {
          if (*ntaptr > 0) {
            (*rinds) = Teuchos::arcp<GlobalOrdinal>(*ntaptr);
          }
        }
      }
    }
    graphData_->rowNumToAlloc_.clear();
    graphData_->indicesAreAllocated_ = true;    
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  void CrsGraph<LocalOrdinal,GlobalOrdinal>::initNumAlloc(const Teuchos::ArrayView<Teuchos_Ordinal> &NumEntriesPerRowToAlloc)
  {
    typename Teuchos::ArrayView<Teuchos_Ordinal>::iterator nnz    = NumEntriesPerRowToAlloc.begin(),
                                                           nnzend = NumEntriesPerRowToAlloc.end();
    typename Teuchos::Array<LocalOrdinal>::iterator ralloc = graphData_->rowNumToAlloc_.begin();
    for (; nnz != nnzend; ++nnz, ++ralloc) {
      if (*nnz > 0) {
        (*ralloc) = *nnz;
      }
      else if (*nnz != 0) {
        TEST_FOR_EXCEPTION(true, std::runtime_error,
            Teuchos::typeName(*this) << "::initNumAlloc(NumEntriesPerRowToAlloc): all entries of NumEntriesPerRowToAlloc must be non-negative.");
      }
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  CrsGraph<LocalOrdinal,GlobalOrdinal>::~CrsGraph()
  {}

  template <class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal CrsGraph<LocalOrdinal,GlobalOrdinal>::numGlobalRows() const
  { return graphData_->rowMap_.getNumGlobalEntries(); }

  template <class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal CrsGraph<LocalOrdinal,GlobalOrdinal>::numGlobalCols() const
  {
    TEST_FOR_EXCEPTION(isFillComplete() != true, std::runtime_error,
      Teuchos::typeName(*this) << ": cannot call numGlobalCols() until fillComplete() is called.");
    return getDomainMap().getNumGlobalEntries();
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsGraph<LocalOrdinal,GlobalOrdinal>::numLocalRows() const
  { return graphData_->rowMap_.getNumMyEntries(); }

  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsGraph<LocalOrdinal,GlobalOrdinal>::numLocalCols() const
  {
    TEST_FOR_EXCEPTION(hasColMap() != true, std::runtime_error,
      Teuchos::typeName(*this) << ": cannot call numLocalCols() without a column map.");
    return graphData_->colMap_.getNumMyEntries();
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsGraph<LocalOrdinal,GlobalOrdinal>::numMyDiagonals() const
  { 
    TEST_FOR_EXCEPTION(!isFillComplete(), std::runtime_error,
      Teuchos::typeName(*this) << ": cannot call numMyDiagonals() until fillComplete() has been called.");
    return graphData_->numLocalDiags_; 
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal CrsGraph<LocalOrdinal,GlobalOrdinal>::numGlobalDiagonals() const
  {
    TEST_FOR_EXCEPTION(!isFillComplete(), std::runtime_error,
      Teuchos::typeName(*this) << ": cannot call numGlobalDiagonals() until fillComplete() has been called.");
    return graphData_->numGlobalDiags_;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  void CrsGraph<LocalOrdinal,GlobalOrdinal>::extractMyRowConstView(LocalOrdinal LocalRow, Teuchos::ArrayView<const LocalOrdinal> &indices) const
  {
    TEST_FOR_EXCEPTION(indicesAreGlobal() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::extractMyRowConstView(): local indices do not exist.");
    TEST_FOR_EXCEPTION(getRowMap().isMyLocalIndex(LocalRow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::extractMyRowConstView(LocalRow,...): LocalRow (== " << LocalRow << ") is not valid on this node.");
    Teuchos_Ordinal rnnz = RNNZ(LocalRow);
    if (indicesAreAllocated() == false || rnnz == 0) {
      indices = Teuchos::ArrayView<const LocalOrdinal>(Teuchos::null);
      return;
    }
    // RNNZ > 0, so there must be something allocated
    if (isStorageOptimized() == true || isStaticProfile() == true) {
      indices = Teuchos::arrayView<const LocalOrdinal>(Teuchos::getRawPtr(graphData_->colLIndsPtrs_[LocalRow]),rnnz);
    }
    else {
      indices = graphData_->colLInds_[LocalRow](0,rnnz);
    }
    return;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  void CrsGraph<LocalOrdinal,GlobalOrdinal>::extractGlobalRowConstView(GlobalOrdinal GlobalRow, Teuchos::ArrayView<const GlobalOrdinal> &indices) const
  {
    TEST_FOR_EXCEPTION(indicesAreLocal() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::extractGlobalRowConstView(): global indices do not exist.");
    Teuchos_Ordinal lrow = getRowMap().getLocalIndex(GlobalRow);
    TEST_FOR_EXCEPTION(lrow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::extractGlobalRowConstView(GlobalRow,...): GlobalRow (== " << GlobalRow << ") does not belong to this node.");
    Teuchos_Ordinal rnnz = RNNZ(lrow);
    if (indicesAreAllocated() == false || rnnz == 0) {
      indices = Teuchos::ArrayView<const GlobalOrdinal>(Teuchos::null);
      return;
    }
    // RNNZ > 0, so there must be something allocated
    if (isStorageOptimized() == true || isStaticProfile() == true) {
      indices = Teuchos::arrayView<const GlobalOrdinal>(Teuchos::getRawPtr(graphData_->colGIndsPtrs_[lrow]),rnnz);
    }
    else {
      indices = graphData_->colGInds_[lrow](0,rnnz);
    }
    return;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  void CrsGraph<LocalOrdinal,GlobalOrdinal>::extractMyRowView(LocalOrdinal LocalRow, Teuchos::ArrayView<LocalOrdinal> &indices)
  { 
    TEST_FOR_EXCEPTION(indicesAreGlobal() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::extractMyRowView(): local indices do not exist.");
    TEST_FOR_EXCEPTION(getRowMap().isMyLocalIndex(LocalRow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::extractMyRowView(LocalRow,...): LocalRow (== " << LocalRow << ") is not valid on this node.");
    Teuchos_Ordinal rnnz = RNNZ(LocalRow);
    if (indicesAreAllocated() == false || rnnz == 0) {
      indices = Teuchos::ArrayView<LocalOrdinal>(Teuchos::null);
      return;
    }
    // RNNZ > 0, so there must be something allocated
    if (isStorageOptimized() == true || isStaticProfile() == true) {
      indices = Teuchos::arrayView<LocalOrdinal>(Teuchos::getRawPtr(graphData_->colLIndsPtrs_[LocalRow]),rnnz);
    }
    else {
      indices = graphData_->colLInds_[LocalRow](0,rnnz);
    }
    return;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  void CrsGraph<LocalOrdinal,GlobalOrdinal>::extractGlobalRowView(GlobalOrdinal GlobalRow, Teuchos::ArrayView<GlobalOrdinal> &indices)
  {
    TEST_FOR_EXCEPTION(indicesAreLocal() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::extractGlobalRowView(): global indices do not exist.");
    Teuchos_Ordinal lrow = getRowMap().getLocalIndex(GlobalRow);
    TEST_FOR_EXCEPTION(lrow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::extractGlobalRowView(GlobalRow,...): GlobalRow (== " << GlobalRow << ") does not belong to this node.");
    Teuchos_Ordinal rnnz = RNNZ(lrow);
    if (indicesAreAllocated() == false || rnnz == 0) {
      indices = Teuchos::ArrayView<GlobalOrdinal>(Teuchos::null);
      return;
    }
    // RNNZ > 0, so there must be something allocated
    if (isStorageOptimized() == true || isStaticProfile() == true) {
      indices = Teuchos::arrayView<GlobalOrdinal>(Teuchos::getRawPtr(graphData_->colGIndsPtrs_[lrow]),rnnz);
    }
    else {
      indices = graphData_->colGInds_[lrow](0,rnnz);
    }
    return;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  void CrsGraph<LocalOrdinal,GlobalOrdinal>::extractMyRowCopy(Teuchos_Ordinal LocalRow, Teuchos::ArrayView<LocalOrdinal> indices, Teuchos_Ordinal& NumIndices) const
  {
    // can only do this if we have local indices
    TEST_FOR_EXCEPTION(indicesAreGlobal() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::extractMyRowCopy(): local indices do not exist.");
    TEST_FOR_EXCEPTION(getRowMap().isMyLocalIndex(LocalRow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::extractMyRowCopy(LocalRow,...): LocalRow (== " << LocalRow << ") is not valid on this node.");
    Teuchos_Ordinal rnnz = RNNZ(LocalRow);
    TEST_FOR_EXCEPTION(indices.size() < rnnz, std::runtime_error,
        Teuchos::typeName(*this) << "::extractMyRowCopy(): specified storage (size==" << indices.size() 
        << ") is not large enough to hold all entries for this row (rnnz=" << rnnz << ").");
    // use one of the view routines to get the proper view, then copy it over
    NumIndices = rnnz;
    Teuchos::ArrayView<const LocalOrdinal> lview;
    extractMyRowConstView(LocalRow,lview);
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(lview.size() != rnnz, std::logic_error,
        typeName(*this) << "::extractMyRowCopy(): Internal logic error. Please contact Tpetra team.");
#endif
    std::copy(lview.begin(),lview.end(),indices.begin());
    return;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  void CrsGraph<LocalOrdinal,GlobalOrdinal>::extractGlobalRowCopy(GlobalOrdinal GlobalRow, Teuchos::ArrayView<GlobalOrdinal> indices, Teuchos_Ordinal &NumIndices) const
  {
    using Teuchos::ArrayView;
    // because the data is copied into user storage, we can always do this
    Teuchos_Ordinal lrow = getRowMap().getLocalIndex(GlobalRow);
    TEST_FOR_EXCEPTION(lrow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::extractGlobalRowCopy(GlobalRow,...): GlobalRow (== " << GlobalRow << ") does not belong to this node.");
    Teuchos_Ordinal rnnz = RNNZ(lrow);
    TEST_FOR_EXCEPTION(indices.size() < rnnz, std::runtime_error,
        Teuchos::typeName(*this) << "::extractGlobalRowCopy(): specified storage (size==" << indices.size() 
        << ") is not large enough to hold all entries for this row (rnnz=" << rnnz << ").");
    // use one of the view routines to get the proper view, then copy it over
    NumIndices = rnnz;
    if (indicesAreLocal()) {
      ArrayView<const LocalOrdinal> lview;
      extractMyRowConstView(lrow,lview);
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION(lview.size() != rnnz, std::logic_error,
          typeName(*this) << "::extractGlobalRowCopy(): Internal logic error. Please contact Tpetra team.");
#endif
      // copy and convert
      typename ArrayView<      GlobalOrdinal>::iterator i=indices.begin();
      typename ArrayView<const LocalOrdinal >::iterator l=lview.begin();
      while (l!=lview.end())
      {
        (*i++) = getColMap().getGlobalIndex(*l++);
      }
    }
    else if (indicesAreGlobal()) {
      ArrayView<const GlobalOrdinal> gview;
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION(gview.size() != rnnz, std::logic_error,
          typeName(*this) << "::extractGlobalRowCopy(): Internal logic error. Please contact Tpetra team.");
#endif
      extractGlobalRowConstView(GlobalRow,gview);
      std::copy(gview.begin(), gview.end(), indices.begin());
    }
    return;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal CrsGraph<LocalOrdinal,GlobalOrdinal>::globalMaxNumRowEntries() const
  {
    TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalMaxRowEntries(): global quantities not computed until fillComplete() has been called.");
    return graphData_->globalMaxNumRowEntries_;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  const Map<LocalOrdinal,GlobalOrdinal> & 
  CrsGraph<LocalOrdinal,GlobalOrdinal>::getRowMap() const
  { return graphData_->rowMap_; }

  template <class LocalOrdinal, class GlobalOrdinal>
  const Map<LocalOrdinal,GlobalOrdinal> & 
  CrsGraph<LocalOrdinal,GlobalOrdinal>::getColMap() const
  { 
    TEST_FOR_EXCEPTION(hasColMap() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getColMap(): graph does not have a column map yet; call fillComplete() first.");
    return graphData_->colMap_; 
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  const Map<LocalOrdinal,GlobalOrdinal> & 
  CrsGraph<LocalOrdinal,GlobalOrdinal>::getDomainMap() const
  {
    TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getDomainMap(): graph does not have a domain map yet; call fillComplete() first.");
    return graphData_->domainMap_; 
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  const Map<LocalOrdinal,GlobalOrdinal> & 
  CrsGraph<LocalOrdinal,GlobalOrdinal>::getRangeMap() const
  {
    TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getRangeMap(): graph does not have a range map yet; call fillComplete() first.");
    return graphData_->rangeMap_; 
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal> >
  CrsGraph<LocalOrdinal,GlobalOrdinal>::getImporter() const
  { 
    TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getImporter(): cannot get importer until fillComplete() has been called.");
    return graphData_->importer_;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal> >
  CrsGraph<LocalOrdinal,GlobalOrdinal>::getExporter() const
  {
    TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getExporter(): cannot get exporter until fillComplete() has been called.");
    return graphData_->exporter_;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal>::hasColMap() const
  { return graphData_->hasColMap_; }

  template <class LocalOrdinal, class GlobalOrdinal>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal>::isStorageOptimized() const
  { return graphData_->storageOptimized_; }

  template <class LocalOrdinal, class GlobalOrdinal>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal>::isStaticProfile() const
  { return graphData_->staticProfile_; }

  template <class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal CrsGraph<LocalOrdinal,GlobalOrdinal>::numGlobalEntries() const
  { return graphData_->numGlobalEntries_; }

  template <class LocalOrdinal, class GlobalOrdinal>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal>::isFillComplete() const
  { return graphData_->fillComplete_; }

  template <class LocalOrdinal, class GlobalOrdinal>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal>::upperTriangular() const
  { return graphData_->upperTriangular_; }

  template <class LocalOrdinal, class GlobalOrdinal>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal>::lowerTriangular() const
  { return graphData_->lowerTriangular_; }

  template <class LocalOrdinal, class GlobalOrdinal>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal>::indicesAreLocal() const
  { return graphData_->indicesAreLocal_; }

  template <class LocalOrdinal, class GlobalOrdinal>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal>::indicesAreGlobal() const
  { return graphData_->indicesAreGlobal_; }

  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsGraph<LocalOrdinal,GlobalOrdinal>::numEntriesForGlobalRow(GlobalOrdinal globalRow) const
  {
    using Teuchos::OrdinalTraits;
    LocalOrdinal rlid = getRowMap().getLocalIndex(globalRow);
    if (rlid == OrdinalTraits<LocalOrdinal>::invalid()) return OrdinalTraits<Teuchos_Ordinal>::invalid();
    return RNNZ(rlid);
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsGraph<LocalOrdinal,GlobalOrdinal>::numEntriesForMyRow(LocalOrdinal localRow) const
  {
    using Teuchos::OrdinalTraits;
    if (!getRowMap().isMyLocalIndex(localRow)) return OrdinalTraits<Teuchos_Ordinal>::invalid();
    return RNNZ(localRow);
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsGraph<LocalOrdinal,GlobalOrdinal>::numAllocatedEntriesForGlobalRow(GlobalOrdinal globalRow) const
  {
    using Teuchos::OrdinalTraits;
    LocalOrdinal rlid = getRowMap().getLocalIndex(globalRow);
    if (rlid == OrdinalTraits<LocalOrdinal>::invalid()) return OrdinalTraits<Teuchos_Ordinal>::invalid();
    return RNumAlloc(rlid);
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsGraph<LocalOrdinal,GlobalOrdinal>::numAllocatedEntriesForMyRow(LocalOrdinal localRow) const
  {
    using Teuchos::OrdinalTraits;
    if (!getRowMap().isMyLocalIndex(localRow)) return OrdinalTraits<Teuchos_Ordinal>::invalid();
    return RNumAlloc(localRow);
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsGraph<LocalOrdinal,GlobalOrdinal>::numMyEntries() const
  { return graphData_->numLocalEntries_; }

  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsGraph<LocalOrdinal,GlobalOrdinal>::myMaxNumRowEntries() const
  { return graphData_->localMaxNumEntries_; }

  template <class LocalOrdinal, class GlobalOrdinal>
  void CrsGraph<LocalOrdinal,GlobalOrdinal>::insertGlobalIndices(GlobalOrdinal grow, const Teuchos::ArrayView<const GlobalOrdinal> &indices)
  {
    // FINISH: add an unchecked, protected version. call it from this method, where this method removes all filtered 
    //         indices, to avoid reallocating/throwing when unnecessary
    using Teuchos::OrdinalTraits;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    TEST_FOR_EXCEPTION(indicesAreLocal() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertGlobalIndices(): graph already contains local indices.");
    if (indicesAreAllocated() == false) {
      bool local = false;
      allocateIndices(local);
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION(indicesAreAllocated() == false, std::logic_error, 
          Teuchos::typeName(*this) << "::insertGlobalIndices(): Internal logic error. Please contact Tpetra team.");
#endif
    }
    graphData_->indicesAreSorted_ = false;
    graphData_->haveGlobalConstants_ = false;
    // insert the indices into the graph
    LocalOrdinal lrow = getRowMap().getLocalIndex(grow);
    if (lrow != OrdinalTraits<LocalOrdinal>::invalid()) {
      // add to allocated space
      Teuchos_Ordinal rowNE = graphData_->rowNumEntries_[lrow],
                      rowNA = RNumAlloc(lrow),
                      toAdd = indices.size();
      if (rowNE+toAdd > rowNA) {
        TEST_FOR_EXCEPTION(isStaticProfile(), std::runtime_error,
            Teuchos::typeName(*this) << "::insertGlobalIndices(): new indices exceed statically allocated graph structure.");
#       if defined(THROW_TPETRA_EFFICIENCY_WARNINGS) || defined(PRINT_TPETRA_EFFICIENCY_WARNINGS)
          std::string err = Teuchos::typeName(*this) + "::insertGlobalIndices(): Pre-allocated space has been exceeded, requiring new allocation. To improve efficiency, suggest larger allocation.";
#         if defined(THROW_TPETRA_EFFICIENCY_WARNINGS)
            TEST_FOR_EXCEPTION(true, std::runtime_error, err);
#         else
            std::cerr << err << std::endl;
#         endif
#       endif
        // increase allocation to necessary amount, copy previous data, reassign ArrayRCP
        ArrayView<const GlobalOrdinal> curInds;
        extractGlobalRowConstView(grow,curInds);
        rowNA = rowNE+toAdd;
        ArrayRCP<GlobalOrdinal> newInds = Teuchos::arcp<GlobalOrdinal>(rowNA);
        std::copy(curInds.begin(), curInds.end(), newInds.begin());
        graphData_->colGInds_[lrow] = newInds;
      }
      typename ArrayRCP<GlobalOrdinal>::iterator dstbeg, dstind;
      dstbeg = getGptr(lrow);
      dstind = dstbeg + rowNE;
      if (hasColMap()) {
        // check the global indices against the column map; only add ones that are defined
        typename ArrayView<const GlobalOrdinal>::iterator srcind = indices.begin();
        while (srcind != indices.end()) 
        {
          if (getColMap().isMyGlobalIndex(*srcind)) {
            (*dstind++) = (*srcind);
          }
          ++srcind;
        }
        graphData_->rowNumEntries_[lrow] = dstind - dstbeg;
      }
      else {
        std::copy( indices.begin(), indices.end(), dstind );
        graphData_->rowNumEntries_[lrow] += toAdd;
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

  template <class LocalOrdinal, class GlobalOrdinal>
  void CrsGraph<LocalOrdinal,GlobalOrdinal>::insertMyIndices(LocalOrdinal lrow, const Teuchos::ArrayView<const LocalOrdinal> &indices)
  {
    // FINISH: add an unchecked, protected version. call it from this method, where this method removes all filtered 
    //         indices, to avoid reallocating/throwing when unnecessary
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(indicesAreGlobal() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertGlobalIndices(): graph already contains local indices.");
    TEST_FOR_EXCEPTION(isFillComplete() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertGlobalIndices(): graph fill already complete; cannot insert new entries.");
    TEST_FOR_EXCEPTION(isStorageOptimized() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertGlobalIndices(): graph storage is already optimized; cannot insert new entires.");
    TEST_FOR_EXCEPTION(hasColMap() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::insertMyIndices(): requires a predefined column map.");
    if (indicesAreAllocated() == false) {
      bool local = true;
      allocateIndices(local);
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION(indicesAreAllocated() == false, std::logic_error, 
          Teuchos::typeName(*this) << "::insertGlobalIndices(): Internal logic error. Please contact Tpetra team.");
#endif
    }
    // insert the indices into the graph
    TEST_FOR_EXCEPTION(getRowMap().isMyLocalIndex(lrow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::insertMyIndices(lrow,...): lrow is not valid for this node.");

    graphData_->indicesAreSorted_ = false;
    graphData_->haveGlobalConstants_ = false;

    // add to allocated space
    Teuchos_Ordinal rowNE = graphData_->rowNumEntries_[lrow],
                    rowNA = RNumAlloc(lrow),
                    toAdd = indices.size();
    if (rowNE+toAdd > rowNA) {
        TEST_FOR_EXCEPTION(isStaticProfile()==true, std::runtime_error,
            Teuchos::typeName(*this) << "::insertMyIndices(): new indices exceed statically allocated graph structure.");
#     if defined(THROW_TPETRA_EFFICIENCY_WARNINGS) || defined(PRINT_TPETRA_EFFICIENCY_WARNINGS)
        std::string err = Teuchos::typeName(*this) + "::insertGlobalIndices(): Pre-allocated space has been exceeded, requiring new allocation. To improve efficiency, suggest larger allocation.";
#       if defined(THROW_TPETRA_EFFICIENCY_WARNINGS)
          TEST_FOR_EXCEPTION(true, std::runtime_error, err);
#       else
          std::cerr << err << std::endl;
#       endif
#     endif
      // increase allocation to necessary amount, copy previous data, reassign ArrayRCP
      ArrayView<const LocalOrdinal> curInds;
      extractMyRowConstView(lrow,curInds);
      rowNA = rowNE+toAdd;             
      ArrayRCP<LocalOrdinal> newInds = Teuchos::arcp<LocalOrdinal>(rowNA);
      std::copy(curInds.begin(), curInds.end(), newInds.begin());
      graphData_->colLInds_[lrow] = newInds;
    }
    // check the local indices against the column map; only add ones that are defined
    typename ArrayView<const LocalOrdinal>::iterator srcind = indices.begin();
    typename ArrayRCP<LocalOrdinal>::iterator dstbeg, dstind;
    dstbeg = getLptr(lrow);
    dstind = dstbeg + rowNE;
    while (srcind != indices.end()) 
    {
      if (getColMap().isMyLocalIndex(*srcind)) {
        (*dstind++) = (*srcind);
      }
      ++srcind;
    }
    graphData_->rowNumEntries_[lrow] = dstind - dstbeg;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP<const Teuchos::Comm<int> > 
  CrsGraph<LocalOrdinal,GlobalOrdinal>::getComm() const
  { return graphData_->comm_; }

  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsGraph<LocalOrdinal,GlobalOrdinal>::getIndexBase() const
  { return graphData_->rowMap_.getIndexBase(); }

  template <class LocalOrdinal, class GlobalOrdinal>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal>::indicesAreSorted() const
  { return graphData_->indicesAreSorted_; }

  template <class LocalOrdinal, class GlobalOrdinal>
  void CrsGraph<LocalOrdinal,GlobalOrdinal>::indicesAreSorted(bool sorted) 
  { graphData_->indicesAreSorted_ = sorted; }

  template <class LocalOrdinal, class GlobalOrdinal>
  void CrsGraph<LocalOrdinal,GlobalOrdinal>::fillComplete(bool OptimizeStorage)
  {
    fillComplete(getRowMap(),getRowMap(),OptimizeStorage);
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal>
  void CrsGraph<LocalOrdinal,GlobalOrdinal>::optimizeStorage()
  {
    TEST_FOR_EXCEPT(true);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal>
  void CrsGraph<LocalOrdinal,GlobalOrdinal>::globalAssemble()
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
    Teuchos_Ordinal MyNonlocals = nonlocals_.size(), 
                    MaxGlobalNonlocals;
    Teuchos::reduceAll<Teuchos_Ordinal>(*getComm(),Teuchos::REDUCE_MAX,MyNonlocals,&MaxGlobalNonlocals);
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
    Teuchos_Ordinal numRecvs = recvIDs.size();

    // we know how many we're sending to already
    // form a contiguous list of all data to be sent
    // track the number of entries for each ID
    Array<pair<GlobalOrdinal,GlobalOrdinal> > IJSendBuffer;
    Array<Teuchos_Ordinal> sendSizes(sendIDs.size(), 0);
    Teuchos_Ordinal numSends = 0;
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
    Array<Teuchos_Ordinal> recvSizes(numRecvs);
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
      Teuchos_Ordinal cur = 0;
      for (Teuchos_Ordinal s=0; s<numSends; ++s) {
        sendBuffers[s] = IJSendBuffer(cur,sendSizes[s]);
        cur += sendSizes[s];
      }
    }
    // perform non-blocking sends
    for (Teuchos_Ordinal s=0; s < numSends ; ++s)
    {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      ArrayRCP<pair<GlobalOrdinal,GlobalOrdinal> > tmparcp = arcp(sendBuffers[s].getRawPtr(),0,sendBuffers[s].size(),false);
      sendRequests.push_back( Teuchos::isend<int,pair<GlobalOrdinal,GlobalOrdinal> >(*getComm(),tmparcp,sendIDs[s]) );
    }
    // calculate amount of storage needed for receives
    // setup pointers for the receives as well
    Teuchos_Ordinal totalRecvSize = std::accumulate(recvSizes.begin(),recvSizes.end(),0);
    Array<pair<GlobalOrdinal,GlobalOrdinal> > IJRecvBuffer(totalRecvSize);
    // from the size info, build the ArrayViews into IJRecvBuffer
    Array<ArrayView<pair<GlobalOrdinal,GlobalOrdinal> > > recvBuffers(numRecvs,Teuchos::null);
    {
      Teuchos_Ordinal cur = 0;
      for (Teuchos_Ordinal r=0; r<numRecvs; ++r) {
        recvBuffers[r] = IJRecvBuffer(cur,recvSizes[r]);
        cur += recvSizes[r];
      }
    }
    // perform non-blocking recvs
    for (Teuchos_Ordinal r=0; r < numRecvs ; ++r)
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
  template <class LocalOrdinal, class GlobalOrdinal>
  void CrsGraph<LocalOrdinal,GlobalOrdinal>::fillComplete(const Map<LocalOrdinal,GlobalOrdinal> &domainMap, const Map<LocalOrdinal,GlobalOrdinal> &rangeMap, bool OptimizeStorage)
  {
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(isFillComplete() == true, std::runtime_error,
      Teuchos::typeName(*this) << "::fillComplete(): fillComplete() has already been called.");

    Teuchos::ArrayView<const GlobalOrdinal> myGlobalEntries = graphData_->rowMap_.getMyGlobalEntries();
    graphData_->domainMap_ = domainMap;
    graphData_->rangeMap_  = rangeMap;

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
    if (graphData_->haveGlobalConstants_ == false) {
      GlobalOrdinal lcl[2], gbl[2];
      lcl[0] = graphData_->numLocalEntries_;
      lcl[1] = graphData_->numLocalDiags_;
      Teuchos::reduceAll<int,GlobalOrdinal>(*getComm(),Teuchos::REDUCE_SUM,2,lcl,gbl);
      graphData_->numGlobalEntries_ = gbl[0]; 
      graphData_->numGlobalDiags_   = gbl[1];
      Teuchos::reduceAll<int,GlobalOrdinal>(*getComm(),Teuchos::REDUCE_MAX,graphData_->localMaxNumEntries_,&graphData_->globalMaxNumRowEntries_);
      graphData_->haveGlobalConstants_ = true;
    }

    // mark transformation as successfully completed
    graphData_->fillComplete_ = true;

    if (OptimizeStorage) optimizeStorage();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal>
  void CrsGraph<LocalOrdinal,GlobalOrdinal>::makeIndicesLocal(const Map<LocalOrdinal,GlobalOrdinal> &domainMap, const Map<LocalOrdinal,GlobalOrdinal> &rangeMap)
  {
    using Teuchos::ArrayRCP;
    computeIndexState(); // Update index state by checking IndicesAreLocal/Global on all nodes
    TEST_FOR_EXCEPTION(indicesAreLocal() && indicesAreGlobal(), std::logic_error,
        Teuchos::typeName(*this) << "::makeIndicesLocal(): indices are marked as both global and local.");

    makeColMap(domainMap, rangeMap); // If user has not prescribed column map, create one from indices

    // Transform indices to local index space
    const Teuchos_Ordinal nlrs = numLocalRows();

    if (indicesAreGlobal() && indicesAreAllocated()) {
      // allocate data for local indices
      if (isStaticProfile()) {
        graphData_->colLIndsPtrs_.resize(numLocalRows()+1);
        if (graphData_->contigColGInds_ != Teuchos::null) {
          graphData_->contigColLInds_ = Teuchos::arcp_reinterpret_cast<LocalOrdinal>(graphData_->contigColGInds_);
          Teuchos_Ordinal sofar = 0, r = 0;
          while (r < numLocalRows()) {
            // RNumAlloc(r) uses colGIndsPtrs_[r+1] and colGIndsPtrs_[r], can't erase till done
            Teuchos_Ordinal numalloc = RNumAlloc(r),
                            numentry = RNNZ(r);
            graphData_->colLIndsPtrs_[r] = graphData_->contigColLInds_.persistingView(sofar,numalloc).begin();
            for (Teuchos_Ordinal j=0; j<numentry; ++j) {
              graphData_->colLIndsPtrs_[r][j] = getColMap().getLocalIndex( graphData_->colGIndsPtrs_[r][j] );
#ifdef HAVE_TPETRA_DEBUG
              TEST_FOR_EXCEPTION(graphData_->colLIndsPtrs_[r][j] == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::logic_error,
                  Teuchos::typeName(*this) << ": Internal error in fillComplete(). Please contact Tpetra team.");
#endif
            }
            graphData_->colGIndsPtrs_[r] = Teuchos::null;
            ++r;
            sofar += numalloc;
          }
          graphData_->colLIndsPtrs_[r] = graphData_->contigColLInds_.end();
#ifdef HAVE_TPETRA_DEBUG
          TEST_FOR_EXCEPTION(sofar != graphData_->contigColGInds_.size(), std::logic_error,
              Teuchos::typeName(*this) << "::makeIndicesLocal(): Internal logic error. Please contact Tpetra team.");
#endif
        }
        graphData_->colGIndsPtrs_.clear();
        graphData_->contigColGInds_ = Teuchos::null;
      }
      else {
        graphData_->colLInds_.resize(nlrs);
        for (Teuchos_Ordinal r=0; r < nlrs; ++r)
        {
          graphData_->colLInds_[r] = Teuchos::arcp_reinterpret_cast<LocalOrdinal>(graphData_->colGInds_[r]);
          typename ArrayRCP<GlobalOrdinal>::const_iterator 
            cindG = graphData_->colGInds_[r].begin(),
            cendG = graphData_->colGInds_[r].begin() + graphData_->rowNumEntries_[r];
          typename ArrayRCP<LocalOrdinal>::iterator cindL = graphData_->colLInds_[r].begin();
          for (; cindG != cendG; ++cindG, ++cindL)
          {
            (*cindL) = graphData_->colMap_.getLocalIndex(*cindG);
#ifdef HAVE_TPETRA_DEBUG
            TEST_FOR_EXCEPTION((*cindL) == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::logic_error,
                Teuchos::typeName(*this) << ": Internal error in fillComplete(). Please contact Tpetra team.");
#endif
          }
          graphData_->colGInds_[r] = Teuchos::null;   // data is owned by colLInds_[r] now
        }
      }
    }
    graphData_->indicesAreLocal_  = true;
    graphData_->indicesAreGlobal_ = false;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal>
  void CrsGraph<LocalOrdinal,GlobalOrdinal>::sortIndices()
  {
    using Teuchos::ArrayView;
    TEST_FOR_EXCEPT(indicesAreGlobal()==true);   // this should be called only after makeIndicesLocal()
    if (indicesAreSorted()) return;
    if (indicesAreAllocated()) {                 // are there any indices?
      const Teuchos_Ordinal nlrs = numLocalRows();
      for (Teuchos_Ordinal r=0; r < nlrs; ++r)
      {
        Teuchos_Ordinal rnnz = RNNZ(r);
        typename Teuchos::ArrayRCP<LocalOrdinal>::iterator rowbeg = getLptr(r);
        std::sort(rowbeg, rowbeg + rnnz);
      }
    }
    graphData_->indicesAreSorted_ = true;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal>
  void CrsGraph<LocalOrdinal,GlobalOrdinal>::makeColMap(const Map<LocalOrdinal,GlobalOrdinal> &domainMap, const Map<LocalOrdinal,GlobalOrdinal> &rangeMap)
  {
    using Teuchos::ArrayRCP;
    using Teuchos::Array;
    using Teuchos::ArrayView;
    typedef Teuchos::OrdinalTraits<GlobalOrdinal> GOT;
    const Teuchos_Ordinal nlrs = numLocalRows();
    if (hasColMap()) return;
    computeIndexState();
    TEST_FOR_EXCEPTION(indicesAreLocal() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::makeColMap(): indices must still be global when making the column map.");

    // ultimate: list of indices for column map
    Array<GlobalOrdinal> myColumns;
    if (indicesAreGlobal() == true) { // otherwise, we have no indices to bother with
      // Construct two lists of columns for which we have a non-zero
      // Local GIDs: Column GIDs which are locally present on the Domain map
      // Remote GIDs: Column GIDs which are not locally present present on the Domain map
      //
      // instead of list of local GIDs, we will use a vector<bool> LocalGID, where LocalGID[LID]
      // indicates that, for LID=DomainMap.LID(GID), GID is local. on modern compilers, vector<bool>
      // should be a partial specialization of vector<> which represents each entry via a single
      // bit, instead of a full byte.
      const LocalOrdinal LINV = Teuchos::OrdinalTraits<LocalOrdinal>::invalid();
      Teuchos_Ordinal numLocalColGIDs = 0, numRemoteColGIDs = 0;
      //
      // intitial: partitioning into local and remote
      std::vector<bool>       GIDisLocal(domainMap.getNumMyEntries());
      std::set<GlobalOrdinal> RemoteGIDSet;
      for (Teuchos_Ordinal r=0; r < nlrs; ++r)
      {
        typename ArrayRCP<GlobalOrdinal>::const_iterator cind, cend;
        cind = getGptr(r);
        cend = cind + graphData_->rowNumEntries_[r];
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
        if (numLocalColGIDs == domainMap.getNumMyEntries()) {
          graphData_->colMap_ = domainMap;
          graphData_->hasColMap_ = true;
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
      if (numLocalColGIDs == domainMap.getNumMyEntries()) {
        std::copy(mge.begin(), mge.end(), LocalColGIDs.begin());
      }
      else {
        Teuchos_Ordinal numlocalagain = 0;
        for (Teuchos_Ordinal i=0; i<domainMap.getNumMyEntries(); ++i) {
          if (GIDisLocal[i]) {
            LocalColGIDs[numlocalagain++] = mge[i];
          }
        }
        TEST_FOR_EXCEPTION(numlocalagain != numLocalColGIDs, std::logic_error,
            Teuchos::typeName(*this) << "::makeColMap(): Internal logic error. Please contact Tpetra team.");
      }
    }
    graphData_->colMap_ = Map<LocalOrdinal,GlobalOrdinal>(GOT::invalid(), myColumns, domainMap.getIndexBase(), domainMap.getComm());
    graphData_->hasColMap_ = true;
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal>
  void CrsGraph<LocalOrdinal,GlobalOrdinal>::removeRedundantIndices() 
  {
    TEST_FOR_EXCEPT(indicesAreGlobal()==true);    // this should be called only after makeIndicesLocal()
    TEST_FOR_EXCEPT(indicesAreSorted()==false);   // this should be called only after sortIndices()
    TEST_FOR_EXCEPT(isStorageOptimized()==true);  // assumptions about available data structures
    const Teuchos_Ordinal nlrs = numLocalRows();
    Teuchos::ArrayView<const GlobalOrdinal> myGlobalEntries = getRowMap().getMyGlobalEntries();
    // reset all local quantities
    graphData_->upperTriangular_ = true;
    graphData_->lowerTriangular_ = true;
    graphData_->localMaxNumEntries_ = 0;
    graphData_->numLocalEntries_    = 0;
    graphData_->numLocalDiags_      = 0;
    // indices are already sorted in each row
    if (indicesAreAllocated()) {
      for (Teuchos_Ordinal r=0; r < nlrs; ++r)
      {
        GlobalOrdinal rgid = myGlobalEntries[r];
        LocalOrdinal rlcid = getColMap().getLocalIndex(rgid);   // determine the local column index for this row, used for delimiting the diagonal
        typename Teuchos::ArrayRCP<LocalOrdinal>::iterator beg, end;
        beg = getLptr(r);
        end = beg + graphData_->rowNumEntries_[r];
        typename Teuchos::ArrayRCP<LocalOrdinal>::iterator newend, cur;
        newend = beg;
        if (beg != end) {
          if (rlcid == (*newend)) ++graphData_->numLocalDiags_;
          for (cur = beg + 1; cur != end; ++cur) {
            // it is unique; add it and test it for the diagonal
            // otherwise, it is a dup and should be ignored
            if (*cur != *newend) {
              ++newend;
              (*newend) = (*cur);
              // is this the diagonal?
              if (rlcid == (*newend)) ++graphData_->numLocalDiags_;
            }
          }
          // because of sorting, smallest column index is (*beg); it indicates upper triangularity
          if ((*beg) < r) graphData_->upperTriangular_ = false;
          // because of sorting, largest column index is (*newend); it indicates lower triangularity
          if (r < (*newend)) graphData_->lowerTriangular_ = false;
          // increment newend so that our range is [beg,newend) instead of [beg,newend]
          ++newend;
        }
        // compute num entries for this row, accumulate into numLocalEntries_, update localMaxNumEntries_
        graphData_->rowNumEntries_[r] = newend - beg;
        graphData_->localMaxNumEntries_ = std::max( graphData_->localMaxNumEntries_, Teuchos::as<Teuchos_Ordinal>(graphData_->rowNumEntries_[r]) );
        graphData_->numLocalEntries_ += graphData_->rowNumEntries_[r];
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal>
  void CrsGraph<LocalOrdinal,GlobalOrdinal>::makeImportExport()
  {
    TEST_FOR_EXCEPT(hasColMap()==false); // must have column map
    // create import, export
    if (!graphData_->domainMap_.isSameAs(graphData_->colMap_)) {
      graphData_->importer_ = Teuchos::rcp( new Import<LocalOrdinal,GlobalOrdinal>(graphData_->domainMap_,graphData_->colMap_) );
    }
    else {
      graphData_->importer_ = Teuchos::null;
    }
    if (!graphData_->rangeMap_.isSameAs(graphData_->rowMap_)) {
      graphData_->exporter_ = Teuchos::rcp( new Export<LocalOrdinal,GlobalOrdinal>(graphData_->rowMap_,graphData_->rangeMap_) );
    }
    else {
      graphData_->exporter_ = Teuchos::null;
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal>
  void CrsGraph<LocalOrdinal,GlobalOrdinal>::computeIndexState()
  {
    char myIndicesAreLocal = 0;
    char myIndicesAreGlobal = 0;
    if(graphData_->indicesAreLocal_) 
      myIndicesAreLocal = 1;
    if(graphData_->indicesAreGlobal_) 
      myIndicesAreGlobal = 1;
    char allIndicesAreLocal;
    char allIndicesAreGlobal;
    Teuchos::reduceAll(*getComm(),Teuchos::REDUCE_MAX,myIndicesAreLocal ,&allIndicesAreLocal );
    Teuchos::reduceAll(*getComm(),Teuchos::REDUCE_MAX,myIndicesAreGlobal,&allIndicesAreGlobal );
    graphData_->indicesAreLocal_  = (allIndicesAreLocal ==1);  // If indices are local on one PE, should be local on all
    graphData_->indicesAreGlobal_ = (allIndicesAreGlobal==1);  // If indices are global on one PE should be local on all
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal>::indicesAreAllocated() const
  { return graphData_->indicesAreAllocated_; }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsGraph<LocalOrdinal,GlobalOrdinal>::RNNZ(Teuchos_Ordinal row) const
  {
    // if storage is optimized, then numalloc == numnz, and indices are local
    if (isStorageOptimized()) {
      return (graphData_->colLIndsPtrs_[row+1] - graphData_->colLIndsPtrs_[row]);
    }
    // otherwise, we are still numnz
    else {
      return graphData_->rowNumEntries_[row];
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsGraph<LocalOrdinal,GlobalOrdinal>::RNumAlloc(Teuchos_Ordinal row) const
  {
    if (indicesAreAllocated() == false) {
      return graphData_->rowNumToAlloc_[row];
    }
    // if static graph or optimize storage, pointers tell us the allocation size
    else if (isStorageOptimized() || isStaticProfile()) {
      if (indicesAreLocal()) {
        return graphData_->colLIndsPtrs_[row+1] - graphData_->colLIndsPtrs_[row];
      }
      else {
        return graphData_->colGIndsPtrs_[row+1] - graphData_->colGIndsPtrs_[row];
      }
    }
    // otherwise, the ArrayRCP knows
    else {
      if (indicesAreLocal()) {
        return graphData_->colLInds_[row].size();
      }
      else {
        return graphData_->colGInds_[row].size();
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal>
  void CrsGraph<LocalOrdinal,GlobalOrdinal>::staticAssertions()
  {
    using Teuchos::OrdinalTraits;
    // Assumption: sizeof(GlobalOrdinal) >= sizeof(LocalOrdinal)
    //    This is so that we can store LocalOrdinals in the memory formerly occupied by GlobalOrdinals
    // Assumption: max(GlobalOrdinal) >= max(Teuchos_Ordinal) >= max(LocalOrdinal)
    //    This is so that we can represent any LocalOrdinal as a Teuchos_Ordinal, and any Teuchos_Ordinal as a GlobalOrdinal, transitively
    Teuchos::CompileTimeAssert<sizeof(GlobalOrdinal) < sizeof(LocalOrdinal)> cta_size;
    (void)cta_size;
    // can't call max() with CompileTimeAssert, because it isn't a constant expression; will need to make this a runtime check
    TEST_FOR_EXCEPTION( OrdinalTraits<GlobalOrdinal>::max() < OrdinalTraits<LocalOrdinal>::max(), std::runtime_error, 
        Teuchos::typeName(*this) << ": Object cannot be allocated with stated template arguments: size assumptions are not valid.");
    TEST_FOR_EXCEPTION( OrdinalTraits<Teuchos_Ordinal>::max() < OrdinalTraits<LocalOrdinal>::max(), std::runtime_error,
        Teuchos::typeName(*this) << ": Object cannot be allocated with stated template arguments: size assumptions are not valid.");
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  typename Teuchos::ArrayRCP<const LocalOrdinal>::iterator 
  CrsGraph<LocalOrdinal,GlobalOrdinal>::getLptr(Teuchos_Ordinal row) const
  {
    if (isStaticProfile() || isFillComplete()) {
      return graphData_->colLIndsPtrs_[row];
    }
    else {
      return graphData_->colLInds_[row].begin();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  typename Teuchos::ArrayRCP<LocalOrdinal>::iterator 
  CrsGraph<LocalOrdinal,GlobalOrdinal>::getLptr(Teuchos_Ordinal row)
  {
    if (isStaticProfile() || isFillComplete()) {
      return graphData_->colLIndsPtrs_[row];
    }
    else {
      return graphData_->colLInds_[row].begin();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  typename Teuchos::ArrayRCP<const GlobalOrdinal>::iterator 
  CrsGraph<LocalOrdinal,GlobalOrdinal>::getGptr(Teuchos_Ordinal row) const
  {
    if (isStaticProfile() || isFillComplete()) {
      return graphData_->colGIndsPtrs_[row];
    }
    else {
      return graphData_->colGInds_[row].begin();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  typename Teuchos::ArrayRCP<GlobalOrdinal>::iterator 
  CrsGraph<LocalOrdinal,GlobalOrdinal>::getGptr(Teuchos_Ordinal row) 
  {
    if (isStaticProfile() || isFillComplete()) {
      return graphData_->colGIndsPtrs_[row];
    }
    else {
      return graphData_->colGInds_[row].begin();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  void
  CrsGraph<LocalOrdinal,GlobalOrdinal>::print(std::ostream &os) const
  {
    // support this operation whether fillComplete() or not
    using std::endl;
    int myImageID = Teuchos::rank(*getComm());
    if (myImageID == 0)
    {
      os << "Tpetra::CrsGraph, label = " << this->getObjectLabel() << endl;
      os << "Number of global rows   = " << numGlobalRows() << endl;
      if (isFillComplete())
      {
        os << "Number of global columns    = " << numGlobalCols() << endl;
        os << "Status = fill complete" << endl;
        os << "Number of global nonzeros   = " << numGlobalEntries() << endl;
        os << "Global max nonzeros per row = " << globalMaxNumRowEntries() << endl;
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
          Teuchos::Array<GlobalOrdinal> indices(myMaxNumRowEntries());
          Teuchos::ArrayView<const GlobalOrdinal> mgis = getRowMap().getMyGlobalEntries();
          os << "% Number of rows on image " << myImageID << " = " << numLocalRows() << endl;
          for (Teuchos_Ordinal r=0; r < numLocalRows(); ++r)
          {                                                        
            GlobalOrdinal globalRow = mgis[r];
            Teuchos_Ordinal rowSize;
            extractGlobalRowCopy(globalRow, indices(), rowSize);
            if (rowSize > 0) {
              os << "Row " << globalRow << ":";
              for (Teuchos_Ordinal j=0; j < rowSize; ++j) {
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
