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

#include <Teuchos_Describable.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_NullIteratorTraits.hpp>

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_CrsGraph.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_RowGraph.hpp"
#include "Tpetra_Util.hpp"

namespace Tpetra 
{

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration
  template <class S, class LO, class GO, class N, class SpMV>
  class CrsMatrix;
#endif

  //! \brief A class for constructing and using sparse compressed index graphs with row access.
  /*! This class is templated on \c LocalOrdinal and \c GlobalOrdinal. If the \c GlobalOrdinal is not specified, then 
   *  it takes the same type as the \c LocalOrdinal.
   */
  template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class CrsGraph : public RowGraph<LocalOrdinal,GlobalOrdinal,Node> {
    template <class S, class LO, class GO, class N, class SpMV>
    friend class CrsMatrix;

    public: 

      //! @name Constructor/Destructor Methods
      //@{ 

      //! Constructor with fixed number of indices per row.
      CrsGraph(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype = DynamicProfile);

      //! Constructor with variable number of indices per row.
      CrsGraph(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const Teuchos::ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, ProfileType pftype = DynamicProfile);

      //! Constructor with fixed number of indices per row and specified column map.
      CrsGraph(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, size_t maxNumEntriesPerRow, ProfileType pftype = DynamicProfile);

      //! Constructor with variable number of indices per row and specified column map.
      CrsGraph(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, const Teuchos::ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, ProfileType pftype = DynamicProfile);

      // !Destructor.
      virtual ~CrsGraph();

      //@}

      //! @name Insertion/Removal Methods
      //@{ 

      //! Submit graph indices, using global IDs.
      void insertGlobalIndices(GlobalOrdinal row, const Teuchos::ArrayView<const GlobalOrdinal> &indices);

      //! Submit graph indices, using local IDs.
      void insertLocalIndices(LocalOrdinal row, const Teuchos::ArrayView<const LocalOrdinal> &indices);

      //@}

      //! @name Transformational Methods
      //@{ 

      //! \brief Communicate non-local contributions to other nodes.
      void globalAssemble();

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

      //! \brief Re-allocate the data into contiguous storage.
      void optimizeStorage();

      //@}

      //! @name Methods implementing RowGraph.
      //@{ 

      //! Returns the communicator.
      const Teuchos::RCP<const Teuchos::Comm<int> > & getComm() const;

      //! Returns the underlying node.
      Teuchos::RCP<Node> getNode() const;

      //! Returns the Map that describes the row distribution in this graph.
      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRowMap() const;

      //! \brief Returns the Map that describes the column distribution in this graph.
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
      global_size_t getGlobalNumRows() const;

      //! \brief Returns the number of global columns in the graph.
      global_size_t getGlobalNumCols() const;

      //! Returns the number of rows owned on the calling node.
      size_t getNodeNumRows() const;

      //! Returns the number of columns connected to the locally owned rows of this graph.
      size_t getNodeNumCols() const;

      //! Returns the index base for global indices for this graph. 
      GlobalOrdinal getIndexBase() const;

      //! Returns the global number of entries in the graph.
      global_size_t getGlobalNumEntries() const;

      //! Returns the local number of entries in the graph.
      size_t getNodeNumEntries() const;

      //! \brief Returns the current number of entries on this node in the specified global row.
      /*! Returns Teuchos::OrdinalTraits<size_t>::invalid() if the specified global row does not belong to this graph. */
      size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const;

      //! Returns the current number of entries on this node in the specified local row.
      /*! Returns Teuchos::OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this graph. */
      size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const;

      //! \brief Returns the total number of indices allocated for the graph, across all rows on this node.
      /*! This is the allocation available to the user. Actual allocation may be larger, for example, after 
          calling fillComplete(), and thus this does not necessarily reflect the memory consumption of the 
          this graph.  */
      size_t getNodeAllocationSize() const;

      //! \brief Returns the current number of allocated entries for this node in the specified global row .
      /*! Throws exception std::runtime_error if the specified global row does not belong to this node. */
      size_t getNumAllocatedEntriesInGlobalRow(GlobalOrdinal globalRow) const;

      //! Returns the current number of allocated entries on this node in the specified local row.
      /*! Throws exception std::runtime_error if the specified local row is not valid for this node. */
      size_t getNumAllocatedEntriesInLocalRow(LocalOrdinal localRow) const;

      //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons. 
      global_size_t getGlobalNumDiags() const;

      //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons. 
      size_t getNodeNumDiags() const;

      //! \brief Returns the maximum number of entries across all rows/columns on all nodes. 
      size_t getGlobalMaxNumRowEntries() const;

      //! \brief Returns the maximum number of entries across all rows/columns on this node. 
      size_t getNodeMaxNumRowEntries() const;

      //! \brief Indicates whether the matrix has a well-defined column map. 
      bool hasColMap() const; 

      //! \brief Indicates whether the graph is lower triangular.
      bool isLowerTriangular() const;

      //! \brief Indicates whether the graph is upper triangular.
      bool isUpperTriangular() const;

      //! \brief If graph indices are in the local range, this function returns true. Otherwise, this function returns false. */
      bool isLocallyIndexed() const;

      //! \brief If graph indices are in the global range, this function returns true. Otherwise, this function returns false. */
      bool isGloballyIndexed() const;

      //! Returns \c true if fillComplete() has been called.
      bool isFillComplete() const;

      //! Extract a list of elements in a specified global row of the graph. Put into pre-allocated storage.
      /*!
        \param LocalRow - (In) Global row number for which indices are desired.
        \param Indices - (Out) Global column indices corresponding to values.
        \param NumIndices - (Out) Number of indices.

         Note: A std::runtime_error exception is thrown indices is not large enough to hold the column indices associated
         with row \c GlobalRow. If \c GlobalRow does not belong to this node, then \c indices is unchanged and \c NumIndices is 
         returned as Teuchos::OrdinalTraits<size_t>::invalid().
       */
      void getGlobalRowCopy(GlobalOrdinal GlobalRow, 
                            const Teuchos::ArrayView<GlobalOrdinal> &Indices, 
                            size_t &NumIndices) const;

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

      //! @name Miscellaneous Query Methods
      //@{

      //! Returns \c true if optimizeStorage() has been called.
      bool isStorageOptimized() const;

      //! Returns \c true if the graph data was allocated in static data structures.
      ProfileType getProfileType() const;

      //@}

      //! @name Overridden from Teuchos::Describable 
      //@{

      /** \brief Return a simple one-line description of this object. */
      std::string description() const;

      /** \brief Print the object with some verbosity level to an FancyOStream object. */
      void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

      //@}


    private:
      // copy constructor disabled
      CrsGraph(const CrsGraph<LocalOrdinal,GlobalOrdinal,Node> &Source);
      // operator= disabled
      CrsGraph<LocalOrdinal,GlobalOrdinal,Node> & operator=(const CrsGraph<LocalOrdinal,GlobalOrdinal,Node> &rhs);
    protected:
      enum AllocateLocalGlobal {
        AllocateLocal,
        AllocateGlobal
      };
      void allocateIndices(AllocateLocalGlobal lg);
      void makeIndicesLocal(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap);
      void makeColMap(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap);
      void computeIndexState();
      void sortIndices();
      void removeRedundantIndices();
      void makeImportExport();
      bool notRedundant() const;
      bool isSorted() const;
      void setSorted(bool sorted);
      bool indicesAreAllocated() const;
      void staticAssertions();
      size_t findMyIndex(size_t row, LocalOrdinal ind) const;
      size_t findGlobalIndex(size_t row, GlobalOrdinal ind) const;
      inline size_t RNNZ(size_t row) const;
      inline size_t RNumAlloc(size_t row) const;
      void checkInternalState() const;
      void clearGlobalConstants();
      void updateLocalAllocation(size_t lrow, size_t allocSize);
      void updateGlobalAllocation(size_t lrow, size_t allocSize);

      //! Get a persisting non-const view of the elements in a specified global row of the graph.
      /*!
        \param GlobalRow - (In) Global row number to get indices.

         Note: If \c GlobalRow does not belong to this node, then returns <tt>Teuchos::null</tt>.

        \pre isGloballyIndexed()==true
       */
      Teuchos::ArrayRCP<GlobalOrdinal> getFullGlobalRowView(GlobalOrdinal GlobalRow);

      //! Get a persisting non-const view of the elements in a specified local row of the graph.
      /*!
        \param LocalRow - (In) Local row number to get indices.

         Note: If \c LocalRow is not valid for this node, then returns <tt>Teuchos::null</tt>.

        \pre isLocallyIndexed()==true
       */
      Teuchos::ArrayRCP<LocalOrdinal> getFullLocalRowView(LocalOrdinal LocalRow);

      // Tpetra support objects
      Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > rowMap_, colMap_, rangeMap_, domainMap_;
      Teuchos::RCP<Import<LocalOrdinal,GlobalOrdinal,Node> > importer_;
      Teuchos::RCP<Export<LocalOrdinal,GlobalOrdinal,Node> > exporter_;

      // local data, stored in a Kokkos::CrsGraph. only initialized after fillComplete()
      Kokkos::CrsGraph<LocalOrdinal,Node> lclGraph_;

      // Local and Global Counts
      global_size_t globalNumEntries_, globalNumDiags_, globalMaxNumRowEntries_;
      size_t          nodeNumEntries_,   nodeNumDiags_,   nodeMaxNumRowEntries_, nodeNumAllocated_;

      // allocate static or dynamic?
      ProfileType pftype_;
      // number of non-zeros to allocate per row; set to Teuchos::null after they are allocated.
      Teuchos::ArrayRCP<const size_t> numAllocPerRow_;
      // number of non-zeros to allocate for all row; either this or numAllocPerRow_ is used, but not both.
      size_t numAllocForAllRows_;

      // number of valid entries in a row; this may be less than the number of allocated per row. Used for static and dynamic allocations.
      // set to Teuchos::null after optimizeStorage(), because is no longer needed (all entries are valid)
      // set to null if nodeNumEntries_ == 0
      Teuchos::ArrayRCP<size_t>       numEntriesPerRow_;

      // graph indices. before allocation, both are Teuchos::null. 
      // after allocation, except during makeIndicesLocal(), one of these is Teuchos::Null.
      // this is a parallel compute buffer, not host memory
      // 1D == StaticAllocation, 2D == DynamicAllocation
      Teuchos::ArrayRCP< LocalOrdinal>                     pbuf_lclInds1D_;
      Teuchos::ArrayRCP<GlobalOrdinal>                     pbuf_gblInds1D_;
      Teuchos::ArrayRCP<Teuchos::ArrayRCP< LocalOrdinal> > pbuf_lclInds2D_;
      Teuchos::ArrayRCP<Teuchos::ArrayRCP<GlobalOrdinal> > pbuf_gblInds2D_;
      // offset to the first entry of each row. length numRows + 1. only used if pftype_ == StaticAllocation, otherwise set to Teuchos::null.
      // also, set to null if nodeNumAllocated_ == 0
      Teuchos::ArrayRCP<size_t>                            pbuf_rowOffsets_;

      bool indicesAreAllocated_,
           indicesAreLocal_,
           indicesAreGlobal_,
           fillComplete_, 
           storageOptimized_,
           lowerTriangular_,
           upperTriangular_,
           indicesAreSorted_,
           haveGlobalConstants_, 
           noRedundancies_;

      // non-local data
      std::map<GlobalOrdinal, std::deque<GlobalOrdinal> > nonlocals_;
  };

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::CrsGraph(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype)
  : rowMap_(rowMap)
  , lclGraph_(rowMap->getNode())
  , nodeNumEntries_(0)
  , nodeNumDiags_(0)
  , nodeMaxNumRowEntries_(0)
  , nodeNumAllocated_(0)
  , pftype_(pftype)
  , numAllocForAllRows_(maxNumEntriesPerRow)
  , indicesAreAllocated_(false)
  , indicesAreLocal_(false)
  , indicesAreGlobal_(false)
  , fillComplete_(false)
  , storageOptimized_(false)
  , lowerTriangular_(false)
  , upperTriangular_(false)
  , indicesAreSorted_(true) 
  , noRedundancies_(false) {
    staticAssertions();
    TEST_FOR_EXCEPTION(maxNumEntriesPerRow > Teuchos::OrdinalTraits<size_t>::max() || (maxNumEntriesPerRow < 1 && maxNumEntriesPerRow != 0), std::runtime_error,
        Teuchos::typeName(*this) << "::CrsGraph(rowMap,maxNumEntriesPerRow): maxNumEntriesPerRow must be non-negative.");
    clearGlobalConstants();
    checkInternalState();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::CrsGraph(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, 
                                                      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, 
                                                      size_t maxNumEntriesPerRow, ProfileType pftype)
  : rowMap_(rowMap)
  , colMap_(colMap)
  , lclGraph_(rowMap->getNode())
  , nodeNumEntries_(0)
  , nodeNumDiags_(0)
  , nodeMaxNumRowEntries_(0)
  , nodeNumAllocated_(0)
  , pftype_(pftype)
  , numAllocForAllRows_(maxNumEntriesPerRow) 
  , indicesAreAllocated_(false)
  , indicesAreLocal_(false)
  , indicesAreGlobal_(false)
  , fillComplete_(false)
  , storageOptimized_(false)
  , lowerTriangular_(false)
  , upperTriangular_(false)
  , indicesAreSorted_(true) 
  , noRedundancies_(false) {
    staticAssertions();
    TEST_FOR_EXCEPTION(maxNumEntriesPerRow > Teuchos::OrdinalTraits<size_t>::max() || (maxNumEntriesPerRow < 1 && maxNumEntriesPerRow != 0), std::runtime_error,
        Teuchos::typeName(*this) << "::CrsGraph(rowMap,colMap,maxNumEntriesPerRow): maxNumEntriesPerRow must be non-negative.");
    clearGlobalConstants();
    checkInternalState();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::CrsGraph(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, 
                                                      const Teuchos::ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, ProfileType pftype)
  : rowMap_(rowMap)
  , lclGraph_(rowMap->getNode())
  , nodeNumEntries_(0)
  , nodeNumDiags_(0)
  , nodeMaxNumRowEntries_(0)
  , nodeNumAllocated_(0)
  , pftype_(pftype)
  , numAllocPerRow_(NumEntriesPerRowToAlloc) 
  , numAllocForAllRows_(0)
  , indicesAreAllocated_(false)
  , indicesAreLocal_(false)
  , indicesAreGlobal_(false)
  , fillComplete_(false)
  , storageOptimized_(false)
  , lowerTriangular_(false)
  , upperTriangular_(false)
  , indicesAreSorted_(true) 
  , noRedundancies_(false) {
    staticAssertions();
    TEST_FOR_EXCEPTION((size_t)NumEntriesPerRowToAlloc.size() != getNodeNumRows(), std::runtime_error,
        Teuchos::typeName(*this) << "::CrsGraph(rowMap,NumEntriesPerRowToAlloc): NumEntriesPerRowToAlloc must have as many entries as specified by rowMap for this node.");
    size_t numMin = Teuchos::OrdinalTraits<size_t>::max(),
           numMax = Teuchos::OrdinalTraits<size_t>::zero();
    for (size_t r=0; r < getNodeNumRows(); ++r) {
      numMin = std::min<size_t>( numMin, NumEntriesPerRowToAlloc[r] );
      numMax = std::max<size_t>( numMax, NumEntriesPerRowToAlloc[r] );
    }
    TEST_FOR_EXCEPTION((numMin < Teuchos::OrdinalTraits<size_t>::one() && numMin != Teuchos::OrdinalTraits<size_t>::zero()) || numMax > Teuchos::OrdinalTraits<size_t>::max(),
        std::runtime_error, Teuchos::typeName(*this) << "::allocateIndices(): Invalid user-specified number of non-zeros per row.");
    clearGlobalConstants();
    checkInternalState();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::CrsGraph(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, const Teuchos::ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, ProfileType pftype)
  : rowMap_(rowMap)
  , colMap_(colMap)
  , lclGraph_(rowMap->getNode())
  , nodeNumEntries_(0)
  , nodeNumDiags_(0)
  , nodeMaxNumRowEntries_(0)
  , nodeNumAllocated_(0)
  , pftype_(pftype)
  , numAllocPerRow_(NumEntriesPerRowToAlloc) 
  , numAllocForAllRows_(0)
  , indicesAreAllocated_(false)
  , indicesAreLocal_(false)
  , indicesAreGlobal_(false)
  , fillComplete_(false)
  , storageOptimized_(false)
  , lowerTriangular_(false)
  , upperTriangular_(false)
  , indicesAreSorted_(true) 
  , noRedundancies_(false) {
    staticAssertions();
    TEST_FOR_EXCEPTION(NumEntriesPerRowToAlloc.size() != getNodeNumRows(), std::runtime_error,
        Teuchos::typeName(*this) << "::CrsGraph(rowMap,colMap,NumEntriesPerRowToAlloc): NumEntriesPerRowToAlloc must have as many entries as specified by rowMap for this node.");
    clearGlobalConstants();
    checkInternalState();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::allocateIndices(AllocateLocalGlobal lorg) {
    using Teuchos::null;
    // this is a protected function, only callable by us. if it was called incorrectly, it is our fault.
    TEST_FOR_EXCEPTION(isLocallyIndexed() && lorg==AllocateGlobal, std::logic_error,
        Teuchos::typeName(*this) << "::allocateIndices(): Internal logic error. Please contact Tpetra team.");
    TEST_FOR_EXCEPTION(isGloballyIndexed() && lorg==AllocateLocal, std::logic_error,
        Teuchos::typeName(*this) << "::allocateIndices(): Internal logic error. Please contact Tpetra team.");
    TEST_FOR_EXCEPTION( indicesAreAllocated() == true, std::logic_error, 
        Teuchos::typeName(*this) << "::allocateIndices(): Internal logic error. Please contact Tpetra team.");
    indicesAreLocal_  = (lorg == AllocateLocal);
    indicesAreGlobal_ = (lorg == AllocateGlobal);
    Teuchos::RCP<Node> node = lclGraph_.getNode();
    const size_t numRows = getNodeNumRows();
    // if we have no row, then we have nothing to do
    if (numRows > 0) {
      if (getProfileType() == StaticProfile) {
        // determine how many entries to allocate and setup offsets into 1D arrays
        pbuf_rowOffsets_ = node->template allocBuffer<size_t>(numRows+1);
        Teuchos::ArrayRCP<size_t> view_offsets = node->template viewBufferNonConst(Kokkos::WriteOnly,numRows+1,pbuf_rowOffsets_);
        if (numAllocPerRow_ != Teuchos::null) {
          // allocate offsets, get host view
          nodeNumAllocated_ = 0;
          for (size_t i=0; i < numRows; ++i) {
            view_offsets[i] = nodeNumAllocated_;
            nodeNumAllocated_ += numAllocPerRow_[i];
          }
          view_offsets[numRows] = nodeNumAllocated_;
        }
        else {
          nodeNumAllocated_ = numAllocForAllRows_ * numRows;
          view_offsets[0] = 0;
          for (size_t i=1; i <= numRows; ++i) {
            view_offsets[i] = view_offsets[i-1] + numAllocForAllRows_;
          }
        }
        // release this view
        view_offsets = Teuchos::null;
        // in hind-sight, we know nodeNumAllocated_ and whether pbuf_rowOffsets_ is needed at all
        if (nodeNumAllocated_ == 0) {
          pbuf_rowOffsets_ = null;
        }
        // allocate the indices
        if (nodeNumAllocated_ > 0) {
          if (lorg == AllocateLocal) {
            pbuf_lclInds1D_ = node->template allocBuffer<LocalOrdinal>(nodeNumAllocated_);
          }
          else {
            pbuf_gblInds1D_ = node->template allocBuffer<GlobalOrdinal>(nodeNumAllocated_);
          }
        }
      }
      else { // getProfileType() == DynamicProfile
        Teuchos::ArrayRCP<const size_t> numalloc = numAllocPerRow_;
        size_t howmany = numAllocForAllRows_;
        if (lorg == AllocateLocal) {
          pbuf_lclInds2D_ = Teuchos::arcp< Teuchos::ArrayRCP<LocalOrdinal> >(numRows);
          for (size_t i=0; i < numRows; ++i) {
            if (numalloc != Teuchos::null) howmany = *numalloc++;
            nodeNumAllocated_ += howmany;
            if (howmany > 0) pbuf_lclInds2D_[i] = node->template allocBuffer<LocalOrdinal>(howmany);
          }
        }
        else { // allocate global indices
          pbuf_gblInds2D_ = Teuchos::arcp< Teuchos::ArrayRCP<GlobalOrdinal> >(numRows);
          for (size_t i=0; i < numRows; ++i) {
            if (numalloc != Teuchos::null) howmany = *numalloc++;
            nodeNumAllocated_ += howmany;
            if (howmany > 0) pbuf_gblInds2D_[i] = node->template allocBuffer<GlobalOrdinal>(howmany);
          }
        }
      }
      if (nodeNumAllocated_ > 0) {
        numEntriesPerRow_ = Teuchos::arcp<size_t>(numRows);
        std::fill(numEntriesPerRow_.begin(), numEntriesPerRow_.end(), 0);
      }
      else {
        // in hind-sight, we know nodeNumAllocated_ and whether pbuf_lclInds2D_ and pbuf_lclInds2D_ are needed at all
        pbuf_lclInds2D_ = Teuchos::null;
        pbuf_gblInds2D_ = Teuchos::null;
      }
    }
    // done with these
    numAllocForAllRows_ = 0;
    numAllocPerRow_     = Teuchos::null;
    indicesAreAllocated_ = true;    
    checkInternalState();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::~CrsGraph()
  {}


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  global_size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getGlobalNumRows() const { 
    return rowMap_->getGlobalNumElements(); 
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  global_size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getGlobalNumCols() const {
    TEST_FOR_EXCEPTION(hasColMap() == false, std::runtime_error,
      Teuchos::typeName(*this) << ": cannot call getGlobalNumCols() without column map.");
    return colMap_->getGlobalNumElements();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNodeNumRows() const { 
    return rowMap_->getNodeNumElements(); 
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNodeNumCols() const {
    TEST_FOR_EXCEPTION(hasColMap() != true, std::runtime_error,
      Teuchos::typeName(*this) << ": cannot call getNodeNumCols() without a column map.");
    return colMap_->getNodeNumElements();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNodeNumDiags() const { 
    TEST_FOR_EXCEPTION(!isFillComplete(), std::runtime_error,
      Teuchos::typeName(*this) << ": cannot call getNodeNumDiags() until fillComplete() has been called.");
    return nodeNumDiags_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  global_size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getGlobalNumDiags() const {
    return globalNumDiags_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<const LocalOrdinal> CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getLocalRowView(LocalOrdinal LocalRow) const {
    TEST_FOR_EXCEPTION(isGloballyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowView(): local indices do not exist.");
    TEST_FOR_EXCEPTION(rowMap_->isNodeLocalElement(LocalRow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowView(LocalRow): LocalRow (== " << LocalRow << ") is not valid on this node.");
    const size_t rnnz = RNNZ(LocalRow);
    Teuchos::ArrayRCP<const LocalOrdinal> ret;
    Teuchos::RCP<Node> node = lclGraph_.getNode();
    if (indicesAreAllocated() == false || rnnz == 0) {
      ret = Teuchos::null;
    }
    else {
      // RNNZ > 0, so there must be something allocated
      if (getProfileType() == StaticProfile) {
        Teuchos::ArrayRCP<const size_t> view_offs = node->template viewBuffer<size_t>(1,pbuf_rowOffsets_+LocalRow);
        ret = node->template viewBuffer<LocalOrdinal>(rnnz, pbuf_lclInds1D_ + view_offs[0]);
      }
      else {  // dynamic profile
        ret = node->template viewBuffer<LocalOrdinal>(rnnz, pbuf_lclInds2D_[LocalRow]);
      }
    }
    return ret;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<const GlobalOrdinal> CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getGlobalRowView(GlobalOrdinal GlobalRow) const {
    TEST_FOR_EXCEPTION(isLocallyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalRowView(): global indices do not exist.");
    const LocalOrdinal lrow = rowMap_->getLocalElement(GlobalRow);
    TEST_FOR_EXCEPTION(lrow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalRowView(GlobalRow): GlobalRow (== " << GlobalRow << ") does not belong to this node.");
    const size_t rnnz = RNNZ(lrow);
    Teuchos::ArrayRCP<const GlobalOrdinal> ret;
    Teuchos::RCP<Node> node = lclGraph_.getNode();
    if (indicesAreAllocated() == false || rnnz == 0) {
      ret = Teuchos::null;
    }
    else {
      // RNNZ > 0, so there must be something allocated
      if (getProfileType() == StaticProfile) {
        Teuchos::ArrayRCP<const size_t> view_offs = node->template viewBuffer<size_t>(1,pbuf_rowOffsets_+lrow);
        ret = node->template viewBuffer<GlobalOrdinal>(rnnz, pbuf_gblInds1D_ + view_offs[0]);
      }
      else {  // dynamic profile
        ret = node->template viewBuffer<GlobalOrdinal>(rnnz, pbuf_gblInds2D_[lrow]);
      }
    }
    return ret;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<LocalOrdinal> CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getFullLocalRowView(LocalOrdinal LocalRow) { 
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(isLocallyIndexed() == false, std::logic_error, 
        Teuchos::typeName(*this) << "::getFullLocalRowView(): Internal logic error. Please contact Tpetra team.");
    TEST_FOR_EXCEPTION(rowMap_->isNodeLocalElement(LocalRow) == false, std::logic_error, 
        Teuchos::typeName(*this) << "::getFullLocalRowView(): Internal logic error. Please contact Tpetra team.");
#endif
    Teuchos::ArrayRCP<LocalOrdinal> ret = Teuchos::null;
    if (indicesAreAllocated() && nodeNumAllocated_ > 0) {
      Teuchos::RCP<Node> node = lclGraph_.getNode();
      const size_t rnnz = RNNZ(LocalRow);
      Kokkos::ReadWriteOption rw = (rnnz == 0 ? Kokkos::WriteOnly : Kokkos::ReadWrite);
      if (getProfileType() == StaticProfile) {
        Teuchos::ArrayRCP<const size_t> offs = node->template viewBuffer<size_t>(2,pbuf_rowOffsets_+LocalRow);
        const size_t rna = offs[1] - offs[0];
        if (rna > 0) {
          ret = node->template viewBufferNonConst<LocalOrdinal>(rw, rna, pbuf_lclInds1D_ + offs[0]);
        }
      }
      else {  // dynamic profile
        const size_t rna = pbuf_lclInds2D_[LocalRow].size();
        if (rna > 0) {
          ret = node->template viewBufferNonConst<LocalOrdinal>(rw, rna, pbuf_lclInds2D_[LocalRow]);
        }
      }
    }
    return ret;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<GlobalOrdinal> CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getFullGlobalRowView(GlobalOrdinal GlobalRow) {
    const LocalOrdinal lrow = rowMap_->getLocalElement(GlobalRow);
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(isGloballyIndexed() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getFullLocalRowView(): Internal logic error. Please contact Tpetra team.");
    TEST_FOR_EXCEPTION(lrow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::logic_error, 
        Teuchos::typeName(*this) << "::getFullLocalRowView(): Internal logic error. Please contact Tpetra team.");
#endif
    Teuchos::ArrayRCP<GlobalOrdinal> ret = Teuchos::null;
    if (indicesAreAllocated() && nodeNumAllocated_ > 0) {
      Teuchos::RCP<Node> node = lclGraph_.getNode();
      const size_t rnnz = RNNZ(lrow);
      Kokkos::ReadWriteOption rw = (rnnz == 0 ? Kokkos::WriteOnly : Kokkos::ReadWrite);
      if (getProfileType() == StaticProfile) {
        Teuchos::ArrayRCP<const size_t> offs = node->template viewBuffer<size_t>(2,pbuf_rowOffsets_+lrow);
        const size_t rna = offs[1] - offs[0];
        if (rna > 0) {
          ret = node->template viewBufferNonConst<GlobalOrdinal>(rw, rna, pbuf_gblInds1D_ + offs[0]);
        }
      }
      else {  // dynamic profile
        const size_t rna = pbuf_gblInds2D_[lrow].size();
        if (rna > 0) {
          ret = node->template viewBufferNonConst<GlobalOrdinal>(rw, rna, pbuf_gblInds2D_[lrow]);
        }
      }
    }
    return ret;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getLocalRowCopy(LocalOrdinal LocalRow, const Teuchos::ArrayView<LocalOrdinal> &indices, size_t &NumIndices) const {
    // can only do this if 
    // * we have local indices: isLocallyIndexed()
    // * we are capable of producing them: isGloballyIndexed() && hasColMap()
    // short circuit if we aren't allocated
    TEST_FOR_EXCEPTION(isGloballyIndexed() == true && hasColMap() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowCopy(): local indices cannot be produced.");
    TEST_FOR_EXCEPTION(rowMap_->isNodeLocalElement(LocalRow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowCopy(LocalRow,...): LocalRow (== " << LocalRow << ") is not valid on this node.");
    NumIndices = RNNZ(LocalRow);
    TEST_FOR_EXCEPTION((size_t)indices.size() < NumIndices, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowCopy(): specified storage (size==" << indices.size() 
        << ") is not large enough to hold all entries for this row (NumIndices == " << NumIndices << ").");
    // use one of the view routines to get the proper view, then copy it over
    if (isLocallyIndexed()) {
      Teuchos::ArrayRCP<const LocalOrdinal> lview = getLocalRowView(LocalRow);
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION((size_t)lview.size() != NumIndices, std::logic_error,
          Teuchos::typeName(*this) << "::getLocalRowCopy(): Internal logic error. Please contact Tpetra team.");
#endif
      std::copy(lview.begin(),lview.end(),indices.begin());
      lview = Teuchos::null;
    }
    else if (isGloballyIndexed()) {
      const GlobalOrdinal grow = rowMap_->getGlobalElement(LocalRow);
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION(grow == Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(), std::logic_error, 
          Teuchos::typeName(*this) << "::getLocalRowCopy(): Internal logic error. Please contact Tpetra team.");
#endif
      Teuchos::ArrayRCP<const GlobalOrdinal> gview = getGlobalRowView(grow);
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION((size_t)gview.size() != NumIndices, std::logic_error,
          Teuchos::typeName(*this) << "::getLocalRowCopy(): Internal logic error. Please contact Tpetra team.");
#endif
      for (size_t j=0; j < NumIndices; ++j) {
        indices[j] = colMap_->getLocalElement(gview[j]);
      }
      gview = Teuchos::null;
    }
    else {
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION( indicesAreAllocated() == true, std::logic_error, 
          Teuchos::typeName(*this) << "::getLocalRowCopy(): Internal logic error. Please contact Tpetra team.");
#endif
      NumIndices = 0;
    }
    return;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getGlobalRowCopy(GlobalOrdinal GlobalRow, const Teuchos::ArrayView<GlobalOrdinal> &indices, size_t &NumIndices) const {
    // we either currently store global indices, or we have a column map with which to transcribe our local indices for the user
    const LocalOrdinal lrow = rowMap_->getLocalElement(GlobalRow);
    TEST_FOR_EXCEPTION(lrow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalRowCopy(GlobalRow,...): GlobalRow (== " << GlobalRow << ") does not belong to this node.");
    NumIndices = RNNZ(lrow);
    TEST_FOR_EXCEPTION((size_t)indices.size() < NumIndices, std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalRowCopy(): specified storage (size==" << indices.size() 
        << ") is not large enough to hold all entries for this row (rnnz == " << NumIndices << ").");
    // use one of the view routines to get the proper view, then copy it over
    if (isLocallyIndexed()) {
      Teuchos::ArrayRCP<const LocalOrdinal> lview = getLocalRowView(lrow);
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION((size_t)lview.size() != NumIndices, std::logic_error,
          Teuchos::typeName(*this) << "::getGlobalRowCopy(): Internal logic error. Please contact Tpetra team.");
#endif
      // copy and convert
      typename Teuchos::ArrayView<    GlobalOrdinal>::iterator dstptr = indices.begin();
      typename Teuchos::ArrayRCP<const LocalOrdinal>::iterator srcptr = lview.begin();
      while ( srcptr!=lview.end() ) {
        (*dstptr++) = colMap_->getGlobalElement(*srcptr++);
      }
      lview = Teuchos::null;
    }
    else if (isGloballyIndexed()) {
      Teuchos::ArrayRCP<const GlobalOrdinal> gview = getGlobalRowView(GlobalRow);
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION((size_t)gview.size() != NumIndices, std::logic_error,
          Teuchos::typeName(*this) << "::getGlobalRowCopy(): Internal logic error. Please contact Tpetra team.");
#endif
      std::copy(gview.begin(), gview.end(), indices.begin());
      gview = Teuchos::null;
    }
    return;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<Node>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNode() const {
    return lclGraph_.getNode();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getRowMap() const {
    return rowMap_; 
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getColMap() const { 
    return colMap_; 
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getDomainMap() const {
    return domainMap_; 
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getRangeMap() const {
    return rangeMap_; 
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> >
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getImporter() const { 
    return importer_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> >
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getExporter() const {
    return exporter_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::hasColMap() const { 
    return (colMap_ != Teuchos::null); 
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::isStorageOptimized() const { 
    return storageOptimized_; 
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ProfileType CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getProfileType() const { 
    return pftype_; 
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  global_size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getGlobalNumEntries() const { 
    return globalNumEntries_; 
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNodeNumEntries() const { 
    return nodeNumEntries_; 
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  global_size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getGlobalMaxNumRowEntries() const {
    return globalMaxNumRowEntries_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNodeMaxNumRowEntries() const { 
    return nodeMaxNumRowEntries_; 
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::isFillComplete() const { 
    return fillComplete_; 
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::isUpperTriangular() const { 
    return upperTriangular_; 
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::isLowerTriangular() const { 
    return lowerTriangular_; 
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::isLocallyIndexed() const { 
    return indicesAreLocal_; 
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::isGloballyIndexed() const { 
    return indicesAreGlobal_; 
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const {
    using Teuchos::OrdinalTraits;
    const LocalOrdinal rlid = rowMap_->getLocalElement(globalRow);
    size_t ret;
    if (rlid == OrdinalTraits<LocalOrdinal>::invalid()) {
      ret = OrdinalTraits<size_t>::invalid();
    }
    else { 
      ret = RNNZ(rlid);
    }
    return ret;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNumEntriesInLocalRow(LocalOrdinal localRow) const {
    using Teuchos::OrdinalTraits;
    size_t ret;
    if (!rowMap_->isNodeLocalElement(localRow)) {
      ret = OrdinalTraits<size_t>::invalid();
    }
    else {
      ret = RNNZ(localRow);
    }
    return ret;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNodeAllocationSize() const {
    return nodeNumAllocated_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNumAllocatedEntriesInGlobalRow(GlobalOrdinal globalRow) const {
    using Teuchos::OrdinalTraits;
    const LocalOrdinal rlid = rowMap_->getLocalElement(globalRow);
    size_t ret;
    if (rlid == OrdinalTraits<LocalOrdinal>::invalid()) {
      ret = OrdinalTraits<size_t>::invalid();
    }
    else {
      ret = RNumAlloc(rlid);
    }
    return ret;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNumAllocatedEntriesInLocalRow(LocalOrdinal localRow) const {
    using Teuchos::OrdinalTraits;
    size_t ret;
    if (!rowMap_->isNodeLocalElement(localRow)) {
      ret = OrdinalTraits<size_t>::invalid();
    }
    else {
      ret = RNumAlloc(localRow);
    }
    return ret; 
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::insertLocalIndices(LocalOrdinal lrow, const Teuchos::ArrayView<const LocalOrdinal> &indices) {
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(isStorageOptimized() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertLocalIndices(): cannot insert new indices after optimizeStorage() has been called.");
    TEST_FOR_EXCEPTION(isGloballyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertLocalIndices(): graph indices are global; use insertGlobalIndices().");
    TEST_FOR_EXCEPTION(hasColMap() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::insertLocalIndices(): cannot insert local indices without a column map.");
    TEST_FOR_EXCEPTION(rowMap_->isNodeLocalElement(lrow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::insertLocalIndices(): row does not belong to this node.");
    if (indicesAreAllocated() == false) {
      allocateIndices(AllocateLocal);
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION(indicesAreAllocated() == false, std::logic_error, 
          Teuchos::typeName(*this) << "::insertLocalIndices(): Internal logic error. Please contact Tpetra team.");
#endif
    }
    //
    indicesAreSorted_ = false;
    noRedundancies_ = false;
    clearGlobalConstants();
    //
    Teuchos::RCP<Node> node = lclGraph_.getNode();
    // add to allocated space
    const size_t rowNE = RNNZ(lrow),
                 toAdd = indices.size(),
                 rowNA = RNumAlloc(lrow);
    if (rowNE+toAdd > rowNA) {
      TEST_FOR_EXCEPTION(getProfileType() == StaticProfile, std::runtime_error,
          Teuchos::typeName(*this) << "::insertLocalIndices(): new indices exceed statically allocated graph structure.");
      TPETRA_EFFICIENCY_WARNING(true, std::runtime_error,
          "::insertLocalIndices(): Pre-allocated space has been exceeded, requiring new allocation. To improve efficiency, suggest larger allocation.");
      // update allocation
      size_t newAlloc = rowNE + toAdd;
      updateLocalAllocation(lrow,newAlloc);
    }
    // get pointers to row allocation
    ArrayRCP<LocalOrdinal> rowview, rowptr;
    rowview = getFullLocalRowView(lrow);
    rowptr = rowview + rowNE;
    // check the local indices against the column map; only add ones that are defined
    typename ArrayView<const LocalOrdinal>::iterator srcind = indices.begin();
    while (srcind != indices.end()) {
      if (colMap_->isNodeLocalElement(*srcind)) {
        (*rowptr++) = (*srcind);
      }
      ++srcind;
    }
    numEntriesPerRow_[lrow] = rowptr - rowview;
    checkInternalState();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::insertGlobalIndices(GlobalOrdinal grow, const Teuchos::ArrayView<const GlobalOrdinal> &indices) {
    using Teuchos::OrdinalTraits;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    TEST_FOR_EXCEPTION(isLocallyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertGlobalIndices(): graph indices are local; use insertLocalIndices().");
    if (indicesAreAllocated() == false) {
      allocateIndices(AllocateGlobal);
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION(indicesAreAllocated() == false, std::logic_error, 
          Teuchos::typeName(*this) << "::insertGlobalIndices(): Internal logic error. Please contact Tpetra team.");
#endif
    }
    // 
    indicesAreSorted_ = false;
    noRedundancies_ = false;
    clearGlobalConstants();
    //
    Teuchos::RCP<Node> node = lclGraph_.getNode();
    const LocalOrdinal lrow = rowMap_->getLocalElement(grow);
    if (lrow != OrdinalTraits<LocalOrdinal>::invalid()) {
      // add to allocated space
      size_t rowNE = RNNZ(lrow),
             toAdd = indices.size(), 
             rowNA = RNumAlloc(lrow);
      if (rowNE+toAdd > rowNA) {
        TEST_FOR_EXCEPTION(getProfileType() == StaticProfile, std::runtime_error,
            Teuchos::typeName(*this) << "::insertGlobalIndices(): new indices exceed statically allocated graph structure.");
        TPETRA_EFFICIENCY_WARNING(true,std::runtime_error,
            "::insertGlobalIndices(): Pre-allocated space has been exceeded, requiring new allocation. To improve efficiency, suggest larger allocation.");
        // update allocation
        size_t newAlloc = rowNE + toAdd;
        updateGlobalAllocation(lrow,newAlloc);
      }
      ArrayRCP<GlobalOrdinal> rowview, rowptr;
      rowview = getFullGlobalRowView(grow);
      rowptr = rowview + rowNE;
      if (hasColMap()) {
        // check the global indices against the column map; only add ones that are defined
        typename ArrayView<const GlobalOrdinal>::iterator srcind = indices.begin();
        while (srcind != indices.end()) {
          if (colMap_->isNodeGlobalElement(*srcind)) {
            (*rowptr++) = (*srcind);
          }
          ++srcind;
        }
        numEntriesPerRow_[lrow] = rowptr - rowview;
      }
      else {
        std::copy( indices.begin(), indices.end(), rowptr );
        numEntriesPerRow_[lrow] += toAdd;
      }
      rowptr = Teuchos::null;
      rowview = Teuchos::null;
    }
    else {
      // a nonlocal
      for (typename Teuchos::ArrayView<const GlobalOrdinal>::iterator i=indices.begin(); i != indices.end(); ++i) {
        nonlocals_[grow].push_back(*i);
      }
    }
    checkInternalState();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const Teuchos::RCP<const Teuchos::Comm<int> > &
  CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getComm() const { 
    return rowMap_->getComm();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  GlobalOrdinal CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getIndexBase() const { 
    return rowMap_->getIndexBase(); 
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::isSorted() const { 
    return indicesAreSorted_; 
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::notRedundant() const { 
    return noRedundancies_; 
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::setSorted(bool sorted) { 
    indicesAreSorted_ = sorted; 
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::indicesAreAllocated() const { 
    return indicesAreAllocated_; 
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::fillComplete(OptimizeOption os) {
    fillComplete(rowMap_,rowMap_,os);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::globalAssemble() {
    using Teuchos::arcp;
    using Teuchos::Array;
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
    {
      size_t MyNonlocals = nonlocals_.size(), MaxGlobalNonlocals;
      Teuchos::reduceAll<int,size_t>(*getComm(),Teuchos::REDUCE_MAX,MyNonlocals,&MaxGlobalNonlocals);
      if (MaxGlobalNonlocals == 0) return;  // no entries to share
    }

    // compute a list of NLRs from nonlocals_ and use it to compute:
    //      IdsAndRows: a vector of (id,row) pairs
    //          NLR2Id: a map from NLR to the Id that owns it
    // globalNeighbors: a global graph of connectivity between images: globalNeighbors(i,j) indicates that j sends to i
    //         sendIDs: a list of all images I send to
    //         recvIDs: a list of all images I receive from (constructed later)
    Array<pair<int,GlobalOrdinal> > IdsAndRows;
    std::map<GlobalOrdinal,int> NLR2Id;
    Teuchos::SerialDenseMatrix<int,char> globalNeighbors;
    Array<int> sendIDs, recvIDs;
    {
      // nonlocals_ contains the entries we are holding for all non-local rows
      // we want a list of the rows for which we have data
      Array<GlobalOrdinal> NLRs;
      std::set<GlobalOrdinal> setOfRows;
      for (NLITER iter = nonlocals_.begin(); iter != nonlocals_.end(); ++iter) {
        setOfRows.insert(iter->first);
      }
      // copy the elements in the set into an Array
      NLRs.resize(setOfRows.size());
      std::copy(setOfRows.begin(), setOfRows.end(), NLRs.begin());

      // get a list of ImageIDs for the non-local rows (NLRs)
      Array<int> NLRIds(NLRs.size());
      {
        LookupStatus stat = rowMap_->getRemoteIndexList(NLRs(),NLRIds());
        char lclerror = ( stat == IDNotPresent ? 1 : 0 );
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
           nlr != NLRs.end(); ++nlr, ++id) {
        NLR2Id[*nlr] = *id;
        localNeighbors[*id] = 1;
        IdsAndRows.push_back(make_pair<int,GlobalOrdinal>(*id,*nlr));
      }
      for (int j=0; j<numImages; ++j) {
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
    for (int j=0; j<numImages; ++j) {
      if (globalNeighbors(myImageID,j)) {
        recvIDs.push_back(j);
      }
    }
    const size_t numRecvs = recvIDs.size();

    // we know how many we're sending to already
    // form a contiguous list of all data to be sent
    // track the number of entries for each ID
    Array<pair<GlobalOrdinal,GlobalOrdinal> > IJSendBuffer;
    Array<size_t> sendSizes(sendIDs.size(), 0);
    size_t numSends = 0;
    for (typename Array<pair<int,GlobalOrdinal> >::const_iterator IdAndRow = IdsAndRows.begin();
         IdAndRow != IdsAndRows.end(); ++IdAndRow) {
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
    for (size_t s=0; s < numSends ; ++s) {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      sendRequests.push_back( Teuchos::isend<int,size_t>(*getComm(),rcp<size_t>(&sendSizes[s],false),sendIDs[s]) );
    }
    // perform non-blocking receives: receive sizes from our senders
    Array<RCP<Teuchos::CommRequest> > recvRequests;
    Array<size_t> recvSizes(numRecvs);
    for (size_t r=0; r < numRecvs; ++r) {
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
    for (size_t r=0; r < numRecvs ; ++r) {
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
    for (typename Array<pair<GlobalOrdinal,GlobalOrdinal> >::const_iterator ij = IJRecvBuffer.begin(); ij != IJRecvBuffer.end(); ++ij) {
      insertGlobalIndices(ij->first, Teuchos::tuple<GlobalOrdinal>(ij->second));
    }
    checkInternalState();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::fillComplete(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, 
                                                               const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap, 
                                                               OptimizeOption os) {
    using Teuchos::Array;
    using Teuchos::ArrayRCP;

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
      global_size_t lcl[2], gbl[2];
      lcl[0] = nodeNumEntries_;
      lcl[1] = nodeNumDiags_;
      Teuchos::reduceAll<int,global_size_t>(*getComm(),Teuchos::REDUCE_SUM,2,lcl,gbl);
      globalNumEntries_ = gbl[0]; 
      globalNumDiags_   = gbl[1];
      Teuchos::reduceAll<int,global_size_t>(*getComm(),Teuchos::REDUCE_MAX,nodeMaxNumRowEntries_,&globalMaxNumRowEntries_);
      haveGlobalConstants_ = true;
    }

    // mark transformation as successfully completed
    fillComplete_ = true;

#ifdef HAVE_TPETRA_DEBUG 
    // check post-conditions
    TEST_FOR_EXCEPTION( !( pbuf_gblInds1D_ == Teuchos::null && pbuf_gblInds2D_ == Teuchos::null ), std::logic_error, 
        Teuchos::typeName(*this) << "::fillComplete(): Internal logic error. Please contact Tpetra team.");
#endif
    if (os == DoOptimizeStorage) optimizeStorage();

    // TEST_FOR_EXCEPT(true); // FINISH fill localGraph_
    checkInternalState();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::makeIndicesLocal(
                                      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, 
                                      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap) {
    using Teuchos::ArrayRCP;
    using Teuchos::NullIteratorTraits;
    // All nodes must be in the same index state.
    // Update index state by checking isLocallyIndexed/Global on all nodes
    computeIndexState(); 
    TEST_FOR_EXCEPTION(isLocallyIndexed() && isGloballyIndexed(), std::logic_error,
        Teuchos::typeName(*this) << "::makeIndicesLocal(): indices are marked as both global and local.");
    // If user has not prescribed column map, create one from indices
    makeColMap(domainMap, rangeMap);
    // Transform indices to local index space
    const size_t nlrs = getNodeNumRows();
    Teuchos::RCP<Node> node = lclGraph_.getNode();
    // 
    if (isGloballyIndexed() && indicesAreAllocated() && nlrs > 0) {
      // allocate data for local indices
      if (nodeNumAllocated_) {
        if (getProfileType() == StaticProfile) {
          ArrayRCP<const size_t> view_offsets = node->template viewBuffer<size_t>(pbuf_rowOffsets_.size(), pbuf_rowOffsets_);
          ArrayRCP<GlobalOrdinal> view_ginds = node->template viewBufferNonConst<GlobalOrdinal>(Kokkos::ReadWrite,pbuf_gblInds1D_.size(), pbuf_gblInds1D_);
          // do the conversion in situ. this must be done from front to back.
          ArrayRCP< LocalOrdinal> view_linds = Teuchos::arcp_reinterpret_cast<LocalOrdinal>(view_ginds);
          for (size_t r=0; r < getNodeNumRows(); ++r) {
            const size_t offset   = view_offsets[r],
                         numentry = numEntriesPerRow_[r];
            for (size_t j=0; j<numentry; ++j) {
              GlobalOrdinal gid = view_ginds[offset + j];
              LocalOrdinal  lid = colMap_->getLocalElement(gid);
              view_linds[offset + j] = lid;
#ifdef HAVE_TPETRA_DEBUG
              TEST_FOR_EXCEPTION(view_linds[offset + j] == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::logic_error,
                  Teuchos::typeName(*this) << ": Internal error in fillComplete(). Please contact Tpetra team.");
#endif
            }
          }
          view_linds = Teuchos::null;
          view_ginds = Teuchos::null;
          // reinterpret the the compute buffer as LocalOrdinal; keep the original size
          pbuf_lclInds1D_ = Teuchos::arcp_reinterpret_cast<LocalOrdinal>(pbuf_gblInds1D_).persistingView(0,nodeNumAllocated_);
          pbuf_gblInds1D_ = Teuchos::null;
        }
        else {  // getProfileType() == DynamicProfile
          pbuf_lclInds2D_ = Teuchos::arcp< Teuchos::ArrayRCP<LocalOrdinal> >(nlrs);
          for (size_t r=0; r < getNodeNumRows(); ++r) {
            if (pbuf_gblInds2D_[r] != Teuchos::null) {
              const size_t rna = pbuf_gblInds2D_[r].size();
              ArrayRCP<GlobalOrdinal> view_ginds = node->template viewBufferNonConst<GlobalOrdinal>(Kokkos::ReadWrite,rna,pbuf_gblInds2D_[r]);
              // do the conversion in situ. this must be done from front to back.
              ArrayRCP< LocalOrdinal> view_linds = Teuchos::arcp_reinterpret_cast<LocalOrdinal>(view_ginds);
              const size_t numentry = numEntriesPerRow_[r];
              for (size_t j=0; j < numentry; ++j) {
                GlobalOrdinal gid = view_ginds[j];
                LocalOrdinal  lid = colMap_->getLocalElement(gid);
                view_linds[j] = lid;
#ifdef HAVE_TPETRA_DEBUG
                TEST_FOR_EXCEPTION(view_linds[j] == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::logic_error,
                    Teuchos::typeName(*this) << ": Internal error in makeIndicesLocal(). Please contact Tpetra team.");
#endif
              }
              view_linds = Teuchos::null;
              view_ginds = Teuchos::null;
              // reinterpret the the compute buffer as LocalOrdinal; keep the original size
              pbuf_lclInds2D_[r] = Teuchos::arcp_reinterpret_cast<LocalOrdinal>(pbuf_gblInds2D_[r]).persistingView(0,rna);
              pbuf_gblInds2D_[r] = Teuchos::null;
            }
          }
          pbuf_gblInds2D_ = Teuchos::null;
        }
      }
    }
    indicesAreLocal_  = true;
    indicesAreGlobal_ = false;
    checkInternalState();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::computeIndexState() {
    char myIndices[2] = {0,0};
    if (indicesAreLocal_)  myIndices[0] = 1;
    if (indicesAreGlobal_) myIndices[1] = 1;
    char allIndices[2];
    Teuchos::reduceAll(*getComm(),Teuchos::REDUCE_MAX,2,myIndices,allIndices);
    indicesAreLocal_  = (allIndices[0]==1);  // If indices are local on one PE, should be local on all
    indicesAreGlobal_ = (allIndices[1]==1);  // If indices are global on one PE should be local on all
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::sortIndices() {
    TEST_FOR_EXCEPT(isGloballyIndexed()==true);   // this should be called only after makeIndicesLocal()
    if (isSorted()) return;
    // are there any indices to sort?
    if (nodeNumAllocated_ > 0) {
      const size_t nlrs = getNodeNumRows();
      for (size_t r=0; r < nlrs; ++r) {
        const size_t rnnz = RNNZ(r);
        Teuchos::ArrayRCP<LocalOrdinal> row_view = getFullLocalRowView(r);
        std::sort(row_view, row_view + rnnz);
      }
    }
    setSorted(true);
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
    // 
    if (hasColMap()) return;
    const size_t nlrs = getNodeNumRows();
    // 
    computeIndexState();
    TEST_FOR_EXCEPTION(isLocallyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::makeColMap(): indices must still be global when making the column map.");
    // ultimate goal: list of indices for column map
    Array<GlobalOrdinal> myColumns;
    // if isGloballyIndexed() == false and isLocallyIndexed() == false, then we have nothing to do
    if (isGloballyIndexed() == true) {
      // Construct two lists of columns for which we have a non-zero
      // Local GIDs: Column GIDs which are locally present on the Domain map
      // Remote GIDs: Column GIDs which are not locally present present on the Domain map
      //
      // instead of list of local GIDs, we will use a set<LocalOrdinal> LocalGID.
      //
      const LocalOrdinal LINV = Teuchos::OrdinalTraits<LocalOrdinal>::invalid();
      size_t numLocalColGIDs = 0, numRemoteColGIDs = 0;
      //
      // intitial: partitioning into local and remote
      Teuchos::Array<char>    GIDisLocal(domainMap->getNodeNumElements(),0);
      std::set<GlobalOrdinal> RemoteGIDSet;
      for (size_t r=0; r < nlrs; ++r) {
        typename ArrayRCP<GlobalOrdinal>::const_iterator cind;
        const size_t rnnz = RNNZ(r);
        if (rnnz > 0) {
          ArrayRCP<GlobalOrdinal> rowgids = getFullGlobalRowView( rowMap_->getGlobalElement(r) ).persistingView(0,rnnz);
          for (cind = rowgids.begin(); cind != rowgids.end(); ++cind) {
            GlobalOrdinal gid = (*cind);
            LocalOrdinal lid = domainMap->getLocalElement(gid);
            if (lid != LINV) {
              char alreadyFound = GIDisLocal[lid];
              if (!alreadyFound) {
                GIDisLocal[lid] = 1;
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
      if (Teuchos::size(*domainMap->getComm()) == 1) {
        TEST_FOR_EXCEPTION(numRemoteColGIDs != 0, std::runtime_error,
            Teuchos::typeName(*this) << "::makeColMap(): Some column IDs are not in the domain map." << std::endl 
            << "Either these column IDs are invalid or the domain map is invalid." << std::endl
            << "Remember, for a rectangular matrix, the domain map must be passed to fillComplete().");
        if (numLocalColGIDs == domainMap->getNodeNumElements()) {
          colMap_ = domainMap;
          checkInternalState();
          return;
        }
      }

      // Now, populate myColumns() with a list of all column GIDs. 
      // Put local GIDs at the front: they correspond to "same" and "permuted" entries between the column map and the domain map
      // Put remote GIDs at the back
      myColumns.resize(numLocalColGIDs + numRemoteColGIDs);
      // get pointers into myColumns for each part
      ArrayView<GlobalOrdinal> LocalColGIDs, RemoteColGIDs;
      LocalColGIDs  = myColumns(0,numLocalColGIDs);
      RemoteColGIDs = myColumns(numLocalColGIDs,numRemoteColGIDs);

      // Copy the Remote GIDs into myColumns
      std::copy(RemoteGIDSet.begin(), RemoteGIDSet.end(), RemoteColGIDs.begin());
      // We will make a list of correspodning node IDs
      Array<int> RemoteImageIDs(numRemoteColGIDs);
      // Lookup the Remote Nodes IDs in the Domain map
      {
        LookupStatus stat = domainMap->getRemoteIndexList(RemoteColGIDs, RemoteImageIDs());
        // This error check is crucial: this tells us that one of the remote indices was not present in the domain map. 
        // This means that the Import object cannot be constructed, because of incongruity between the column map and domain map.
        //   * The user has made a mistake in the column indices
        //   * The user has made a mistake w.r.t. the domain map 
        // Same error message as above for serial case.
        char missingID_lcl = (stat == IDNotPresent ? 1 : 0);
        char missingID_gbl;
        Teuchos::reduceAll<int,char>(*getComm(),Teuchos::REDUCE_MAX,missingID_lcl,&missingID_gbl);
        TEST_FOR_EXCEPTION(missingID_gbl == 1, std::runtime_error,
            Teuchos::typeName(*this) << "::makeColMap(): Some column IDs are not in the domain map." << std::endl 
            << "Either these column IDs are invalid or the domain map is invalid." << std::endl
            << "Common cause: for a rectangular matrix, the domain map must be passed to fillComplete().");
      }
      // Sort External column indices so that all columns coming from a given remote processor are contiguous
      // This obliviates the need for the Distributor associated with the Import from having to reorder data.
      sort2(RemoteImageIDs.begin(), RemoteImageIDs.end(), RemoteColGIDs.begin());

      // Copy the Local GIDs into myColumns. Two cases:
      // (1) If the number of Local column GIDs is the same as the number of Local domain GIDs, we
      //     can simply read the domain GIDs into the front part of ColIndices (see logic above from the serial short circuit case)
      // (2) We step through the GIDs of the DomainMap, checking to see if each domain GID is a column GID.
      //     we want to do this to maintain a consistent ordering of GIDs between the columns and the domain.
      ArrayView<const GlobalOrdinal> mge = domainMap->getNodeElementList();
      if (numLocalColGIDs == domainMap->getNodeNumElements()) {
        std::copy(mge.begin(), mge.end(), LocalColGIDs.begin());
      }
      else {
        size_t numlocalagain = 0;
        for (size_t i=0; i < domainMap->getNodeNumElements(); ++i) {
          if (GIDisLocal[i]) {
            LocalColGIDs[numlocalagain++] = mge[i];
          }
        }
        TEST_FOR_EXCEPTION(numlocalagain != numLocalColGIDs, std::logic_error,
            Teuchos::typeName(*this) << "::makeColMap(): Internal logic error. Please contact Tpetra team.");
      }
    }
    colMap_ = Teuchos::rcp(new Map<LocalOrdinal,GlobalOrdinal,Node>(GOT::invalid(), myColumns, domainMap->getIndexBase(), domainMap->getComm()) );
    checkInternalState();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::removeRedundantIndices() {
    TEST_FOR_EXCEPT( isGloballyIndexed() != false );      // this should be called only after makeIndicesLocal()
    TEST_FOR_EXCEPT( isSorted() != true );                // this should be called only after sortIndices()
    TEST_FOR_EXCEPT( isStorageOptimized() == true );      // assumptions about available data structures
    if ( notRedundant() ) return;
    const size_t nlrs = getNodeNumRows();
    Teuchos::ArrayView<const GlobalOrdinal> myGlobalEntries = rowMap_->getNodeElementList();
    // reset all local quantities
    upperTriangular_ = true;
    lowerTriangular_ = true;
    nodeMaxNumRowEntries_ = 0;
    nodeNumEntries_       = 0;
    nodeNumDiags_         = 0;
    // indices are already sorted in each row
    if (nodeNumAllocated_) {
      for (size_t r=0; r < nlrs; ++r) {
        GlobalOrdinal rgid = myGlobalEntries[r];
        // determine the local column index for this row, used for delimiting the diagonal
        const LocalOrdinal rlcid = colMap_->getLocalElement(rgid);   
        Teuchos::ArrayRCP<LocalOrdinal> rview = getFullLocalRowView(r);
        typename Teuchos::ArrayRCP<LocalOrdinal>::iterator beg, end;
        beg = rview.begin();
        end = beg + RNNZ(r);
        typename Teuchos::ArrayRCP<LocalOrdinal>::iterator newend, cur;
        newend = beg;
        if (beg != end) {
          if (rlcid == (*newend)) ++nodeNumDiags_;
          for (cur = beg + 1; cur != end; ++cur) {
            // it is unique; add it and test it for the diagonal
            // otherwise, it is a dup and should be ignored
            if (*cur != *newend) {
              ++newend;
              (*newend) = (*cur);
              // is this the diagonal?
              if (rlcid == (*newend)) ++nodeNumDiags_;
            }
          }
          // because of sorting, smallest column index is (*beg); it indicates upper triangularity
          if (Teuchos::as<size_t>(*beg) < r) upperTriangular_ = false;
          // because of sorting, largest column index is (*newend); it indicates lower triangularity
          if (r < Teuchos::as<size_t>(*newend)) lowerTriangular_ = false;
          // increment newend so that our range is [beg,newend) instead of [beg,newend]
          ++newend;
        }
        // compute num entries for this row, accumulate into nodeNumEntries_, update nodeMaxNumRowEntries_
        numEntriesPerRow_[r] = newend - beg;
        nodeMaxNumRowEntries_ = std::max( nodeMaxNumRowEntries_, numEntriesPerRow_[r] );
        nodeNumEntries_ += numEntriesPerRow_[r];
      }
    }
    noRedundancies_ = true;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::makeImportExport() {
    TEST_FOR_EXCEPT(hasColMap()==false); // must have column map
    // create import, export
    if (!domainMap_->isSameAs(*colMap_)) {
      importer_ = Teuchos::rcp( new Import<LocalOrdinal,GlobalOrdinal,Node>(domainMap_,colMap_) );
    }
    else {
      importer_ = Teuchos::null;
    }
    if (!rangeMap_->isSameAs(*rowMap_)) {
      exporter_ = Teuchos::rcp( new Export<LocalOrdinal,GlobalOrdinal,Node>(rowMap_,rangeMap_) );
    }
    else {
      exporter_ = Teuchos::null;
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::optimizeStorage() {
    using Teuchos::ArrayRCP;
    // optimizeStorage will perform two functions:
    // 1) create a single allocation of memory
    // 2) pack data in that allocation
    // if getProfileType() == StaticProfile, then 1) has already been done
    // 
    // post-condition:
    //   getProfileType() == StaticProfile
    //   numEntriesPerRow_ == Teuchos::null; size and allocation are the same, compute by subtracting offsets
    if (isStorageOptimized() == true) return;

    if (indicesAreAllocated() == false) {
      // won't ever allocate now
      indicesAreAllocated_ = true;
      nodeNumAllocated_ = 0;
      numAllocPerRow_ = Teuchos::null;
      numAllocForAllRows_ = 0;
    }

    TEST_FOR_EXCEPTION(isFillComplete() == false || isSorted() == false || notRedundant() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::optimizeStorage(): fillComplete() must be called before optimizeStorage().");

    Teuchos::RCP<Node> node = lclGraph_.getNode();

    // 1) allocate single memory block
    const size_t nlrs = getNodeNumRows();
    if (nlrs > 0 && nodeNumAllocated_ > 0) {
      if (getProfileType() == DynamicProfile) {
        pbuf_rowOffsets_ = node->template allocBuffer<size_t>(nlrs+1);
        ArrayRCP<size_t> view_offsets = node->template viewBufferNonConst(Kokkos::WriteOnly,nlrs+1,pbuf_rowOffsets_);
        if (nodeNumEntries_ > 0) {
          pbuf_lclInds1D_ = node->template allocBuffer<LocalOrdinal>(nodeNumEntries_);
          ArrayRCP<LocalOrdinal> curptr = pbuf_lclInds1D_;
          size_t sofar = 0;
          for (size_t r=0; r<nlrs; ++r) {
            const size_t rne = numEntriesPerRow_[r];
            if (rne > 0) {
              node->template copyBuffers<LocalOrdinal>(rne, pbuf_lclInds2D_[r], curptr);
            }
            view_offsets[r] = sofar;
            curptr += rne;
            sofar += rne;
          }
          view_offsets[nlrs] = sofar;
#ifdef HAVE_TPETRA_DEBUG
          TEST_FOR_EXCEPTION( nodeNumEntries_ != sofar, std::logic_error, 
              Teuchos::typeName(*this) << "::optimizeStorage(): Internal Tpetra logic error. Please contact Tpetra team.");
#endif
        }
        else {
          std::fill(view_offsets.begin(), view_offsets.end(), 0);
        }
        view_offsets = Teuchos::null;
        pbuf_lclInds2D_ = Teuchos::null;
      }
      else {
        // storage is already allocated; just need to pack
        if (nodeNumEntries_ > 0) {
          ArrayRCP<size_t> view_offsets = node->template viewBufferNonConst(Kokkos::WriteOnly,nlrs+1,pbuf_rowOffsets_);
          ArrayRCP<LocalOrdinal> curptr = pbuf_lclInds1D_,
                                 oldptr = pbuf_lclInds1D_;
          size_t sofar = 0;
          for (size_t r=0; r<nlrs; ++r) {
            const size_t rne = numEntriesPerRow_[r],
                          na = view_offsets[r+1] - view_offsets[r];
            if (curptr != oldptr) {
              node->template copyBuffers<LocalOrdinal>(rne, oldptr, curptr);
              view_offsets[r] = sofar;
            }
            sofar += rne;
            curptr += rne;
            oldptr += na;
          }
          view_offsets[nlrs] = sofar;
          view_offsets = Teuchos::null;
#ifdef HAVE_TPETRA_DEBUG
          TEST_FOR_EXCEPTION( nodeNumEntries_ != sofar, std::logic_error, 
              Teuchos::typeName(*this) << "::optimizeStorage(): Internal Tpetra logic error. Please contact Tpetra team.");
#endif
        }
      }
      nodeNumAllocated_ = nodeNumEntries_;
    }
    numEntriesPerRow_ = Teuchos::null;
    if (nodeNumAllocated_ == 0) {
      pbuf_rowOffsets_ = Teuchos::null;
    }
    storageOptimized_ = true;
    pftype_ = StaticProfile;
    checkInternalState();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::RNNZ(size_t row) const {
    size_t rnnz;
    // if storage is optimized, then numalloc == numnz, and indices are local
    if (isStorageOptimized()) {
      rnnz = RNumAlloc(row);
    }
    else if (nodeNumAllocated_ == 0) {
      rnnz = 0;
    }
    else {
      rnnz = numEntriesPerRow_[row];
    }
    return rnnz;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::RNumAlloc(size_t row) const {
    size_t numalloc;
    if (indicesAreAllocated() == false) {
      if (numAllocPerRow_ == Teuchos::null) {
        numalloc = numAllocForAllRows_;
      }
      else {
        numalloc = numAllocPerRow_[row];
      }
    }
    else if (nodeNumAllocated_ == 0) {
      return 0;
    }
    // if static graph or optimize storage, offsets tell us the allocation size
    else if (getProfileType() == StaticProfile) {
      Teuchos::RCP<Node> node = lclGraph_.getNode();
      Teuchos::ArrayRCP<const size_t> offs = node->template viewBuffer<size_t>(2,pbuf_rowOffsets_+row);
      numalloc = offs[1] - offs[0];
      offs = Teuchos::null;
    }
    // otherwise, the ArrayRCP knows
    else {
      if (isLocallyIndexed()) {
        numalloc = pbuf_lclInds2D_[row].size();
      }
      else {
        numalloc = pbuf_gblInds2D_[row].size();
      }
    }
    return numalloc;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::staticAssertions() {
    using Teuchos::OrdinalTraits;
    // Assumption: sizeof(GlobalOrdinal) >= sizeof(LocalOrdinal)
    //    This is so that we can store LocalOrdinals in the memory formerly occupied by GlobalOrdinals
    // Assumption: max(GlobalOrdinal) >= max(LocalOrdinal)  and  max(size_t) >= max(LocalOrdinal)
    //    This is so that we can represent any LocalOrdinal as a size_t, and any LocalOrdinal as a GlobalOrdinal
    Teuchos::CompileTimeAssert<sizeof(GlobalOrdinal) < sizeof(LocalOrdinal)> cta_size;
    (void)cta_size;
    // can't call max() with CompileTimeAssert, because it isn't a constant expression; will need to make this a runtime check
    const char * err = ": Object cannot be allocated with stated template arguments: size assumptions are not valid.";
    TEST_FOR_EXCEPTION( (size_t)OrdinalTraits<LocalOrdinal>::max() > OrdinalTraits<size_t>::max(),          std::runtime_error, Teuchos::typeName(*this) << err);
    TEST_FOR_EXCEPTION( OrdinalTraits<LocalOrdinal>::max() > OrdinalTraits<GlobalOrdinal>::max(),   std::runtime_error, Teuchos::typeName(*this) << err);
    TEST_FOR_EXCEPTION( (size_t)OrdinalTraits<GlobalOrdinal>::max() > OrdinalTraits<global_size_t>::max(),  std::runtime_error, Teuchos::typeName(*this) << err);
    TEST_FOR_EXCEPTION( OrdinalTraits<size_t>::max() > OrdinalTraits<global_size_t>::max(),         std::runtime_error, Teuchos::typeName(*this) << err);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::updateLocalAllocation(size_t lrow, size_t allocSize) {
    using Teuchos::ArrayRCP;
    Teuchos::RCP<Node> node = lclGraph_.getNode();
    const size_t curNA = RNumAlloc(lrow),
                  rnnz = RNNZ(lrow);
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPT( rowMap_->isNodeLocalElement(lrow) == false );
    TEST_FOR_EXCEPT( allocSize < curNA );
    TEST_FOR_EXCEPT( isGloballyIndexed() );
    TEST_FOR_EXCEPT( allocSize == 0 );
#endif
    if (pbuf_lclInds2D_ == Teuchos::null) {
      pbuf_lclInds2D_ = Teuchos::arcp< ArrayRCP<LocalOrdinal> >(getNodeNumRows());
    }
    ArrayRCP<LocalOrdinal> old_row, new_row;
    old_row = pbuf_lclInds2D_[lrow];
    new_row = node->template allocBuffer<LocalOrdinal>(allocSize);
    if (rnnz) {
      node->template copyBuffers<LocalOrdinal>(rnnz,old_row,new_row);
    }
    old_row = Teuchos::null;
    pbuf_lclInds2D_[lrow] = new_row;
    nodeNumAllocated_ += (allocSize - curNA);
    if (numEntriesPerRow_ == Teuchos::null) {
      numEntriesPerRow_ = Teuchos::arcp<size_t>( getNodeNumRows() );
      std::fill(numEntriesPerRow_.begin(), numEntriesPerRow_.end(), 0);
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::updateGlobalAllocation(size_t lrow, size_t allocSize) {
    using Teuchos::ArrayRCP;
    Teuchos::RCP<Node> node = lclGraph_.getNode();
    const size_t curNA = RNumAlloc(lrow),
                  rnnz = RNNZ(lrow);
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPT( rowMap_->isNodeLocalElement(lrow) == false );
    TEST_FOR_EXCEPT( allocSize < curNA );
    TEST_FOR_EXCEPT( isLocallyIndexed() );
    TEST_FOR_EXCEPT( allocSize == 0 );
#endif
    if (pbuf_gblInds2D_ == Teuchos::null) {
      pbuf_gblInds2D_ = Teuchos::arcp< ArrayRCP<GlobalOrdinal> >(getNodeNumRows());
    }
    ArrayRCP<GlobalOrdinal> old_row, new_row;
    old_row = pbuf_gblInds2D_[lrow];
    new_row = node->template allocBuffer<GlobalOrdinal>(allocSize);
    if (rnnz) {
      node->template copyBuffers<GlobalOrdinal>(rnnz,old_row,new_row);
    }
    old_row = Teuchos::null;
    pbuf_gblInds2D_[lrow] = new_row;
    nodeNumAllocated_ += (allocSize - curNA);
    if (numEntriesPerRow_ == Teuchos::null) {
      numEntriesPerRow_ = Teuchos::arcp<size_t>( getNodeNumRows() );
      std::fill(numEntriesPerRow_.begin(), numEntriesPerRow_.end(), 0);
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::findMyIndex(size_t row, LocalOrdinal ind) const {
    typedef typename Teuchos::ArrayRCP<const LocalOrdinal>::iterator IT;
    bool found = true;
    const size_t nE = RNNZ(row);
    Teuchos::ArrayRCP<const LocalOrdinal> rowview = getLocalRowView(row);
    IT rptr, locptr;
    rptr = rowview.begin();
    if (isSorted()) {
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
    size_t ret;
    if (!found) {
      ret = Teuchos::OrdinalTraits<size_t>::invalid();
    }
    else {
      ret = (locptr - rptr);
    }
    locptr = Teuchos::NullIteratorTraits<IT>::getNull();
    rptr = Teuchos::NullIteratorTraits<IT>::getNull();
    rowview = Teuchos::null;
    return ret;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::findGlobalIndex(size_t row, GlobalOrdinal ind) const {
    typedef typename Teuchos::ArrayRCP<const GlobalOrdinal>::iterator IT;
    bool found = true;
    const size_t nE = RNNZ(row);
    Teuchos::ArrayRCP<const GlobalOrdinal> rowview = getGlobalRowView(row);
    IT rptr, locptr;
    rptr = rowview.begin();
    if (isSorted()) {
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
    size_t ret;
    if (!found) {
      ret = Teuchos::OrdinalTraits<size_t>::invalid();
    }
    else {
      ret = (locptr - rptr);
    }
    locptr = Teuchos::NullIteratorTraits<IT>::getNull();
    rptr = Teuchos::NullIteratorTraits<IT>::getNull();
    rowview = Teuchos::null;
    return ret;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::clearGlobalConstants() {
    globalNumEntries_ = Teuchos::OrdinalTraits<global_size_t>::invalid();
    globalNumDiags_ = Teuchos::OrdinalTraits<global_size_t>::invalid();
    globalMaxNumRowEntries_ = Teuchos::OrdinalTraits<global_size_t>::invalid();
    haveGlobalConstants_ = false;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::checkInternalState() const {
#ifdef HAVE_TPETRA_DEBUG
    Teuchos::RCP<Node> node = lclGraph_.getNode();
    const global_size_t gsti = Teuchos::OrdinalTraits<global_size_t>::invalid();
    using Teuchos::null;
    std::string err = Teuchos::typeName(*this) + "::optimizeStorage(): Internal logic error. Please contact Tpetra team.";
    // check the internal state of this data structure
    // this is called by numerous state-changing methods, in a debug build, to ensure that the object 
    // always remains in a valid state
    TEST_FOR_EXCEPTION( rowMap_ == null, std::logic_error, err );
    TEST_FOR_EXCEPTION( hasColMap() == (colMap_ == null), std::logic_error, err );
    TEST_FOR_EXCEPTION( isFillComplete() == true && (colMap_ == null || rangeMap_ == null || domainMap_ == null), std::logic_error, err );
    TEST_FOR_EXCEPTION( isStorageOptimized() == true && nodeNumAllocated_ != nodeNumEntries_, std::logic_error, err );
    TEST_FOR_EXCEPTION( haveGlobalConstants_ == false && ( globalNumEntries_ != gsti || globalNumDiags_ != gsti || globalMaxNumRowEntries_ != gsti ), std::logic_error, err ); 
    TEST_FOR_EXCEPTION( haveGlobalConstants_ == true && ( globalNumEntries_ < nodeNumEntries_ || globalNumDiags_ < nodeNumDiags_ || globalMaxNumRowEntries_ < nodeMaxNumRowEntries_ ),
                        std::logic_error, err );
    TEST_FOR_EXCEPTION( nodeNumAllocated_ != 0 && indicesAreAllocated_ == false, std::logic_error, err );
    TEST_FOR_EXCEPTION( nodeNumAllocated_ != 0 && (numAllocPerRow_ != null || numAllocForAllRows_ != 0), std::logic_error, err );
    TEST_FOR_EXCEPTION( isStorageOptimized() && pftype_ == DynamicProfile, std::logic_error, err );
    TEST_FOR_EXCEPTION( pftype_ == DynamicProfile && nodeNumAllocated_ != 0 && pbuf_lclInds2D_ == null && pbuf_gblInds2D_ == null, std::logic_error, err );
    TEST_FOR_EXCEPTION( pftype_ == DynamicProfile && (pbuf_lclInds1D_ != null || pbuf_gblInds1D_ != null), std::logic_error, err );
    TEST_FOR_EXCEPTION( pftype_ == StaticProfile && nodeNumAllocated_ != 0 && pbuf_lclInds1D_ == null && pbuf_gblInds1D_ == null, std::logic_error, err );
    TEST_FOR_EXCEPTION( pftype_ == StaticProfile && (pbuf_lclInds2D_ != null || pbuf_gblInds2D_ != null), std::logic_error, err );
    TEST_FOR_EXCEPTION( pftype_ == DynamicProfile && pbuf_rowOffsets_ != null, std::logic_error, err );
    TEST_FOR_EXCEPTION( pftype_ == StaticProfile && nodeNumAllocated_ != 0 && pbuf_rowOffsets_ == null, std::logic_error, err );
    TEST_FOR_EXCEPTION( indicesAreAllocated_ == false && (pbuf_rowOffsets_ != null || numEntriesPerRow_ != null||
                                                          pbuf_lclInds1D_ != null || pbuf_lclInds2D_ != null ||
                                                          pbuf_gblInds1D_ != null || pbuf_gblInds2D_ != null), std::logic_error, err );
    TEST_FOR_EXCEPTION( nodeNumAllocated_ == 0 && (pbuf_rowOffsets_ != null || numEntriesPerRow_ != null), std::logic_error, err );
    TEST_FOR_EXCEPTION( nodeNumAllocated_ != 0 && storageOptimized_ == true  && numEntriesPerRow_ != null, std::logic_error, err );
    TEST_FOR_EXCEPTION( nodeNumAllocated_ != 0 && storageOptimized_ == false && numEntriesPerRow_ == null, std::logic_error, err );
    TEST_FOR_EXCEPTION( indicesAreLocal_ == true && (pbuf_gblInds1D_ != null || pbuf_gblInds2D_ != null), std::logic_error, err );
    TEST_FOR_EXCEPTION( indicesAreGlobal_ == true && (pbuf_lclInds1D_ != null || pbuf_lclInds2D_ != null), std::logic_error, err );
    TEST_FOR_EXCEPTION( indicesAreLocal_ == true && nodeNumAllocated_ != 0 && pbuf_lclInds1D_ == null && pbuf_lclInds2D_ == null, std::logic_error, err );
    TEST_FOR_EXCEPTION( indicesAreGlobal_ == true && nodeNumAllocated_ != 0 && pbuf_gblInds1D_ == null && pbuf_gblInds2D_ == null, std::logic_error, err );
    TEST_FOR_EXCEPTION( indicesAreAllocated_ == false && (nodeNumAllocated_ != 0 || nodeNumEntries_ != 0), std::logic_error, err );
    size_t actualNumAllocated = 0;
    if (pftype_ == DynamicProfile) {
      if (isGloballyIndexed() && pbuf_gblInds2D_ != Teuchos::null) {
        for (size_t r = 0; r < getNodeNumRows(); ++r) {
          actualNumAllocated += pbuf_gblInds2D_[r].size();
        }
      }
      else if (isLocallyIndexed() && pbuf_lclInds2D_ != Teuchos::null) {
        for (size_t r = 0; r < getNodeNumRows(); ++r) {
          actualNumAllocated += pbuf_lclInds2D_[r].size();
        }
      }
    }
    else { // pftype_ == StaticProfile)
      if (pbuf_rowOffsets_ != Teuchos::null) {
        Teuchos::ArrayRCP<const size_t> last_offset = node->template viewBuffer<size_t>(1,pbuf_rowOffsets_+getNodeNumRows());
        actualNumAllocated = last_offset[0];
        last_offset = Teuchos::null;
      }
      else {
        actualNumAllocated = 0;
      }
      TEST_FOR_EXCEPTION( storageOptimized_ == false && isLocallyIndexed() == true && (size_t)pbuf_lclInds1D_.size() != actualNumAllocated, std::logic_error, err );
      TEST_FOR_EXCEPTION(                              isGloballyIndexed() == true && (size_t)pbuf_gblInds1D_.size() != actualNumAllocated, std::logic_error, err );
    }
    TEST_FOR_EXCEPTION(actualNumAllocated != nodeNumAllocated_, std::logic_error, err );
#endif
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::description() const {
    std::ostringstream oss;
    oss << Teuchos::Describable::description();
    if (isFillComplete()) {
      oss << "{status = fill complete"
          << ", global rows = " << getGlobalNumRows()
          << ", global cols = " << getGlobalNumCols()
          << ", global num entries = " << getGlobalNumEntries()
          << "}";
    }
    else {
      oss << "{status = fill not complete"
          << ", global rows = " << getGlobalNumRows()
          << "}";
    }
    return oss.str();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
    using std::endl;
    using std::setw;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    Teuchos::EVerbosityLevel vl = verbLevel;
    if (vl == VERB_DEFAULT) vl = VERB_LOW;
    Teuchos::RCP<const Teuchos::Comm<int> > comm = this->getComm();
    const int myImageID = comm->getRank(),
              numImages = comm->getSize();
    size_t width = 1;
    for (size_t dec=10; dec<getGlobalNumRows(); dec *= 10) {
      ++width;
    }
    width = std::max<size_t>(width,11) + 2;
    Teuchos::OSTab tab(out);
    //    none: print nothing
    //     low: print O(1) info from node 0
    //  medium: print O(P) info, num entries per node
    //    high: print O(N) info, num entries per row
    // extreme: print O(NNZ) info: print graph indices
    // 
    // for medium and higher, print constituent objects at specified verbLevel
    if (vl != VERB_NONE) {
      if (myImageID == 0) out << this->description() << std::endl; 
      // O(1) globals, minus what was already printed by description()
      if (isFillComplete() && myImageID == 0) {
        out << "Global number of diagonals = " << globalNumDiags_ << std::endl;
        out << "Global max number of entries = " << globalMaxNumRowEntries_ << std::endl;
      }
      // constituent objects
      if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
        if (myImageID == 0) out << "\nRow map: " << std::endl;
        rowMap_->describe(out,vl);
        if (colMap_ != Teuchos::null) {
          if (myImageID == 0) out << "\nColumn map: " << std::endl;
          colMap_->describe(out,vl);
        }
        if (domainMap_ != Teuchos::null) {
          if (myImageID == 0) out << "\nDomain map: " << std::endl;
          domainMap_->describe(out,vl);
        }
        if (rangeMap_ != Teuchos::null) {
          if (myImageID == 0) out << "\nRange map: " << std::endl;
          rangeMap_->describe(out,vl);
        }
      }
      // O(P) data
      if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
        for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
          if (myImageID == imageCtr) {
            out << "Node ID = " << imageCtr << std::endl
                << "Node number of entries = " << nodeNumEntries_ << std::endl
                << "Node number of diagonals = " << nodeNumDiags_ << std::endl
                << "Node max number of entries = " << nodeMaxNumRowEntries_ << std::endl
                << "Node number of allocated entries = " << nodeNumAllocated_ << std::endl;
          }
          comm->barrier();
          comm->barrier();
          comm->barrier();
        }
      }
      // O(N) and O(NNZ) data
      if (vl == VERB_HIGH || vl == VERB_EXTREME) {
        for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
          if (myImageID == imageCtr) {
            out << std::setw(width) << "Node ID"
                << std::setw(width) << "Global Row" 
                << std::setw(width) << "Num Entries";
            if (vl == VERB_EXTREME) {
              out << "Entries";
            }
            out << std::endl;
            for (size_t r=0; r < getNodeNumRows(); ++r) {
              const size_t nE = RNNZ(r);
              GlobalOrdinal gid = rowMap_->getGlobalElement(r);
              out << std::setw(width) << myImageID 
                  << std::setw(width) << gid
                  << std::setw(width) << nE;
              if (vl == VERB_EXTREME) {
                if (isGloballyIndexed()) {
                  Teuchos::ArrayRCP<const GlobalOrdinal> rowview = getGlobalRowView(gid);
                  for (size_t j=0; j < nE; ++j) out << rowview[j] << " ";
                }
                else if (isLocallyIndexed()) {
                  Teuchos::ArrayRCP<const LocalOrdinal> rowview = getLocalRowView(r);
                  for (size_t j=0; j < nE; ++j) out << colMap_->getGlobalElement(rowview[j]) << " ";
                }
              }
              out << std::endl;
            }
          }
          comm->barrier();
          comm->barrier();
          comm->barrier();
        }
      }
    }
  }

} // namespace Tpetra

#endif
