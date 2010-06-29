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

#ifndef TPETRA_CRSGRAPH_DECL_HPP
#define TPETRA_CRSGRAPH_DECL_HPP

// TODO: filter column indices first in insertLocalIndices()
// TODO: filter column indices first in insertGlobalIndices()

#include <Teuchos_Describable.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_NullIteratorTraits.hpp>
#include <Teuchos_CompileTimeAssert.hpp>

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_CrsGraph.hpp>
#include <Kokkos_NodeHelpers.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_RowGraph.hpp"
#include "Tpetra_DistObject.hpp"
#include "Tpetra_Util.hpp"

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration
  template <class S, class LO, class GO, class N, class SpMatOps>
  class CrsMatrix;
#endif

  struct RowInfo {
    size_t allocSize;
    size_t numEntries;
    size_t offset1D;
  };

  //! \brief A class for constructing and using sparse compressed index graphs with row access.
  /*! This class is templated on \c LocalOrdinal and \c GlobalOrdinal. If the \c GlobalOrdinal is not specified, then 
   *  it takes the same type as the \c LocalOrdinal.
   */
  template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class CrsGraph : public RowGraph<LocalOrdinal,GlobalOrdinal,Node>,
                   public DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Node> {
    template <class S, class LO, class GO, class N, class SpMatOps>
    friend class CrsMatrix;

    public: 
      typedef LocalOrdinal  local_ordinal_type;
      typedef GlobalOrdinal global_ordinal_type;
      typedef Node          node_type;

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

      //! Remove graph indices from local row.
      void removeLocalIndices(LocalOrdinal row);

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
          this graph.  

          This quantity is computed during the actual allocation. Therefore, if <tt>indicesAreAllocated() == false</tt>,
          this method returns <tt>Teuchos::OrdinalTraits<size_t>::invalid()</tt>.
      */
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

      //! @name Methods implementing Tpetra::DistObject
      //@{

      bool checkSizes(const DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Node>& source);

      void copyAndPermute(const DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Node> & source,
                          size_t numSameIDs,
                          const Teuchos::ArrayView<const LocalOrdinal> &permuteToLIDs,
                          const Teuchos::ArrayView<const LocalOrdinal> &permuteFromLIDs);

      void packAndPrepare(const DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Node> & source,
                          const Teuchos::ArrayView<const LocalOrdinal> &exportLIDs,
                          Teuchos::Array<GlobalOrdinal> &exports,
                          const Teuchos::ArrayView<size_t> & numPacketsPerLID,
                          size_t& constantNumPackets,
                          Distributor &distor);

      void unpackAndCombine(const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
                            const Teuchos::ArrayView<const GlobalOrdinal> &imports,
                            const Teuchos::ArrayView<size_t> &numPacketsPerLID,
                            size_t constantNumPackets,
                            Distributor &distor,
                            CombineMode CM);
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
      size_t findLocalIndex(size_t row, LocalOrdinal ind, const Teuchos::ArrayRCP<const LocalOrdinal> &alreadyHaveAView) const;
      size_t findGlobalIndex(size_t row, GlobalOrdinal ind, const Teuchos::ArrayRCP<const GlobalOrdinal> &alreadyHaveAView) const;
      inline size_t RNNZ(size_t row) const;
      inline size_t RNumAlloc(size_t row) const;
      void checkInternalState() const;
      void clearGlobalConstants();
      void updateLocalAllocation(size_t lrow, size_t allocSize);
      void updateGlobalAllocation(size_t lrow, size_t allocSize);
      void fillLocalGraph();

      //! \brief Get the sizes associated with the allocated rows.
      /*! This is used by the row view routines and others. It computes the size and offset information
          for a particular row. It is designed to do this with minimum overhead. No checking is done except in a debug build.

        \param myRow      - (In) \c size_t specifying the local row.
        \returns \c RowInfo struct specifying the size of the allocation for the specified row, the number of entries, and the 
                 offset into 1D allocation, if <tt>getProfileType() == StaticProfile</tt>.
      */
      RowInfo getRowInfo(size_t myRow) const;

      //! \brief Get a persisting const view of the elements in a specified local row of the graph, along with other details.
      /*! This protected method is used internally for almost all access to the graph elements. It is designed to provide the information 
          needed by CrsGraph and CrsMatrix with as little overhead as possible. No checking is done except in a debug build.

        \param myRow      - (In) \c size_t specifying the local row.
        \param indices    - (Out) persisting, const view of the local indices. <tt>indices.size()</tt> specifies the size of the allocation.

        \returns Returns row info; see getRowInfo().

        \pre isGloballyIndexed()==false
       */
      RowInfo getFullLocalView(size_t myRow, Teuchos::ArrayRCP<const LocalOrdinal> &indices) const;

      //! \brief Get a persisting non-const view of the elements in a specified local row of the graph, along with other details.
      /*! This protected method is used internally for almost all access to the graph elements. It is designed to provide the information 
          needed by CrsGraph and CrsMatrix with as little overhead as possible. No checking is done except in a debug build.

        \param myRow      - (In) \c size_t specifying the local row.
        \param indices    - (Out) persisting, non-const view of the local indices. <tt>indices.size()</tt> specifies the size of the allocation.

        \returns Returns row info; see getRowInfo().

        \pre isGloballyIndexed()==false
       */
      RowInfo getFullLocalViewNonConst(size_t myRow, Teuchos::ArrayRCP<LocalOrdinal> &indices);

      //! \brief Get a persisting const view of the elements in a specified global row of the graph, along with other details.
      /*! This protected method is used internally for almost all access to the graph elements. It is designed to provide the information 
          needed by CrsGraph and CrsMatrix with as little overhead as possible. No checking is done except in a debug build.

        \param myRow      - (In) \c size_t specifying the local row.
        \param indices    - (Out) persisting, const view of the local indices. <tt>indices.size()</tt> specifies the size of the allocation.

        \returns Returns row info; see getRowInfo().

        \pre isLocallyIndexed()==false
       */
      RowInfo getFullGlobalView(size_t myRow, Teuchos::ArrayRCP<const GlobalOrdinal> &indices) const;

      //! \brief Get a persisting non-const view of the elements in a specified local row of the graph, along with other details.
      /*! This protected method is used internally for almost all access to the graph elements. It is designed to provide the information 
          needed by CrsGraph and CrsMatrix with as little overhead as possible. No checking is done except in a debug build.

        \param myRow      - (In) \c size_t specifying the local row.
        \param indices    - (Out) persisting, non-const view of the local indices. <tt>indices.size()</tt> specifies the size of the allocation.

        \returns Returns row info; see getRowInfo().

        \pre isLocallyIndexed()==false
       */
      RowInfo getFullGlobalViewNonConst(size_t myRow, Teuchos::ArrayRCP<GlobalOrdinal> &indices);

      void insertLocalIndicesViaView(size_t myRow, const Teuchos::ArrayView<const LocalOrdinal> &indices, const Teuchos::ArrayRCP<LocalOrdinal> &inds_view);
      void insertGlobalIndicesViaView(size_t myRow, const Teuchos::ArrayView<const GlobalOrdinal> &indices, const Teuchos::ArrayRCP<GlobalOrdinal> &inds_view);

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
      // after allocation, except during makeIndicesLocal(), one of these is Teuchos::null.
      // this is a parallel compute buffer, not host memory
      // 1D == StaticAllocation, 2D == DynamicAllocation
      Teuchos::ArrayRCP< LocalOrdinal>                     pbuf_lclInds1D_, view_lclInds1D_;
      Teuchos::ArrayRCP<GlobalOrdinal>                     pbuf_gblInds1D_, view_gblInds1D_;
      Teuchos::ArrayRCP<Teuchos::ArrayRCP< LocalOrdinal> > pbuf_lclInds2D_, view_lclInds2D_;
      Teuchos::ArrayRCP<Teuchos::ArrayRCP<GlobalOrdinal> > pbuf_gblInds2D_, view_gblInds2D_;
      // offset to the first entry of each row. length numRows + 1. only used if pftype_ == StaticAllocation, otherwise set to Teuchos::null.
      // also, set to null if nodeNumAllocated_ == 0
      Teuchos::ArrayRCP<      size_t>                      pbuf_rowOffsets_, view_rowOffsets_;

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

  }; // class CrsGrap

} // namespace Tpetra

#endif
