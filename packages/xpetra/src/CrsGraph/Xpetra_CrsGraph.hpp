// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef XPETRA_CRSGRAPH_HPP
#define XPETRA_CRSGRAPH_HPP

#include <Teuchos_ParameterList.hpp>

#include <Teuchos_Describable.hpp>
#include <Kokkos_DefaultNode.hpp>
#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DistObject.hpp"
#include "Xpetra_Exceptions.hpp"

#include "Xpetra_Map.hpp"

#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
#ifdef HAVE_XPETRA_TPETRA
#include <Kokkos_StaticCrsGraph.hpp>
#endif
#endif

namespace Xpetra {

  using Teuchos::ParameterList;

  struct RowInfo {
    size_t localRow;
    size_t allocSize;
    size_t numEntries;
    size_t offset1D;
  };

  enum ELocalGlobal {
    LocalIndices,
    GlobalIndices
  };

  template <class LocalOrdinal,
            class GlobalOrdinal,
            class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class CrsGraph
    : /*public RowGraph<>,*/ public DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Node>
  {
  public:
    typedef LocalOrdinal local_ordinal_type;
    typedef GlobalOrdinal global_ordinal_type;
    typedef Node node_type;

    //! @name Constructor/Destructor Methods
    //@{

    //! Destructor.
    virtual ~CrsGraph() { }

   //@}

    //! @name Insertion/Removal Methods
    //@{

    //! Insert global indices into the graph.
    virtual void insertGlobalIndices(GlobalOrdinal globalRow, const ArrayView< const GlobalOrdinal > &indices)= 0;

    //! Insert local indices into the graph.
    virtual void insertLocalIndices(const LocalOrdinal localRow, const ArrayView< const LocalOrdinal > &indices)= 0;

    //! Remove all graph indices from the specified local row.
    virtual void removeLocalIndices(LocalOrdinal localRow)= 0;

    //@}

    //! @name Transformational Methods
    //@{

    //! Signal that data entry is complete, specifying domain and range maps.
    virtual void fillComplete(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &domainMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rangeMap, const RCP< ParameterList > &params=null)= 0;

    //! Signal that data entry is complete.
    virtual void fillComplete(const RCP< ParameterList > &params=null)= 0;

    //@}

    //! @name Methods implementing RowGraph.
    //@{

    //! Returns the communicator.
    virtual RCP< const Comm< int > > getComm() const = 0;

    //! Returns the Map that describes the row distribution in this graph.
    virtual RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > getRowMap() const = 0;

    //! Returns the Map that describes the column distribution in this graph.
    virtual RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > getColMap() const = 0;

    //! Returns the Map associated with the domain of this graph.
    virtual RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > getDomainMap() const = 0;

    //! Returns the Map associated with the domain of this graph.
    virtual RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > getRangeMap() const = 0;

    //! Returns the importer associated with this graph.
    virtual RCP< const Import< LocalOrdinal, GlobalOrdinal, Node > > getImporter() const = 0;

    //! Returns the exporter associated with this graph.
    virtual RCP< const Export< LocalOrdinal, GlobalOrdinal, Node > > getExporter() const = 0;

    //! Returns the number of global rows in the graph.
    virtual global_size_t getGlobalNumRows() const = 0;

    //! Returns the number of global columns in the graph.
    virtual global_size_t getGlobalNumCols() const = 0;

    //! Returns the number of graph rows owned on the calling node.
    virtual size_t getNodeNumRows() const = 0;

    //! Returns the number of columns connected to the locally owned rows of this graph.
    virtual size_t getNodeNumCols() const = 0;

    //! Returns the index base for global indices for this graph.
    virtual GlobalOrdinal getIndexBase() const = 0;

    //! Returns the global number of entries in the graph.
    virtual global_size_t getGlobalNumEntries() const = 0;

    //! Returns the local number of entries in the graph.
    virtual size_t getNodeNumEntries() const = 0;

    //! Returns the current number of entries on this node in the specified global row.
    virtual size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const = 0;

    //! Returns the current number of entries on this node in the specified local row.
    virtual size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const = 0;

    //! Returns the current number of allocated entries for this node in the specified global row .
    virtual size_t getNumAllocatedEntriesInGlobalRow(GlobalOrdinal globalRow) const = 0;

    //! Returns the current number of allocated entries on this node in the specified local row.
    virtual size_t getNumAllocatedEntriesInLocalRow(LocalOrdinal localRow) const = 0;

    //! Maximum number of entries in all rows over all processes.
    virtual size_t getGlobalMaxNumRowEntries() const = 0;

    //! Maximum number of entries in all rows owned by the calling process.
    virtual size_t getNodeMaxNumRowEntries() const = 0;

    //! Whether the graph has a column Map.
    virtual bool hasColMap() const = 0;

    //! Whether column indices are stored using local indices on the calling process.
    virtual bool isLocallyIndexed() const = 0;

    //! Whether column indices are stored using global indices on the calling process.
    virtual bool isGloballyIndexed() const = 0;

    //! Whether fillComplete() has been called and the graph is in compute mode.
    virtual bool isFillComplete() const = 0;

    //! Returns true if storage has been optimized.
    virtual bool isStorageOptimized() const = 0;

    //! Return a const, nonpersisting view of global indices in the given row.
    virtual void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView< const GlobalOrdinal > &Indices) const = 0;

    //! Return a const, nonpersisting view of local indices in the given row.
    virtual void getLocalRowView(LocalOrdinal LocalRow, ArrayView< const LocalOrdinal > &indices) const = 0;

    //! Force the computation of global constants if we don't have them
    virtual void computeGlobalConstants() =0;

    //@}

    //! @name Tpetra-specific routines
    //@{
#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
#ifdef HAVE_XPETRA_TPETRA
    typedef typename node_type::execution_space execution_space;
    typedef typename node_type::device_type device_type;
    typedef Kokkos::StaticCrsGraph<LocalOrdinal, Kokkos::LayoutLeft, device_type, void, size_t> local_graph_type;

    /// \brief Get the local graph.
    ///
    /// \warning THIS IS AN EXPERT MODE FUNCTION.  THIS IS AN
    ///   IMPLEMENTATION DETAIL.  DO NOT CALL THIS FUNCTION!!!
    ///
    /// This is only a valid representation of the local graph if the
    /// (global) graph is fill complete.
    virtual local_graph_type getLocalGraph () const = 0;
#else
#ifdef __GNUC__
#warning "Xpetra Kokkos interface for CrsMatrix is enabled (HAVE_XPETRA_KOKKOS_REFACTOR) but Tpetra is disabled. The Kokkos interface needs Tpetra to be enabled, too."
#endif
#endif
#endif

    //@}

    //! @name Overridden from Teuchos::Describable
    //@{

    //! Return a simple one-line description of this object.
    virtual std::string description() const = 0;

    //! Print the object with some verbosity level to an FancyOStream object.
    virtual void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const = 0;

    //@}

  }; // CrsGraph class

} // Xpetra namespace

#define XPETRA_CRSGRAPH_SHORT
#endif // XPETRA_CRSGRAPH_HPP
