// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
// @HEADER

#ifndef TPETRA_ROWGRAPH_HPP
#define TPETRA_ROWGRAPH_HPP

#include <Teuchos_Describable.hpp>
#include <Kokkos_DefaultNode.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"

namespace Tpetra 
{
  
  //! \brief A pure virtual interface for row-partitioned graphs.
  /*!
     This class is templated on \c LocalOrdinal, \c GlobalOrdinal and \c Node.
     The \c LocalOrdinal type, if omitted, defaults to \c int. 
     The \c GlobalOrdinal type defaults to the \c LocalOrdinal type.
     The \c Node type defaults to the default node in Kokkos.
   */
  template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class RowGraph : virtual public Teuchos::Describable {
    public: 
      typedef LocalOrdinal  local_ordinal_type;
      typedef GlobalOrdinal global_ordinal_type;
      typedef Node          node_type;

      //! @name Destructor Method
      //@{ 

      //! Destructor.
      virtual ~RowGraph();

      //@}

      //! @name Graph Query Methods
      //@{ 

      //! Returns the communicator.
      virtual const Teuchos::RCP<const Teuchos::Comm<int> > & getComm() const = 0;

      //! Returns the underlying node.
      virtual Teuchos::RCP<Node> getNode() const = 0;

      //! Returns the Map that describes the row distribution in this graph.
      virtual const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRowMap() const = 0;

      //! \brief Returns the Map that describes the column distribution in this graph.
      virtual const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getColMap() const = 0;

      //! Returns the Map associated with the domain of this graph.
      virtual const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getDomainMap() const = 0;

      //! Returns the Map associated with the domain of this graph.
      virtual const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRangeMap() const = 0;

      //! Returns the importer associated with this graph.
      virtual Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > getImporter() const = 0;

      //! Returns the exporter associated with this graph.
      virtual Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > getExporter() const = 0;

      //! Returns the number of global rows in the graph.
      virtual global_size_t getGlobalNumRows() const = 0;

      //! \brief Returns the number of global columns in the graph.
      virtual global_size_t getGlobalNumCols() const = 0;

      //! Returns the number of rows owned on the calling node.
      virtual size_t getNodeNumRows() const = 0;

      //! Returns the number of columns connected to the locally owned rows of this graph.
      virtual size_t getNodeNumCols() const = 0;

      //! Returns the index base for global indices for this graph. 
      virtual GlobalOrdinal getIndexBase() const = 0;

      //! Returns the global number of entries in the graph.
      virtual global_size_t getGlobalNumEntries() const = 0;

      //! Returns the local number of entries in the graph.
      virtual size_t getNodeNumEntries() const = 0;

      //! \brief Returns the current number of entries on this node in the specified global row.
      /*! Returns Teuchos::OrdinalTraits<size_t>::invalid() if the specified global row does not belong to this graph. */
      virtual size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const = 0;

      //! Returns the current number of entries on this node in the specified local row.
      /*! Returns Teuchos::OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this graph. */
      virtual size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const = 0;

      //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons. 
      virtual global_size_t getGlobalNumDiags() const = 0;

      //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons. 
      virtual size_t getNodeNumDiags() const = 0;

      //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
      virtual size_t getGlobalMaxNumRowEntries() const = 0;

      //! \brief Returns the maximum number of entries across all rows/columns on this node.
      virtual size_t getNodeMaxNumRowEntries() const = 0;

      //! \brief Indicates whether the graph has a well-defined column map. 
      virtual bool hasColMap() const = 0;

      //! \brief Indicates whether the graph is lower triangular.
      virtual bool isLowerTriangular() const = 0;

      //! \brief Indicates whether the graph is upper triangular.
      virtual bool isUpperTriangular() const = 0;

      //! \brief If graph indices are in the local range, this function returns true. Otherwise, this function returns false. */
      virtual bool isLocallyIndexed() const = 0;

      //! \brief If graph indices are in the global range, this function returns true. Otherwise, this function returns false. */
      virtual bool isGloballyIndexed() const = 0;

      //! Returns \c true if fillComplete() has been called.
      virtual bool isFillComplete() const = 0;

      //@}

      //! @name Extraction Methods
      //@{ 
          
      //! Extract a list of entries in a specified global row of the graph. Put into pre-allocated storage.
      /*!
        \param LocalRow - (In) Global row number for which indices are desired.
        \param Indices - (Out) Global column indices corresponding to values.
        \param NumIndices - (Out) Number of indices.

         Note: A std::runtime_error exception is thrown if \c Indices is not large enough to hold the column indices associated
         with row \c GlobalRow. If \c GlobalRow does not belong to this node, then \c Indices is unchanged and \c NumIndices is 
         returned as Teuchos::OrdinalTraits<size_t>::invalid().
       */
      virtual void getGlobalRowCopy(GlobalOrdinal GlobalRow,
                                    const Teuchos::ArrayView<GlobalOrdinal> &Indices,
                                    size_t &NumIndices) const = 0;

      //! Extract a list of entries in a specified local row of the graph. Put into storage allocated by calling routine.
      /*!
        \param LocalRow - (In) Local row number for which indices are desired.
        \param Indices - (Out) Local column indices corresponding to values.
        \param NumIndices - (Out) Number of indices.

         Note: A std::runtime_error exception is thrown if \c Indices is not large enough to hold the column indices associated
         with row \c LocalRow. If \c LocalRow is not valid for this node, then \c Indices is unchanged and \c NumIndices is 
         returned as Teuchos::OrdinalTraits<size_t>::invalid().
       */
      virtual void getLocalRowCopy(LocalOrdinal LocalRow, 
                                   const Teuchos::ArrayView<LocalOrdinal> &Indices,
                                   size_t &NumIndices) const = 0;

      //@}

  }; // class RowGraph

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  RowGraph<LocalOrdinal,GlobalOrdinal,Node>::~RowGraph() {
  }


} // namespace Tpetra

#endif
