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

#ifndef TPETRA_ROWGRAPH_HPP
#define TPETRA_ROWGRAPH_HPP

#include <Teuchos_Describable.hpp>
#include <Teuchos_ArrayView.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Map.hpp"

namespace Tpetra 
{
  
  //! \brief A pure virtual interface for row-partitioned graphs.
  /*!
     This class is templated on \c LocalOrdinal and \c GlobalOrdinal. 
     The \c GlobalOrdinal type, if omitted, defaults to the \c LocalOrdinal type.
   */
  template<class LocalOrdinal, class GlobalOrdinal = LocalOrdinal>
  class RowGraph : public Teuchos::Describable
  {
    public: 

      //! @name Graph Query Methods
      //@{ 

      //! Returns the communicator.
      virtual Teuchos::RCP<const Teuchos::Comm<int> > getComm() const = 0;

      //! Returns the Map that describes the row distribution in this graph.
      virtual const Map<LocalOrdinal,GlobalOrdinal> & getRowMap() const = 0;

      //! \brief Returns the Map that describes the column distribution in this graph.
      /*! */
      virtual const Map<LocalOrdinal,GlobalOrdinal> & getColMap() const = 0;

      //! Returns the Map associated with the domain of this graph.
      virtual const Map<LocalOrdinal,GlobalOrdinal> & getDomainMap() const = 0;

      //! Returns the Map associated with the domain of this graph.
      virtual const Map<LocalOrdinal,GlobalOrdinal> & getRangeMap() const = 0;

      //! Returns the importer associated with this graph.
      virtual Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal> > getImporter() const = 0;

      //! Returns the exporter associated with this graph.
      virtual Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal> > getExporter() const = 0;

      //! Returns the number of global rows in the graph.
      virtual GlobalOrdinal numGlobalRows() const = 0;

      //! \brief Returns the number of global columns in the graph.
      /*! */
      virtual GlobalOrdinal numGlobalCols() const = 0;

      //! Returns the number of rows owned on the calling node.
      virtual Teuchos_Ordinal numLocalRows() const = 0;

      //! Returns the index base for global indices for this graph. 
      virtual Teuchos_Ordinal getIndexBase() const = 0;

      //! Returns the number of columns connected to the locally owned rows of this graph.
      /*! */
      virtual Teuchos_Ordinal numLocalCols() const = 0;

      //! Returns the global number of entries in the graph.
      virtual GlobalOrdinal numGlobalEntries() const = 0;

      //! Returns the local number of entries in the graph.
      virtual Teuchos_Ordinal numMyEntries() const = 0;

      //! \brief Returns the current number of graph entries on this node in the specified global row .
      /*! Returns Teuchos::OrdinalTraits<Teuchos_Ordinal>::invalid() if the specified global row does not belong to this graph. */
      virtual Teuchos_Ordinal numEntriesForGlobalRow(GlobalOrdinal globalRow) const = 0;

      //! Returns the current number of graph entries on this node in the specified local row.
      /*! Returns Teuchos::OrdinalTraits<Teuchos_Ordinal>::invalid() if the specified local row is not valid for this graph. */
      virtual Teuchos_Ordinal numEntriesForMyRow(LocalOrdinal localRow) const = 0;

      //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons. 
      /*! */
      virtual GlobalOrdinal numGlobalDiagonals() const = 0;

      //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons. 
      /*! */
      virtual Teuchos_Ordinal numMyDiagonals() const = 0;

      //! \brief Returns the maximum number of entries across all rows/columns on all nodes. 
      /*! */
      virtual GlobalOrdinal globalMaxNumRowEntries() const = 0;

      //! \brief Returns the maximum number of entries across all rows/columns on this node. 
      /*! */
      virtual Teuchos_Ordinal myMaxNumRowEntries() const = 0;

      //! \brief Indicates whether the graph has a well-defined column map. 
      /*! */
      virtual bool hasColMap() const = 0;

      //! \brief Indicates whether the graph is lower triangular.
      virtual bool lowerTriangular() const = 0;

      //! \brief Indicates whether the graph is upper triangular.
      virtual bool upperTriangular() const = 0;

      //! \brief If graph indices have been transformed to local, this function returns true. Otherwise, it returns false. */
      /*! 
      */
      virtual bool indicesAreLocal() const = 0;

      //! \brief If graph indices are global, this function returns true. Otherwise, it returns false. */
      /*! 
      */
      virtual bool indicesAreGlobal() const = 0;

      //! Returns \c true if fillComplete() has been called.
      virtual bool isFillComplete() const = 0;

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
      virtual void extractGlobalRowCopy(GlobalOrdinal GlobalRow, Teuchos::ArrayView<GlobalOrdinal> indices, Teuchos_Ordinal &NumIndices) const = 0;

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
      virtual void extractMyRowCopy(LocalOrdinal LocalRow, Teuchos::ArrayView<LocalOrdinal> indices, Teuchos_Ordinal& NumIndices) const = 0;

      //! Get a non-persisting view of the elements in a specified global row of the graph.
      /*!
        \param GlobalRow - (In) Global row number to get indices.
        \param Indices - (Out) Global column indices corresponding to values.

         Note: If \c GlobalRow does not belong to this node, then \c indices is unchanged and \c NumIndices is 
         returned as Teuchos::OrdinalTraits<Teuchos_Ordinal>::invalid().

        \pre indicesAreLocal()==false
       */
      virtual void extractGlobalRowConstView(GlobalOrdinal GlobalRow, Teuchos::ArrayView<const GlobalOrdinal> &indices) const = 0;

      //! Get a view of the elements in a specified local row of the graph.
      /*!
        \param LocalRow - (In) Local row number to get indices.
        \param Indices - (Out) Local column indices corresponding to values.

         Note: If \c LocalRow is not valid for this node, then \c indices is unchanged and \c NumIndices is 
         returned as Teuchos::OrdinalTraits<Teuchos_Ordinal>::invalid().

        \pre indicesAreLocal()==true
       */
      virtual void extractMyRowConstView(LocalOrdinal LocalRow, Teuchos::ArrayView<const LocalOrdinal> &indices) const = 0;

      //@}

      //! @name I/O Methods
      //@{ 
      
      //! Prints the matrix on the specified stream. This is very verbose.
      void print(std::ostream& os) const;

      // @}
  };

} // namespace Tpetra

#endif
