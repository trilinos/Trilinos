/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_DETAILS_OVERLAPPINGROWGRAPH_DECL_HPP
#define IFPACK2_DETAILS_OVERLAPPINGROWGRAPH_DECL_HPP

#include <Ifpack2_ConfigDefs.hpp>
#include <Tpetra_RowGraph.hpp>
#include <Tpetra_Import_decl.hpp>
#include <Tpetra_Export_decl.hpp>

namespace Ifpack2 {
namespace Details {

/// \class OverlappingRowGraph
/// \brief Sparse graph (Tpetra::RowGraph subclass) with ghost rows.
/// \tparam GraphType Tpetra::RowGraph or Tpetra::CrsGraph specialization.
///
/// This class is meant to be created by and used with
/// OverlappingRowMatrix.  It is the subclass of Tpetra::RowGraph
/// returned by that class' getGraph() method.
///
/// \warning This class is an implementation detail of Ifpack2.  Users
///   should not rely on its interface.
template<class GraphType>
class OverlappingRowGraph :
    virtual public Tpetra::RowGraph<typename GraphType::local_ordinal_type,
                                    typename GraphType::global_ordinal_type,
                                    typename GraphType::node_type> {
public:
  //! \name Typedefs
  //@{
  typedef typename GraphType::local_ordinal_type local_ordinal_type;
  typedef typename GraphType::global_ordinal_type global_ordinal_type;
  typedef typename GraphType::node_type node_type;

  typedef Tpetra::Export<local_ordinal_type, global_ordinal_type, node_type> export_type;
  typedef Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type> import_type;
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;
  typedef Tpetra::RowGraph<local_ordinal_type, global_ordinal_type, node_type> row_graph_type;
  //@}
  //! \name Constructors and destructor
  //@{

  /// Standard constructor.
  ///
  /// \param nonoverlappingGraph [in] Graph containing only the
  ///   nonoverlapping rows.
  /// \param overlappingGraph [in] Graph containing only the
  ///   overlapping rows.
  /// \param rowMap [in] Map which is the union of the row Maps of the
  ///   two graphs.
  /// \param colMap [in] Column Map.
  /// \param numGlobalRows [in] Global number of rows.
  /// \param numGlobalCols [in] Global number of columns.
  /// \param numGlobalNonzeros [in] Global number of entries.
  /// \param maxNumEntries [in] Max number of entries in any row owned
  ///   by this process.
  /// \param nonoverlappingImporter [in] Import from the
  ///   nonoverlapping graph's row Map, to the row Map representing
  ///   the entire graph (nonoverlapping + overlapping).
  /// \param overlappingImporter [in] Import from the nonoverlapping
  ///   graph's row Map, to the overlapping graph's row Map.
  OverlappingRowGraph (const Teuchos::RCP<const row_graph_type>& nonoverlappingGraph,
                       const Teuchos::RCP<const row_graph_type>& overlappingGraph,
                       const Teuchos::RCP<const map_type>& rowMap,
                       const Teuchos::RCP<const map_type>& colMap,
                       const Tpetra::global_size_t numGlobalRows,
                       const Tpetra::global_size_t numGlobalCols,
                       const Tpetra::global_size_t numGlobalNonzeros,
                       const size_t maxNumEntries,
                       const Teuchos::RCP<const import_type>& nonoverlappingImporter,
                       const Teuchos::RCP<const import_type>& overlappingImporter);
  //! Destructor
  virtual ~OverlappingRowGraph ();

  //@}
  //! @name Matrix query methods
  //@{

  //! The communicator over which the graph is distributed.
  virtual Teuchos::RCP<const Teuchos::Comm<int> > getComm () const;

  //! The graph's Node instance.
  virtual Teuchos::RCP<node_type> getNode () const;

  //! The Map that describes the distribution of rows over processes.
  virtual Teuchos::RCP<const map_type> getRowMap () const;

  //! The Map that describes the distribution of columns over processes.
  virtual Teuchos::RCP<const map_type> getColMap () const;

  /// \brief The Map that describes the domain of this graph.
  ///
  /// The domain is the distribution of valid input vectors of
  /// apply(), for a matrix whose graph is <tt>*this</tt>.
  virtual Teuchos::RCP<const map_type> getDomainMap () const;

  /// \brief The Map that describes the range of this graph.
  ///
  /// The range is the distribution of valid output vectors of
  /// apply(), for a matrix whose graph is <tt>*this</tt>.
  virtual Teuchos::RCP<const map_type> getRangeMap () const;

  //! Import object (from domain Map to column Map).
  virtual Teuchos::RCP<const import_type> getImporter () const;

  //! Export object (from row Map to range Map).
  virtual Teuchos::RCP<const export_type> getExporter () const;

  //! The global number of rows in this graph.
  virtual global_size_t getGlobalNumRows () const;

  //! The global number of columns in this graph.
  virtual global_size_t getGlobalNumCols () const;

  //! The number of rows owned by the calling process.
  virtual size_t getNodeNumRows () const;

  /// \brief The number of columns owned by the calling process.
  ///
  /// This is the number of columns needed to apply the forward
  /// operator on the calling process, that is, the number of elements
  /// listed in the column Map on the calling process.
  virtual size_t getNodeNumCols () const;

  //! The index base for global indices for this graph.
  virtual global_ordinal_type getIndexBase () const;

  //! The global number of entries in this graph.
  virtual global_size_t getGlobalNumEntries () const;

  //! The number of entries in this graph owned by the calling process.
  virtual size_t getNodeNumEntries () const;

  /// \brief The number of entries in the given global row that are
  ///   owned by the calling process.
  ///
  /// \param globalRow [in] Global index of the row.
  ///
  /// \return Teuchos::OrdinalTraits<size_t>::invalid() if the
  ///   specified row (either in the input graph or in the overlap
  ///   graph) is not owned by the calling process, else the number of
  ///   entries in that row that are owned by the calling process.
  virtual size_t getNumEntriesInGlobalRow (global_ordinal_type globalRow) const;

  /// \brief The number of entries in the given local row that are
  ///   owned by the calling process.
  ///
  /// \param globalRow [in] Local index of the row.
  ///
  /// \return Teuchos::OrdinalTraits<size_t>::invalid() if the
  ///   specified row (either in the input graph or in the overlap
  ///   graph) is not owned by the calling process, else the number of
  ///   entries in that row that are owned by the calling process.
  virtual size_t getNumEntriesInLocalRow (local_ordinal_type localRow) const;

  //! The global number of diagonal entries.
  virtual global_size_t getGlobalNumDiags () const;

  //! The number of diagonal entries owned by the calling process.
  virtual size_t getNodeNumDiags () const;

  //! The maximum number of entries in any row on any process.
  virtual size_t getGlobalMaxNumRowEntries () const;

  //! The maximum number of entries in any row on the calling process.
  virtual size_t getNodeMaxNumRowEntries() const;

  //! Whether this graph has a column Map.
  virtual bool hasColMap() const;

  //! Whether this graph is lower triangular.
  virtual bool isLowerTriangular() const;

  //! Whether this graph is upper triangular.
  virtual bool isUpperTriangular() const;

  //! Whether this graph is locally indexed.
  virtual bool isLocallyIndexed () const;

  //! Whether this graph is globally indexed.
  virtual bool isGloballyIndexed () const;

  //! \c true if fillComplete() has been called, else \c false.
  virtual bool isFillComplete() const;

  //@}
  //! @name Extraction methods
  //@{

  /// \brief Copy out a list of column indices in the given global row
  ///   that are owned by the calling process.
  ///
  /// \param globalRow [in] Global index of the row.
  /// \param indices [out] Global column indices in that row that are
  ///   owned by the calling process.
  /// \param numIndices [out] Number of indices returned in \c indices.
  ///
  /// This method throws std::runtime_error if \c indices is not large
  /// enough to hold the column indices in row \c globalRow.  If row
  /// \c globalRow does not belong to this process, then \c indices is
  /// not modified and \c numIndices is set to
  /// Teuchos::OrdinalTraits<size_t>::invalid() on output.
  virtual void
  getGlobalRowCopy (global_ordinal_type globalRow,
                    const Teuchos::ArrayView<global_ordinal_type>& indices,
                    size_t& numIndices) const;

  /// \brief Copy out a list of local column indices in the given
  ///   local row that are owned by the calling process.
  ///
  /// This method only makes sense to call if the graph has a column
  /// Map, since the local column indices are computed with respect to
  /// the column Map.
  ///
  /// \param localRow [in] Local index of the row.
  /// \param indices [out] Local column indices in that row that are
  ///   owned by the calling process.
  /// \param numIndices [out] Number of indices returned in \c indices.
  ///
  /// This method throws std::runtime_error if \c indices is not large
  /// enough to hold the column indices in row \c localRow.  If row
  /// <tt>localRow</tt> does not belong to this process, then
  /// <tt>indices</tt> is not modified and \c numIndices is set to
  /// Teuchos::OrdinalTraits<size_t>::invalid() on output.
  virtual void
  getLocalRowCopy (local_ordinal_type localRow,
                   const Teuchos::ArrayView<local_ordinal_type>& indices,
                   size_t& numIndices) const;
  //@}
private:
  //! \name Internal data
  //@{
  Teuchos::RCP<const row_graph_type> nonoverlappingGraph_;
  Teuchos::RCP<const row_graph_type> overlappingGraph_;
  Teuchos::RCP<const map_type> rowMap_;
  Teuchos::RCP<const map_type> colMap_;
  const Tpetra::global_size_t numGlobalRows_;
  const Tpetra::global_size_t numGlobalCols_;
  const Tpetra::global_size_t numGlobalNonzeros_;
  const size_t maxNumEntries_;
  Teuchos::RCP<const import_type> nonoverlappingImporter_;
  Teuchos::RCP<const import_type> overlappingImporter_;
  //@}
};

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_OVERLAPPINGROWGRAPH_DECL_HPP
