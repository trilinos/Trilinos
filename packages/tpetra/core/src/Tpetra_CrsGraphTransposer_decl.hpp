
// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_CRSGRAPHTRANSPOSER_DECL_HPP
#define TPETRA_CRSGRAPHTRANSPOSER_DECL_HPP

/// \file Tpetra_CrsGraphTransposer_decl.hpp
///
/// Declaration of Tpetra::CrsGraphTransposer.

#include "Tpetra_CrsGraphTransposer_fwd.hpp"
#include "Tpetra_CrsGraph_fwd.hpp"
#include "Tpetra_Map_fwd.hpp"
#include "Teuchos_RCP.hpp"
#include <string>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Teuchos {
  // Forward declaration of ParameterList
  class ParameterList;
} // namespace Teuchos
#endif // DOXYGEN_SHOULD_SKIP_THIS

namespace Tpetra {

/// \class CrsGraphTransposer
/// \brief Construct and (optionally) redistribute the explicitly
///   stored transpose of a CrsGraph.
///
/// This class is based on the EpetraExt version.  It first transposes
/// the graph to an intermediate version with overlapping row map.
/// That graph is then converted to a final version whose row map is
/// "unique", i.e., a row is wholly owned by one process.
///
/// This class takes the same template parameters as CrsGraph.
template<class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
class CrsGraphTransposer {
public:
  //! @name Typedefs
  //@{
  typedef LocalOrdinal local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef Node node_type;

  typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
  typedef CrsGraph<LocalOrdinal, GlobalOrdinal, Node> crs_graph_type;

  //@}
  //! @name Constructors
  //@{

  //! Constructor that takes the graph to transpose.
  CrsGraphTransposer (const Teuchos::RCP<const crs_graph_type>& origGraph,const std::string & label = std::string());

  //@}
  //! @name Methods for computing the explicit transpose.
  //@{

  //! Compute and return graph+graph^T of the graph given to the constructor.
  Teuchos::RCP<crs_graph_type> symmetrize(const Teuchos::RCP<Teuchos::ParameterList> &params=Teuchos::null);

  //! Compute and return the transpose of the graph given to the constructor.
  Teuchos::RCP<crs_graph_type> createTranspose(const Teuchos::RCP<Teuchos::ParameterList> &params=Teuchos::null);

  /// \brief Compute and return the transpose of the graph given to the constructor.
  ///
  /// In this call, we (potentially) leave the graph with an
  /// overlapping row Map.  This is a perfectly valid graph, but
  /// won't work correctly with some routines in Ifpack or Muelu.
  ///
  /// \warning This routine leaves overlapping rows.  Unless you're
  /// sure that's OK, call createTranspose() instead.
  Teuchos::RCP<crs_graph_type> createTransposeLocal(const Teuchos::RCP<Teuchos::ParameterList> &params=Teuchos::null);

private:
  //! The original graph to be transposed.
  Teuchos::RCP<const crs_graph_type> origGraph_;

  //! Label for timers
  std::string label_;
};

} // namespace Tpetra

#endif /* TPETRA_CRSGRAPHTRANSPOSER_DECL_HPP */
