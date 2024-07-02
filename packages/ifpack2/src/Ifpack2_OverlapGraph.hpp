// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_OVERLAPGRAPH_HPP
#define IFPACK2_OVERLAPGRAPH_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Import.hpp"
#include "Teuchos_RCP.hpp"
#include "Ifpack2_CreateOverlapGraph.hpp"

namespace Teuchos {
  class ParameterList;
}

namespace Ifpack2 {

/// \class OverlapGraph
/// \brief Construct an overlapped graph from a given nonoverlapping graph.
///
/// FIXME (mfh 26 Sep 2013) It would be more appropriate for this
/// class to have a single template parameter, namely a specialization
/// of Tpetra::CrsGraph or Tpetra::RowGraph.  That would avoid issues
/// with the optional fourth template parameter of Tpetra::CrsGraph.
///
/// FIXME (mfh 26 Sep 2013) I seem to have found this class in a
/// half-implemented state.  I did not add this class and I make no
/// promises that it works.  I heartily encourage developers to make
/// use of the Ifpack2::Details namespace for classes that they are
/// not yet ready to make public.

template<class LocalOrdinal = typename Tpetra::CrsGraph<>::local_ordinal_type,
         class GlobalOrdinal = typename Tpetra::CrsGraph<LocalOrdinal>::global_ordinal_type,
         class Node = typename Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal>::node_type>
class OverlapGraph : public Teuchos::Describable {
public:
  //! The Tpetra::CrsGraph specialization that this class uses.
  typedef Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> graph_type;

  /// \brief Constructor that takes a graph and the level of overlap.
  ///
  /// \param UserMatrixGraph_in [in] The input graph.  We assume that
  ///   its row Map is nonoverlapping.
  ///
  /// \param OverlapLevel_in [in] The level of overlap; zero means none.
  OverlapGraph (const Teuchos::RCP<const graph_type>& UserMatrixGraph_in,
                int OverlapLevel_in);

  //! Copy constructor.
  OverlapGraph (const OverlapGraph<LocalOrdinal,GlobalOrdinal,Node>& Source);

  //! Destructor (virtual for memory safety of derived classes).
  virtual ~OverlapGraph () {}

  //! Return the overlap graph.
  const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>&
  getOverlapGraph () const { return *OverlapGraph_; }

  //! Return the overlap graph's row Map.
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>&
  getOverlapRowMap () const {return *OverlapRowMap_; }

  //! Return the Import object.
  const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node>&
  getOverlapImporter () const { return *OverlapImporter_; }

  /// \brief Return the level of overlap used to create this graph.
  ///
  /// The graph created by this class uses a recursive definition of
  /// overlap.  Level one overlap is created by copying all
  /// off-process rows that are reached to be at least one column of
  /// the rows that are on processor.  Level two overlap is the same
  /// process used on the level one graph, and so on.
  int OverlapLevel () const { return OverlapLevel_; }
  //@}

protected:
  Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > OverlapGraph_;
  Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > UserMatrixGraph_;
  Teuchos::RCP<Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > OverlapRowMap_;
  Teuchos::RCP<Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > OverlapImporter_;
  int OverlapLevel_;
  bool IsOverlapped_;
};

template<class LocalOrdinal, class GlobalOrdinal, class Node>
OverlapGraph<LocalOrdinal,GlobalOrdinal,Node>::
OverlapGraph (const Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >& UserMatrixGraph_in,
              int OverlapLevel_in)
  : UserMatrixGraph_ (UserMatrixGraph_in),
    OverlapLevel_ (OverlapLevel_in),
    IsOverlapped_ (OverlapLevel_in > 0 && UserMatrixGraph_in->getDomainMap ()->isDistributed ())
{
  OverlapGraph_ = createOverlapGraph (UserMatrixGraph_, OverlapLevel_);
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
OverlapGraph<LocalOrdinal,GlobalOrdinal,Node>::
OverlapGraph (const OverlapGraph<LocalOrdinal,GlobalOrdinal,Node>& Source)
  : UserMatrixGraph_ (Source.UserMatrixGraph_),
    OverlapRowMap_ (Source.OverlapRowMap_),
    OverlapLevel_ (Source.OverlapLevel_),
    IsOverlapped_ (Source.IsOverlapped_)
{
  using Teuchos::rcp;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> map_type;

  if (IsOverlapped_) {
    if (! OverlapGraph_.is_null ()) {
      OverlapGraph_ = rcp (new graph_type (*OverlapGraph_));
    }
    if (! OverlapRowMap_.is_null ()) {
      OverlapRowMap_ = rcp (new map_type (*OverlapRowMap_));
    }
  }
}

}//namespace Ifpack2

#endif // IFPACK2_OVERLAPGRAPH_HPP
