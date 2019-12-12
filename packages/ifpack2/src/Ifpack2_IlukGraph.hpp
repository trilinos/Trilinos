/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
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

/// \file Ifpack2_IlukGraph.hpp
/// \brief Declaration and definition of IlukGraph
///
/// \warning This header file is an implementation detail of RILUK.
///   Its contents may change or it may go away at any time.

#ifndef IFPACK2_ILUK_GRAPH_HPP
#define IFPACK2_ILUK_GRAPH_HPP

#include <algorithm>
#include <vector>

#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_Import.hpp>
#include <Ifpack2_CreateOverlapGraph.hpp>
#include <Ifpack2_Parameters.hpp>

namespace Ifpack2 {

/// \class IlukGraph
/// \tparam GraphType A specialization of Tpetra::CrsGraph (tested) or
///   Tpetra::RowGraph (not tested).
/// \brief Construct a level filled graph for use in computing an
///   ILU(k) incomplete factorization.
///
/// \warning This class is an implementation detail of RILUK.  Its
///   interface may change or it may go away at any time.
///
/// Ifpack2::IlukGraph enables the construction of matrix graphs using
/// level-fill algorithms.  The only function required for
/// construction is a getGlobalRowView capability, i.e., the graph
/// that is passed in to the constructor must implement the
/// Tpetra::RowGraph interface defined in Tpetra_RowGraph.hpp.
///
/// \section Ifpack2_IlukGraph_cnstr Constructing IlukGraph objects
///
/// Constructing an IlukGraph is a two step process.  First, call the
/// constructor, passing in a CrsGraph object and an integer
/// indicating the desired level of fill.  Then, call the initialize()
/// method to complete the process.  This naturally matches the
/// three-stage initialization of Preconditioner objects, and in
/// particular of RILUK.
///
/// It is worth noting that an IlukGraph object creates two
/// Tpetra::CrsGraph objects containing L and U, the graphs for the
/// lower and upper triangular parts of the ILU(k) graph.  Thus, it is
/// possible to manually insert and delete graph entries in L and U
/// via the Tpetra::CrsGraph InsertIndices and RemoveIndices
/// functions.  However, in this case FillComplete must be called
/// before the graph is used for subsequent operations.
template<class GraphType>
class IlukGraph : public Teuchos::Describable {
public:
  typedef typename GraphType::local_ordinal_type local_ordinal_type;
  typedef typename GraphType::global_ordinal_type global_ordinal_type;
  typedef typename GraphType::node_type node_type;

  //! Tpetra::RowGraph specialization used by this class.
  typedef Tpetra::RowGraph<local_ordinal_type,
                           global_ordinal_type,
                           node_type> row_graph_type;
  //! Tpetra::CrsGraph specialization used by this class.
  typedef Tpetra::CrsGraph<local_ordinal_type,
                           global_ordinal_type,
                           node_type> crs_graph_type;

  /// \brief Constructor.
  ///
  /// Create a IlukGraph object using the input graph and specified
  /// level of fill.
  ///
  /// \param G [in] An existing graph.
  /// \param levelFill [in] The level of fill to compute; the k of ILU(k).
  /// \param levelOverlap [in] The level of overlap between subdomains.
  ///
  /// \note Actual construction occurs in initialize().
  IlukGraph (const Teuchos::RCP<const GraphType>& G,
             const int levelFill,
             const int levelOverlap);

  //! IlukGraph Destructor
  virtual ~IlukGraph ();

  /// \brief Set parameters.
  ///
  /// This method recognizes two parameter names: Level_fill and
  /// Level_overlap.  Both are case insensitive, and in both cases the
  /// parameter must have type int.
  void setParameters (const Teuchos::ParameterList& parameterlist);

  /// \brief Set up the graph structure of the L and U factors.
  ///
  /// This method is called "initialize" by analogy with
  /// Preconditioner, where initialize() computes the symbolic
  /// (incomplete) factorization, and compute() computes the
  /// corresponding numeric factorization.  IlukGraph is just a graph,
  /// so it can only compute a symbolic factorization (i.e., the graph
  /// structure of the factorization).  Hence, it implements
  /// initialize(), but not compute().  RILUK calls IlukGraph's
  /// initialize() method in its own initialize() method, as one would
  /// expect.
  void initialize ();

  //! The level of fill used to construct this graph.
  int getLevelFill () const { return LevelFill_; }

  //! The level of overlap used to construct this graph.
  int getLevelOverlap () const { return LevelOverlap_; }

  //! Returns the graph of lower triangle of the ILU(k) graph as a Tpetra::CrsGraph.
  Teuchos::RCP<crs_graph_type> getL_Graph () const {
    return L_Graph_;
  }

  //! Returns the graph of upper triangle of the ILU(k) graph as a Tpetra::CrsGraph.
  Teuchos::RCP<crs_graph_type> getU_Graph () const {
    return U_Graph_;
  }

  //! Returns the the overlapped graph.
  Teuchos::RCP<const crs_graph_type> getOverlapGraph () const {
    return OverlapGraph_;
  }

  //! Returns the global number of diagonals in the ILU(k) graph.
  size_t getNumGlobalDiagonals() const { return NumGlobalDiagonals_; }

private:
  typedef typename GraphType::map_type map_type;

  /// \brief Copy constructor (UNIMPLEMENTED; DO NOT USE).
  ///
  /// This copy constructor is declared private and unimplemented, in
  /// order to forbid its use syntactically.  If you decide that you
  /// need to implement this method, it should do deep copies of all
  /// internal graphs.  It may do a shallow copy of the input graph
  /// (which is not modified).  Also, you must implement operator= in
  /// a similar way.
  IlukGraph (const IlukGraph<GraphType>&);

  /// \brief Assignment operator (UNIMPLEMENTED; DO NOT USE).
  ///
  /// This assignment operator is declared private and unimplemented,
  /// in order to forbid its use syntactically.  If you decide that
  /// you need to implement this method, it should do deep copies of
  /// all internal graphs.  It may do a shallow copy of the input
  /// graph (which is not modified).  Also, you must implement the
  /// copy constructor in a similar way.
  IlukGraph& operator= (const IlukGraph<GraphType>&);

  //! Construct the overlap graph.
  void constructOverlapGraph();

  Teuchos::RCP<const GraphType> Graph_;
  Teuchos::RCP<const crs_graph_type> OverlapGraph_;
  int LevelFill_;
  int LevelOverlap_;
  Teuchos::RCP<crs_graph_type> L_Graph_;
  Teuchos::RCP<crs_graph_type> U_Graph_;
  size_t NumMyDiagonals_;
  size_t NumGlobalDiagonals_;
};


template<class GraphType>
IlukGraph<GraphType>::
IlukGraph (const Teuchos::RCP<const GraphType>& G,
           const int levelFill,
           const int levelOverlap)
  : Graph_ (G),
    LevelFill_ (levelFill),
    LevelOverlap_ (levelOverlap),
    NumMyDiagonals_ (0),
    NumGlobalDiagonals_ (0)
{}


template<class GraphType>
IlukGraph<GraphType>::~IlukGraph()
{}


template<class GraphType>
void IlukGraph<GraphType>::
setParameters (const Teuchos::ParameterList& parameterlist)
{
  getParameter (parameterlist, "iluk level-of-fill", LevelFill_);
  getParameter (parameterlist, "iluk level-of-overlap", LevelOverlap_);
}


template<class GraphType>
void IlukGraph<GraphType>::constructOverlapGraph () {
  // FIXME (mfh 22 Dec 2013) This won't do if we want
  // RILUK::initialize() to do the right thing (that is,
  // unconditionally recompute the "symbolic factorization").
  if (OverlapGraph_ == Teuchos::null) {
    OverlapGraph_ = createOverlapGraph<GraphType> (Graph_, LevelOverlap_);
  }
}


template<class GraphType>
void IlukGraph<GraphType>::initialize()
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_SUM;
  using Teuchos::reduceAll;

  size_t NumIn, NumL, NumU;
  bool DiagFound;

  constructOverlapGraph();

  // Get Maximum Row length
  const int MaxNumIndices = OverlapGraph_->getNodeMaxNumRowEntries ();

  // FIXME (mfh 23 Dec 2013) Use size_t or whatever
  // getNodeNumElements() returns, instead of ptrdiff_t.
  const int NumMyRows = OverlapGraph_->getRowMap ()->getNodeNumElements ();

  // Heuristic to get the maximum number of entries per row.
  const int MaxNumEntriesPerRow = (LevelFill_ == 0)
                                ? MaxNumIndices
                                : MaxNumIndices + 5*LevelFill_;
  L_Graph_ = rcp (new crs_graph_type (OverlapGraph_->getRowMap (),
                                      OverlapGraph_->getRowMap (), MaxNumEntriesPerRow));
  U_Graph_ = rcp (new crs_graph_type (OverlapGraph_->getRowMap (),
                                      OverlapGraph_->getRowMap (), MaxNumEntriesPerRow));

  Array<local_ordinal_type> L (MaxNumIndices);
  Array<local_ordinal_type> U (MaxNumIndices);

  // First we copy the user's graph into L and U, regardless of fill level

  NumMyDiagonals_ = 0;

  for (int i = 0; i< NumMyRows; ++i) {
    ArrayView<const local_ordinal_type> my_indices;
    OverlapGraph_->getLocalRowView (i, my_indices);

    // Split into L and U (we don't assume that indices are ordered).

    NumL = 0;
    NumU = 0;
    DiagFound = false;
    NumIn = my_indices.size();

    for (size_t j = 0; j < NumIn; ++j) {
      const local_ordinal_type k = my_indices[j];

      if (k<NumMyRows) { // Ignore column elements that are not in the square matrix

        if (k==i) {
          DiagFound = true;
        }
        else if (k < i) {
          L[NumL] = k;
          NumL++;
        }
        else {
          U[NumU] = k;
          NumU++;
        }
      }
    }

    // Check in things for this row of L and U

    if (DiagFound) {
      ++NumMyDiagonals_;
    }
    if (NumL) {
      ArrayView<const local_ordinal_type> L_view = L.view (0, NumL);
      L_Graph_->insertLocalIndices (i, L_view);
    }
    if (NumU) {
      ArrayView<const local_ordinal_type> U_view = U.view (0, NumU);
      U_Graph_->insertLocalIndices (i, U_view);
    }
  }

  if (LevelFill_ > 0) {
    // Complete Fill steps
    RCP<const map_type> L_DomainMap = OverlapGraph_->getRowMap ();
    RCP<const map_type> L_RangeMap = Graph_->getRangeMap ();
    RCP<const map_type> U_DomainMap = Graph_->getDomainMap ();
    RCP<const map_type> U_RangeMap = OverlapGraph_->getRowMap ();
    RCP<Teuchos::ParameterList> params = Teuchos::parameterList ();
    params->set ("Optimize Storage",false);
    L_Graph_->fillComplete (L_DomainMap, L_RangeMap, params);
    U_Graph_->fillComplete (U_DomainMap, U_RangeMap, params);
    L_Graph_->resumeFill ();
    U_Graph_->resumeFill ();

    // At this point L_Graph and U_Graph are filled with the pattern of input graph,
    // sorted and have redundant indices (if any) removed.  Indices are zero based.
    // LevelFill is greater than zero, so continue...

    int MaxRC = NumMyRows;
    std::vector<std::vector<int> > Levels(MaxRC);
    std::vector<int> LinkList(MaxRC);
    std::vector<int> CurrentLevel(MaxRC);
    Array<local_ordinal_type> CurrentRow (MaxRC + 1);
    std::vector<int> LevelsRowU(MaxRC);

    for (int i = 0; i < NumMyRows; ++i) {
      int First, Next;

      // copy column indices of row into workspace and sort them

      size_t LenL = L_Graph_->getNumEntriesInLocalRow(i);
      size_t LenU = U_Graph_->getNumEntriesInLocalRow(i);
      size_t Len = LenL + LenU + 1;

      CurrentRow.resize(Len);

      L_Graph_->getLocalRowCopy(i, CurrentRow(), LenL);      // Get L Indices
      CurrentRow[LenL] = i;                                     // Put in Diagonal
      if (LenU > 0) {
        ArrayView<local_ordinal_type> URowView = CurrentRow.view (LenL+1, LenU);
        // Get U Indices
        U_Graph_->getLocalRowCopy (i, URowView, LenU);
      }

      // Construct linked list for current row

      for (size_t j=0; j<Len-1; j++) {
        LinkList[CurrentRow[j]] = CurrentRow[j+1];
        CurrentLevel[CurrentRow[j]] = 0;
      }

      LinkList[CurrentRow[Len-1]] = NumMyRows;
      CurrentLevel[CurrentRow[Len-1]] = 0;

      // Merge List with rows in U

      First = CurrentRow[0];
      Next = First;
      while (Next < i) {
        int PrevInList = Next;
        int NextInList = LinkList[Next];
        int RowU = Next;
        // Get Indices for this row of U
        ArrayView<const local_ordinal_type> IndicesU;
        U_Graph_->getLocalRowView (RowU, IndicesU);
        // FIXME (mfh 23 Dec 2013) size() returns ptrdiff_t, not int.
        int LengthRowU = IndicesU.size ();

        int ii;

        // Scan RowU

        for (ii = 0; ii < LengthRowU; /*nop*/) {
          int CurInList = IndicesU[ii];
          if (CurInList < NextInList) {
            // new fill-in
            int NewLevel = CurrentLevel[RowU] + Levels[RowU][ii+1] + 1;
            if (NewLevel <= LevelFill_) {
              LinkList[PrevInList]  = CurInList;
              LinkList[CurInList] = NextInList;
              PrevInList = CurInList;
              CurrentLevel[CurInList] = NewLevel;
            }
            ii++;
          }
          else if (CurInList == NextInList) {
            PrevInList = NextInList;
            NextInList = LinkList[PrevInList];
            int NewLevel = CurrentLevel[RowU] + Levels[RowU][ii+1] + 1;
            CurrentLevel[CurInList] = std::min (CurrentLevel[CurInList], NewLevel);
            ii++;
          }
          else { // (CurInList > NextInList)
            PrevInList = NextInList;
            NextInList = LinkList[PrevInList];
          }
        }
        Next = LinkList[Next];
      }

      // Put pattern into L and U

      CurrentRow.resize (0);

      Next = First;

      // Lower

      while (Next < i) {
        CurrentRow.push_back (Next);
        Next = LinkList[Next];
      }

      // FIXME (mfh 23 Dec 2013) It's not clear to me that
      // removeLocalIndices works like people expect it to work.  In
      // particular, it does not actually change the column Map.
      L_Graph_->removeLocalIndices (i); // Delete current set of Indices
      if (CurrentRow.size() > 0) {
        L_Graph_->insertLocalIndices (i, CurrentRow ());
      }

      // Diagonal

      TEUCHOS_TEST_FOR_EXCEPTION(
        Next != i, std::runtime_error,
        "Ifpack2::IlukGraph::initialize: FATAL: U has zero diagonal")

      LevelsRowU[0] = CurrentLevel[Next];
      Next = LinkList[Next];

      // Upper

      CurrentRow.resize (0);
      LenU = 0;

      while (Next < NumMyRows) {
        LevelsRowU[LenU+1] = CurrentLevel[Next];
        CurrentRow.push_back (Next);
        ++LenU;
        Next = LinkList[Next];
      }

      // FIXME (mfh 23 Dec 2013) It's not clear to me that
      // removeLocalIndices works like people expect it to work.  In
      // particular, it does not actually change the column Map.

      U_Graph_->removeLocalIndices (i); // Delete current set of Indices
      if (LenU > 0) {
        U_Graph_->insertLocalIndices (i, CurrentRow ());
      }

      // Allocate and fill Level info for this row
      Levels[i] = std::vector<int> (LenU+1);
      for (size_t jj=0; jj<LenU+1; jj++) {
        Levels[i][jj] = LevelsRowU[jj];
      }
    }
  }

  // Complete Fill steps
  RCP<const map_type> L_DomainMap = OverlapGraph_->getRowMap ();
  RCP<const map_type> L_RangeMap  = Graph_->getRangeMap ();
  RCP<const map_type> U_DomainMap = Graph_->getDomainMap ();
  RCP<const map_type> U_RangeMap  = OverlapGraph_->getRowMap ();
  L_Graph_->fillComplete (L_DomainMap, L_RangeMap);//DoOptimizeStorage is default here...
  U_Graph_->fillComplete (U_DomainMap, U_RangeMap);//DoOptimizeStorage is default here...

  reduceAll<int, size_t> (* (L_DomainMap->getComm ()), REDUCE_SUM, 1,
                          &NumMyDiagonals_, &NumGlobalDiagonals_);
}

} // namespace Ifpack2

#endif /* IFPACK2_ILUK_GRAPH_HPP */
