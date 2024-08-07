// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Ifpack2_IlukGraph.hpp
/// \brief Declaration and definition of IlukGraph
///
/// \warning This header file is an implementation detail of RILUK.
///   Its contents may change or it may go away at any time.

#ifndef IFPACK2_ILUK_GRAPH_HPP
#define IFPACK2_ILUK_GRAPH_HPP

#include <algorithm>
#include <vector>

#include "KokkosSparse_spiluk.hpp"

#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_Details_WrappedDualView.hpp>
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
template<class GraphType, class KKHandleType>
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



  typedef typename crs_graph_type::nonconst_global_inds_host_view_type nonconst_global_inds_host_view_type;
  typedef typename crs_graph_type::nonconst_local_inds_host_view_type nonconst_local_inds_host_view_type;
  typedef typename crs_graph_type::global_inds_host_view_type global_inds_host_view_type;
  typedef typename crs_graph_type::local_inds_host_view_type  local_inds_host_view_type;

  /// \brief Constructor.
  ///
  /// Create a IlukGraph object using the input graph and specified
  /// level of fill.
  ///
  /// \param G [in] An existing graph.
  /// \param levelFill [in] The level of fill to compute; the k of ILU(k).
  /// \param levelOverlap [in] The level of overlap between subdomains.
  /// \param overalloc [in] The estimated number of nonzeros per row in the resulting matrices is (maxNodeNumRowEntries of the input * overalloc^levelFill).  Must be greater than 1.  Smaller values are more conservative with memory, but may require recomputation if the estimate is too low.  Default value is two.
  ///
  /// \note Actual construction occurs in initialize().
  IlukGraph (const Teuchos::RCP<const GraphType>& G,
             const int levelFill,
             const int levelOverlap,
             const double overalloc = 2.);

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
  void initialize();
  void initialize(const Teuchos::RCP<KKHandleType>& KernelHandle);

  //! The level of fill used to construct this graph.
  int getLevelFill () const { return LevelFill_; }

  //! The level of overlap used to construct this graph.
  int getLevelOverlap () const { return LevelOverlap_; }

  //! Returns the original graph given
  Teuchos::RCP<const GraphType> getA_Graph () const {
    return Graph_;
  }

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
  IlukGraph (const IlukGraph<GraphType, KKHandleType>&);

  /// \brief Assignment operator (UNIMPLEMENTED; DO NOT USE).
  ///
  /// This assignment operator is declared private and unimplemented,
  /// in order to forbid its use syntactically.  If you decide that
  /// you need to implement this method, it should do deep copies of
  /// all internal graphs.  It may do a shallow copy of the input
  /// graph (which is not modified).  Also, you must implement the
  /// copy constructor in a similar way.
  IlukGraph& operator= (const IlukGraph<GraphType, KKHandleType>&);

  //! Construct the overlap graph.
  void constructOverlapGraph();

  Teuchos::RCP<const GraphType> Graph_;
  Teuchos::RCP<const crs_graph_type> OverlapGraph_;
  int LevelFill_;
  int LevelOverlap_;
  const double Overalloc_;
  Teuchos::RCP<crs_graph_type> L_Graph_;
  Teuchos::RCP<crs_graph_type> U_Graph_;
  size_t NumMyDiagonals_;
  size_t NumGlobalDiagonals_;
};


template<class GraphType, class KKHandleType>
IlukGraph<GraphType, KKHandleType>::
IlukGraph (const Teuchos::RCP<const GraphType>& G,
           const int levelFill,
           const int levelOverlap,
           const double overalloc)
  : Graph_ (G),
    LevelFill_ (levelFill),
    LevelOverlap_ (levelOverlap),
    Overalloc_ (overalloc),
    NumMyDiagonals_ (0),
    NumGlobalDiagonals_ (0)
{
  TEUCHOS_TEST_FOR_EXCEPTION(Overalloc_ <= 1., std::runtime_error,
    "Ifpack2::IlukGraph: FATAL: overalloc must be greater than 1.")
}


template<class GraphType, class KKHandleType>
IlukGraph<GraphType, KKHandleType>::~IlukGraph()
{}


template<class GraphType, class KKHandleType>
void IlukGraph<GraphType, KKHandleType>::
setParameters (const Teuchos::ParameterList& parameterlist)
{
  getParameter (parameterlist, "iluk level-of-fill", LevelFill_);
  getParameter (parameterlist, "iluk level-of-overlap", LevelOverlap_);
}


template<class GraphType, class KKHandleType>
void IlukGraph<GraphType, KKHandleType>::constructOverlapGraph () {
  // FIXME (mfh 22 Dec 2013) This won't do if we want
  // RILUK::initialize() to do the right thing (that is,
  // unconditionally recompute the "symbolic factorization").
  if (OverlapGraph_ == Teuchos::null) {
    OverlapGraph_ = createOverlapGraph<GraphType> (Graph_, LevelOverlap_);
  }
}


template<class GraphType, class KKHandleType>
void IlukGraph<GraphType, KKHandleType>::initialize()
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
  const int MaxNumIndices = OverlapGraph_->getLocalMaxNumRowEntries ();

  // FIXME (mfh 23 Dec 2013) Use size_t or whatever
  // getLocalNumElements() returns, instead of ptrdiff_t.
  const int NumMyRows = OverlapGraph_->getRowMap ()->getLocalNumElements ();

  using device_type = typename node_type::device_type;
  using execution_space = typename device_type::execution_space;
  using dual_view_type = Kokkos::DualView<size_t*,device_type>;
  dual_view_type numEntPerRow_dv("numEntPerRow",NumMyRows);
  Tpetra::Details::WrappedDualView<dual_view_type> numEntPerRow(numEntPerRow_dv);

  const auto overalloc = Overalloc_;
  const auto levelfill = LevelFill_;
  {
    // Scoping for the  localOverlapGraph access
    auto numEntPerRow_d = numEntPerRow.getDeviceView(Tpetra::Access::OverwriteAll);
    auto localOverlapGraph = OverlapGraph_->getLocalGraphDevice();
    Kokkos::parallel_for("CountOverlapGraphRowEntries",
                         Kokkos::RangePolicy<execution_space>(0, NumMyRows),
                         KOKKOS_LAMBDA(const int i)
                         {
                           // Heuristic to get the maximum number of entries per row.
                           int RowMaxNumIndices = localOverlapGraph.rowConst(i).length;
                           numEntPerRow_d(i) = (levelfill == 0) ? RowMaxNumIndices  // No additional storage needed
                             : Kokkos::ceil(static_cast<double>(RowMaxNumIndices)
                                    * Kokkos::pow(overalloc, levelfill));
                         });

  };

  bool insertError;  // No error found yet while inserting entries
  do {
    insertError = false;
    Teuchos::ArrayView<const size_t> a_numEntPerRow(numEntPerRow.getHostView(Tpetra::Access::ReadOnly).data(),NumMyRows);
    L_Graph_ = rcp (new crs_graph_type (OverlapGraph_->getRowMap (),
                                        OverlapGraph_->getRowMap (),
                                        a_numEntPerRow));
    U_Graph_ = rcp (new crs_graph_type (OverlapGraph_->getRowMap (),
                                        OverlapGraph_->getRowMap (),
                                        a_numEntPerRow));

    Array<local_ordinal_type> L (MaxNumIndices);
    Array<local_ordinal_type> U (MaxNumIndices);

    // First we copy the user's graph into L and U, regardless of fill level

    NumMyDiagonals_ = 0;

    for (int i = 0; i< NumMyRows; ++i) {
      local_inds_host_view_type my_indices;
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
        L_Graph_->insertLocalIndices (i, NumL, L.data());
      }
      if (NumU) {
        U_Graph_->insertLocalIndices (i, NumU, U.data());
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

      try {
        for (int i = 0; i < NumMyRows; ++i) {
          int First, Next;

          // copy column indices of row into workspace and sort them

          size_t LenL = L_Graph_->getNumEntriesInLocalRow(i);
          size_t LenU = U_Graph_->getNumEntriesInLocalRow(i);
          size_t Len = LenL + LenU + 1;
          CurrentRow.resize(Len);
          nonconst_local_inds_host_view_type CurrentRow_view(CurrentRow.data(),CurrentRow.size());
          L_Graph_->getLocalRowCopy(i, CurrentRow_view, LenL);  // Get L Indices
          CurrentRow[LenL] = i;                              // Put in Diagonal
          if (LenU > 0) {
            ArrayView<local_ordinal_type> URowView = CurrentRow.view (LenL+1,LenU);
            nonconst_local_inds_host_view_type URowView_v(URowView.data(),URowView.size());

            // Get U Indices
            U_Graph_->getLocalRowCopy (i, URowView_v, LenU);
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
            local_inds_host_view_type IndicesU;
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
                CurrentLevel[CurInList] = std::min (CurrentLevel[CurInList],
                                                    NewLevel);
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
          CurrentRow.resize(0);

          Next = First;

          // Lower
          while (Next < i) {
            CurrentRow.push_back(Next);
            Next = LinkList[Next];
          }

          // FIXME (mfh 23 Dec 2013) It's not clear to me that
          // removeLocalIndices works like people expect it to work.  In
          // particular, it does not actually change the column Map.
          L_Graph_->removeLocalIndices (i); // Delete current set of Indices
          if (CurrentRow.size() > 0) {
            L_Graph_->insertLocalIndices (i, CurrentRow.size(),CurrentRow.data());
          }

          // Diagonal

          TEUCHOS_TEST_FOR_EXCEPTION(
            Next != i, std::runtime_error,
            "Ifpack2::IlukGraph::initialize: FATAL: U has zero diagonal")

          LevelsRowU[0] = CurrentLevel[Next];
          Next = LinkList[Next];

          // Upper
          CurrentRow.resize(0);
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
            U_Graph_->insertLocalIndices (i, CurrentRow.size(),CurrentRow.data());
          }

          // Allocate and fill Level info for this row
          Levels[i] = std::vector<int> (LenU+1);
          for (size_t jj=0; jj<LenU+1; jj++) {
            Levels[i][jj] = LevelsRowU[jj];
          }
        }
      }
      catch (std::runtime_error &e) {
        insertError = true;
        auto numEntPerRow_d = numEntPerRow.getDeviceView(Tpetra::Access::OverwriteAll);
        Kokkos::parallel_for("CountOverlapGraphRowEntries",
          Kokkos::RangePolicy<execution_space>(0, NumMyRows),
          KOKKOS_LAMBDA(const int i)
          {
            const auto numRowEnt = numEntPerRow_d(i);
            numEntPerRow_d(i) = ceil(static_cast<double>((numRowEnt != 0 ? numRowEnt : 1)) * overalloc);
          });
      }
      const int localInsertError = insertError ? 1 : 0;
      int globalInsertError = 0;
      reduceAll (* (OverlapGraph_->getRowMap ()->getComm ()), REDUCE_SUM, 1,
                 &localInsertError, &globalInsertError);
      insertError = globalInsertError > 0;
    }
  } while (insertError);  // do until all insertions complete successfully

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


template<class GraphType, class KKHandleType>
void IlukGraph<GraphType, KKHandleType>::initialize(const Teuchos::RCP<KKHandleType>& KernelHandle)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_SUM;
  using Teuchos::reduceAll;

  typedef typename crs_graph_type::local_graph_device_type local_graph_device_type;
  typedef typename local_graph_device_type::size_type      size_type;
  typedef typename local_graph_device_type::data_type      data_type;
  typedef typename local_graph_device_type::array_layout   array_layout;
  typedef typename local_graph_device_type::device_type    device_type;

  typedef typename Kokkos::View<size_type*, array_layout, device_type> lno_row_view_t;
  typedef typename Kokkos::View<data_type*, array_layout, device_type> lno_nonzero_view_t;

  constructOverlapGraph();

  // FIXME (mfh 23 Dec 2013) Use size_t or whatever
  // getLocalNumElements() returns, instead of ptrdiff_t.
  const int NumMyRows = OverlapGraph_->getRowMap()->getLocalNumElements();
  auto localOverlapGraph = OverlapGraph_->getLocalGraphDevice();

  if (KernelHandle->get_spiluk_handle()->get_nrows() < static_cast<size_type>(NumMyRows)) {
    KernelHandle->get_spiluk_handle()->reset_handle(NumMyRows,
                                                    KernelHandle->get_spiluk_handle()->get_nnzL(),
                                                    KernelHandle->get_spiluk_handle()->get_nnzU());
  }

  lno_row_view_t     L_row_map("L_row_map", NumMyRows + 1);
  lno_nonzero_view_t L_entries("L_entries", KernelHandle->get_spiluk_handle()->get_nnzL());
  lno_row_view_t     U_row_map("U_row_map", NumMyRows + 1);
  lno_nonzero_view_t U_entries("U_entries", KernelHandle->get_spiluk_handle()->get_nnzU());

  bool symbolicError;
  do {
    symbolicError = false;
    try {
      KokkosSparse::Experimental::spiluk_symbolic( KernelHandle.getRawPtr(), LevelFill_,
                                                   localOverlapGraph.row_map, localOverlapGraph.entries,
                                                   L_row_map, L_entries, U_row_map, U_entries );
    }
    catch (std::runtime_error &e) {
      symbolicError = true;
      data_type nnzL = static_cast<data_type>(Overalloc_)*L_entries.extent(0);
      data_type nnzU = static_cast<data_type>(Overalloc_)*U_entries.extent(0);
      KernelHandle->get_spiluk_handle()->reset_handle(NumMyRows, nnzL, nnzU);
      Kokkos::resize(L_entries, KernelHandle->get_spiluk_handle()->get_nnzL());
      Kokkos::resize(U_entries, KernelHandle->get_spiluk_handle()->get_nnzU());
    }
    const int localSymbolicError = symbolicError ? 1 : 0;
    int globalSymbolicError = 0;
    reduceAll (* (OverlapGraph_->getRowMap ()->getComm ()), REDUCE_SUM, 1,
                   &localSymbolicError, &globalSymbolicError);
    symbolicError = globalSymbolicError > 0;
  } while (symbolicError);

  Kokkos::resize(L_entries, KernelHandle->get_spiluk_handle()->get_nnzL());
  Kokkos::resize(U_entries, KernelHandle->get_spiluk_handle()->get_nnzU());

  RCP<Teuchos::ParameterList> params = Teuchos::parameterList ();
  params->set ("Optimize Storage",false);

  L_Graph_ = rcp (new crs_graph_type (OverlapGraph_->getRowMap (),
                                      OverlapGraph_->getRowMap (),
                                      L_row_map, L_entries));
  U_Graph_ = rcp (new crs_graph_type (OverlapGraph_->getRowMap (),
                                      OverlapGraph_->getRowMap (),
                                      U_row_map, U_entries));

  RCP<const map_type> L_DomainMap = OverlapGraph_->getRowMap ();
  RCP<const map_type> L_RangeMap  = Graph_->getRangeMap ();
  RCP<const map_type> U_DomainMap = Graph_->getDomainMap ();
  RCP<const map_type> U_RangeMap  = OverlapGraph_->getRowMap ();
  L_Graph_->fillComplete (L_DomainMap, L_RangeMap, params);
  U_Graph_->fillComplete (U_DomainMap, U_RangeMap, params);
}

} // namespace Ifpack2

#endif /* IFPACK2_ILUK_GRAPH_HPP */
