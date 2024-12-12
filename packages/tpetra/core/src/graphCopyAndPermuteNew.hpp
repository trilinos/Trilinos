template <class LocalOrdinal, class GlobalOrdinal, class Node>
void
CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
copyAndPermuteNew(
                  const row_graph_type& srcRowGraph,
                  row_graph_type& tgtRowGraph,
                  const size_t numSameIDs,
                  const Kokkos::DualView<const local_ordinal_type*,
                  buffer_device_type>& permuteToLIDs,
                  const Kokkos::DualView<const local_ordinal_type*,
                  buffer_device_type>& permuteFromLIDs,
                  const CombineMode CM)
{
  using std::endl;
  using LO = local_ordinal_type;
  using GO = global_ordinal_type;
  using this_CRS_type = CrsGraph<LO, GO, node_type>;
  const char tfecfFuncName[] = "copyAndPermuteNew: ";
  const bool verbose = verbose_;

  Details::ProfilingRegion regionCAP("Tpetra::CrsGraph::copyAndPermuteNew");
  INCL_EXP(double capTime = Teuchos::Time::wallTime());
  std::unique_ptr<std::string> prefix;
  if (verbose) {
    prefix = this->createPrefix("CrsGraph", "copyAndPermuteNew");
    std::ostringstream os;
    os << *prefix << endl;
    std::cerr << os.str ();
  }

  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
    (permuteToLIDs.extent (0) != permuteFromLIDs.extent (0),
     std::runtime_error, "permuteToLIDs.extent(0) = "
     << permuteToLIDs.extent (0) << " != permuteFromLIDs.extent(0) = "
     << permuteFromLIDs.extent (0) << ".");

  if (verbose) {
    std::ostringstream os;
    os << *prefix << "Compute padding" << endl;
    std::cerr << os.str ();
  }

  // std::cout << "here 0" << std::endl;

  using crs_graph_type = CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
  const crs_graph_type *srcCrsGraphPtr = dynamic_cast<const crs_graph_type *>(&srcRowGraph);
  if (!srcCrsGraphPtr) {
    std::cout << "srk error srcGraph type= " << typeid(srcRowGraph).name() << std::endl;
    std::terminate();
  }
  const crs_graph_type& srcCrsGraph = *srcCrsGraphPtr;

  crs_graph_type *tgtCrsGraphPtr = dynamic_cast<crs_graph_type *>(&tgtRowGraph);
  if (!tgtCrsGraphPtr) {
    std::cout << "srk error tgtGraph type= " << typeid(tgtRowGraph).name() << std::endl;
    std::terminate();
  }
  crs_graph_type& tgtCrsGraph = *tgtCrsGraphPtr;
  // std::cout << "here 1" << std::endl;

  INCL_EXP(double padTime = Teuchos::Time::wallTime());
  auto padding = tgtCrsGraph.computeCrsPadding(srcRowGraph, numSameIDs,
                                               permuteToLIDs, permuteFromLIDs, verbose);

  INCL_EXP(if (IN_EVAL_J) Timers["capsg_G_pad"].first += -padTime + Teuchos::Time::wallTime());
  INCL_EXP(double apadTime = Teuchos::Time::wallTime());
  tgtCrsGraph.applyCrsPadding(*padding, verbose);
  INCL_EXP(if (IN_EVAL_J) Timers["capsg_G_apad"].first += -apadTime + Teuchos::Time::wallTime());

  const map_type& srcRowMap = *(srcRowGraph.getRowMap());
  const map_type& tgtRowMap = *(tgtRowGraph.getRowMap());
  const bool src_filled = srcRowGraph.isFillComplete();
  nonconst_global_inds_host_view_type row_copy;
  LO myid = 0;

  // std::cout << "here 2" << std::endl;

  //
  // "Copy" part of "copy and permute."
  //
  if (src_filled || srcCrsGraphPtr == nullptr) {
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "src_filled || srcCrsGraph == nullptr" << endl;
      std::cerr << os.str ();
    }
    // If the source graph is fill complete, we can't use view mode,
    // because the data might be stored in a different format not
    // compatible with the expectations of view mode.  Also, if the
    // source graph is not a CrsGraph, we can't use view mode,
    // because RowGraph only provides copy mode access to the data.
    INCL_EXP(double time = Teuchos::Time::wallTime());
#if 0
    // std::cout << "here 3" << std::endl;

    for (size_t i = 0; i < numSameIDs; ++i, ++myid) {
      const GO gid = srcRowMap.getGlobalElement (myid);
      size_t row_length = srcRowGraph.getNumEntriesInGlobalRow (gid);
      Kokkos::resize(row_copy,row_length);
      size_t check_row_length = 0;
      srcRowGraph.getGlobalRowCopy (gid, row_copy, check_row_length);
      tgtCrsGraph.insertGlobalIndices (gid, row_length, row_copy.data());
    }
    // std::cout << "here 4" << std::endl;

#else


#include "inner.hpp"
    
#endif
    INCL_EXP(if (IN_EVAL_J) Timers["capsg_G_1"].first += -time + Teuchos::Time::wallTime());
  } else {
    INCL_EXP(double time = Teuchos::Time::wallTime());
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "! src_filled && srcCrsGraph != nullptr" << endl;
      std::cerr << os.str ();
    }
    for (size_t i = 0; i < numSameIDs; ++i, ++myid) {
      const GO gid = srcRowMap.getGlobalElement (myid);
      global_inds_host_view_type row;
      srcCrsGraph.getGlobalRowView (gid, row);
      tgtCrsGraph.insertGlobalIndices (gid, row.extent(0), row.data());
    }
    INCL_EXP(if (IN_EVAL_J) Timers["capsg_G_2"].first += -time + Teuchos::Time::wallTime());
  }

  //
  // "Permute" part of "copy and permute."
  //
  auto permuteToLIDs_h = permuteToLIDs.view_host ();
  auto permuteFromLIDs_h = permuteFromLIDs.view_host ();

  if (src_filled || srcCrsGraphPtr == nullptr) {
    INCL_EXP(double time = Teuchos::Time::wallTime());
    for (LO i = 0; i < static_cast<LO> (permuteToLIDs_h.extent (0)); ++i) {
      const GO mygid = tgtRowMap.getGlobalElement (permuteToLIDs_h[i]);
      const GO srcgid = srcRowMap.getGlobalElement (permuteFromLIDs_h[i]);
      size_t row_length = srcRowGraph.getNumEntriesInGlobalRow (srcgid);
      Kokkos::resize(row_copy,row_length);
      size_t check_row_length = 0;
      srcRowGraph.getGlobalRowCopy (srcgid, row_copy, check_row_length);
      tgtCrsGraph.insertGlobalIndices (mygid, row_length, row_copy.data());
    }
    INCL_EXP(if (IN_EVAL_J) Timers["capsg_G_3"].first += -time + Teuchos::Time::wallTime());
  } else {
    INCL_EXP(double time = Teuchos::Time::wallTime());
    for (LO i = 0; i < static_cast<LO> (permuteToLIDs_h.extent (0)); ++i) {
      const GO mygid = tgtRowMap.getGlobalElement (permuteToLIDs_h[i]);
      const GO srcgid = srcRowMap.getGlobalElement (permuteFromLIDs_h[i]);
      global_inds_host_view_type row;
      srcCrsGraph.getGlobalRowView (srcgid, row);
      tgtCrsGraph.insertGlobalIndices (mygid, row.extent(0), row.data());
    }
    INCL_EXP(if (IN_EVAL_J) Timers["capsg_G_4"].first += -time + Teuchos::Time::wallTime());
  }

  if (verbose) {
    std::ostringstream os;
    os << *prefix << "Done" << endl;
    std::cerr << os.str ();
  }
  INCL_EXP(if (IN_EVAL_J) Timers["capsg_G"].first += -capTime + Teuchos::Time::wallTime());
}

