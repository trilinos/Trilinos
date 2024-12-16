template <class LocalOrdinal, class GlobalOrdinal, class Node>
void
CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
insertGlobalIndicesDevice(const CrsGraph<LocalOrdinal,GlobalOrdinal,Node>& srcCrsGraph,
                          CrsGraph<LocalOrdinal,GlobalOrdinal,Node>& tgtCrsGraph,
                          const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& permuteToLIDs,
                          const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& permuteFromLIDs,
                          LocalOrdinal loopEnd)
{
  using crs_graph_type = CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
  using LO = LocalOrdinal;
  using GO = GlobalOrdinal;
  typedef typename crs_graph_type::global_inds_device_view_type::non_const_value_type global_inds_device_value_t;
  typedef typename crs_graph_type::local_graph_device_type k_local_graph_device_type;
  typedef typename Node::execution_space exec_space;
  typedef Kokkos::RangePolicy<exec_space, LO> range_type;

  const LocalOrdinal LINV = Teuchos::OrdinalTraits<LocalOrdinal>::invalid ();
  const GlobalOrdinal GINV = Teuchos::OrdinalTraits<GlobalOrdinal>::invalid ();

  const k_local_graph_device_type & srcGraphDevice = srcCrsGraph.getLocalGraphDevice();
  const k_local_graph_device_type & tgtGraphDevice = tgtCrsGraph.getLocalGraphDevice();

  using local_map_type = typename crs_graph_type::map_type::local_map_type;
  local_map_type srcRowMapLocal = srcCrsGraph.getRowMap()->getLocalMap();
  local_map_type srcColMapLocal = srcCrsGraph.getColMap()->getLocalMap();
  local_map_type tgtRowMapLocal = tgtCrsGraph.getRowMap()->getLocalMap();

  auto tgtLocalRowPtrsDevice = tgtCrsGraph.getRowPtrsUnpackedDevice();
  auto tgtGlobalColInds = tgtCrsGraph.gblInds_wdv.getDeviceView(Access::ReadWrite);
  auto srcLocalRowPtrsDevice = srcCrsGraph.getLocalRowPtrsDevice();
  auto srcLocalColIndsDevice = srcCrsGraph.lclIndsUnpacked_wdv.getDeviceView(Access::ReadOnly);

  typename crs_graph_type::num_row_entries_type::non_const_type h_numRowEnt = tgtCrsGraph.k_numRowEntries_;

  auto k_numRowEnt = Kokkos::create_mirror_view_and_copy (device_type (), h_numRowEnt);

  const bool sorted = false;

  bool hasMap = permuteFromLIDs.extent(0) > 0;
  auto permuteToLIDs_d = permuteToLIDs.view_device ();
  auto permuteFromLIDs_d = permuteFromLIDs.view_device ();

#ifdef PANZER_DO_CHECK_INNER_HPP
#undef PANZER_DO_CHECK_INNER_HPP
#endif
#define PANZER_DO_CHECK_INNER_HPP 0
#if PANZER_DO_CHECK_INNER_HPP
#define CHECK(a,i) do {                                                 \
    if ((int)(i) >= (int)a.extent(0)) {                                 \
      printf("ERROR: i= %d a= %s e= %d", (int)(i), #a, (int)a.extent(0)); \
      Kokkos::abort("inding error");                                               \
    } } while(0)
#else
#define CHECK(a,i) do { } while(0)
#endif

#ifdef PANZER_INNER_ABORT
#undef PANZER_INNER_ABORT
#endif

#define PANZER_INNER_ABORT(lin) do {            \
    printf("ERROR: line= %d", lin);        \
    Kokkos::abort("error");                         \
  } while(0)

  Kokkos::parallel_for("Tpetra_CrsGraph::copyAndPermuteNew2",
                       range_type (0, loopEnd),
                       KOKKOS_LAMBDA(const LO sourceLID)
                       {
                         auto srcLid = sourceLID;
                         auto tgtLid = sourceLID;
                         if (hasMap) {
                           srcLid = permuteFromLIDs_d(srcLid);
                           tgtLid = permuteToLIDs_d(tgtLid);
                         }
                         auto srcGid = srcRowMapLocal.getGlobalElement(srcLid);
                         if (srcGid == GINV) PANZER_INNER_ABORT(__LINE__);
                         auto tgtGid = tgtRowMapLocal.getGlobalElement(tgtLid);

                         auto tgtLocalRow = tgtRowMapLocal.getLocalElement(tgtGid);
                         if (tgtLocalRow == LINV) PANZER_INNER_ABORT(__LINE__);
                         if (tgtLocalRow != tgtLid) PANZER_INNER_ABORT(__LINE__);
                         CHECK(k_numRowEnt, tgtLocalRow);
                         auto tgtNumEntries = k_numRowEnt(tgtLocalRow);

                         // FIXME no auto use
                         CHECK(srcLocalRowPtrsDevice, srcLid);
                         auto start     = srcLocalRowPtrsDevice(srcLid);
                         CHECK(srcLocalRowPtrsDevice, srcLid + 1);
                         auto end       = srcLocalRowPtrsDevice(srcLid + 1);
                         auto rowLength = (end - start);

                         //KOKKOS_ASSERT(rowLength <= max_row_entries);

                         CHECK(tgtLocalRowPtrsDevice, tgtLocalRow);
                         auto tstart      = tgtLocalRowPtrsDevice(tgtLocalRow);
                         auto tend        = tstart + tgtNumEntries;
                         CHECK(tgtLocalRowPtrsDevice, tgtLocalRow + 1);
                         auto tend1       = tgtLocalRowPtrsDevice(tgtLocalRow + 1);

                         const size_t num_avail = (tend1 < tend) ? size_t (0) : tend1 - tend;
                         size_t num_inserted = 0;

                         CHECK(tgtGlobalColInds, tstart);
                         global_inds_device_value_t *tgtGlobalColIndsPtr = tgtGlobalColInds.data();

                         size_t hint=0;
                         for (size_t j = 0; j < rowLength; j++) {
                           CHECK(srcLocalColIndsDevice, start + j);
                           auto ci = srcLocalColIndsDevice(start + j);
                           GO gi = srcColMapLocal.getGlobalElement(ci);
                           if (gi == GINV) PANZER_INNER_ABORT(__LINE__);
                           auto numInTgtRow = (tend - tstart);

                           const size_t offset =
                             KokkosSparse::findRelOffset (tgtGlobalColIndsPtr+tstart,
                                                          numInTgtRow,
                                                          gi, hint, sorted);

                           if (offset == numInTgtRow) {
                             if (num_inserted >= num_avail) { // not enough room
                               //return Teuchos::OrdinalTraits<size_t>::invalid();
                               Kokkos::abort("num_avail");
                             }
                             //Kokkos::atomic_store (&tgtRowVals[offset], newVals);
                             CHECK(tgtGlobalColInds, tstart + offset);
                             tgtGlobalColIndsPtr[tstart + offset] = gi;
                             ++tend;
                             hint = offset + 1;
                             ++num_inserted;
                           }
                         }
                         CHECK(k_numRowEnt, tgtLocalRow);
                         k_numRowEnt(tgtLocalRow) += num_inserted;

                         return size_t(0);
                       });

  Kokkos::fence("here 10");
  Kokkos::deep_copy(tgtCrsGraph.k_numRowEntries_, k_numRowEnt);
  tgtCrsGraph.setLocallyModified();
}


template <class LocalOrdinal, class GlobalOrdinal, class Node>
void
CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
copyAndPermuteNew(const row_graph_type& srcRowGraph,
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

  using crs_graph_type = CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
  const crs_graph_type *srcCrsGraphPtr = dynamic_cast<const crs_graph_type *>(&srcRowGraph);
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(!srcCrsGraphPtr, std::runtime_error,
                                        "error srcGraph type= " << typeid(srcRowGraph).name());
  const crs_graph_type& srcCrsGraph = *srcCrsGraphPtr;

  crs_graph_type *tgtCrsGraphPtr = dynamic_cast<crs_graph_type *>(&tgtRowGraph);
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(!srcCrsGraphPtr, std::runtime_error,
                                        "error tgtGraph type= " << typeid(tgtRowGraph).name());

  crs_graph_type& tgtCrsGraph = *tgtCrsGraphPtr;

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

  //
  // "Copy" part of "copy and permute."
  //
  LO numSameIDs_as_LID = static_cast<LO>(numSameIDs);
  using LidMapType = std::function<LocalOrdinal(const LocalOrdinal lid)> ;

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

    Kokkos::DualView<const local_ordinal_type*, buffer_device_type> noPermute;
    insertGlobalIndicesDevice(srcCrsGraph, tgtCrsGraph,
                              noPermute, noPermute,
                              numSameIDs_as_LID);
    
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
  auto permuteToLIDs_d = permuteToLIDs.view_device ();
  auto permuteFromLIDs_d = permuteFromLIDs.view_device ();

  if (src_filled || srcCrsGraphPtr == nullptr) {
    INCL_EXP(double time = Teuchos::Time::wallTime());
#if 0
    for (LO i = 0; i < static_cast<LO> (permuteToLIDs_h.extent (0)); ++i) {
      const GO mygid = tgtRowMap.getGlobalElement (permuteToLIDs_h[i]);
      const GO srcgid = srcRowMap.getGlobalElement (permuteFromLIDs_h[i]);
      size_t row_length = srcRowGraph.getNumEntriesInGlobalRow (srcgid);
      Kokkos::resize(row_copy,row_length);
      size_t check_row_length = 0;
      srcRowGraph.getGlobalRowCopy (srcgid, row_copy, check_row_length);
      tgtCrsGraph.insertGlobalIndices (mygid, row_length, row_copy.data());
    }
#else
    insertGlobalIndicesDevice(srcCrsGraph, tgtCrsGraph,
                              permuteToLIDs, permuteFromLIDs,  // note reversed arg order, tgt, then src 
                              static_cast<LO> (permuteToLIDs_h.extent (0)));
#endif
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

