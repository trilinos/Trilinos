  // not yet a member function of CrsMatrix
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  copyAndPermuteStaticGraph_new(
                                const RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& srcMat,
                                RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tgtMat,
                                const size_t numSameIDs,
                                const LocalOrdinal permuteToLIDs[],
                                const LocalOrdinal permuteFromLIDs[],
                                const size_t numPermutes)
  {
    using Details::ProfilingRegion;
    using Teuchos::Array;
    //using Teuchos::ArrayView;
    using std::endl;
    using LO = LocalOrdinal;
    using GO = GlobalOrdinal;

    using impl_scalar_type = typename Kokkos::ArithTraits<Scalar>::val_type;
    typedef typename Kokkos::View<impl_scalar_type*, typename Node::device_type>::non_const_type  nonconst_values_device_view_type;
    typedef typename Kokkos::View<GlobalOrdinal *, typename Node::device_type>::non_const_type nonconst_global_inds_device_view_type;
    
    const LocalOrdinal LINV = Teuchos::OrdinalTraits<LocalOrdinal>::invalid ();

    const char tfecfFuncName[] = "copyAndPermuteStaticGraph";
    ProfilingRegion regionCAP
      ("Tpetra::CrsMatrix::copyAndPermuteStaticGraph");

    // const bool debug = Details::Behavior::debug("CrsGraph");
    // const bool verbose = Details::Behavior::verbose("CrsGraph");

    using crs_matrix_type = CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

    const crs_matrix_type *srcMatCrsPtr = dynamic_cast<const crs_matrix_type *>(&srcMat);
    if (!srcMatCrsPtr) {
      std::cout << "srk error srcMat type= " << typeid(srcMat).name() << std::endl;
      std::terminate();
    }
    const crs_matrix_type& srcMatCrs = *srcMatCrsPtr;

    crs_matrix_type *tgtMatCrsPtr = dynamic_cast<crs_matrix_type *>(&tgtMat);
    if (!tgtMatCrsPtr) {
      std::cout << "srk error tgtMat type= " << typeid(tgtMat).name() << std::endl;
      std::terminate();
    }
    crs_matrix_type& tgtMatCrs = *tgtMatCrsPtr;

    std::string prefix = tfecfFuncName;
    // const char* const prefix_raw = prefix.c_str();

    const bool sourceIsLocallyIndexed = srcMat.isLocallyIndexed ();
    //const bool targetIsLocallyIndexed = tgtMat.isLocallyIndexed ();
    //
    // Copy the first numSame row from source to target (this matrix).
    // This involves copying rows corresponding to LIDs [0, numSame-1].
    //
    const auto& srcRowMap = * (srcMat.getRowMap ());
    auto comm = srcRowMap.getComm();


    const LO numSameIDs_as_LID = static_cast<LO> (numSameIDs);

    if (sourceIsLocallyIndexed) {

      typedef typename crs_matrix_type::local_matrix_device_type k_local_matrix_device_type;
      //typedef typename k_local_matrix_device_type::StaticCrsGraphType k_graph_t;
      // typedef typename k_graph_t::row_map_type::non_const_type k_row_map_t;
      // typedef typename k_graph_t::entries_type::non_const_type k_nnz_t;
      // typedef typename k_local_matrix_device_type::values_type::non_const_type k_scalar_view_t;

      const k_local_matrix_device_type & srcMatDevice = srcMatCrs.getLocalMatrixDevice();
      const k_local_matrix_device_type & tgtMatDevice = tgtMatCrs.getLocalMatrixDevice();
      const k_local_matrix_device_type * srcMatDevicePtr = &srcMatDevice;
      const k_local_matrix_device_type * tgtMatDevicePtr = &tgtMatDevice;

      // auto lclIndsUnpacked_device = srcMat.getLocalIndicesDevice ();

      // auto nr = srcMatDevice.graph.numRows();
      // // if ((size_t)numSameIDs_as_LID >= (size_t)nr) {
      // //   std::cout << "numSameIDs_as_LID= " << numSameIDs_as_LID << " nr= " << nr << std::endl;
      // // }
      // TEUCHOS_ASSERT((size_t)numSameIDs_as_LID <= (size_t)nr);

#define PR1(a) std::cout << "[srk] " << #a << "= " << a << std::endl
      
      typename crs_matrix_type::row_ptrs_device_view_type tgtLocalRowPtrsDevice     = tgtMatCrs.getLocalRowPtrsDevice();
      typename crs_matrix_type::local_inds_device_view_type tgtLocalColIndsDevice   = tgtMatCrs.getLocalIndicesDevice();
      typename crs_matrix_type::row_ptrs_host_view_type srcLocalRowPtrsHost         = srcMatCrs.getLocalRowPtrsHost();
      typename crs_matrix_type::row_ptrs_device_view_type srcLocalRowPtrsDevice     = srcMatCrs.getLocalRowPtrsDevice();
      typename crs_matrix_type::local_inds_device_view_type srcLocalColIndsDevice   = srcMatCrs.getLocalIndicesDevice();

      nonconst_global_inds_device_view_type srowInfo(Kokkos::ViewAllocateWithoutInitializing("srowInfo"), numSameIDs_as_LID);

      printf("here fence 0 numSameIDs_as_LID= %ld\n", numSameIDs_as_LID);
      Kokkos::fence("srk0");
      printf("here fence 1\n");

      typedef typename Node::execution_space exec_space;
      typedef Kokkos::RangePolicy<exec_space, LO> range_type;

      size_t mre=0;
      for (LO sourceLID=0; sourceLID < numSameIDs_as_LID; sourceLID++) {
           auto start = srcLocalRowPtrsHost(sourceLID);
           auto end = srcLocalRowPtrsHost(sourceLID+1);
           size_t rowLength = static_cast<size_t>(end - start);
           printf("sourceLID= %d start= %d end= %d rowLength= %ld\n", sourceLID, start, end, rowLength);
           if (rowLength > mre) mre = rowLength;
      }
      printf("here b4 row_map, max_row_entries=%ld\n", mre);  // prints 33
      Kokkos::parallel_for
        ("Tpetra_CrsMatrix::copyAndPermuteStaticGraph",
         range_type (0, numSameIDs_as_LID),
         KOKKOS_LAMBDA(const LO sourceLID)
         {
           auto start = srcMatDevice.graph.row_map(sourceLID);  // always print 0
           auto end = srcMatDevice.graph.row_map(sourceLID+1);  // these print correctly
           size_t rowLength = static_cast<size_t>(end - start);
           printf("0 k_sourceLID= %d start= %d end= %d rowLength= %ld\n", sourceLID, start, end, rowLength);
           //printf("k_sourceLID= %d\n", sourceLID);
           //srowInfo(sourceLID) = rowLength;
         });  // kokkos parallel_for

      printf("here fence 2.0\n");
      Kokkos::fence("srk00");
      printf("here fence 2\n");

      printf("here b4 srowInfo, max_row_entries=%ld\n", mre);
      Kokkos::parallel_for
        ("Tpetra_CrsMatrix::copyAndPermuteStaticGraph",
         range_type (0, numSameIDs_as_LID),
         KOKKOS_LAMBDA(const LO sourceLID)
         {
           auto start = srcLocalRowPtrsDevice(sourceLID);
           auto end = srcLocalRowPtrsDevice(sourceLID+1);
           size_t rowLength = static_cast<size_t>(end - start);
           printf("1 k_sourceLID= %d start= %d end= %d rowLength= %ld\n", sourceLID, start, end, rowLength);
           //printf("k_sourceLID= %d\n", sourceLID);
           //srowInfo(sourceLID) = rowLength;
         });  // kokkos parallel_for

      printf("here fence 2.0\n");
      Kokkos::fence("srk00");
      printf("here fence 2\n");

      size_t max_row_entries = 0;
      Kokkos::parallel_reduce ("Tpetra_CrsMatrix_capsg_get_max_nc", range_type (0, numSameIDs_as_LID),
                               KOKKOS_LAMBDA(const LO sourceLID, size_t& gmax) {
                                 auto start = srcLocalRowPtrsDevice(sourceLID);
                                 auto end = srcLocalRowPtrsDevice(sourceLID+1);
                                 size_t ct = static_cast<size_t>(end - start);

                                 if (ct > gmax) gmax = ct;
                               }, max_row_entries);

      printf("here 0-af-pr: max_row_entries= %ld mre= %ld\n", max_row_entries, mre);
      max_row_entries = mre;

      Kokkos::fence("srk1");
      printf("here 0-af-fence:\n");

      auto local_map = srcMat.getRowMap()->getLocalMap();
      auto local_map_ptr = &local_map;
      auto local_col_map = srcMat.getColMap()->getLocalMap();
      auto local_col_map_ptr = &local_col_map;

      nonconst_global_inds_device_view_type rowInds(Kokkos::ViewAllocateWithoutInitializing("srk_rowInds"), max_row_entries);
      nonconst_values_device_view_type rowVals(Kokkos::ViewAllocateWithoutInitializing("srk_rowVals"), max_row_entries);

      bool tgtMatIsSorted = tgtMatCrs.getCrsGraph()->isSorted();

      using local_map_type = typename crs_matrix_type::map_type::local_map_type;
      //using crs_graph_type = CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;

      local_map_type tgt_local_map = tgtMatCrs.getRowMap()->getLocalMap();
      local_map_type tgt_local_col_map = tgtMatCrs.getColMap()->getLocalMap();
      auto tgt_local_map_ptr = &tgt_local_map;
      auto tgt_local_col_map_ptr = &tgt_local_col_map;

      auto my_replaceGlobalValuesImpl
        = KOKKOS_LAMBDA(
                        const bool sorted, const bool atomic, size_t hint[], 
                        const size_t numInTgtRow,   const LO tgtColInds[],       impl_scalar_type tgtRowVals[], 
                        const size_t numToReplace,  const GO inds[],       const impl_scalar_type newVals[]
                        ) -> LO
        {
         LO numValid = 0; // number of valid input column indices

         if (atomic) {
           for (LO j = 0; j < (LO)numToReplace; ++j) {
             const LO lclColInd = tgt_local_col_map_ptr->getLocalElement(inds[j]);
             if (lclColInd != LINV) {
               const size_t offset =
                 KokkosSparse::findRelOffset (tgtColInds, numInTgtRow,
                                              lclColInd, hint[0], sorted);
               if (offset != numInTgtRow) {
                 Kokkos::atomic_store (&tgtRowVals[offset], newVals[j]);
                 hint[0] = offset + 1;
                 numValid++;
               }
             }
           }
         } else {
           for (LO j = 0; j < (LO)numToReplace; ++j) {
             const LO lclColInd = tgt_local_col_map_ptr->getLocalElement (inds[j]);
             if (lclColInd != LINV) {
               const size_t offset =
                 KokkosSparse::findRelOffset (tgtColInds, numInTgtRow,
                                              lclColInd, hint[0], sorted);
               if (offset != numInTgtRow) {
                 tgtRowVals[offset] = newVals[j];
                 hint[0] = offset + 1;
                 numValid++;
               }
             }
           }
         }
         return numValid;
        };      

      printf("here 0: R: %d  row_map.size= %d numSameIDs_as_LID= %d\n", comm->getRank(), srcMatDevicePtr->graph.row_map.extent(0), numSameIDs_as_LID);

      printf("here 001: %d\n", comm->getRank());
      Kokkos::fence("srk01");
      printf("here 002: %d\n", comm->getRank());

      const GO rl0 =
        Tpetra::Details::getEntryOnHost(srowInfo, 0);

      printf("here 01: %d %ld\n", comm->getRank(), rl0);

      auto vals = srcMatCrs.getLocalValuesDevice (Access::ReadOnly);
      auto tvals = tgtMatCrs.getLocalValuesDevice (Access::ReadWrite);

      Kokkos::fence("srk01");

      Kokkos::parallel_for
        ("Tpetra_CrsMatrix::copyAndPermuteStaticGraph",
         range_type (0, numSameIDs_as_LID),
         KOKKOS_LAMBDA(const LO sourceLID)
         {
           //printf("sourceLID= %d\n", sourceLID);

           auto start = srcLocalRowPtrsDevice(sourceLID);
           auto end = srcLocalRowPtrsDevice(sourceLID+1);
           size_t rowLength = static_cast<size_t>(end - start);
           // srowInfo(sourceLID) = rowLength;
           [[maybe_unused]] size_t checkRowLength = 0;

           const size_t numEntries = rowLength;
           //auto dev_row_info = srcMatDevicePtr->row(sourceLID);

           //KOKKOS_ASSERT(dev_row_info.length == rowLength);

           checkRowLength = numEntries; // first side effect

#ifdef COMP
#  undef COMP
#endif
#define COMP(a,b) do { if (int(a) != int(b)) { std::cout << "error: " << #a << "= " << a << " " << #b << "= " << b << " line= " << __LINE__ << std::endl; std::terminate(); } } while(0)

            for (size_t j = 0; j < rowLength; j++) {
           //   //auto ci = dev_row_info.colidx(j);
              auto ci = srcLocalColIndsDevice(start + j);
              auto gi = local_col_map_ptr->getGlobalElement(ci);
           //   rowInds(j) = gi;
           //   rowVals(j) = vals(start + j);
           }

           // auto tgt_dev_row_info = tgtMatDevicePtr->row(sourceLID);
           // LO *tgtColInds = &tgt_dev_row_info.colidx(0);
           // Scalar *tgtRowVals = &tgt_dev_row_info.value(0);
           // size_t numInTgtRow = tgt_dev_row_info.length;
           auto tstart = tgtLocalRowPtrsDevice(sourceLID);
           auto tend = tgtLocalRowPtrsDevice(sourceLID + 1);
           size_t numInTgtRow = static_cast<size_t>(tend - tstart);
           Scalar *tgtRowVals = &tvals(tstart);
           const LO *tgtColInds = &tgtLocalColIndsDevice(tstart);

           size_t hint=0;
           // my_replaceGlobalValuesImpl(tgtMatIsSorted, false, &hint, //  tgt_local_col_map,
           //                            numInTgtRow, tgtColInds, tgtRowVals, 
           //                            rowLength,  rowInds.data(), rowVals.data());

         });  // kokkos parallel_for
      

      printf("here 02: %d\n", comm->getRank());

      Kokkos::fence("srk02");

      {
        bool tsd = tgtMatCrs.get_values_unpacked_wdv().need_sync_device();
        bool tsh = tgtMatCrs.get_values_unpacked_wdv().need_sync_host();
        bool ssd = srcMatCrs.get_values_unpacked_wdv().need_sync_device();
        bool ssh = srcMatCrs.get_values_unpacked_wdv().need_sync_host();

        PR1(tsd);
        PR1(tsh);
        PR1(ssd);
        PR1(ssh);
      }      

      // if (tsd) tgtMatCrs.get_values_unpacked_wdv().sync_device();
      // if (tsh) tgtMatCrs.get_values_unpacked_wdv().sync_host();
      // if (ssd) srcMatCrs.get_values_unpacked_wdv().sync_device();
      // if (ssh) srcMatCrs.get_values_unpacked_wdv().sync_host();
      // auto tgtSyncView = tgtMatCrs.get_values_unpacked_wdv().getHostView(Access::ReadOnly);
      // const RowInfo rowInfo = tgtMatCrs.getCrsGraph()->getRowInfo(0);
      // auto t2 = tgtMatCrs.getCrsGraph()->getLocalIndsViewHost(rowInfo);
      // std::cout << "[srk] tgtSyncView= " << tgtSyncView.extent(0) << " t2= " << t2.extent(0) << std::endl;

      Kokkos::fence("srk2");

      {
        bool tsd = tgtMatCrs.get_values_unpacked_wdv().need_sync_device();
        bool tsh = tgtMatCrs.get_values_unpacked_wdv().need_sync_host();
        bool ssd = srcMatCrs.get_values_unpacked_wdv().need_sync_device();
        bool ssh = srcMatCrs.get_values_unpacked_wdv().need_sync_host();

        PR1(tsd);
        PR1(tsh);
        PR1(ssd);
        PR1(ssh);
      }
      
      printf("here 1: %d\n", comm->getRank());


    } else {
      for (LO sourceLID = 0; sourceLID < numSameIDs_as_LID; ++sourceLID) {
        // Global ID for the current row index in the source matrix.
        // The first numSameIDs GIDs in the two input lists are the
        // same, so sourceGID == targetGID in this case.
        const GO sourceGID = srcRowMap.getGlobalElement (sourceLID);
        const GO targetGID = sourceGID;

        Teuchos::ArrayView<const GO> rowIndsConstView;
        Teuchos::ArrayView<const Scalar> rowValsConstView;

        typename crs_matrix_type::global_inds_host_view_type rowIndsView;
        typename crs_matrix_type::values_host_view_type rowValsView;
        srcMat.getGlobalRowView(sourceGID, rowIndsView, rowValsView);
        // KDDKDD UVM TEMPORARY:  refactor combineGlobalValues to take
        // KDDKDD UVM TEMPORARY:  Kokkos::View instead of ArrayView
        // KDDKDD UVM TEMPORARY:  For now, wrap the view in ArrayViews
        // KDDKDD UVM TEMPORARY:  Should be safe because we hold the KokkosViews
        rowIndsConstView = Teuchos::ArrayView<const GO> (  // BAD BAD BAD
                                                         rowIndsView.data(), rowIndsView.extent(0),
                                                         Teuchos::RCP_DISABLE_NODE_LOOKUP);
        rowValsConstView = Teuchos::ArrayView<const Scalar> (  // BAD BAD BAD
                                                             reinterpret_cast<const Scalar*>(rowValsView.data()), rowValsView.extent(0),
                                                             Teuchos::RCP_DISABLE_NODE_LOOKUP);
        // KDDKDD UVM TEMPORARY:  Add replace, sum, transform methods with
        // KDDKDD UVM TEMPORARY:  KokkosView interface

        // Applying a permutation to a matrix with a static graph
        // means REPLACE-ing entries.
        tgtMatCrs.replaceGlobalValues(targetGID, rowIndsConstView,
                                      rowValsConstView);
      }
    }

    //
    // "Permute" part of "copy and permute."
    //

    typename crs_matrix_type::nonconst_global_inds_host_view_type rowInds;
    typename crs_matrix_type::nonconst_values_host_view_type rowVals;

    const auto& tgtRowMap = * (tgtMat.getRowMap ());
    for (size_t p = 0; p < numPermutes; ++p) {
      const GO sourceGID = srcRowMap.getGlobalElement (permuteFromLIDs[p]);
      const GO targetGID = tgtRowMap.getGlobalElement (permuteToLIDs[p]);

      Teuchos::ArrayView<const GO> rowIndsConstView;
      Teuchos::ArrayView<const Scalar> rowValsConstView;

      if (sourceIsLocallyIndexed) {
        const size_t rowLength = srcMat.getNumEntriesInGlobalRow (sourceGID);
        if (rowLength > static_cast<size_t> (rowInds.size ())) {
          Kokkos::resize(rowInds,rowLength);
          Kokkos::resize(rowVals,rowLength);
        }
        // Resizing invalidates an Array's views, so we must make new
        // ones, even if rowLength hasn't changed.
        typename crs_matrix_type::nonconst_global_inds_host_view_type rowIndsView = Kokkos::subview(rowInds,std::make_pair((size_t)0, rowLength));
        typename crs_matrix_type::nonconst_values_host_view_type rowValsView = Kokkos::subview(rowVals,std::make_pair((size_t)0, rowLength));

        // The source matrix is locally indexed, so we have to get a
        // copy.  Really it's the GIDs that have to be copied (because
        // they have to be converted from LIDs).
        size_t checkRowLength = 0;
        srcMat.getGlobalRowCopy(sourceGID, rowIndsView,
                                rowValsView, checkRowLength);

        // KDDKDD UVM TEMPORARY:  refactor combineGlobalValues to take
        // KDDKDD UVM TEMPORARY:  Kokkos::View instead of ArrayView
        // KDDKDD UVM TEMPORARY:  For now, wrap the view in ArrayViews
        // KDDKDD UVM TEMPORARY:  Should be safe because we hold the KokkosViews
        rowIndsConstView = Teuchos::ArrayView<const GO> (  // BAD BAD BAD
                                                         rowIndsView.data(), rowIndsView.extent(0),
                                                         Teuchos::RCP_DISABLE_NODE_LOOKUP);
        rowValsConstView = Teuchos::ArrayView<const Scalar> (  // BAD BAD BAD
                                                             reinterpret_cast<const Scalar*>(rowValsView.data()), rowValsView.extent(0),
                                                             Teuchos::RCP_DISABLE_NODE_LOOKUP);
        // KDDKDD UVM TEMPORARY:  Add replace, sum, transform methods with
        // KDDKDD UVM TEMPORARY:  KokkosView interface
      }
      else {
        typename crs_matrix_type::global_inds_host_view_type rowIndsView;
        typename crs_matrix_type::values_host_view_type rowValsView;
        srcMat.getGlobalRowView(sourceGID, rowIndsView, rowValsView);
        // KDDKDD UVM TEMPORARY:  refactor combineGlobalValues to take
        // KDDKDD UVM TEMPORARY:  Kokkos::View instead of ArrayView
        // KDDKDD UVM TEMPORARY:  For now, wrap the view in ArrayViews
        // KDDKDD UVM TEMPORARY:  Should be safe because we hold the KokkosViews
        rowIndsConstView = Teuchos::ArrayView<const GO> (  // BAD BAD BAD
                                                         rowIndsView.data(), rowIndsView.extent(0),
                                                         Teuchos::RCP_DISABLE_NODE_LOOKUP);
        rowValsConstView = Teuchos::ArrayView<const Scalar> (  // BAD BAD BAD
                                                             reinterpret_cast<const Scalar*>(rowValsView.data()), rowValsView.extent(0),
                                                             Teuchos::RCP_DISABLE_NODE_LOOKUP);
        // KDDKDD UVM TEMPORARY:  Add replace, sum, transform methods with
        // KDDKDD UVM TEMPORARY:  KokkosView interface
      }

      tgtMatCrs.replaceGlobalValues(targetGID, rowIndsConstView,
                                    rowValsConstView);
    }

  }
