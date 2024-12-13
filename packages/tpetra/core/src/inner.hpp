// const GO gid = srcRowMapLocal.getGlobalElement (sourceLID);
//       size_t row_length = srcRowGraph.getNumEntriesInGlobalRow (gid);
//       Kokkos::resize(row_copy,row_length);
//       size_t check_row_length = 0;
//       srcRowGraph.getGlobalRowCopy (gid, row_copy, check_row_length);
//       tgtCrsGraph.insertGlobalIndices (gid, row_length, row_copy.data());

using global_inds_device_value_t = GlobalOrdinal;
using row_ptrs_device_value_t = size_t;
typedef typename crs_graph_type::local_graph_device_type k_local_graph_device_type;
typedef typename Node::execution_space exec_space;
typedef Kokkos::RangePolicy<exec_space, LO> range_type;

const LocalOrdinal LINV = Teuchos::OrdinalTraits<LocalOrdinal>::invalid ();
const GlobalOrdinal GINV = Teuchos::OrdinalTraits<GlobalOrdinal>::invalid ();

typedef typename crs_graph_type::local_graph_device_type k_local_graph_device_type;
typedef typename Node::execution_space exec_space;
typedef Kokkos::RangePolicy<exec_space, LO> range_type;
typedef typename Kokkos::GraphRowViewConst<local_graph_device_type> graph_row_view_const_type;

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

typedef typename Node::execution_space exec_space;
typedef Kokkos::RangePolicy<exec_space, LO> range_type;
typename num_row_entries_type::non_const_type h_numRowEnt = tgtCrsGraph.k_numRowEntries_;

auto k_numRowEnt = Kokkos::create_mirror_view_and_copy (device_type (), h_numRowEnt);
LO numSameIDs_as_LID = static_cast<LO>(numSameIDs);
const bool sorted = false;

#define CHECK(a,i) do {                                                 \
    if (i >= a.extent(0) || i < 0) {                                    \
      char buf[100];                                                    \
      sprintf(buf,"ERROR: i= %d a= %s e= %d", i, #a, a.extent(0));      \
      Kokkos::abort(buf);                                               \
    } } while(0)

#define AB(lin) do {                            \
    char buf[100];                              \
    sprintf(buf,"ERROR: line= %d", lin);        \
    Kokkos::abort(buf);                         \
  } while(0)

Kokkos::parallel_for("Tpetra_CrsGraph::copyAndPermuteNew2",
                     range_type (0, numSameIDs_as_LID),
                     KOKKOS_LAMBDA(const LO sourceLID)
                     {
                       auto srcGid = srcRowMapLocal.getGlobalElement(sourceLID);
                       auto tgtLocalRow = tgtRowMapLocal.getLocalElement(srcGid);
                       if (tgtLocalRow == LINV) {
                         AB(__LINE__);
                       }
                       CHECK(k_numRowEnt, tgtLocalRow);
                       auto tgtNumEntries = k_numRowEnt(tgtLocalRow);

                       // FIXME no auto use
                       CHECK(srcLocalRowPtrsDevice, sourceLID);
                       auto start     = srcLocalRowPtrsDevice(sourceLID);
                       CHECK(srcLocalRowPtrsDevice, sourceLID+1);
                       auto end       = srcLocalRowPtrsDevice(sourceLID+1);
                       auto rowLength = (end - start);

                       //KOKKOS_ASSERT(rowLength <= max_row_entries);

                       CHECK(tgtLocalRowPtrsDevice, tgtLocalRow);
                       auto tstart      = tgtLocalRowPtrsDevice(tgtLocalRow);
                       auto tend        = tstart + tgtNumEntries;
                       CHECK(tgtLocalRowPtrsDevice, tgtLocalRow + 1);
                       auto tend1       = tgtLocalRowPtrsDevice(tgtLocalRow + 1);

                       const size_t num_avail = (tend1 < tend) ? size_t (0) : tend1 - tend;
                       const size_t num_new_indices = rowLength;
                       size_t num_inserted = 0;

                       CHECK(tgtGlobalColInds, tstart);
                       //global_inds_device_value_t *tgtColInds = tgtGlobalColInds.data()+tstart;
                       global_inds_device_value_t *tgtGlobalColIndsPtr = tgtGlobalColInds.data();

                       size_t hint=0;
                       for (LO j = 0; j < rowLength; j++) {
                         auto ci = srcLocalColIndsDevice(start + j);
                         GO gi = srcColMapLocal.getGlobalElement(ci);
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
                           tgtGlobalColIndsPtr[tstart + offset] = gi;
                           ++tend;
                           hint = offset + 1;
                           ++num_inserted;
                         }
                       }
                       k_numRowEnt(tgtLocalRow) += num_inserted;

                       return size_t(0);
                     });

std::cout << "here 10 host size= " << tgtCrsGraph.k_numRowEntries_.extent(0) << " dev: " << k_numRowEnt.extent(0) << std::endl;

Kokkos::fence("here 10");
  std::cout << "here 10.1" << std::endl;

Kokkos::deep_copy(tgtCrsGraph.k_numRowEntries_, k_numRowEnt);
  std::cout << "here 11" << std::endl;

tgtCrsGraph.setLocallyModified();
  std::cout << "here 12" << std::endl;
