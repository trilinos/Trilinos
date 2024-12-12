// const GO gid = srcRowMapLocal.getGlobalElement (sourceLID);
//       size_t row_length = srcRowGraph.getNumEntriesInGlobalRow (gid);
//       Kokkos::resize(row_copy,row_length);
//       size_t check_row_length = 0;
//       srcRowGraph.getGlobalRowCopy (gid, row_copy, check_row_length);
//       tgtCrsGraph.insertGlobalIndices (gid, row_length, row_copy.data());

// using RowInfoViewType = Kokkos::View<graph_row_view_const_type *, Kokkos::LayoutRight, device_type>;
// RowInfoViewType tgtRowInfoView("RowInfoView", numSameIDs);
// kokkos::parallel_for("Tpetra_CrsGraph::copyAndPermuteNew1",
//                      range_type (0, numSameIDs_as_LID),
//                      KOKKOS_LAMBDA(const LO sourceLID) {
//                        tgtRowInfoView(sourceLID) = tgtGraphDevice.rowConst(sourceLID);
//                      }
//                      );

  std::cout << "here 0" << std::endl;

using global_inds_device_value_t = GlobalOrdinal;
using row_ptrs_device_value_t = size_t;
typedef typename crs_graph_type::local_graph_device_type k_local_graph_device_type;
typedef typename Node::execution_space exec_space;
typedef Kokkos::RangePolicy<exec_space, LO> range_type;

const LocalOrdinal LINV = Teuchos::OrdinalTraits<LocalOrdinal>::invalid ();
const GlobalOrdinal GINV = Teuchos::OrdinalTraits<GlobalOrdinal>::invalid ();

// typedef typename crs_graph_type::global_inds_device_view_type::non_const_value_type global_inds_device_value_t; 
// typedef typename crs_graph_type::row_ptrs_device_view_type::non_const_value_type row_ptrs_device_value_t;
typedef typename crs_graph_type::local_graph_device_type k_local_graph_device_type;
typedef typename Node::execution_space exec_space;
typedef Kokkos::RangePolicy<exec_space, LO> range_type;
typedef typename Kokkos::GraphRowViewConst<local_graph_device_type> graph_row_view_const_type;

const k_local_graph_device_type & srcGraphDevice = srcCrsGraph.getLocalGraphDevice();
const k_local_graph_device_type & tgtGraphDevice = tgtCrsGraph.getLocalGraphDevice();
  std::cout << "here 1" << std::endl;

using local_map_type = typename crs_graph_type::map_type::local_map_type;

local_map_type srcRowMapLocal = srcCrsGraph.getRowMap()->getLocalMap();
  std::cout << "here 2" << std::endl;
local_map_type srcColMapLocal = srcCrsGraph.getColMap()->getLocalMap();
  std::cout << "here 3" << std::endl;

local_map_type tgtRowMapLocal = tgtCrsGraph.getRowMap()->getLocalMap();
  std::cout << "here 4" << std::endl;
// std::cout << "here 4.1 " << tgtCrsGraph.getColMap() << std::endl;
// local_map_type tgtColMapLocal = tgtCrsGraph.getColMap()->getLocalMap();
//   std::cout << "here 5" << std::endl;

auto tgtLocalRowPtrsDevice = tgtCrsGraph.getLocalRowPtrsDevice();
auto tgtLocalColIndsDevice = tgtCrsGraph.getLocalIndicesDevice();
auto tgtLocalColIndsDeviceNonConst = tgtCrsGraph.lclIndsUnpacked_wdv.getDeviceView(Access::ReadWrite);
auto tgtGlobalColInds = tgtCrsGraph.gblInds_wdv.getDeviceView(Access::ReadWrite);
//auto srcLocalRowPtrsHost   = srcCrsGraph.getLocalRowPtrsHost();
auto srcLocalRowPtrsDevice = srcCrsGraph.getLocalRowPtrsDevice();
auto srcLocalColIndsDevice = srcCrsGraph.getLocalIndicesDevice();
  std::cout << "here 7" << std::endl;

typedef typename Node::execution_space exec_space;
typedef Kokkos::RangePolicy<exec_space, LO> range_type;
typename num_row_entries_type::non_const_type h_numRowEnt = tgtCrsGraph.k_numRowEntries_;
std::cout << "here 8: " << h_numRowEnt.extent(0) << " numSameIDs= " << numSameIDs << std::endl;

auto k_numRowEnt = Kokkos::create_mirror_view_and_copy (device_type (), h_numRowEnt);

  std::cout << "here 9" << std::endl;

LO numSameIDs_as_LID = static_cast<LO>(numSameIDs);
const bool sorted = false;

Kokkos::parallel_for("Tpetra_CrsGraph::copyAndPermuteNew2",
                     range_type (0, numSameIDs_as_LID),
                     KOKKOS_LAMBDA(const LO sourceLID)
                     {
                       auto srcGid = srcRowMapLocal.getGlobalElement(sourceLID);
                       auto tgtLocalRow = tgtRowMapLocal.getLocalElement(srcGid);
                       auto tgtRowInfo = tgtGraphDevice.rowConst(tgtLocalRow);
                       auto tgtNumEntries = tgtRowInfo.length;

                       auto start     = srcLocalRowPtrsDevice(sourceLID);
                       auto end       = srcLocalRowPtrsDevice(sourceLID+1);
                       auto rowLength = (end - start);

                       //KOKKOS_ASSERT(rowLength <= max_row_entries);

                       auto tstart      = tgtLocalRowPtrsDevice(tgtLocalRow);
                       auto tend        = tstart + tgtNumEntries;
                       auto tend1       = tgtLocalRowPtrsDevice(tgtLocalRow + 1);

                       const size_t num_avail = (tend1 < tend) ? size_t (0) : tend1 - tend;
                       const size_t num_new_indices = rowLength;
                       size_t num_inserted = 0;

                       //local_inds_device_value_t *tgtColInds = tgtLocalColIndsDeviceNonConst.data()+tstart;
                       global_inds_device_value_t *tgtColInds = tgtGlobalColInds.data()+tstart;

                       size_t hint=0;
                       for (LO j = 0; j < rowLength; j++) {
                         auto ci = srcLocalColIndsDevice(start + j);
                         GO gi = srcColMapLocal.getGlobalElement(ci);
                         //const auto lclColInd = tgtColMapLocal.getLocalElement(gi);

                         auto numInTgtRow = (tend - tstart);

                         const size_t offset =
                           KokkosSparse::findRelOffset (tgtColInds, numInTgtRow,
                                                        gi, hint, sorted);

                         if (offset == numInTgtRow) {
                           if (num_inserted >= num_avail) { // not enough room
                             return Teuchos::OrdinalTraits<size_t>::invalid();
                           }
                           //Kokkos::atomic_store (&tgtRowVals[offset], newVals);
                           tgtColInds[tend] = gi;
                           ++tend;
                           hint = offset + 1;
                           ++num_inserted;
                         }
                         k_numRowEnt(tgtLocalRow) += num_inserted;

                       }             
                       return size_t(0);
                     });
  std::cout << "here 10" << std::endl;

Kokkos::deep_copy(exec_space(), tgtCrsGraph.k_numRowEntries_, k_numRowEnt);
  std::cout << "here 11" << std::endl;

tgtCrsGraph.setLocallyModified();
  std::cout << "here 12" << std::endl;
