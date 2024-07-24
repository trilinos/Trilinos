// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GDSW_DEF_HPP
#define GDSW_DEF_HPP
#include "GDSW_Proxy_decl.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "Tpetra_Details_PackTraits.hpp"

using Teuchos::RCP;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraFunctions<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
importSquareMatrix(RCP<const CrsMatrix> inputMatrix, 
                   RCP<const Map> outputRowMap, 
                   RCP<CrsMatrix> & outputMatrix)
{
    RCP<Tpetra::Distributor> distributor;
    std::vector<GO> ownedRowGIDs;
    std::vector<LO> localRowsSend, localRowsSendBegin, localRowsRecv, localRowsRecvBegin;
    // ownedRowGIDs = rows of inputMatrix contained in outputRowMap
    // localRowsSend[indicesSend] = localRows of sending process i corresponding to
    //                              localRowsRecv[indicesRecv], where
    //               indicesSend = localRowsSendBegin[i]:localRowsSendBegin[i+1]-1
    //               indicesRecv = localRowsRecvBegin[i]:localRowsRecvBegin[i+1]-1
    constructDistributor(inputMatrix, outputRowMap, distributor, ownedRowGIDs,
                         localRowsSend, localRowsSendBegin,
                         localRowsRecv, localRowsRecvBegin);
    std::vector<GO> targetMapGIDs;
    std::vector<LO> targetMapGIDsBegin;
    // targetMapGIDs[indices] = globalIDs of outputRowMap for the i'th process
    //                          receiving matrix data
    //         indices = targetMapGIDsBegin[i]:targetMapGIDsBegin[i+1]-1;
    // Note: length of targetMapGIDsBegin is the number of processes receiving 
    //       matrix data plus 1
    communicateRowMap(outputRowMap, distributor, targetMapGIDs, targetMapGIDsBegin);
    communicateMatrixData(inputMatrix, outputRowMap, distributor, targetMapGIDs, 
                          targetMapGIDsBegin, ownedRowGIDs,
                          localRowsSend, localRowsSendBegin, localRowsRecv,
                          localRowsRecvBegin, outputMatrix);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraFunctions<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
importSquareMatrixFromImporter(RCP<const CrsMatrix> inputMatrix, 
                               RCP<const Import> importer,
                               RCP<CrsMatrix> & outputMatrix)
{
  RCP<const Map> outputRowMap = importer->getTargetMap();
  RCP<const Map> sourceMap = importer->getSourceMap();
  RCP<Tpetra::Distributor> distributor = Teuchos::rcpFromRef(importer->getDistributor());
  Teuchos::ArrayView<const LO> localRowsSend = importer->getExportLIDs();
  Teuchos::ArrayView<const LO> localRowsRecv = importer->getRemoteLIDs();


  // Scansum to get the begin arrays
  size_t numSends = distributor->getNumSends();
  size_t numRecvs = distributor->getNumReceives();
  Teuchos::ArrayView<const size_t> lengthsTo = distributor->getLengthsTo();
  Teuchos::ArrayView<const size_t> lengthsFrom = distributor->getLengthsFrom();
  std::vector<LO> localRowsSendBegin(numSends+1);
  std::vector<LO> localRowsRecvBegin(numRecvs+1);
  for (size_t i=0; i<numSends; i++)
    localRowsSendBegin[i+1] = localRowsSendBegin[i]  + lengthsTo[i];
  for (size_t i=0; i<numRecvs; i++)
    localRowsRecvBegin[i+1] = localRowsRecvBegin[i]  + lengthsFrom[i];


  // Combine sames and permutes for ownedRowGIDs
  // QUESTION: Does iSM handle permutes correctly?
  size_t numSames = importer->getNumSameIDs();
  size_t numPermutes = importer->getNumPermuteIDs();
  Teuchos::ArrayView<const LO> permutes = importer->getPermuteFromLIDs();
  std::vector<GO> ownedRowGIDs(numSames + numPermutes);
  for(size_t i=0; i <numSames; i++)
    ownedRowGIDs[i] = sourceMap->getGlobalElement(i);
  for(size_t i=0; i <numPermutes; i++)
    ownedRowGIDs[i] = numSames + sourceMap->getGlobalElement(permutes[i]);
  


  //    std::vector<LO> localRowsSend, localRowsSendBegin, localRowsRecv, localRowsRecvBegin;
    // ownedRowGIDs = rows of inputMatrix contained in outputRowMap
    // localRowsSend[indicesSend] = localRows of sending process i corresponding to
    //                              localRowsRecv[indicesRecv], where
    //               indicesSend = localRowsSendBegin[i]:localRowsSendBegin[i+1]-1
    //               indicesRecv = localRowsRecvBegin[i]:localRowsRecvBegin[i+1]-1
    //    constructDistributor(inputMatrix, outputRowMap, distributor, ownedRowGIDs,
    //                         localRowsSend, localRowsSendBegin,
    //                         localRowsRecv, localRowsRecvBegin);
    std::vector<GO> targetMapGIDs;
    std::vector<LO> targetMapGIDsBegin;
    // targetMapGIDs[indices] = globalIDs of outputRowMap for the i'th process
    //                          receiving matrix data
    //         indices = targetMapGIDsBegin[i]:targetMapGIDsBegin[i+1]-1;
    // Note: length of targetMapGIDsBegin is the number of processes receiving 
    //       matrix data plus 1
    communicateRowMap(outputRowMap, distributor, targetMapGIDs, targetMapGIDsBegin);
    communicateMatrixData(inputMatrix, outputRowMap, distributor, targetMapGIDs, 
                          targetMapGIDsBegin, ownedRowGIDs,
                          localRowsSend, localRowsSendBegin, localRowsRecv,
                          localRowsRecvBegin, outputMatrix);
}

/****************************************************************************************************/
/****************************************************************************************************/
/****************************************************************************************************/

#ifdef AINT_WORKING_YET
template<class OutputOffsetsViewType,
         class CountsViewType,
         class InputOffsetsViewType,
         class InputColIndicesViewType,
         class InputColMapViewType,
         class InputLocalRowIndicesViewType,
         class InputLocalRowPidsViewType,
         class AllowedColumnsViewType,
         class AllowedColumnsBeginViewType,
         const bool debug =
#ifdef HAVE_TPETRA_DEBUG
         true
#else
         false
#endif // HAVE_TPETRA_DEBUG
         >
class NumPacketsAndOffsetsFunctorRestricted {
public:
  typedef typename OutputOffsetsViewType::non_const_value_type output_offset_type;
  typedef typename CountsViewType::non_const_value_type count_type;
  typedef typename InputOffsetsViewType::data_type input_offset_data_type;
  typedef typename InputLocalRowIndicesViewType::non_const_value_type local_row_index_type;
  typedef typename InputLocalRowPidsViewType::non_const_value_type local_row_pid_type;  
  // output Views drive where execution happens.
  typedef typename OutputOffsetsViewType::device_type device_type;
  static_assert (std::is_same<typename CountsViewType::device_type::execution_space,
                   typename device_type::execution_space>::value,
                 "OutputOffsetsViewType and CountsViewType must have the same execution space.");
  static_assert (Kokkos::is_view<OutputOffsetsViewType>::value,
                 "OutputOffsetsViewType must be a Kokkos::View.");
  static_assert (std::is_same<typename OutputOffsetsViewType::value_type, output_offset_type>::value,
                 "OutputOffsetsViewType must be a nonconst Kokkos::View.");
  static_assert (std::is_integral<output_offset_type>::value,
                 "The type of each entry of OutputOffsetsViewType must be a built-in integer type.");
  static_assert (Kokkos::is_view<CountsViewType>::value,
                 "CountsViewType must be a Kokkos::View.");
  static_assert (std::is_same<typename CountsViewType::value_type, output_offset_type>::value,
                 "CountsViewType must be a nonconst Kokkos::View.");
  static_assert (std::is_integral<count_type>::value,
                 "The type of each entry of CountsViewType must be a built-in integer type.");
  static_assert (Kokkos::is_view<InputOffsetsViewType>::value,
                 "InputOffsetsViewType must be a Kokkos::View.");
  static_assert (std::is_integral<input_offset_data_type>::value,
                 "The type of each entry of InputOffsetsViewType must be a built-in integer type.");
  static_assert (Kokkos::is_view<InputLocalRowIndicesViewType>::value,
                 "InputLocalRowIndicesViewType must be a Kokkos::View.");
  static_assert (std::is_integral<local_row_index_type>::value,
                 "The type of each entry of InputLocalRowIndicesViewType must be a built-in integer type.");

  NumPacketsAndOffsetsFunctor (const OutputOffsetsViewType& outputOffsets,
                               const CountsViewType& counts,
                               const InputOffsetsViewType& rowOffsets,
                               const InputColIndicesViewType & entries,
                               const InputColMapViewType & colMap,
                               const InputLocalRowIndicesViewType& lclRowInds,
                               const InputLocalRowPidsViewType& lclRowPids,
                               const AllowedColumnsViewType &allowedColumns,
                               const AllowedColumnsBeginViewType &allowedColumnsBegin,
                               const count_type sizeOfLclCount,
                               const count_type sizeOfGblColInd,
                               const count_type sizeOfPid,
                               const count_type sizeOfValue) :
    outputOffsets_ (outputOffsets),
    counts_ (counts),
    rowOffsets_ (rowOffsets),
    entries_ (entries),
    colMap_ (colMap),
    lclRowInds_ (lclRowInds),
    lclRowPids_ (lclRowPids),
    allowedColumns_ (allowedColumns),
    allowedColumnsBegin_ (allowedColumnsBegin),
    sizeOfLclCount_ (sizeOfLclCount),
    sizeOfGblColInd_ (sizeOfGblColInd),
    sizeOfPid_ (sizeOfPid),
    sizeOfValue_ (sizeOfValue),
    error_ ("error") // don't forget this, or you'll get segfaults!
  {
    if (debug) {
      const size_t numRowsToPack = static_cast<size_t> (lclRowInds_.extent (0));

      if (numRowsToPack != static_cast<size_t> (counts_.extent (0))) {
        std::ostringstream os;
        os << "lclRowInds.extent(0) = " << numRowsToPack
           << " != counts.extent(0) = " << counts_.extent (0)
           << ".";
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str ());
      }
      if (static_cast<size_t> (numRowsToPack + 1) !=
          static_cast<size_t> (outputOffsets_.extent (0))) {
        std::ostringstream os;
        os << "lclRowInds.extent(0) + 1 = " << (numRowsToPack + 1)
           << " != outputOffsets.extent(0) = " << outputOffsets_.extent (0)
           << ".";
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str ());
      }
    }
  }

  KOKKOS_INLINE_FUNCTION void
  operator() (const local_row_index_type& curInd,
              output_offset_type& update,
              const bool final) const
  {
    if (debug) {
      if (curInd < static_cast<local_row_index_type> (0)) {
        error_ () = 1;
        return;
      }
    }

    if (final) {
      if (debug) {
        if (curInd >= static_cast<local_row_index_type> (outputOffsets_.extent (0))) {
          error_ () = 2;
          return;
        }
      }
      outputOffsets_(curInd) = update;
    }

    if (curInd < static_cast<local_row_index_type> (counts_.extent (0))) {
      const auto lclRow = lclRowInds_(curInd);
      if (static_cast<size_t> (lclRow + 1) >= static_cast<size_t> (rowOffsets_.extent (0)) ||
          static_cast<local_row_index_type> (lclRow) < static_cast<local_row_index_type> (0)) {
        error_ () = 3;
        return;
      }
      // count_type could differ from the type of each row offset.
      // For example, row offsets might each be 64 bits, but if their
      // difference always fits in 32 bits, we may then safely use a
      // 32-bit count_type.
      const count_type count;

      // Here we need to check each entry in lclRow to see if it is something
      // we can send based on the allowedColumns list
      for(auto j=rowOffsets_[lclRow], j<rowOffsets_[lclRow+1]; j++) {
        auto l_cid = entries_[j];
        auto l_gid = columnMap.getGlobalElement(l_cid);



      }
      



        static_cast<count_type> (rowOffsets_(lclRow+1) - rowOffsets_(lclRow));

      // We pack first the number of entries in the row, then that
      // many global column indices, then that many pids (if any),
      // then that many values.  However, if the number of entries in
      // the row is zero, we pack nothing.
      const count_type numBytes = (count == 0) ?
        static_cast<count_type> (0) :
        sizeOfLclCount_ + count * (sizeOfGblColInd_ +
                                   (lclRowPids_.size() > 0 ? sizeOfPid_ : 0) +
                                   sizeOfValue_);

      if (final) {
        counts_(curInd) = numBytes;
      }
      update += numBytes;
    }
  }

  // mfh 31 May 2017: Don't need init or join.  If you have join, MUST
  // have join both with and without volatile!  Otherwise intrawarp
  // joins are really slow on GPUs.

  //! Host function for getting the error.
  int getError () const {
    typedef typename device_type::execution_space execution_space;
    auto error_h = Kokkos::create_mirror_view (error_);
    // DEEP_COPY REVIEW - DEVICE-TO-HOSTMIRROR
    Kokkos::deep_copy (execution_space(), error_h, error_);
    return error_h ();
  }

private:
  OutputOffsetsViewType outputOffsets_;
  CountsViewType counts_;
  typename InputOffsetsViewType::const_type rowOffsets_;
  typename InputColIndicesViewType::const_type entries_;
  typename InputColMapViewType::const_type colMap_,  
  typename InputLocalRowIndicesViewType::const_type lclRowInds_;
  typename InputLocalRowPidsViewType::const_type lclRowPids_;
  typename AllowedColumnsViewType::const_type &allowedColumns_;
  typename AllowedColumnsBeginViewType::const_type &allowedColumnsBegin_,

  count_type sizeOfLclCount_;
  count_type sizeOfGblColInd_;
  count_type sizeOfPid_;
  count_type sizeOfValue_;
  Kokkos::View<int, device_type> error_;
};



template<class OutputOffsetsViewType,
         class CountsViewType,
         class InputOffsetsViewType,
         class InputColIndicesViewType,
         class InputColMapViewType,
         class InputLocalRowIndicesViewType,
         class InputLocalRowPidsViewType,
         class AllowedColumnsViewType,
         class AllowedColumnsBeginViewType>
typename CountsViewType::non_const_value_type
computeNumPacketsAndOffsetsRestricted (const OutputOffsetsViewType& outputOffsets,
                                       const CountsViewType& counts,
                                       const InputOffsetsViewType& rowOffsets,
                                       const InputColIndicesViewType & entries,
                                       const InputColMapViewType & colMap,
                                       const InputLocalRowIndicesViewType& lclRowInds,
                                       const InputLocalRowPidsViewType& lclRowPids,
                                       const AllowedColumnsViewType &allowedColumns,
                                       const AllowedColumnsBeginViewType &allowedColumnsBegin,
                                       const typename CountsViewType::non_const_value_type sizeOfLclCount,
                                       const typename CountsViewType::non_const_value_type sizeOfGblColInd,
                                       const typename CountsViewType::non_const_value_type sizeOfPid,
                                       const typename CountsViewType::non_const_value_type sizeOfValue)
{
  typedef NumPacketsAndOffsetsFunctor<OutputOffsetsViewType,
    CountsViewType, typename InputOffsetsViewType::const_type,
    typename InputLocalRowIndicesViewType::const_type,
    typename InputLocalRowPidsViewType::const_type> functor_type;
  typedef typename CountsViewType::non_const_value_type count_type;
  typedef typename OutputOffsetsViewType::size_type size_type;
  typedef typename OutputOffsetsViewType::execution_space execution_space;
  typedef typename functor_type::local_row_index_type LO;
  typedef Kokkos::RangePolicy<execution_space, LO> range_type;
  const char prefix[] = "computeNumPacketsAndOffsets: ";

  count_type count = 0;
  const count_type numRowsToPack = lclRowInds.extent (0);

  if (numRowsToPack == 0) {
    return count;
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION
      (rowOffsets.extent (0) <= static_cast<size_type> (1),
       std::invalid_argument, prefix << "There is at least one row to pack, "
       "but the matrix has no rows.  lclRowInds.extent(0) = " <<
       numRowsToPack << ", but rowOffsets.extent(0) = " <<
       rowOffsets.extent (0) << " <= 1.");
    TEUCHOS_TEST_FOR_EXCEPTION
      (outputOffsets.extent (0) !=
       static_cast<size_type> (numRowsToPack + 1), std::invalid_argument,
       prefix << "Output dimension does not match number of rows to pack.  "
       << "outputOffsets.extent(0) = " << outputOffsets.extent (0)
       << " != lclRowInds.extent(0) + 1 = "
       << static_cast<size_type> (numRowsToPack + 1) << ".");
    TEUCHOS_TEST_FOR_EXCEPTION
      (counts.extent (0) != numRowsToPack, std::invalid_argument,
       prefix << "counts.extent(0) = " << counts.extent (0)
       << " != numRowsToPack = " << numRowsToPack << ".");

    functor_type f (outputOffsets, counts, rowOffsets,
                    entries,colMap,
                    lclRowInds, lclRowPids, allowedColumns,
                    allowedColumnsBegin, sizeOfLclCount,
                    sizeOfGblColInd, sizeOfPid, sizeOfValue);
    Kokkos::parallel_scan (range_type (0, numRowsToPack + 1), f);

    // At least in debug mode, this functor checks for errors.
    const int errCode = f.getError ();
    TEUCHOS_TEST_FOR_EXCEPTION
      (errCode != 0, std::runtime_error, prefix << "parallel_scan error code "
       << errCode << " != 0.");

    // Get last entry of outputOffsets, which is the sum of the entries
    // of counts.  Don't assume UVM.
    using Tpetra::Details::getEntryOnHost;
    return static_cast<count_type> (getEntryOnHost (outputOffsets,
                                                    numRowsToPack));
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node,typename BufferDeviceType>
void
packCrsMatrixRestricted (const CrsMatrix<ST, LO, GO, NT>& sourceMatrix,
                         Kokkos::DualView<char*, BufferDeviceType>& exports,
                         const Kokkos::View<size_t*, BufferDeviceType>& num_packets_per_lid,
                         const Kokkos::View<const LO*, BufferDeviceType>& export_lids,
                         const Kokkos::View<const int*, typename NT::device_type>& export_pids,
                         const Kokkos::View<const GO*, BufferDeviceType> &restricted_col_gids;
                         const Kokkos::View<const size_t*, BufferDeviceType> &restricted_col_gids_begin,
                         size_t& constant_num_packets,
                         const bool pack_pids)
{
  using Kokkos::View;
  typedef BufferDeviceType DT;
  typedef Kokkos::DualView<char*, BufferDeviceType> exports_view_type;
  auto local_matrix = sourceMatrix.getLocalMatrixDevice();
  auto local_col_map = sourceMatrix.getColMap ()->getLocalMap();


  // Setting this to zero tells the caller to expect a possibly
  // different ("nonconstant") number of packets per local index
  // (i.e., a possibly different number of entries per row).
  constant_num_packets = 0;

  const size_t num_export_lids =
    static_cast<size_t> (export_lids.extent (0));
  TEUCHOS_TEST_FOR_EXCEPTION
    (num_export_lids !=
     static_cast<size_t> (num_packets_per_lid.extent (0)),
     std::invalid_argument, prefix << "num_export_lids.extent(0) = "
     << num_export_lids << " != num_packets_per_lid.extent(0) = "
     << num_packets_per_lid.extent (0) << ".");
  if (num_export_lids != 0) {
    TEUCHOS_TEST_FOR_EXCEPTION
      (num_packets_per_lid.data () == NULL, std::invalid_argument,
       prefix << "num_export_lids = "<< num_export_lids << " != 0, but "
       "num_packets_per_lid.data() = "
       << num_packets_per_lid.data () << " == NULL.");
  }

  // Useful information
  const size_t num_bytes_per_lid = PackTraits<LO>::packValueCount (LO (0));
  const size_t num_bytes_per_gid = PackTraits<GO>::packValueCount (GO (0));
  const size_t num_bytes_per_pid = PackTraits<int>::packValueCount (int (0));

  size_t num_bytes_per_value = 0;
  if (PackTraits<ST>::compileTimeSize) {
    // Assume ST is default constructible; packValueCount wants an instance.
    num_bytes_per_value = PackTraits<ST>::packValueCount (ST ());
  }
  else {
    // Since the packed data come from the source matrix, we can use
    // the source matrix to get the number of bytes per Scalar value
    // stored in the matrix.  This assumes that all Scalar values in
    // the source matrix require the same number of bytes.  If the
    // source matrix has no entries on the calling process, then we
    // hope that some process does have some idea how big a Scalar
    // value is.  Of course, if no processes have any entries, then no
    // values should be packed (though this does assume that in our
    // packing scheme, rows with zero entries take zero bytes).
    size_t num_bytes_per_value_l = 0;
    if (local_matrix.values.extent(0) > 0) {
      const ST& val = local_matrix.values(0);
      num_bytes_per_value_l = PackTraits<ST>::packValueCount (val);
    }
    using Teuchos::reduceAll;
    reduceAll<int, size_t> (* (inputMatrix.getComm ()),
                            Teuchos::REDUCE_MAX,
                            num_bytes_per_value_l,
                            Teuchos::outArg (num_bytes_per_value));
  }
  if (num_export_lids == 0) {
    exports = exports_view_type ("exports", 0);
    return;
  }

  // Array of offsets into the pack buffer.
  Kokkos::View<size_t*, DT> offsets ("offsets", num_export_lids + 1);
  
  // Compute number of packets per LID (row to send), as well as
  // corresponding offsets (the prefix sum of the packet counts).
  const size_t count =
    computeNumPacketsAndOffsetsRestricted (offsets, num_packets_per_lid,
                                           local_matrix.graph.row_map, 
                                           local_matrix.graph.entries,
                                           local_col_map,
                                           export_lids,
                                           export_pids,
                                           restricted_col_gids,
                                           restricted_col_gids_begin,
                                           num_bytes_per_lid, num_bytes_per_gid,
                                           num_bytes_per_pid, num_bytes_per_value);

  // Resize the output pack buffer if needed.
  if (count > static_cast<size_t> (exports.extent (0))) {
    exports = exports_view_type ("exports", count);
    if (debug) {
      std::ostringstream os;
      os << "*** exports resized to " << count << std::endl;
      std::cerr << os.str ();
    }
  }
  if (debug) {
    std::ostringstream os;
    os << "*** count: " << count << ", exports.extent(0): "
       << exports.extent (0) << std::endl;
    std::cerr << os.str ();
  }

  // If exports has nonzero length at this point, then the matrix has
  // at least one entry to pack.  Thus, if packing process ranks, we
  // had better have at least one process rank to pack.
  TEUCHOS_TEST_FOR_EXCEPTION
    (pack_pids && exports.extent (0) != 0 &&
     export_pids.extent (0) == 0, std::invalid_argument, prefix <<
     "pack_pids is true, and exports.extent(0) = " <<
     exports.extent (0)  << " != 0, meaning that we need to pack at least "
     "one matrix entry, but export_pids.extent(0) = 0.");

  typedef typename std::decay<decltype (local_matrix)>::type
    local_matrix_device_type;
  typedef typename std::decay<decltype (local_col_map)>::type
    local_map_type;

  exports.modify_device ();
  auto exports_d = exports.view_device ();
  do_pack<local_matrix_device_type, local_map_type, DT>
    (local_matrix, local_col_map, exports_d, num_packets_per_lid,
     export_lids, export_pids, offsets, num_bytes_per_value,
     pack_pids);
  // If we got this far, we succeeded.
}

#endif



template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraFunctions<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
importSquareMatrixFromImporter4(RCP<const CrsMatrix> inputMatrix, 
                               RCP<const Import> importer,
                                RCP<CrsMatrix> & outputMatrix) {  
  // Approach:  Modify the packing to only pack the rows we need, then use the IAFC comm & unpack pat

  // This will need to get Kokkos-ified later, but for now this is fine
  /*
  std::vector<GO> restrictedColGIDs;
  std::vector<LO> restrictedColGIDsBegin;
  communicateRowMap(outputRowMap, distributor, targetColGIDs, targetColGIDsBegin);
  */

  
  throw std::runtime_error("Implementation not yet complete");



}


/****************************************************************************************************/
/****************************************************************************************************/
/****************************************************************************************************/



template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraFunctions<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
importSquareMatrixFromImporter3(RCP<const CrsMatrix> inputMatrix, 
                               RCP<const Import> importer,
                               RCP<CrsMatrix> & outputMatrix)
{
    std::vector<GO> targetMapGIDs;
    std::vector<LO> targetMapGIDsBegin;
    // targetMapGIDs[indices] = globalIDs of outputRowMap for the i'th process
    //                          receiving matrix data
    //         indices = targetMapGIDsBegin[i]:targetMapGIDsBegin[i+1]-1;
    // Note: length of targetMapGIDsBegin is the number of processes receiving 
    //       matrix data plus 1
    communicateRowMap( importer->getTargetMap(),Teuchos::rcpFromRef(importer->getDistributor()) , targetMapGIDs, targetMapGIDsBegin);
    communicateMatrixData3(inputMatrix, importer, targetMapGIDs, 
                           targetMapGIDsBegin, outputMatrix);
}




template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraFunctions<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
communicateMatrixData3(RCP<const CrsMatrix> inputMatrix, 
                       RCP<const Tpetra::Import<LO,GO,NT> > importer,
                      const std::vector<GO> & targetMapGIDs, 
                      const std::vector<LO> & targetMapGIDsBegin,
                      RCP<CrsMatrix> & outputMatrix)
{
  RCP<const Map> rowMap = importer->getTargetMap();
  RCP<const Map> sourceMap = importer->getSourceMap();
  RCP<Tpetra::Distributor> distributor = Teuchos::rcpFromRef(importer->getDistributor());
  Teuchos::ArrayView<const LO> localRowsSend = importer->getExportLIDs();
  Teuchos::ArrayView<const LO> localRowsRecv = importer->getRemoteLIDs();
  Teuchos::ArrayView<const size_t> localRowsSendSize = distributor->getLengthsTo();
  Teuchos::ArrayView<const size_t> localRowsRecvSize = distributor->getLengthsFrom();

  // Combine sames and permutes for ownedRowGIDs
  // QUESTION: Does iSM handle permutes correctly?
  size_t numSames = importer->getNumSameIDs();
  size_t numPermutes = importer->getNumPermuteIDs();
  Teuchos::ArrayView<const LO> permutes = importer->getPermuteFromLIDs();
  std::vector<GO> ownedRowGIDs(numSames + numPermutes);
  for(size_t i=0; i <numSames; i++)
    ownedRowGIDs[i] = sourceMap->getGlobalElement(i);
  for(size_t i=0; i <numPermutes; i++)
    ownedRowGIDs[i] = numSames + sourceMap->getGlobalElement(permutes[i]);
  


  const size_t numSends = distributor->getNumSends();
  const size_t numRecvs = distributor->getNumReceives();
  TEUCHOS_TEST_FOR_EXCEPTION(numSends != (size_t)localRowsSendSize.size(),std::runtime_error,
                             "invalid size of localRowsSendSize");
  TEUCHOS_TEST_FOR_EXCEPTION(numRecvs != (size_t)localRowsRecvSize.size(),std::runtime_error,
                             "invalid size of localRowsRecvSize");

    // CMS: For each send, figure out how many nonzero entry columns of the local inputMatrix exist correspond
    // to rows on that process.
    const size_t numRowsSendTotal = localRowsSend.size();
    std::vector<size_t> count(numRowsSendTotal, 0);
    IndicesViewT indices;
    ValuesViewT values;
    RCP<const Map> columnMap = inputMatrix->getColMap();
    auto IGO = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    auto ILO = Teuchos::OrdinalTraits<LO>::invalid();
    RCP<const Teuchos::Comm<int> > serialComm = RCP( new Teuchos::SerialComm<int>() );
    for (size_t i=0, send_buff_start = 0; i<numSends; i++) {
        const GO* targetGIDs = &targetMapGIDs[targetMapGIDsBegin[i]];
        const size_t numTargetGIDs = targetMapGIDsBegin[i+1] - targetMapGIDsBegin[i];
        RCP<const Map> targetMap = 
            RCP( new Map(IGO, Teuchos::ArrayView<const GO>(targetGIDs, numTargetGIDs),
                         0, serialComm) );
        for (LO jj=0; jj<(LO)localRowsSendSize[i]; jj++) {
          LO j = send_buff_start + jj;
          const LO localRow = localRowsSend[j];
          inputMatrix->getLocalRowView(localRow, indices, values);
          const int sz = indices.size();
          for (int k=0; k<sz; k++) {
            const GO colGlobalID = columnMap->getGlobalElement(indices[k]);
            LO targetIndex = targetMap->getLocalElement(colGlobalID);
            if (targetIndex != ILO) count[j]++;
          }
        }
        send_buff_start +=localRowsSendSize[i];
    }
    /*
      std::cout << "local rows and send counts for myPID = " << myPID << std::endl;
      Teuchos::ArrayView<const int> procsTo = distributor->getProcsTo();
      for (size_t i=0; i<numSends; i++) {
      std::cout << "--- for sending to proc " << procsTo[i] << " ---" << std::endl;
      for (LO j=localRowsSendBegin[i]; j<localRowsSendBegin[i+1]; j++) {
      std::cout << localRowsSend[j] << " " << count[j] << std::endl;
      }
      }
    */
    // CMS: Pack the nonzero entries whose columns correspond to rows on that process
    // only.
    size_t sumCount = 0;
    for (size_t i=0; i<count.size(); i++) sumCount += count[i];
    std::vector<LO> localColIDs(sumCount);
    std::vector<SC> matrixValues(sumCount);
    sumCount = 0;
    for (size_t i=0, send_buff_start = 0; i<numSends; i++) {
        const GO* targetGIDs = &targetMapGIDs[targetMapGIDsBegin[i]];
        const size_t numTargetGIDs = targetMapGIDsBegin[i+1] - targetMapGIDsBegin[i];
        RCP<const Map> targetMap = 
            RCP( new Map(IGO, Teuchos::ArrayView<const GO>(targetGIDs, numTargetGIDs),
                         0, serialComm) );
        for (LO jj=0; jj<(LO)localRowsSendSize[i]; jj++) {
          LO j = send_buff_start + jj;
          const LO localRow = localRowsSend[j];
          inputMatrix->getLocalRowView(localRow, indices, values);
          const int sz = indices.size();
          for (int k=0; k<sz; k++) {
            const GO colGlobalID = columnMap->getGlobalElement(indices[k]);
            LO targetIndex = targetMap->getLocalElement(colGlobalID);
            if (targetIndex != ILO) {
              localColIDs[sumCount] = targetIndex;
              matrixValues[sumCount++] = values[k];
            }
          }
        }
        send_buff_start +=localRowsSendSize[i];
    }
    // CMS: Compute sizes and communicate the count data
    const size_t numRowsRecv = localRowsRecv.size();
    std::vector<size_t> recvRowNumNonzeros(numRowsRecv);
    distributor->doPostsAndWaits(Kokkos::View<const size_t*, Kokkos::HostSpace>(count.data(), count.size()),
                                 localRowsSendSize,
                                 Kokkos::View<size_t*, Kokkos::HostSpace>(recvRowNumNonzeros.data(), recvRowNumNonzeros.size()),
                                 localRowsRecvSize);

    /*
      const int myPID = rowMap->getComm()->getRank();
      Teuchos::ArrayView<const int> procsFrom = distributor->getProcsFrom();
      for (size_t i=0; i<numRecvs; i++) {
      std::cout << "row globalIDs and counts for proc " << myPID << " received from proc "
      << procsFrom[i] << std::endl;
      for (LO j=localRowsRecvBegin[i]; j<localRowsRecvBegin[i+1]; j++) {
      std::cout << rowMap->getGlobalElement(localRowsRecv[j]) << " "
      << recvRowNumNonzeros[j] << std::endl;
      }
      }
    */ 
    std::vector<size_t> sourceSize(numSends);
    for (size_t i=0,send_buff_start=0; i<numSends; i++) {
        size_t procCount(0);
        for (LO jj=0; jj<(LO)localRowsSendSize[i]; jj++) {
          LO j = send_buff_start + jj;
          procCount += count[j];
        }
        sourceSize[i] = procCount;
        send_buff_start +=localRowsSendSize[i];
    }
    std::vector<size_t> targetSize(numRecvs);
    size_t numTerms(0);
    for (size_t i=0,recv_buff_start=0; i<numRecvs; i++) {
        size_t procCount(0);
        for (LO jj=0; jj<(LO)localRowsRecvSize[i]; jj++) {
          LO j = recv_buff_start + jj;
          procCount += recvRowNumNonzeros[j];
        }
        targetSize[i] = procCount;
        numTerms += procCount;
        recv_buff_start +=localRowsRecvSize[i];
    }

    // CMS: Communicate the nonzeros
    std::vector<LO> columnsRecv(numTerms);
    std::vector<SC> valuesRecv(numTerms);
    distributor->doPostsAndWaits(Kokkos::View<const LO*, Kokkos::HostSpace>(localColIDs.data(), localColIDs.size()),
                                 Teuchos::ArrayView<const size_t>(sourceSize),
                                 Kokkos::View<LO*, Kokkos::HostSpace>(columnsRecv.data(), columnsRecv.size()),
                                 Teuchos::ArrayView<const size_t>(targetSize));
    
    using KSX = typename Kokkos::ArithTraits<SC>::val_type;
    const KSX* matrixValues_K = reinterpret_cast<const KSX*>(matrixValues.data());
    KSX* valuesRecv_K = reinterpret_cast<KSX*>(valuesRecv.data());
    const size_t sizeSend = matrixValues.size();
    const size_t sizeRecv = valuesRecv.size();
    distributor->doPostsAndWaits(Kokkos::View<const KSX*, Kokkos::HostSpace>(matrixValues_K, sizeSend),
                                 Teuchos::ArrayView<const size_t>(sourceSize),
                                 Kokkos::View<KSX*, Kokkos::HostSpace>(valuesRecv_K, sizeRecv),
                                 Teuchos::ArrayView<const size_t>(targetSize));
    RCP<const Map> rowMap1to1 = inputMatrix->getRowMap();
    const size_t numRows = rowMap->getLocalNumElements();
    std::vector<size_t> rowCount(numRows, 0);
    // rowCounts for owned rows
    RCP<const Map> colMapSource = inputMatrix->getColMap();
    for (size_t i=0; i<ownedRowGIDs.size(); i++) {
        const LO localRowSource = rowMap1to1->getLocalElement(ownedRowGIDs[i]);
        const LO localRowTarget = rowMap->getLocalElement(ownedRowGIDs[i]);
        TEUCHOS_TEST_FOR_EXCEPTION(localRowSource == ILO, std::runtime_error,"globalID not found in rowMap1to1");
        TEUCHOS_TEST_FOR_EXCEPTION(localRowTarget == ILO, std::runtime_error,"globalID not found in rowMap");
        inputMatrix->getLocalRowView(localRowSource, indices, values);
        const int sz = indices.size();
        for (int j=0; j<sz; j++) {
            const GO colGID = colMapSource->getGlobalElement(indices[j]);
            LO localColTarget = rowMap->getLocalElement(colGID);
            if (localColTarget != ILO) rowCount[localRowTarget]++;
        }
    }
    // rowCounts for received rows
    for (LO i=0; i<localRowsRecv.size(); i++) {
        const LO localRow = localRowsRecv[i];
        rowCount[localRow] = recvRowNumNonzeros[i];
    }
    /*
      const int myPID = rowMap->getComm()->getRank();
      std::cout << "rowGIDs and numEntries for myPID = " << myPID << std::endl;
      for (size_t i=0; i<numRows; i++) {
      std::cout << rowMap->getGlobalElement(i) << ": " << rowCount[i] << std::endl;
      }
    */
    // CMS: Unpack the nonzeros
    outputMatrix = 
        RCP( new CrsMatrix(rowMap, rowMap, Teuchos::ArrayView<const size_t>(rowCount)) );
    std::vector<LO> indicesVec(numRows);
    std::vector<SC> valuesVec(numRows);
    // insert owned constributions
    for (size_t i=0; i<ownedRowGIDs.size(); i++) {
        indicesVec.resize(0);
        valuesVec.resize(0);
        const LO localRowSource = rowMap1to1->getLocalElement(ownedRowGIDs[i]);
        const LO localRowTarget = rowMap->getLocalElement(ownedRowGIDs[i]);
        inputMatrix->getLocalRowView(localRowSource, indices, values);
        const int sz = indices.size();
        for (int j=0; j<sz; j++) {
            const GO colGID = colMapSource->getGlobalElement(indices[j]);
            LO localColTarget = rowMap->getLocalElement(colGID);
            if (localColTarget != ILO) {
                indicesVec.push_back(localColTarget);
                valuesVec.push_back(values[j]);
            }
        }
        outputMatrix->insertLocalValues(localRowTarget, 
                                        Teuchos::ArrayView<const LO>(indicesVec),
                                        Teuchos::ArrayView<const SC>(valuesVec));
    }
    // insert received contributions
    numTerms = 0;
    for (LO j=0; j<localRowsRecv.size(); j++) {
        const LO localRow = localRowsRecv[j];
        indicesVec.resize(recvRowNumNonzeros[j]);
        valuesVec.resize(recvRowNumNonzeros[j]);
        for (size_t k=0; k<recvRowNumNonzeros[j]; k++) {
            indicesVec[k] = columnsRecv[numTerms];
            valuesVec[k] = valuesRecv[numTerms++];
        }
        outputMatrix->insertLocalValues(localRow, 
                                        Teuchos::ArrayView<const LO>(indicesVec),
                                        Teuchos::ArrayView<const SC>(valuesVec));
    }
    outputMatrix->fillComplete(rowMap1to1, rowMap1to1);
    /*
      const int myPID = rowMap->getComm()->getRank();
      std::cout << "matrix entries for proc (global indices) = " << myPID << std::endl;
      for (size_t i=0; i<numRows; i++) {
      std::cout << rowMap->getGlobalElement(i) << ": ";
      outputMatrix->getLocalRowView(i, indices, values);
      for (int j=0; j<indices.size(); j++) {
      std::cout << rowMap->getGlobalElement(indices[j]) << "("
      << values[j] << ") ";
      }
      std::cout << std::endl;
      }
    */
}



/****************************************************************************************************/
/****************************************************************************************************/
/****************************************************************************************************/



template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraFunctions<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
importSquareMatrixFromImporter2(RCP<const CrsMatrix> inputMatrix, 
                               RCP<const Import> importer,
                               RCP<CrsMatrix> & outputMatrix)
{
  RCP<const Map> outputRowMap = importer->getTargetMap();
  RCP<const Map> sourceMap = importer->getSourceMap();
  RCP<Tpetra::Distributor> distributor = Teuchos::rcpFromRef(importer->getDistributor());
  Teuchos::ArrayView<const LO> localRowsSend = importer->getExportLIDs();
  Teuchos::ArrayView<const LO> localRowsRecv = importer->getRemoteLIDs();


  // Scansum to get the begin arrays
  //  size_t numSends = distributor->getNumSends();
  //  size_t numRecvs = distributor->getNumReceives();
  Teuchos::ArrayView<const size_t> lengthsTo = distributor->getLengthsTo();
  Teuchos::ArrayView<const size_t> lengthsFrom = distributor->getLengthsFrom();

  // Combine sames and permutes for ownedRowGIDs
  // QUESTION: Does iSM handle permutes correctly?
  size_t numSames = importer->getNumSameIDs();
  size_t numPermutes = importer->getNumPermuteIDs();
  Teuchos::ArrayView<const LO> permutes = importer->getPermuteFromLIDs();
  std::vector<GO> ownedRowGIDs(numSames + numPermutes);
  for(size_t i=0; i <numSames; i++)
    ownedRowGIDs[i] = sourceMap->getGlobalElement(i);
  for(size_t i=0; i <numPermutes; i++)
    ownedRowGIDs[i] = numSames + sourceMap->getGlobalElement(permutes[i]);
  


  //    std::vector<LO> localRowsSend, localRowsSendBegin, localRowsRecv, localRowsRecvBegin;
    // ownedRowGIDs = rows of inputMatrix contained in outputRowMap
    // localRowsSend[indicesSend] = localRows of sending process i corresponding to
    //                              localRowsRecv[indicesRecv], where
    //               indicesSend = localRowsSendBegin[i]:localRowsSendBegin[i+1]-1
    //               indicesRecv = localRowsRecvBegin[i]:localRowsRecvBegin[i+1]-1
    //    constructDistributor(inputMatrix, outputRowMap, distributor, ownedRowGIDs,
    //                         localRowsSend, localRowsSendBegin,
    //                         localRowsRecv, localRowsRecvBegin);
    std::vector<GO> targetMapGIDs;
    std::vector<LO> targetMapGIDsBegin;
    // targetMapGIDs[indices] = globalIDs of outputRowMap for the i'th process
    //                          receiving matrix data
    //         indices = targetMapGIDsBegin[i]:targetMapGIDsBegin[i+1]-1;
    // Note: length of targetMapGIDsBegin is the number of processes receiving 
    //       matrix data plus 1
    communicateRowMap(outputRowMap, distributor, targetMapGIDs, targetMapGIDsBegin);
    communicateMatrixData2(inputMatrix, outputRowMap, distributor, targetMapGIDs, 
                          targetMapGIDsBegin, ownedRowGIDs,
                          localRowsSend, lengthsTo, localRowsRecv,
                           lengthsFrom, outputMatrix);
}





template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
template<class LOVector, class STVector>
void TpetraFunctions<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
communicateMatrixData2(RCP<const CrsMatrix> inputMatrix, 
                      RCP<const Map> rowMap,
                      RCP<Tpetra::Distributor> distributor,
                      const std::vector<GO> & targetMapGIDs, 
                      const std::vector<LO> & targetMapGIDsBegin,
                      const std::vector<GO> & ownedRowGIDs,                       
                      const LOVector & localRowsSend,
                      const STVector & localRowsSendSize,
                      const LOVector & localRowsRecv,
                      const STVector & localRowsRecvSize,
                      RCP<CrsMatrix> & outputMatrix)
{
    const size_t numSends = distributor->getNumSends();
    const size_t numRecvs = distributor->getNumReceives();
    TEUCHOS_TEST_FOR_EXCEPTION(numSends != (size_t)localRowsSendSize.size(),std::runtime_error,
                    "invalid size of localRowsSendSize");
    TEUCHOS_TEST_FOR_EXCEPTION(numRecvs != (size_t)localRowsRecvSize.size(),std::runtime_error,
                    "invalid size of localRowsRecvSize");

    // CMS: For each send, figure out how many nonzero entry columns of the local inputMatrix exist correspond
    // to rows on that process.
    const size_t numRowsSendTotal = localRowsSend.size();
    std::vector<size_t> count(numRowsSendTotal, 0);
    IndicesViewT indices;
    ValuesViewT values;
    RCP<const Map> columnMap = inputMatrix->getColMap();
    auto IGO = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    auto ILO = Teuchos::OrdinalTraits<LO>::invalid();
    RCP<const Teuchos::Comm<int> > serialComm = RCP( new Teuchos::SerialComm<int>() );
    for (size_t i=0, send_buff_start = 0; i<numSends; i++) {
        const GO* targetGIDs = &targetMapGIDs[targetMapGIDsBegin[i]];
        const size_t numTargetGIDs = targetMapGIDsBegin[i+1] - targetMapGIDsBegin[i];
        RCP<const Map> targetMap = 
            RCP( new Map(IGO, Teuchos::ArrayView<const GO>(targetGIDs, numTargetGIDs),
                         0, serialComm) );
        for (LO jj=0; jj<(LO)localRowsSendSize[i]; jj++) {
          LO j = send_buff_start + jj;
          const LO localRow = localRowsSend[j];
          inputMatrix->getLocalRowView(localRow, indices, values);
          const int sz = indices.size();
          for (int k=0; k<sz; k++) {
            const GO colGlobalID = columnMap->getGlobalElement(indices[k]);
            LO targetIndex = targetMap->getLocalElement(colGlobalID);
            if (targetIndex != ILO) count[j]++;
          }
        }
        send_buff_start +=localRowsSendSize[i];
    }
    /*
      std::cout << "local rows and send counts for myPID = " << myPID << std::endl;
      Teuchos::ArrayView<const int> procsTo = distributor->getProcsTo();
      for (size_t i=0; i<numSends; i++) {
      std::cout << "--- for sending to proc " << procsTo[i] << " ---" << std::endl;
      for (LO j=localRowsSendBegin[i]; j<localRowsSendBegin[i+1]; j++) {
      std::cout << localRowsSend[j] << " " << count[j] << std::endl;
      }
      }
    */
    // CMS: Pack the nonzero entries whose columns correspond to rows on that process
    // only.
    size_t sumCount = 0;
    for (size_t i=0; i<count.size(); i++) sumCount += count[i];
    std::vector<LO> localColIDs(sumCount);
    std::vector<SC> matrixValues(sumCount);
    sumCount = 0;
    for (size_t i=0, send_buff_start = 0; i<numSends; i++) {
        const GO* targetGIDs = &targetMapGIDs[targetMapGIDsBegin[i]];
        const size_t numTargetGIDs = targetMapGIDsBegin[i+1] - targetMapGIDsBegin[i];
        RCP<const Map> targetMap = 
            RCP( new Map(IGO, Teuchos::ArrayView<const GO>(targetGIDs, numTargetGIDs),
                         0, serialComm) );
        for (LO jj=0; jj<(LO)localRowsSendSize[i]; jj++) {
          LO j = send_buff_start + jj;
          const LO localRow = localRowsSend[j];
          inputMatrix->getLocalRowView(localRow, indices, values);
          const int sz = indices.size();
          for (int k=0; k<sz; k++) {
            const GO colGlobalID = columnMap->getGlobalElement(indices[k]);
            LO targetIndex = targetMap->getLocalElement(colGlobalID);
            if (targetIndex != ILO) {
              localColIDs[sumCount] = targetIndex;
              matrixValues[sumCount++] = values[k];
            }
          }
        }
        send_buff_start +=localRowsSendSize[i];
    }
    // CMS: Compute sizes and communicate the count data
    const size_t numRowsRecv = localRowsRecv.size();
    std::vector<size_t> recvRowNumNonzeros(numRowsRecv);
    distributor->doPostsAndWaits(Kokkos::View<const size_t*, Kokkos::HostSpace>(count.data(), count.size()),
                                 localRowsSendSize,
                                 Kokkos::View<size_t*, Kokkos::HostSpace>(recvRowNumNonzeros.data(), recvRowNumNonzeros.size()),
                                 localRowsRecvSize);

    /*
      const int myPID = rowMap->getComm()->getRank();
      Teuchos::ArrayView<const int> procsFrom = distributor->getProcsFrom();
      for (size_t i=0; i<numRecvs; i++) {
      std::cout << "row globalIDs and counts for proc " << myPID << " received from proc "
      << procsFrom[i] << std::endl;
      for (LO j=localRowsRecvBegin[i]; j<localRowsRecvBegin[i+1]; j++) {
      std::cout << rowMap->getGlobalElement(localRowsRecv[j]) << " "
      << recvRowNumNonzeros[j] << std::endl;
      }
      }
    */ 
    std::vector<size_t> sourceSize(numSends);
    for (size_t i=0,send_buff_start=0; i<numSends; i++) {
        size_t procCount(0);
        for (LO jj=0; jj<(LO)localRowsSendSize[i]; jj++) {
          LO j = send_buff_start + jj;
          procCount += count[j];
        }
        sourceSize[i] = procCount;
        send_buff_start +=localRowsSendSize[i];
    }
    std::vector<size_t> targetSize(numRecvs);
    size_t numTerms(0);
    for (size_t i=0,recv_buff_start=0; i<numRecvs; i++) {
        size_t procCount(0);
        for (LO jj=0; jj<(LO)localRowsRecvSize[i]; jj++) {
          LO j = recv_buff_start + jj;
          procCount += recvRowNumNonzeros[j];
        }
        targetSize[i] = procCount;
        numTerms += procCount;
        recv_buff_start +=localRowsRecvSize[i];
    }

    // CMS: Communicate the nonzeros
    std::vector<LO> columnsRecv(numTerms);
    std::vector<SC> valuesRecv(numTerms);
    distributor->doPostsAndWaits(Kokkos::View<const LO*, Kokkos::HostSpace>(localColIDs.data(), localColIDs.size()),
                                 Teuchos::ArrayView<const size_t>(sourceSize),
                                 Kokkos::View<LO*, Kokkos::HostSpace>(columnsRecv.data(), columnsRecv.size()),
                                 Teuchos::ArrayView<const size_t>(targetSize));
    
    using KSX = typename Kokkos::ArithTraits<SC>::val_type;
    const KSX* matrixValues_K = reinterpret_cast<const KSX*>(matrixValues.data());
    KSX* valuesRecv_K = reinterpret_cast<KSX*>(valuesRecv.data());
    const size_t sizeSend = matrixValues.size();
    const size_t sizeRecv = valuesRecv.size();
    distributor->doPostsAndWaits(Kokkos::View<const KSX*, Kokkos::HostSpace>(matrixValues_K, sizeSend),
                                 Teuchos::ArrayView<const size_t>(sourceSize),
                                 Kokkos::View<KSX*, Kokkos::HostSpace>(valuesRecv_K, sizeRecv),
                                 Teuchos::ArrayView<const size_t>(targetSize));
    RCP<const Map> rowMap1to1 = inputMatrix->getRowMap();
    const size_t numRows = rowMap->getLocalNumElements();
    std::vector<size_t> rowCount(numRows, 0);
    // rowCounts for owned rows
    RCP<const Map> colMapSource = inputMatrix->getColMap();
    for (size_t i=0; i<ownedRowGIDs.size(); i++) {
        const LO localRowSource = rowMap1to1->getLocalElement(ownedRowGIDs[i]);
        const LO localRowTarget = rowMap->getLocalElement(ownedRowGIDs[i]);
        TEUCHOS_TEST_FOR_EXCEPTION(localRowSource == ILO, std::runtime_error,"globalID not found in rowMap1to1");
        TEUCHOS_TEST_FOR_EXCEPTION(localRowTarget == ILO, std::runtime_error,"globalID not found in rowMap");
        inputMatrix->getLocalRowView(localRowSource, indices, values);
        const int sz = indices.size();
        for (int j=0; j<sz; j++) {
            const GO colGID = colMapSource->getGlobalElement(indices[j]);
            LO localColTarget = rowMap->getLocalElement(colGID);
            if (localColTarget != ILO) rowCount[localRowTarget]++;
        }
    }
    // rowCounts for received rows
    for (LO i=0; i<localRowsRecv.size(); i++) {
        const LO localRow = localRowsRecv[i];
        rowCount[localRow] = recvRowNumNonzeros[i];
    }
    /*
      const int myPID = rowMap->getComm()->getRank();
      std::cout << "rowGIDs and numEntries for myPID = " << myPID << std::endl;
      for (size_t i=0; i<numRows; i++) {
      std::cout << rowMap->getGlobalElement(i) << ": " << rowCount[i] << std::endl;
      }
    */
    // CMS: Unpack the nonzeros
    outputMatrix = 
        RCP( new CrsMatrix(rowMap, rowMap, Teuchos::ArrayView<const size_t>(rowCount)) );
    std::vector<LO> indicesVec(numRows);
    std::vector<SC> valuesVec(numRows);
    // insert owned constributions
    for (size_t i=0; i<ownedRowGIDs.size(); i++) {
        indicesVec.resize(0);
        valuesVec.resize(0);
        const LO localRowSource = rowMap1to1->getLocalElement(ownedRowGIDs[i]);
        const LO localRowTarget = rowMap->getLocalElement(ownedRowGIDs[i]);
        inputMatrix->getLocalRowView(localRowSource, indices, values);
        const int sz = indices.size();
        for (int j=0; j<sz; j++) {
            const GO colGID = colMapSource->getGlobalElement(indices[j]);
            LO localColTarget = rowMap->getLocalElement(colGID);
            if (localColTarget != ILO) {
                indicesVec.push_back(localColTarget);
                valuesVec.push_back(values[j]);
            }
        }
        outputMatrix->insertLocalValues(localRowTarget, 
                                        Teuchos::ArrayView<const LO>(indicesVec),
                                        Teuchos::ArrayView<const SC>(valuesVec));
    }
    // insert received contributions
    numTerms = 0;
    for (LO j=0; j<localRowsRecv.size(); j++) {
        const LO localRow = localRowsRecv[j];
        indicesVec.resize(recvRowNumNonzeros[j]);
        valuesVec.resize(recvRowNumNonzeros[j]);
        for (size_t k=0; k<recvRowNumNonzeros[j]; k++) {
            indicesVec[k] = columnsRecv[numTerms];
            valuesVec[k] = valuesRecv[numTerms++];
        }
        outputMatrix->insertLocalValues(localRow, 
                                        Teuchos::ArrayView<const LO>(indicesVec),
                                        Teuchos::ArrayView<const SC>(valuesVec));
    }
    outputMatrix->fillComplete(rowMap1to1, rowMap1to1);
    /*
      const int myPID = rowMap->getComm()->getRank();
      std::cout << "matrix entries for proc (global indices) = " << myPID << std::endl;
      for (size_t i=0; i<numRows; i++) {
      std::cout << rowMap->getGlobalElement(i) << ": ";
      outputMatrix->getLocalRowView(i, indices, values);
      for (int j=0; j<indices.size(); j++) {
      std::cout << rowMap->getGlobalElement(indices[j]) << "("
      << values[j] << ") ";
      }
      std::cout << std::endl;
      }
    */
}



/****************************************************************************************************/
/****************************************************************************************************/
/****************************************************************************************************/


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraFunctions<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
communicateRowMap(RCP<const Map> rowMap, 
                  RCP<Tpetra::Distributor> distributor, 
                  std::vector<GO> & rowMapGIDs, 
                  std::vector<LO> & rowMapGIDsBegin)
{
    size_t numSends = distributor->getNumSends();
    size_t numRecvs = distributor->getNumReceives();
    size_t numRows = rowMap->getLocalNumElements();
    std::vector<size_t> targetValues(numSends);
    std::vector<size_t> targetSizes(numSends, 1);
    std::vector<size_t> sourceValues(numRecvs, numRows);
    std::vector<size_t> sourceSizes(numRecvs, 1); 

    // CMS: Reverse sends the # of rows on each proc to all neighbors (I think)
    distributor->doReversePostsAndWaits(Kokkos::View<const size_t*, Kokkos::HostSpace>(sourceValues.data(), sourceValues.size()),
                                        1,
                                        Kokkos::View<size_t*, Kokkos::HostSpace>(targetValues.data(), targetValues.size()));
    // CMS: Compute the total number of rows on reverse neighbors of this rank
    int numTerms(0);
    for (size_t i=0; i<targetValues.size(); i++) numTerms += targetValues[i];
    rowMapGIDs.resize(numTerms);

    // CMS: For each recv, record each of my GIDs
    std::vector<GO> globalIDsSource(numRecvs*numRows);
    numTerms = 0;
    for (size_t i=0; i<numRecvs; i++) {
        for (size_t j=0; j<numRows; j++) {
            globalIDsSource[numTerms++] = rowMap->getGlobalElement(j);
        }
    }

    // CMS: Reverse all of the GIDs owned by the neighboring proc
    distributor->doReversePostsAndWaits(Kokkos::View<const GO*, Kokkos::HostSpace>(globalIDsSource.data(), globalIDsSource.size()),
                                        Teuchos::ArrayView<const size_t>(sourceValues),
                                        Kokkos::View<GO*, Kokkos::HostSpace>(rowMapGIDs.data(), rowMapGIDs.size()),
                                        Teuchos::ArrayView<const size_t>(targetValues));

    // CMS: Note the row beginnings of each reverse neighboring row
    rowMapGIDsBegin.resize(numSends+1, 0);
    for (size_t i=0; i<numSends; i++) {
        rowMapGIDsBegin[i+1] = rowMapGIDsBegin[i] + targetValues[i];
    }
    /*
      const int myPID = rowMap->getComm()->getRank();
      std::cout << "map globalIDs for myPID = " << myPID << std::endl;
      for (size_t i=0; i<numSends; i++) {
      for (LO j=rowMapGIDsBegin[i]; j<rowMapGIDsBegin[i+1]; j++) {
      std::cout << rowMapGIDs[j] << " ";
      }
      std::cout << std::endl;
      }
    */
}






template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
template<class LOVector>
void TpetraFunctions<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
communicateMatrixData(RCP<const CrsMatrix> inputMatrix, 
                      RCP<const Map> rowMap,
                      RCP<Tpetra::Distributor> distributor,
                      const std::vector<GO> & targetMapGIDs, 
                      const std::vector<LO> & targetMapGIDsBegin,
                      const std::vector<GO> & ownedRowGIDs,
                      const LOVector & localRowsSend,
                      const std::vector<LO> & localRowsSendBegin,
                      const LOVector & localRowsRecv,
                      const std::vector<LO> & localRowsRecvBegin,
                      RCP<CrsMatrix> & outputMatrix)
{
    const size_t numSends = distributor->getNumSends();
    const size_t numRecvs = distributor->getNumReceives();
    TEUCHOS_TEST_FOR_EXCEPTION(numSends != localRowsSendBegin.size()-1,std::runtime_error,
                    "invalid size of localRowsSendBegin");
    TEUCHOS_TEST_FOR_EXCEPTION(numRecvs != localRowsRecvBegin.size()-1,std::runtime_error,
                    "invalid size of localRowsRecvBegin");



    const size_t numRowsSendTotal = localRowsSend.size();
    std::vector<size_t> count(numRowsSendTotal, 0);
    IndicesViewT indices;
    ValuesViewT values;
    RCP<const Map> columnMap = inputMatrix->getColMap();
    auto IGO = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    auto ILO = Teuchos::OrdinalTraits<LO>::invalid();
    RCP<const Teuchos::Comm<int> > serialComm = RCP( new Teuchos::SerialComm<int>() );
    for (size_t i=0; i<numSends; i++) {
        const GO* targetGIDs = &targetMapGIDs[targetMapGIDsBegin[i]];
        const size_t numTargetGIDs = targetMapGIDsBegin[i+1] - targetMapGIDsBegin[i];
        RCP<const Map> targetMap = 
            RCP( new Map(IGO, Teuchos::ArrayView<const GO>(targetGIDs, numTargetGIDs),
                         0, serialComm) );
        for (LO j=localRowsSendBegin[i]; j<localRowsSendBegin[i+1]; j++) {
            const LO localRow = localRowsSend[j];
            inputMatrix->getLocalRowView(localRow, indices, values);
            const int sz = indices.size();
            for (int k=0; k<sz; k++) {
                const GO colGlobalID = columnMap->getGlobalElement(indices[k]);
                LO targetIndex = targetMap->getLocalElement(colGlobalID);
                if (targetIndex != ILO) count[j]++;
            }

        }

    }
    /*
      std::cout << "local rows and send counts for myPID = " << myPID << std::endl;
      Teuchos::ArrayView<const int> procsTo = distributor->getProcsTo();
      for (size_t i=0; i<numSends; i++) {
      std::cout << "--- for sending to proc " << procsTo[i] << " ---" << std::endl;
      for (LO j=localRowsSendBegin[i]; j<localRowsSendBegin[i+1]; j++) {
      std::cout << localRowsSend[j] << " " << count[j] << std::endl;
      }
      }
    */


    size_t sumCount = 0;
    for (size_t i=0; i<count.size(); i++) sumCount += count[i];
    std::vector<LO> localColIDs(sumCount);
    std::vector<SC> matrixValues(sumCount);
    sumCount = 0;
    for (size_t i=0; i<numSends; i++) {
        const GO* targetGIDs = &targetMapGIDs[targetMapGIDsBegin[i]];
        const size_t numTargetGIDs = targetMapGIDsBegin[i+1] - targetMapGIDsBegin[i];
        RCP<const Map> targetMap = 
            RCP( new Map(IGO, Teuchos::ArrayView<const GO>(targetGIDs, numTargetGIDs),
                         0, serialComm) );
        for (LO j=localRowsSendBegin[i]; j<localRowsSendBegin[i+1]; j++) {
            const LO localRow = localRowsSend[j];
            inputMatrix->getLocalRowView(localRow, indices, values);
            const int sz = indices.size();
            for (int k=0; k<sz; k++) {
                const GO colGlobalID = columnMap->getGlobalElement(indices[k]);
                LO targetIndex = targetMap->getLocalElement(colGlobalID);
                if (targetIndex != ILO) {
                    localColIDs[sumCount] = targetIndex;
                    matrixValues[sumCount++] = values[k];
                }
            }

        }

    }
    const size_t numRowsRecv = localRowsRecvBegin[numRecvs];
    std::vector<size_t> recvSize(numRecvs), recvRowNumNonzeros(numRowsRecv);
    for (size_t i=0; i<numRecvs; i++)  {
        recvSize[i] = localRowsRecvBegin[i+1] - localRowsRecvBegin[i];
    }
    std::vector<size_t> sendSize(numSends);
    for (size_t i=0; i<numSends; i++) {
        sendSize[i] = localRowsSendBegin[i+1] - localRowsSendBegin[i];
    }
    distributor->doPostsAndWaits(Kokkos::View<const size_t*, Kokkos::HostSpace>(count.data(), count.size()),
                                 Teuchos::ArrayView<const size_t>(sendSize),
                                 Kokkos::View<size_t*, Kokkos::HostSpace>(recvRowNumNonzeros.data(), recvRowNumNonzeros.size()),
                                 Teuchos::ArrayView<const size_t>(recvSize));

    /*
      const int myPID = rowMap->getComm()->getRank();
      Teuchos::ArrayView<const int> procsFrom = distributor->getProcsFrom();
      for (size_t i=0; i<numRecvs; i++) {
      std::cout << "row globalIDs and counts for proc " << myPID << " received from proc "
      << procsFrom[i] << std::endl;
      for (LO j=localRowsRecvBegin[i]; j<localRowsRecvBegin[i+1]; j++) {
      std::cout << rowMap->getGlobalElement(localRowsRecv[j]) << " "
      << recvRowNumNonzeros[j] << std::endl;
      }
      }
    */ 
    std::vector<size_t> sourceSize(numSends);
    for (size_t i=0; i<numSends; i++) {
        size_t procCount(0);
        for (LO j=localRowsSendBegin[i]; j<localRowsSendBegin[i+1]; j++) {
            procCount += count[j];

        }
        sourceSize[i] = procCount;

    }
    std::vector<size_t> targetSize(numRecvs);
    size_t numTerms(0);
    for (size_t i=0; i<numRecvs; i++) {
        size_t procCount(0);
        for (LO j=localRowsRecvBegin[i]; j<localRowsRecvBegin[i+1]; j++) {
            procCount += recvRowNumNonzeros[j];

        }
        targetSize[i] = procCount;
        numTerms += procCount;

    }


    std::vector<LO> columnsRecv(numTerms);
    std::vector<SC> valuesRecv(numTerms);
    distributor->doPostsAndWaits(Kokkos::View<const LO*, Kokkos::HostSpace>(localColIDs.data(), localColIDs.size()),
                                 Teuchos::ArrayView<const size_t>(sourceSize),
                                 Kokkos::View<LO*, Kokkos::HostSpace>(columnsRecv.data(), columnsRecv.size()),
                                 Teuchos::ArrayView<const size_t>(targetSize));
    
    using KSX = typename Kokkos::ArithTraits<SC>::val_type;
    const KSX* matrixValues_K = reinterpret_cast<const KSX*>(matrixValues.data());
    KSX* valuesRecv_K = reinterpret_cast<KSX*>(valuesRecv.data());
    const size_t sizeSend = matrixValues.size();
    const size_t sizeRecv = valuesRecv.size();
    distributor->doPostsAndWaits(Kokkos::View<const KSX*, Kokkos::HostSpace>(matrixValues_K, sizeSend),
                                 Teuchos::ArrayView<const size_t>(sourceSize),
                                 Kokkos::View<KSX*, Kokkos::HostSpace>(valuesRecv_K, sizeRecv),
                                 Teuchos::ArrayView<const size_t>(targetSize));
    RCP<const Map> rowMap1to1 = inputMatrix->getRowMap();
    const size_t numRows = rowMap->getLocalNumElements();
    std::vector<size_t> rowCount(numRows, 0);
    // rowCounts for owned rows
    RCP<const Map> colMapSource = inputMatrix->getColMap();
    for (size_t i=0; i<ownedRowGIDs.size(); i++) {
        const LO localRowSource = rowMap1to1->getLocalElement(ownedRowGIDs[i]);
        const LO localRowTarget = rowMap->getLocalElement(ownedRowGIDs[i]);
        TEUCHOS_TEST_FOR_EXCEPTION(localRowSource == ILO, std::runtime_error,"globalID not found in rowMap1to1");
        TEUCHOS_TEST_FOR_EXCEPTION(localRowTarget == ILO, std::runtime_error,"globalID not found in rowMap");
        inputMatrix->getLocalRowView(localRowSource, indices, values);
        const int sz = indices.size();
        for (int j=0; j<sz; j++) {
            const GO colGID = colMapSource->getGlobalElement(indices[j]);
            LO localColTarget = rowMap->getLocalElement(colGID);
            if (localColTarget != ILO) rowCount[localRowTarget]++;
        }
    }
    // rowCounts for received rows
    for (LO i=0; i<localRowsRecvBegin[numRecvs]; i++) {
        const LO localRow = localRowsRecv[i];
        rowCount[localRow] = recvRowNumNonzeros[i];
    }
    /*
      const int myPID = rowMap->getComm()->getRank();
      std::cout << "rowGIDs and numEntries for myPID = " << myPID << std::endl;
      for (size_t i=0; i<numRows; i++) {
      std::cout << rowMap->getGlobalElement(i) << ": " << rowCount[i] << std::endl;
      }
    */

    outputMatrix = 
        RCP( new CrsMatrix(rowMap, rowMap, Teuchos::ArrayView<const size_t>(rowCount)) );
    std::vector<LO> indicesVec(numRows);
    std::vector<SC> valuesVec(numRows);
    // insert owned constributions
    for (size_t i=0; i<ownedRowGIDs.size(); i++) {
        indicesVec.resize(0);
        valuesVec.resize(0);
        const LO localRowSource = rowMap1to1->getLocalElement(ownedRowGIDs[i]);
        const LO localRowTarget = rowMap->getLocalElement(ownedRowGIDs[i]);
        inputMatrix->getLocalRowView(localRowSource, indices, values);
        const int sz = indices.size();
        for (int j=0; j<sz; j++) {
            const GO colGID = colMapSource->getGlobalElement(indices[j]);
            LO localColTarget = rowMap->getLocalElement(colGID);
            if (localColTarget != ILO) {
                indicesVec.push_back(localColTarget);
                valuesVec.push_back(values[j]);
            }
        }
        outputMatrix->insertLocalValues(localRowTarget, 
                                        Teuchos::ArrayView<const LO>(indicesVec),
                                        Teuchos::ArrayView<const SC>(valuesVec));
    }
    // insert received contributions
    numTerms = 0;
    for (LO j=0; j<localRowsRecvBegin[numRecvs]; j++) {
        const LO localRow = localRowsRecv[j];
        indicesVec.resize(recvRowNumNonzeros[j]);
        valuesVec.resize(recvRowNumNonzeros[j]);
        for (size_t k=0; k<recvRowNumNonzeros[j]; k++) {
            indicesVec[k] = columnsRecv[numTerms];
            valuesVec[k] = valuesRecv[numTerms++];
        }
        outputMatrix->insertLocalValues(localRow, 
                                        Teuchos::ArrayView<const LO>(indicesVec),
                                        Teuchos::ArrayView<const SC>(valuesVec));
    }
    outputMatrix->fillComplete(rowMap1to1, rowMap1to1);
    /*
      const int myPID = rowMap->getComm()->getRank();
      std::cout << "matrix entries for proc (global indices) = " << myPID << std::endl;
      for (size_t i=0; i<numRows; i++) {
      std::cout << rowMap->getGlobalElement(i) << ": ";
      outputMatrix->getLocalRowView(i, indices, values);
      for (int j=0; j<indices.size(); j++) {
      std::cout << rowMap->getGlobalElement(indices[j]) << "("
      << values[j] << ") ";
      }
      std::cout << std::endl;
      }
    */
}









template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraFunctions<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
constructDistributor(RCP<const CrsMatrix> inputMatrix, 
                     RCP<const Map> rowMap, 
                     RCP<Tpetra::Distributor> & distributor,
                     std::vector<GO> & ownedRowGIDs,
                     std::vector<LO> & localRowsSend,
                     std::vector<LO> & localRowsSendBegin,
                     std::vector<LO> & localRowsRecv,
                     std::vector<LO> & localRowsRecvBegin)
{
    const LO numRows = rowMap->getLocalNumElements();
    std::vector<GO> globalIDs(numRows);
    for (LO i=0; i<numRows; i++) globalIDs[i] = rowMap->getGlobalElement(i);
    std::vector<int> remotePIDs(numRows);
    std::vector<LO> remoteLocalRows(numRows);
    RCP<const Map> rowMap1to1 = inputMatrix->getRowMap();
    rowMap1to1->getRemoteIndexList(Teuchos::ArrayView<GO>(globalIDs),
                                   Teuchos::ArrayView<int>(remotePIDs),
                                   Teuchos::ArrayView<LO>(remoteLocalRows));
    std::vector<int> offProcessorPIDs(numRows);
    const int myPID = rowMap->getComm()->getRank();
    size_t numOffProcessorRows(0), numOnProcessorRows(0);
    ownedRowGIDs.resize(numRows);
    for (LO i=0; i<numRows; i++) {
        if (remotePIDs[i] != myPID) {
            globalIDs[numOffProcessorRows] = globalIDs[i];
            remotePIDs[numOffProcessorRows] = remotePIDs[i];
            remoteLocalRows[numOffProcessorRows++] = remoteLocalRows[i];
        }
        else {
            ownedRowGIDs[numOnProcessorRows++] = globalIDs[i];
        }
    }
    remotePIDs.resize(numOffProcessorRows);
    remoteLocalRows.resize(numOffProcessorRows);
    ownedRowGIDs.resize(numOnProcessorRows);
    std::vector<int> recvPIDs;
    getUniqueEntries(remotePIDs, recvPIDs);
    size_t numRecvs = recvPIDs.size();
    std::map<int,int> offProcessorMap;
    for (size_t i=0; i<numRecvs; i++) {
        offProcessorMap.emplace(recvPIDs[i], i);
    }
    std::vector<GO> recvGIDs(numRecvs);
    std::vector<size_t> count(numRecvs, 0);
    for (size_t i=0; i<numOffProcessorRows; i++) {
        auto iter = offProcessorMap.find(remotePIDs[i]);
        recvGIDs[iter->second] = globalIDs[i];
        count[iter->second]++;
    }
    localRowsRecvBegin.resize(numRecvs+1, 0);
    for (size_t i=0; i<numRecvs; i++) {
        localRowsRecvBegin[i+1] = localRowsRecvBegin[i] + count[i];
        count[i] = 0;
    }
    localRowsRecv.resize(numOffProcessorRows);
    for (size_t i=0; i<numOffProcessorRows; i++) {
        auto iter = offProcessorMap.find(remotePIDs[i]);
        const int index = localRowsRecvBegin[iter->second] + count[iter->second];
        localRowsRecv[index] = remoteLocalRows[i];
        count[iter->second]++;
    }
    /*
      for (size_t i=0; i<numRecvs; i++) {
      std::cout << "proc " << myPID << " local rows to be received from proc " << recvPIDs[i]
      << std::endl;
      for (int j=localRowsRecvBegin[i]; j<localRowsRecvBegin[i+1]; j++) {
      std::cout << localRowsRecv[j] << " ";
      }
      std::cout << std::endl;
      }
    */
    Teuchos::Array<GO> sendGIDs;
    Teuchos::Array<int> sendPIDs;
    distributor = RCP( new Tpetra::Distributor(rowMap->getComm()) );
    distributor->createFromRecvs(Teuchos::ArrayView<const GO>(recvGIDs),
                                 Teuchos::ArrayView<const int>(recvPIDs),
                                 sendGIDs, sendPIDs);
    TEUCHOS_TEST_FOR_EXCEPTION(distributor->hasSelfMessage() == true,std::runtime_error,
                               "distributor hasSelfMessage error");
    TEUCHOS_TEST_FOR_EXCEPTION(distributor->getNumReceives() != numRecvs,std::runtime_error,
                               "inconsistent numRecvs");
    /*
      Teuchos::ArrayView<const int> procsTo = distributor->getProcsTo();
      Teuchos::ArrayView<const int> procsFrom = distributor->getProcsFrom();
      Teuchos::ArrayView<const size_t> lengthsFrom = distributor->getLengthsFrom();
      Teuchos::ArrayView<const size_t> lengthsTo = distributor->getLengthsTo();
      std::cout << "procsTo for myPID = " << myPID << ": ";
      for (int i=0; i<procsTo.size(); i++) std::cout << procsTo[i] << " ";
      std::cout << std::endl;
      std::cout << "procsFrom for myPID = " << myPID << ": ";
      for (int i=0; i<procsFrom.size(); i++) std::cout << procsFrom[i] << " ";
      std::cout << std::endl;
      std::cout << "lengthsFrom for myPID = " << myPID << ": ";
      for (int i=0; i<lengthsFrom.size(); i++) std::cout << lengthsFrom[i] << " ";
      std::cout << std::endl;
      std::cout << "lengthsTo for myPID = " << myPID << ": ";
      for (int i=0; i<lengthsTo.size(); i++) std::cout << lengthsTo[i] << " ";
      std::cout << std::endl;
    */
    size_t numSends = distributor->getNumSends();
    std::vector<size_t> targetValues(numSends);
    std::vector<size_t> targetSizes(numSends, 1);
    std::vector<size_t> sourceSizes(numRecvs, 1); 
    distributor->doReversePostsAndWaits(Kokkos::View<const size_t*, Kokkos::HostSpace>(count.data(), count.size()),
                                        1,
                                        Kokkos::View<size_t*, Kokkos::HostSpace>(targetValues.data(), targetValues.size()));
    localRowsSendBegin.resize(numSends+1, 0);
    for (size_t i=0; i<numSends; i++) {
        localRowsSendBegin[i+1] = localRowsSendBegin[i] + targetValues[i];
    }
    int numTerms = localRowsSendBegin[numSends];
    localRowsSend.resize(numTerms);
    distributor->doReversePostsAndWaits(Kokkos::View<const int*, Kokkos::HostSpace>(localRowsRecv.data(), localRowsRecv.size()),
                                        Teuchos::ArrayView<const size_t>(count),
                                        Kokkos::View<int*, Kokkos::HostSpace>(localRowsSend.data(), localRowsSend.size()),
                                        Teuchos::ArrayView<const size_t>(targetValues));
    /*
      Teuchos::ArrayView<const int> procsTo = distributor->getProcsTo();
      for (size_t i=0; i<numSends; i++) {
      std::cout << "proc " << myPID << " globalIDs for rows of matrix to be sent to proc " 
      << procsTo[i] << std::endl;
      for (int j=localRowsSendBegin[i]; j<localRowsSendBegin[i+1]; j++) {
      std::cout << rowMap1to1->getGlobalElement(localRowsSend[j]) << " ";
      }
      std::cout << std::endl;
      }
    */
    // switch localRowsRecv to on-processor rather than off-processor localRows
    for (size_t i=0; i<numRecvs; i++) count[i] = 0;
    for (size_t i=0; i<numOffProcessorRows; i++) {
        auto iter = offProcessorMap.find(remotePIDs[i]);
        const int index = localRowsRecvBegin[iter->second] + count[iter->second];
        localRowsRecv[index] = rowMap->getLocalElement(globalIDs[i]);
        count[iter->second]++;
    }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraFunctions<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getUniqueEntries(const std::vector<int> & vector, 
                 std::vector<int> & vectorUnique)
{
    vectorUnique = vector;
    std::sort(vectorUnique.begin(), vectorUnique.end());
    auto iter = std::unique(vectorUnique.begin(), vectorUnique.end());
    vectorUnique.erase(iter, vectorUnique.end());
}


#endif
