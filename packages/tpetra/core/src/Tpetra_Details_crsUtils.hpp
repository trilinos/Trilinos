// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_CRSUTILS_HPP
#define TPETRA_DETAILS_CRSUTILS_HPP
#include <numeric>
#include <type_traits>

#include "TpetraCore_config.h"
#include "Kokkos_Core.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_CrsPadding.hpp"
#include "Tpetra_Details_WrappedDualView.hpp"
#include <iostream>
#include <memory>
#include <unordered_map>

/// \file Tpetra_Details_crsUtils.hpp
/// \brief Functions for manipulating CRS arrays
/// \warning This file, and its contents, are implementation details
///   of Tpetra.  The file itself or its contents may disappear or
///   change at any time.

namespace Tpetra {
namespace Details {

namespace impl {

template<class ViewType>
ViewType
make_uninitialized_view(
  const std::string& name,
  const size_t size,
  const bool verbose,
  const std::string* const prefix)
{
  if (verbose) {
    std::ostringstream os;
    os << *prefix << "Allocate Kokkos::View " << name
       << ": " << size << std::endl;
    std::cerr << os.str();
  }
  using Kokkos::view_alloc;
  using Kokkos::WithoutInitializing;
  return ViewType(view_alloc(name, WithoutInitializing), size);
}

template<class ViewType>
ViewType
make_initialized_view(
  const std::string& name,
  const size_t size,
  const bool verbose,
  const std::string* const prefix)
{
  if (verbose) {
    std::ostringstream os;
    os << *prefix << "Allocate & initialize Kokkos::View "
       << name << ": " << size << std::endl;
    std::cerr << os.str();
  }
  return ViewType(name, size);
}

template<class OutViewType, class InViewType>
void
assign_to_view(OutViewType& out,
               const InViewType& in,
               const char viewName[],
               const bool verbose,
               const std::string* const prefix)
{
  if (verbose) {
    std::ostringstream os;
    os << *prefix << "Assign to Kokkos::View " << viewName
       << ": Old size: " << out.extent(0)
       << ", New size: " << in.extent(0) << std::endl;
    std::cerr << os.str();
  }
  out = in;
}

template<class MemorySpace, class ViewType>
auto create_mirror_view(
  const MemorySpace& memSpace,
  const ViewType& view,
  const bool verbose,
  const std::string* const prefix) ->
  decltype(Kokkos::create_mirror_view(memSpace, view))
{
  if (verbose) {
    std::ostringstream os;
    os << *prefix << "Create mirror view: "
       << "view.extent(0): " << view.extent(0) << std::endl;
    std::cerr << os.str();
  }
  return Kokkos::create_mirror_view(memSpace, view);
}

enum class PadCrsAction {
  INDICES_ONLY,
  INDICES_AND_VALUES
};

/// \brief Implementation of padCrsArrays
///
/// \param row_ptr_beg [in] Offset to beginning of each row.
/// \param row_ptr_end [in] Offset to end of each row.
///
/// Each row lclRow has row_ptr_end[lclRow] - row_ptr_beg[lclRow]
/// entries.  Offsets row_ptr_end[lclRow] to
/// row_ptr_beg[lclRow+1] - 1 (inclusive) are extra space.
template<class RowPtr, class Indices, class Values, class Padding>
void
pad_crs_arrays(
  const PadCrsAction action,
  const RowPtr& row_ptr_beg,
  const RowPtr& row_ptr_end,
  Indices& indices_wdv,
  Values& values_wdv,
  const Padding& padding,
  const int my_rank,
  const bool verbose)
{
  using execution_space = typename Indices::t_dev::execution_space;
  using Kokkos::view_alloc;
  using Kokkos::WithoutInitializing;
  using std::endl;
  std::unique_ptr<std::string> prefix;

  const size_t maxNumToPrint = verbose ?
    Behavior::verbosePrintCountThreshold() : size_t(0);
  if (verbose) {
    std::ostringstream os;
    os << "Proc " << my_rank << ": Tpetra::...::pad_crs_arrays: ";
    prefix = std::unique_ptr<std::string>(new std::string(os.str()));
    os << "Start" << endl;
    std::cerr << os.str();
  }
  Kokkos::HostSpace hostSpace;

  if (verbose) {
    std::ostringstream os;
    os << *prefix << "On input: ";
    auto row_ptr_beg_h =
      Kokkos::create_mirror_view(hostSpace, row_ptr_beg);
    // DEEP_COPY REVIEW - NOT TESTED
    Kokkos::deep_copy(row_ptr_beg_h, row_ptr_beg);
    verbosePrintArray(os, row_ptr_beg_h, "row_ptr_beg before scan",
                      maxNumToPrint);
    os << ", ";
    auto row_ptr_end_h =
      Kokkos::create_mirror_view(hostSpace, row_ptr_end);
    // DEEP_COPY REVIEW - NOT TESTED
    Kokkos::deep_copy(row_ptr_end_h, row_ptr_end);
    verbosePrintArray(os, row_ptr_end_h, "row_ptr_end before scan",
                      maxNumToPrint);
    os << ", indices.extent(0): " << indices_wdv.extent(0)
       << ", values.extent(0): " << values_wdv.extent(0)
       << ", padding: ";
    padding.print(os);
    os << endl;
    std::cerr << os.str();
  }

  if (row_ptr_beg.size() == 0) {
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done; local matrix has no rows" << endl;
      std::cerr << os.str();
    }
    return; // nothing to do
  }

  const size_t lclNumRows(row_ptr_beg.size() - 1);
  RowPtr newAllocPerRow =
    make_uninitialized_view<RowPtr>("newAllocPerRow", lclNumRows,
                                    verbose, prefix.get());
  if (verbose) {
    std::ostringstream os;
    os << *prefix << "Fill newAllocPerRow & compute increase" << endl;
    std::cerr << os.str();
  }
  size_t increase = 0;
  {
    // Must do on host because padding uses std::map
    execution_space exec_space_instance = execution_space();
    auto row_ptr_end_h = create_mirror_view(
      hostSpace, row_ptr_end, verbose, prefix.get());
    // DEEP_COPY REVIEW - DEVICE-TO-HOSTMIRROR
    Kokkos::deep_copy(exec_space_instance, row_ptr_end_h, row_ptr_end);
    auto row_ptr_beg_h = create_mirror_view(
      hostSpace, row_ptr_beg, verbose, prefix.get());
    // DEEP_COPY REVIEW - DEVICE-TO-HOSTMIRROR
    Kokkos::deep_copy(exec_space_instance, row_ptr_beg_h, row_ptr_beg);

    // lbv 03/15/23: The execution space deep_copy does an asynchronous
    // copy so we really want to fence that space before touching the
    // host data as it is not guarenteed to have arrived by the time we
    // start the parallel_reduce below which might use a different
    // execution space, see:
    // https://kokkos.github.io/kokkos-core-wiki/API/core/view/deep_copy.html#semantics
    exec_space_instance.fence();

    auto newAllocPerRow_h = create_mirror_view(
      hostSpace, newAllocPerRow, verbose, prefix.get());
    using host_range_type = Kokkos::RangePolicy<
      Kokkos::DefaultHostExecutionSpace, size_t>;
    Kokkos::parallel_reduce
      ("Tpetra::CrsGraph: Compute new allocation size per row",
       host_range_type(0, lclNumRows),
       [&] (const size_t lclRowInd, size_t& lclIncrease) {
         const size_t start = row_ptr_beg_h[lclRowInd];
         const size_t end   = row_ptr_beg_h[lclRowInd+1];
         TEUCHOS_ASSERT( end >= start );
         const size_t oldAllocSize = end - start;
         const size_t oldNumEnt = row_ptr_end_h[lclRowInd] - start;
         TEUCHOS_ASSERT( oldNumEnt <= oldAllocSize );

         // This is not a pack routine.  Do not shrink!  Shrinking now
         // to fit the number of entries would ignore users' hint for
         // the max number of entries in each row.  Also, CrsPadding
         // only counts entries and ignores any available free space.

         auto result = padding.get_result(lclRowInd);
         const size_t newNumEnt = oldNumEnt + result.numInSrcNotInTgt;
         if (newNumEnt > oldAllocSize) {
           lclIncrease += (newNumEnt - oldAllocSize);
           newAllocPerRow_h[lclRowInd] = newNumEnt;
         }
         else {
           newAllocPerRow_h[lclRowInd] = oldAllocSize;
         }
       }, increase);

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "increase: " << increase << ", ";
      verbosePrintArray(os, newAllocPerRow_h, "newAllocPerRow",
                        maxNumToPrint);
      os << endl;
      std::cerr << os.str();
    }

    if (increase == 0) {
      return;
    }
    // DEEP_COPY REVIEW - HOSTMIRROR-TO-DEVICE
    Kokkos::deep_copy(execution_space(), newAllocPerRow, newAllocPerRow_h);
  }

  using inds_value_type = 
        typename Indices::t_dev::non_const_value_type;
  using vals_value_type = typename Values::t_dev::non_const_value_type;

  {
    auto indices_old = indices_wdv.getDeviceView(Access::ReadOnly);
    const size_t newIndsSize = size_t(indices_old.size()) + increase;
    auto indices_new = make_uninitialized_view<typename Indices::t_dev>(
      "Tpetra::CrsGraph column indices", newIndsSize, verbose,
      prefix.get());

    typename Values::t_dev values_new;
    auto values_old = values_wdv.getDeviceView(Access::ReadOnly);
    if (action == PadCrsAction::INDICES_AND_VALUES) {
      const size_t newValsSize = newIndsSize;
      // NOTE (mfh 10 Feb 2020) If we don't initialize values_new here,
      // then the CrsMatrix tests fail.
      values_new = make_initialized_view<typename Values::t_dev>(
        "Tpetra::CrsMatrix values", newValsSize, verbose, prefix.get());
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Repack" << endl;
      std::cerr << os.str();
    }

    using range_type = Kokkos::RangePolicy<execution_space, size_t>;
    Kokkos::parallel_scan(
      "Tpetra::CrsGraph or CrsMatrix repack",
      range_type(size_t(0), size_t(lclNumRows+1)),
      KOKKOS_LAMBDA (const size_t lclRow, size_t& newRowBeg,
                    const bool finalPass)
      {
        // row_ptr_beg    has lclNumRows + 1 entries.
        // row_ptr_end    has lclNumRows     entries.
        // newAllocPerRow has lclNumRows     entries.
        const size_t row_beg = row_ptr_beg[lclRow];
        const size_t row_end =
          lclRow < lclNumRows ? row_ptr_end[lclRow] : row_beg;
        const size_t numEnt = row_end - row_beg;
        const size_t newRowAllocSize =
          lclRow < lclNumRows ? newAllocPerRow[lclRow] : size_t(0);
        if (finalPass) {
          if (lclRow < lclNumRows) {
            const Kokkos::pair<size_t, size_t> oldRange(
              row_beg, row_beg + numEnt);
            const Kokkos::pair<size_t, size_t> newRange(
              newRowBeg, newRowBeg + numEnt);
            auto oldColInds = Kokkos::subview(indices_old, oldRange);
            auto newColInds = Kokkos::subview(indices_new, newRange);
            // memcpy works fine on device; the next step is to
            // introduce two-level parallelism and use team copy.
            memcpy(newColInds.data(), oldColInds.data(),
                  numEnt * sizeof(inds_value_type));
            if (action == PadCrsAction::INDICES_AND_VALUES) {
              auto oldVals = 
                  Kokkos::subview(values_old, oldRange);
              auto newVals = Kokkos::subview(values_new, newRange);
              memcpy((void*) newVals.data(), oldVals.data(),
                    numEnt * sizeof(vals_value_type));
            }
          }
          // It's the final pass, so we can modify these arrays.
          row_ptr_beg[lclRow] = newRowBeg;
          if (lclRow < lclNumRows) {
            row_ptr_end[lclRow] = newRowBeg + numEnt;
          }
        }
        newRowBeg += newRowAllocSize;
      });

    if (verbose) 
    {
      std::ostringstream os;

      os << *prefix;
      auto row_ptr_beg_h =
        Kokkos::create_mirror_view(hostSpace, row_ptr_beg);
      // DEEP_COPY REVIEW - NOT TESTED
      Kokkos::deep_copy(row_ptr_beg_h, row_ptr_beg);
      verbosePrintArray(os, row_ptr_beg_h, "row_ptr_beg after scan",
                        maxNumToPrint);
      os << endl;

      os << *prefix;
      auto row_ptr_end_h =
        Kokkos::create_mirror_view(hostSpace, row_ptr_end);
      // DEEP_COPY REVIEW - NOT TESTED
      Kokkos::deep_copy(row_ptr_end_h, row_ptr_end);
      verbosePrintArray(os, row_ptr_end_h, "row_ptr_end after scan",
                        maxNumToPrint);
      os << endl;

      std::cout << os.str();
    }

    indices_wdv = Indices(indices_new);
    values_wdv = Values(values_new);
  }
  
  if (verbose) {
    auto indices_h = indices_wdv.getHostView(Access::ReadOnly);
    auto values_h = values_wdv.getHostView(Access::ReadOnly);
    std::ostringstream os;
    os << "On output: ";
    verbosePrintArray(os, indices_h, "indices", maxNumToPrint);
    os << ", ";
    verbosePrintArray(os, values_h, "values", maxNumToPrint);
    os << ", padding: ";
    padding.print(os);
    os << endl;
  }

  if (verbose) {
    std::ostringstream os;
    os << *prefix << "Done" << endl;
    std::cerr << os.str();
  }
}

/// \brief Implementation of insertCrsIndices
template <class Pointers, class InOutIndices, class InIndices, class IndexMap>
size_t
insert_crs_indices(
    typename Pointers::value_type const row,
    Pointers const& row_ptrs,
    InOutIndices& cur_indices,
    size_t& num_assigned,
    InIndices const& new_indices,
    IndexMap&& map,
    std::function<void(size_t const, size_t const, size_t const)> cb)
{
  if (new_indices.size() == 0) {
    return 0;
  }

  if (cur_indices.size() == 0) {
    // No room to insert new indices
    return Teuchos::OrdinalTraits<size_t>::invalid();
  }

  using offset_type = typename std::decay<decltype (row_ptrs[0])>::type;
  using ordinal_type = typename std::decay<decltype (cur_indices[0])>::type;

  const offset_type start = row_ptrs[row];
  offset_type end = start + static_cast<offset_type> (num_assigned);
  const size_t num_avail = (row_ptrs[row + 1] < end) ? size_t (0) :
    row_ptrs[row + 1] - end;
  const size_t num_new_indices = static_cast<size_t> (new_indices.size ());
  size_t num_inserted = 0;

  size_t numIndicesLookup = num_assigned + num_new_indices;

  // Threshold determined from test/Utils/insertCrsIndicesThreshold.cpp
  const size_t useLookUpTableThreshold = 400; 

  if (numIndicesLookup <= useLookUpTableThreshold || num_new_indices == 1) {
    // For rows with few nonzeros, can use a serial search to find duplicates
    // Or if inserting only one index, serial search is as fast as anything else
    for (size_t k = 0; k < num_new_indices; ++k) {
      const ordinal_type idx = std::forward<IndexMap>(map)(new_indices[k]);
      offset_type row_offset = start;
      for (; row_offset < end; ++row_offset) {
        if (idx == cur_indices[row_offset]) {
          break;
        }
      }
  
      if (row_offset == end) {
        if (num_inserted >= num_avail) { // not enough room
          return Teuchos::OrdinalTraits<size_t>::invalid();
        }
        // This index is not yet in indices
        cur_indices[end++] = idx;
        num_inserted++;
      }
      if (cb) {
        cb(k, start, row_offset - start);
      }
    }
  }
  else {
    // For rows with many nonzeros, use a lookup table to find duplicates
    std::unordered_map<ordinal_type, offset_type> idxLookup(numIndicesLookup);

    // Put existing indices into the lookup table
    for (size_t k = 0; k < num_assigned; k++) {
      idxLookup[cur_indices[start+k]] = start+k;
    }

    // Check for new indices in table; insert if not there yet
    for (size_t k = 0; k < num_new_indices; k++) {
      const ordinal_type idx = std::forward<IndexMap>(map)(new_indices[k]);
      offset_type row_offset;

      auto it = idxLookup.find(idx);
      if (it == idxLookup.end()) {
        if (num_inserted >= num_avail) { // not enough room
          return Teuchos::OrdinalTraits<size_t>::invalid();
        }
        // index not found; insert it
        row_offset = end;
        cur_indices[end++] = idx;
        idxLookup[idx] = row_offset;
        num_inserted++;
      }
      else {
        // index found; note its position
        row_offset = it->second;
      }
      if (cb) {
        cb(k, start, row_offset - start);
      }
    }
  }
  num_assigned += num_inserted;
  return num_inserted;
}

/// \brief Implementation of findCrsIndices
template <class Pointers, class Indices1, class Indices2, class IndexMap, class Callback>
size_t
find_crs_indices(
    typename Pointers::value_type const row,
    Pointers const& row_ptrs,
    const size_t curNumEntries,
    Indices1 const& cur_indices,
    Indices2 const& new_indices,
    IndexMap&& map,
    Callback&& cb)
{
  if (new_indices.size() == 0)
    return 0;

  using ordinal = 
        typename std::remove_const<typename Indices1::value_type>::type;
  auto invalid_ordinal = Teuchos::OrdinalTraits<ordinal>::invalid();

  const size_t start = static_cast<size_t> (row_ptrs[row]);
  const size_t end = start + curNumEntries;
  size_t num_found = 0;
  for (size_t k = 0; k < new_indices.size(); k++)
  {
    auto row_offset = start;
    auto idx = std::forward<IndexMap>(map)(new_indices[k]);
    if (idx == invalid_ordinal)
      continue;
    for (; row_offset < end; row_offset++)
    {
      if (idx == cur_indices[row_offset])
      {
        std::forward<Callback>(cb)(k, start, row_offset - start);
        num_found++;
      }
    }
  }
  return num_found;
}

} // namespace impl


/// \brief Determine if the row pointers and indices arrays need to be resized
///   to accommodate new entries. If they do need to be resized, resize the
///   indices arrays and shift the existing contents to accommodate new entries.
///   Modify values in the row pointers array to point to the newly shifted
///   locations in the indices arrays.
///
///   This routine is called to resize/shift the CRS arrays before attempting to
///   insert new values if the number of new values exceeds the amount of free
///   space in the CRS arrays.
///
/// \param [in/out] rowPtrBeg - rowPtrBeg[i] points to the first
///        column index (in the indices array) of row i.
/// \param [in/out] rowPtrEnd - rowPtrEnd[i] points to the last
///        column index (in the indices array) of row i.
/// \param [in/out] indices - array containing columns indices of nonzeros in
///        CRS representation.
///
template<class RowPtr, class Indices, class Padding>
void
padCrsArrays(
    const RowPtr& rowPtrBeg,
    const RowPtr& rowPtrEnd,
    Indices& indices_wdv,
    const Padding& padding,
    const int my_rank,
    const bool verbose)
{
  using impl::pad_crs_arrays;
  // send empty values array
  Indices values_null; 
  pad_crs_arrays<RowPtr, Indices, Indices, Padding>( 
    impl::PadCrsAction::INDICES_ONLY, rowPtrBeg, rowPtrEnd,
    indices_wdv, values_null, padding, my_rank, verbose);
}

template<class RowPtr, class Indices, class Values, class Padding>
void
padCrsArrays(
    const RowPtr& rowPtrBeg,
    const RowPtr& rowPtrEnd,
    Indices& indices_wdv,
    Values& values_wdv,
    const Padding& padding,
    const int my_rank,
    const bool verbose)
{
  using impl::pad_crs_arrays;
  pad_crs_arrays<RowPtr, Indices, Values, Padding>(
    impl::PadCrsAction::INDICES_AND_VALUES, rowPtrBeg, rowPtrEnd,
    indices_wdv, values_wdv, padding, my_rank, verbose);
}

/// \brief Insert new indices in to current list of indices
///
/// \param row [in] The row in which to insert
/// \param rowPtrs [in] "Pointers" to beginning of each row
/// \param curIndices [in/out] The current indices
/// \param numAssigned [in/out] The number of currently assigned indices in row \c row
/// \param newIndices [in] The indices to insert
/// \param map [in] An optional function mapping newIndices[k] to its actual index
/// \param cb [in] An optional callback function called on every insertion at the local
///     index and the offset in to the inserted location
/// \return numInserted The number of indices inserted. If there is not
///     capacity in curIndices for newIndices, return -1;
///
/// \bf Notes
/// \c curIndices is the current list of CRS indices. it is not assumed to be sorted, but
/// entries are unique. For each \c newIndices[k], we look to see if the index exists in
/// \c cur_indices. If it does, we do not insert it (no repeats). If it does not exist, we
/// first check to make sure there is capacity in \c curIndices and if there is we insert
/// it at the end.
///
/// The actual value of \c newIndices[k] that is inserted is the value returned from \c
/// map(newIndices[k]). If an identity map is provided, \c newIndices[k] is directly
/// inserted. However, any other map can be provided. For instance, for a locally indexed
/// graph on which \c insertGlobalIndices is called, the \c curIndices array can be a
/// view of the graph's local indices, the \c newIndices array are the new *global*
/// indices, and \c map is the graph's column map to convert global indices to local.
/// If this function is called through the overload below without the \c map
/// argument, the identity map is provided.
///
/// The optional function \c cb is called on every valid index. \c cb is sent the
/// current loop index \c k, \c rowPtrs[k] (the start of the row), and the relative
/// offset from \c start in to the \c curIndices array for \c newIndices[k]. This
/// function could, for example, be used by \c CrsMatrix to fill the values array during
/// \c sumInto*Values or \c replace*Values; Eg, \c CrsMatrix::sumIntoLocalValues
/// might have the following:
///
/// <code>
/// CrsMatrix::sumIntoLocalValues(LO row, array<LO> cols, array<S> vals)
/// {
///   this->graph_->insertLocalValues(row, cols,
///       [&](size_t const k, size_t const start, size_t const offset){
///           this->values_[start+offset] += vals[k]; });
/// }
/// </code>
///
template <class Pointers, class InOutIndices, class InIndices>
size_t
insertCrsIndices(
    typename Pointers::value_type const row,
    Pointers const& rowPtrs,
    InOutIndices& curIndices,
    size_t& numAssigned,
    InIndices const& newIndices,
    std::function<void(const size_t, const size_t, const size_t)> cb =
        std::function<void(const size_t, const size_t, const size_t)>())
{
  static_assert(std::is_same<typename std::remove_const<typename InOutIndices::value_type>::type,
                             typename std::remove_const<typename InIndices::value_type>::type>::value,
    "Expected views to have same value type");

  // Provide a unit map for the more general insert_indices
  using ordinal = typename InOutIndices::value_type;
  auto numInserted = impl::insert_crs_indices(row, rowPtrs, curIndices,
    numAssigned, newIndices, [](ordinal const idx) { return idx; }, cb);
  return numInserted;
}

template <class Pointers, class InOutIndices, class InIndices>
size_t
insertCrsIndices(
    typename Pointers::value_type const row,
    Pointers const& rowPtrs,
    InOutIndices& curIndices,
    size_t& numAssigned,
    InIndices const& newIndices,
    std::function<typename InOutIndices::value_type(const typename InIndices::value_type)> map,
    std::function<void(const size_t, const size_t, const size_t)> cb =
      std::function<void(const size_t, const size_t, const size_t)>())
{
  auto numInserted = impl::insert_crs_indices(row, rowPtrs, curIndices,
    numAssigned, newIndices, map, cb);
  return numInserted;
}


/// \brief Finds offsets in to current list of indices
///
/// \param row [in] The row in which to insert
/// \param rowPtrs [in] "Pointers" to beginning of each row
/// \param curIndices [in] The current indices
/// \param numAssigned [in] The number of currently assigned indices in row \c row
/// \param newIndices [in] The indices to insert
/// \param cb [in] An optional function called on every insertion at the local
///     index and the offset in to the inserted location
/// \return numFound The number of indices found.
///
/// \bf Notes
/// \c curIndices is the current list of CRS indices. it is not assumed to be sorted, but
/// entries are unique. For each \c newIndices[k], we look to see if the index exists in
/// \c curIndices. If it does, we do not insert it (no repeats). If it does not exist, we
/// first check to make sure there is capacity in \c curIndices and if there is we insert
/// it at the end.
///
/// The actual value of \c newIndices[k] that is inserted is the value returned from \c
/// map(newIndices[k]). If an identity map is provided, \c newIndices[k] is directly
/// inserted. However, any other map can be provided. For instance, for a locally indexed
/// graph on which \c insertGlobalIndices is called, the \c curIndices array can be a
/// view of the graph's local indices, the \c newIndices array are the new *global*
/// indices, and \c map is the graph's column map to convert global indices to local.
/// If this function is called through the overload below without the \c map
/// argument, the identity map is provided.
///
/// The function \c cb is called on every valid index.
///
template <class Pointers, class Indices1, class Indices2, class Callback>
size_t
findCrsIndices(
    typename Pointers::value_type const row,
    Pointers const& rowPtrs,
    const size_t curNumEntries,
    Indices1 const& curIndices,
    Indices2 const& newIndices,
    Callback&& cb)
{
  static_assert(std::is_same<typename std::remove_const<typename Indices1::value_type>::type,
                             typename std::remove_const<typename Indices2::value_type>::type>::value,
    "Expected views to have same value type");
  // Provide a unit map for the more general find_crs_indices
  using ordinal = typename Indices2::value_type;
  auto numFound = impl::find_crs_indices(row, rowPtrs, curNumEntries, curIndices, newIndices,
    [=](ordinal ind){ return ind; }, cb);
  return numFound;
}

template <class Pointers, class Indices1, class Indices2, class IndexMap, class Callback>
size_t
findCrsIndices(
    typename Pointers::value_type const row,
    Pointers const& rowPtrs,
    const size_t curNumEntries,
    Indices1 const& curIndices,
    Indices2 const& newIndices,
    IndexMap&& map,
    Callback&& cb)
{
  return impl::find_crs_indices(row, rowPtrs, curNumEntries, curIndices, newIndices, map, cb);
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_CRSUTILS_HPP
