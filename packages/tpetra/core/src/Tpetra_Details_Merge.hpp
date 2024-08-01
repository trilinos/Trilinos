// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_MERGE_HPP
#define TPETRA_DETAILS_MERGE_HPP

#include "TpetraCore_config.h"
#include "Teuchos_TestForException.hpp"
#include <algorithm> // std::sort
#include <utility> // std::pair, std::make_pair
#include <stdexcept>

namespace Tpetra {
namespace Details {

/// \brief Count the number of column indices that can be merged into
///   the current row, assuming that both the current row's indices
///   and the input indices are unsorted.
///
/// Neither the current row's entries, nor the input, are sorted.
/// Return the number of input entries that can be merged into the
/// current row.  Don't actually merge them.  'numCurInds' corresponds
/// to 'midPos' in mergeUnsortedIndices.
///
/// The current indices are NOT allowed to have repeats, but the input
/// indices ARE allowed to have repeats.  (The whole point of these
/// methods is to keep the current entries without repeats -- "merged
/// in.")  Repeats in the input are counted separately with respect to
/// merges.
///
/// The unsorted case is bad for asymptotics, but the asymptotics only
/// show up with dense or nearly dense rows, which are bad for other
/// reasons.
template<class OrdinalType, class IndexType>
IndexType
countMergeUnsortedIndices (const OrdinalType curInds[],
                           const IndexType numCurInds,
                           const OrdinalType inputInds[],
                           const IndexType numInputInds)
{
  IndexType mergeCount = 0;

  if (numCurInds <= numInputInds) {
    // More input than current entries, so iterate linearly over
    // input and scan current entries repeatedly.
    for (IndexType inPos = 0; inPos < numInputInds; ++inPos) {
      const OrdinalType inVal = inputInds[inPos];
      for (IndexType curPos = 0; curPos < numCurInds; ++curPos) {
        if (curInds[curPos] == inVal) {
          ++mergeCount;
        }
      }
    }
  }
  else { // numCurInds > numInputInds
    // More current entries than input, so iterate linearly over
    // current entries and scan input repeatedly.
    for (IndexType curPos = 0; curPos < numCurInds; ++curPos) {
      const OrdinalType curVal = curInds[curPos];
      for (IndexType inPos = 0; inPos < numInputInds; ++inPos) {
        if (inputInds[inPos] == curVal) {
          ++mergeCount;
        }
      }
    }
  }

#ifdef HAVE_TPETRA_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION
    (mergeCount > numInputInds, std::logic_error, "mergeCount = " <<
     mergeCount << " > numInputInds = " << numInputInds << ".");
#endif // HAVE_TPETRA_DEBUG
  return mergeCount;
}

/// \brief Count the number of column indices that can be merged into
///   the current row, assuming that both the current row's indices
///   and the input indices are sorted.
///
/// Both the current row's entries and the input are sorted.
/// Return the number of input entries that can be merged into the
/// current row.  Don't actually merge them.  'numCurInds'
/// corresponds to 'midPos' in mergeSortedIndices.
///
/// The current indices are NOT allowed to have repeats, but the input
/// indices ARE allowed to have repeats.  (The whole point of these
/// methods is to keep the current entries without repeats -- "merged
/// in.")  Repeats in the input are counted separately with respect to
/// merges.
///
/// The sorted case is good for asymptotics, but imposes an order
/// on the entries of each row.  Sometimes users don't want that.
template<class OrdinalType, class IndexType>
IndexType
countMergeSortedIndices (const OrdinalType curInds[],
                         const IndexType numCurInds,
                         const OrdinalType inputInds[],
                         const IndexType numInputInds)
{
  // Only count possible merges; don't merge yet.  If the row
  // doesn't have enough space, we want to return without side
  // effects.
  IndexType curPos = 0;
  IndexType inPos = 0;
  IndexType mergeCount = 0;
  while (inPos < numInputInds && curPos < numCurInds) {
    const OrdinalType inVal = inputInds[inPos];
    const OrdinalType curVal = curInds[curPos];

    if (curVal == inVal) { // can merge
      ++mergeCount;
      ++inPos; // go on to next input
    } else if (curVal < inVal) {
      ++curPos; // go on to next row entry
    } else { // curVal > inVal
      ++inPos; // can't merge it ever, since row entries sorted
    }
  }

#ifdef HAVE_TPETRA_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION
    (inPos > numInputInds, std::logic_error, "inPos = " << inPos <<
     " > numInputInds = " << numInputInds << ".");
  TEUCHOS_TEST_FOR_EXCEPTION
    (curPos > numCurInds, std::logic_error, "curPos = " << curPos <<
     " > numCurInds = " << numCurInds << ".");
  TEUCHOS_TEST_FOR_EXCEPTION
    (mergeCount > numInputInds, std::logic_error, "mergeCount = " <<
     mergeCount << " > numInputInds = " << numInputInds << ".");
#endif // HAVE_TPETRA_DEBUG

  // At this point, 2 situations are possible:
  //
  // 1. inPos == numInputInds: We looked at all inputs.  Some
  //    (mergeCount of them) could have been merged.
  // 2. inPos < numInputInds: We didn't get to look at all inputs.
  //    Since the inputs are sorted, we know that those inputs we
  //    didn't examine weren't mergeable.
  //
  // Either way, mergeCount gives the number of mergeable inputs.
  return mergeCount;
}


/// \brief Attempt to merge the input indices into the current row's
///   column indices, assuming that both the current row's indices and
///   the input indices are sorted.
///
/// Both the current row's entries and the input are sorted.  If and
/// only if the current row has enough space for the input (after
/// merging), merge the input with the current row.
///
/// Assume that both curInds and inputInds are sorted.
/// Current indices: curInds[0 ..  midPos-1].
/// Extra space at end: curInds[midPos .. endPos-1]
/// Input indices to merge in: inputInds[0 .. numInputInds].
/// Any of those could be empty.
///
/// If the merge succeeded, return true and the new number of entries
/// in the row.  Else, return false and the new number of entries in
/// the row required to fit the input.
///
/// The sorted case is good for asymptotics, but imposes an order on
/// the entries of each row.  Sometimes users don't want that.
template<class OrdinalType, class IndexType>
std::pair<bool, IndexType>
mergeSortedIndices (OrdinalType curInds[],
                    const IndexType midPos,
                    const IndexType endPos,
                    const OrdinalType inputInds[],
                    const IndexType numInputInds)
{
  // Optimize for the following cases, in decreasing order of
  // optimization concern:
  //
  //   a. Current row has allocated space but no entries
  //   b. All input indices already in the graph
  //
  // If the row has insufficient space for a merge, don't do
  // anything!  Just return an upper bound on the number of extra
  // entries required to fit everything.  This imposes extra cost,
  // but correctly supports the count, allocate, fill, compute
  // pattern.  (If some entries were merged in and others weren't,
  // how would you know which got merged in?  CrsGraph insert is
  // idempotent, but CrsMatrix insert does a += on the value and
  // is therefore not idempotent.)
  if (midPos == 0) {
    // Current row has no entries, but may have preallocated space.
    if (endPos >= numInputInds) {
      // Sufficient space for new entries; copy directly.
      for (IndexType k = 0; k < numInputInds; ++k) {
        curInds[k] = inputInds[k];
      }
      std::sort (curInds, curInds + numInputInds);
      return std::make_pair (true, numInputInds);
    }
    else { // not enough space
      return std::make_pair (false, numInputInds);
    }
  }
  else { // current row contains indices, requiring merge
    // Only count possible merges; don't merge yet.  If the row
    // doesn't have enough space, we want to return without side
    // effects.
    const IndexType mergeCount =
      countMergeSortedIndices<OrdinalType, IndexType> (curInds, midPos,
                                                       inputInds,
                                                       numInputInds);
    const IndexType extraSpaceNeeded = numInputInds - mergeCount;
    const IndexType newRowLen = midPos + extraSpaceNeeded;
    if (newRowLen > endPos) {
      return std::make_pair (false, newRowLen);
    }
    else { // we have enough space; merge in
      IndexType curPos = 0;
      IndexType inPos = 0;
      IndexType newPos = midPos;
      while (inPos < numInputInds && curPos < midPos) {
        const OrdinalType inVal = inputInds[inPos];
        const OrdinalType curVal = curInds[curPos];

        if (curVal == inVal) { // can merge
          ++inPos; // merge and go on to next input
        } else if (curVal < inVal) {
          ++curPos; // go on to next row entry
        } else { // curVal > inVal
          // The input doesn't exist in the row.
          // Copy it to the end; we'll sort it in later.
          curInds[newPos] = inVal;
          ++newPos;
          ++inPos; // move on to next input
        }
      }

      // If any inputs remain, and the current row has space for them,
      // then copy them in.  We'll sort them later.
      for (; inPos < numInputInds && newPos < newRowLen; ++inPos, ++newPos) {
        curInds[newPos] = inputInds[inPos];
      }

#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION
        (newPos != newRowLen, std::logic_error, "mergeSortedIndices: newPos = "
         << newPos << " != newRowLen = " << newRowLen << " = " << midPos <<
         " + " << extraSpaceNeeded << ".  Please report this bug to the Tpetra "
         "developers.");
#endif // HAVE_TPETRA_DEBUG

      if (newPos != midPos) { // new entries at end; sort them in
        // FIXME (mfh 03 Jan 2016) Rather than sorting, it would
        // be faster (linear time) just to iterate backwards
        // through the current entries, pushing them over to make
        // room for unmerged input.  However, I'm not so worried
        // about the asymptotics here, because dense rows in a
        // sparse matrix are ungood anyway.
        std::sort (curInds, curInds + newPos);
      }
      return std::make_pair (true, newPos);
    }
  }
}


/// \brief Attempt to merge the input indices into the current row's
///   column indices, assuming that both the current row's indices and
///   the input indices are unsorted.
///
/// Neither the current row's entries nor the input are sorted.  If
/// and only if the current row has enough space for the input (after
/// merging), merge the input with the current row.
///
/// Assume that neither curInds nor inputInds are sorted.
/// Current indices: curInds[0 ..  midPos-1].
/// Extra space at end: curInds[midPos .. endPos-1]
/// Input indices to merge in: inputInds[0 .. numInputInds].
/// Any of those could be empty.
///
/// If the merge succeeded, return true and the new number of entries
/// in the row.  Else, return false and the new number of entries in
/// the row required to fit the input.
///
/// The unsorted case is bad for asymptotics, but the asymptotics only
/// show up with dense or nearly dense rows, which are bad for other
/// reasons.
template<class OrdinalType, class IndexType>
std::pair<bool, IndexType>
mergeUnsortedIndices (OrdinalType curInds[],
                      const IndexType midPos,
                      const IndexType endPos,
                      const OrdinalType inputInds[],
                      const IndexType numInputInds)
{
  // Optimize for the following cases, in decreasing order of
  // optimization concern:
  //
  //   a. Current row has allocated space but no entries
  //   b. All input indices already in the graph
  //
  // If the row has insufficient space for a merge, don't do
  // anything!  Just return an upper bound on the number of extra
  // entries required to fit everything.  This imposes extra cost,
  // but correctly supports the count, allocate, fill, compute
  // pattern.  (If some entries were merged in and others weren't,
  // how would you know which got merged in?  CrsGraph insert is
  // idempotent, but CrsMatrix insert does a += on the value and
  // is therefore not idempotent.)
  if (midPos == 0) {
    // Current row has no entries, but may have preallocated space.
    if (endPos >= numInputInds) {
      // Sufficient space for new entries; copy directly.
      for (IndexType k = 0; k < numInputInds; ++k) {
        curInds[k] = inputInds[k];
      }
      return std::make_pair (true, numInputInds);
    }
    else { // not enough space
      return std::make_pair (false, numInputInds);
    }
  }
  else { // current row contains indices, requiring merge
    // Only count possible merges; don't merge yet.  If the row
    // doesn't have enough space, we want to return without side
    // effects.
    const IndexType mergeCount =
      countMergeUnsortedIndices<OrdinalType, IndexType> (curInds, midPos,
                                                         inputInds,
                                                         numInputInds);
    const IndexType extraSpaceNeeded = numInputInds - mergeCount;
    const IndexType newRowLen = midPos + extraSpaceNeeded;
    if (newRowLen > endPos) {
      return std::make_pair (false, newRowLen);
    }
    else { // we have enough space; merge in
      // Iterate linearly over input.  Scan current entries
      // repeatedly.  Add new entries at end.
      IndexType newPos = midPos;
      for (IndexType inPos = 0; inPos < numInputInds; ++inPos) {
        const OrdinalType inVal = inputInds[inPos];
        bool merged = false;
        for (IndexType curPos = 0; curPos < midPos; ++curPos) {
          if (curInds[curPos] == inVal) {
            merged = true;
          }
        }
        if (! merged) {
          curInds[newPos] = inVal;
          ++newPos;
        }
      }
      return std::make_pair (true, newPos);
    }
  }
}

/// \brief Attempt to merge the input indices and values into the
///   current row's column indices and corresponding values, assuming
///   that both the current row's indices and the input indices are
///   unsorted.
///
/// Neither the current row's entries nor the input are sorted.  If
/// and only if the current row has enough space for the input (after
/// merging), merge the input with the current row.
///
/// Assume that neither curInds nor inputInds are sorted.
/// Current indices: curInds[0 .. midPos-1].
/// Current values: curVals[0 .. midPos-1].
/// Extra space for indices at end: curInds[midPos .. endPos-1].
/// Extra space for values at end: curVals[midPos .. endPos-1].
/// Input indices to merge in: inputInds[0 .. numInputInds].
/// Input values to merge in: inputVals[0 .. numInputInds].
///
/// If the merge succeeded, return true and the new number of entries
/// in the row.  Else, return false and the new number of entries in
/// the row required to fit the input.
///
/// The unsorted case is bad for asymptotics, but the asymptotics only
/// show up with dense or nearly dense rows, which are bad for other
/// reasons.
template<class OrdinalType, class ValueType, class IndexType>
std::pair<bool, IndexType>
mergeUnsortedIndicesAndValues (OrdinalType curInds[],
                               ValueType curVals[],
                               const IndexType midPos,
                               const IndexType endPos,
                               const OrdinalType inputInds[],
                               const ValueType inputVals[],
                               const IndexType numInputInds)
{
  // Optimize for the following cases, in decreasing order of
  // optimization concern:
  //
  //   a. Current row has allocated space but no entries
  //   b. All input indices already in the graph
  //
  // If the row has insufficient space for a merge, don't do
  // anything!  Just return an upper bound on the number of extra
  // entries required to fit everything.  This imposes extra cost,
  // but correctly supports the count, allocate, fill, compute
  // pattern.  (If some entries were merged in and others weren't,
  // how would you know which got merged in?  CrsGraph insert is
  // idempotent, but CrsMatrix insert does a += on the value and
  // is therefore not idempotent.)
  if (midPos == 0) {
    // Current row has no entries, but may have preallocated space.
    if (endPos >= numInputInds) {
      // Sufficient space for new entries; copy directly.
      for (IndexType k = 0; k < numInputInds; ++k) {
        curInds[k] = inputInds[k];
        curVals[k] = inputVals[k];
      }
      return std::make_pair (true, numInputInds);
    }
    else { // not enough space
      return std::make_pair (false, numInputInds);
    }
  }
  else { // current row contains indices, requiring merge
    // Only count possible merges; don't merge yet.  If the row
    // doesn't have enough space, we want to return without side
    // effects.
    const IndexType mergeCount =
      countMergeUnsortedIndices<OrdinalType, IndexType> (curInds, midPos,
                                                         inputInds,
                                                         numInputInds);
    const IndexType extraSpaceNeeded = numInputInds - mergeCount;
    const IndexType newRowLen = midPos + extraSpaceNeeded;
    if (newRowLen > endPos) {
      return std::make_pair (false, newRowLen);
    }
    else { // we have enough space; merge in
      // Iterate linearly over input.  Scan current entries
      // repeatedly.  Add new entries at end.
      IndexType newPos = midPos;
      for (IndexType inPos = 0; inPos < numInputInds; ++inPos) {
        const OrdinalType inInd = inputInds[inPos];
        bool merged = false;
        for (IndexType curPos = 0; curPos < midPos; ++curPos) {
          if (curInds[curPos] == inInd) {
            merged = true;
            curVals[curPos] += inputVals[inPos];
          }
        }
        if (! merged) {
          curInds[newPos] = inInd;
          curVals[newPos] = inputVals[inPos];
          ++newPos;
        }
      }
      return std::make_pair (true, newPos);
    }
  }
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_MERGE_HPP
