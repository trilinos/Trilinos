// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_COMPUTEOFFSETS_HPP
#define TPETRA_DETAILS_COMPUTEOFFSETS_HPP

/// \file Tpetra_Details_computeOffsets.hpp
/// \brief Declare and define the functions
///   Tpetra::Details::computeOffsetsFromCounts and
///   Tpetra::computeOffsetsFromConstantCount.  These functions are
///   implementation details of Tpetra (in particular, of
///   FixedHashTable, CrsGraph, and CrsMatrix).

#include "TpetraCore_config.h"
#include "Tpetra_Details_getEntryOnHost.hpp"
#include <limits>
#include <type_traits>

namespace Tpetra {
namespace Details {

//
// Implementation details for computeOffsetsFromCounts (see below).
// Users should skip over this anonymous namespace.
//
namespace { // (anonymous)

/// \brief Parallel scan functor for computing offsets from counts.
///
/// \warning This is NOT for users.  It is an implementation detail of
///   the computeOffsetsFromCounts function (see below), which you
///   should call instead.
///
/// OffsetType must be able to store the sum of all the entries of the
/// \c counts input array.  This functor makes no attempt to check for
/// overflow in this sum.
template<class OffsetType,
         class CountType,
         class SizeType>
class ComputeOffsetsFromCounts {
public:
  static_assert (std::is_integral<OffsetType>::value,
                 "OffsetType must be a built-in integer.");
  static_assert (std::is_integral<CountType>::value,
                 "CountType must be a built-in integer.");
  static_assert (std::is_integral<SizeType>::value,
                 "SizeType must be a built-in integer.");

  using offsets_view_type =
    Kokkos::View<OffsetType*, Kokkos::AnonymousSpace>;
  using counts_view_type =
    Kokkos::View<const CountType*, Kokkos::AnonymousSpace>;

  /// \brief Constructor
  ///
  /// \param offsets [out] (Preallocated) offsets; one entry longer
  ///   than \c counts
  /// \param counts [in] View of bucket counts
  ComputeOffsetsFromCounts (const offsets_view_type& offsets,
                            const counts_view_type& counts) :
    offsets_ (offsets),
    counts_ (counts),
    size_ (counts.extent (0))
  {}

  //! Reduction operator.
  KOKKOS_INLINE_FUNCTION void
  operator () (const SizeType i, OffsetType& update,
               const bool finalPass) const
  {
    const auto curVal = (i < size_) ? counts_[i] : OffsetType ();
    if (finalPass) {
      offsets_[i] = update;
    }
    update += (i < size_) ? curVal : OffsetType ();
  }

  template<class ExecutionSpace>
  static OffsetType
  run (const ExecutionSpace& execSpace,
       const offsets_view_type& offsets,
       const counts_view_type& counts)
  {
    const SizeType numCounts (counts.extent (0));
    using range_type = Kokkos::RangePolicy<ExecutionSpace, SizeType>;
    range_type range (execSpace, 0, numCounts + SizeType (1));
    using functor_type =
      ComputeOffsetsFromCounts<OffsetType, CountType, SizeType>;
    functor_type functor (offsets, counts);
    OffsetType total (0);
    const char funcName[] = "Tpetra::Details::computeOffsetsFromCounts";
    Kokkos::parallel_scan (funcName, range, functor, total);
    return total;
  }

private:
  //! Offsets (output argument)
  offsets_view_type offsets_;
  //! Bucket counts (input argument).
  counts_view_type counts_;
  //! Number of entries in counts_.
  SizeType size_;
};

/// \brief Parallel scan functor for computing offsets from a constant
///   count.
///
/// \warning This is NOT for users.  It is an implementation detail of
///   the computeOffsetsFromConstantCount function (see below), which
///   you should call instead.
///
/// OffsetType must be able to store <tt>ptr.extent(0) * count</tt>.
/// This functor makes no attempt to check for overflow in this sum.
template<class OffsetType,
         class CountType,
         class SizeType>
class ComputeOffsetsFromConstantCount {
public:
  static_assert (std::is_integral<OffsetType>::value,
                 "OffsetType must be a built-in integer.");
  static_assert (std::is_integral<CountType>::value,
                 "CountType must be a built-in integer.");
  static_assert (std::is_integral<SizeType>::value,
                 "SizeType must be a built-in integer.");

  using offsets_view_type =
    Kokkos::View<OffsetType*, Kokkos::AnonymousSpace>;

  /// \brief Constructor
  ///
  /// \param offsets [out] (Preallocated) offsets; one entry longer
  ///   than \c counts
  /// \param count [in] The constant count
  ComputeOffsetsFromConstantCount (const offsets_view_type& offsets,
                                   const CountType count) :
    offsets_ (offsets),
    count_ (count)
  {}

  //! Reduction operator.
  KOKKOS_INLINE_FUNCTION void
  operator () (const SizeType i) const
  {
    offsets_[i] = count_*i;
  }

  template<class ExecutionSpace>
  static OffsetType
  run (const ExecutionSpace& execSpace,
       const offsets_view_type& offsets,
       const CountType count)
  {
    const SizeType numOffsets (offsets.extent (0));
    if(numOffsets == SizeType(0))
    {
      // Special case that is possible with zero rows
      return 0;
    }
    using range_type = Kokkos::RangePolicy<ExecutionSpace, SizeType>;
    range_type range (execSpace, 0, numOffsets);
    using functor_type =
      ComputeOffsetsFromConstantCount<OffsetType, CountType, SizeType>;
    functor_type functor (offsets, count);
    const OffsetType total = (numOffsets - 1) * count;
    const char funcName[] =
      "Tpetra::Details::computeOffsetsFromConstantCount";
    Kokkos::parallel_for (funcName, range, functor);
    return total;
  }

private:
  //! Offsets (output argument)
  offsets_view_type offsets_;
  //! "Count" input argument
  CountType count_;
};

} // namespace (anonymous)

/// \brief Compute offsets from counts
///
/// Compute offsets from counts via prefix sum:
///
/// ptr[i+1] = \sum_{j=0}^{i} counts[j]
///
/// Thus, ptr[i+1] - ptr[i] = counts[i], so that ptr[i+1] = ptr[i] +
/// counts[i].  If we stored counts[i] in ptr[i+1] on input, then the
/// formula is ptr[i+1] += ptr[i].
///
/// \return Sum of all counts; last entry of \c ptr.
///
/// \tparam ExecutionSpace Kokkos execution space instance on which to
///   run.
/// \tparam OffsetsViewType Type of the Kokkos::View specialization
///   used to store the offsets; the output array of this function.
/// \tparam CountsViewType Type of the Kokkos::View specialization
///   used to store the counts; the input array of this function.
/// \tparam SizeType The parallel loop index type; a built-in integer
///   type.  Defaults to the type of the input View's dimension.  You
///   may use a shorter type to improve performance.
///
/// The type of each entry of the \c ptr array must be able to store
/// the sum of all the entries of \c counts.  This functor makes no
/// attempt to check for overflow in this sum.
template<class ExecutionSpace,
         class OffsetsViewType,
         class CountsViewType,
         class SizeType = typename OffsetsViewType::size_type>
typename OffsetsViewType::non_const_value_type
computeOffsetsFromCounts (const ExecutionSpace& execSpace,
                          const OffsetsViewType& ptr,
                          const CountsViewType& counts)
{
  static_assert (Kokkos::is_execution_space<ExecutionSpace>::value,
                 "ExecutionSpace must be a Kokkos execution space.");
  static_assert (Kokkos::is_view<OffsetsViewType>::value,
                 "OffsetsViewType (the type of ptr) must be a Kokkos::View.");
  static_assert (Kokkos::is_view<CountsViewType>::value,
                 "CountsViewType (the type of counts) must be a Kokkos::View.");
  static_assert (std::is_same<typename OffsetsViewType::value_type,
                   typename OffsetsViewType::non_const_value_type>::value,
                 "OffsetsViewType (the type of ptr) must be a nonconst Kokkos::View.");
  static_assert (static_cast<int> (OffsetsViewType::rank) == 1,
                 "OffsetsViewType (the type of ptr) must be a rank-1 Kokkos::View.");
  static_assert (static_cast<int> (CountsViewType::rank) == 1,
                 "CountsViewType (the type of counts) must be a rank-1 Kokkos::View.");

  using offset_type = typename OffsetsViewType::non_const_value_type;
  static_assert (std::is_integral<offset_type>::value,
                 "The entries of ptr must be built-in integers.");
  using count_type = typename CountsViewType::non_const_value_type;
  static_assert (std::is_integral<count_type>::value,
                 "The entries of counts must be built-in integers.");
  static_assert (std::is_integral<SizeType>::value,
                 "SizeType must be a built-in integer type.");

  const char funcName[] = "Tpetra::Details::computeOffsetsFromCounts";

  const auto numOffsets = ptr.size ();
  const auto numCounts = counts.size ();
  offset_type total (0);

  if (numOffsets != 0) {
    TEUCHOS_TEST_FOR_EXCEPTION
      (numCounts >= numOffsets, std::invalid_argument, funcName <<
       ": counts.size() = " << numCounts << " >= ptr.size() = " <<
       numOffsets << ".");

    using Kokkos::AnonymousSpace;
    using Kokkos::View;
    View<offset_type*, AnonymousSpace> ptr_a = ptr;
    View<const count_type*, AnonymousSpace> counts_a;

    using offsets_device_type = typename OffsetsViewType::device_type;
    using counts_copy_type = View<count_type*, offsets_device_type>;
    counts_copy_type counts_copy;

    using offsets_memory_space =
      typename offsets_device_type::memory_space;
    using counts_memory_space = typename CountsViewType::memory_space;
    constexpr bool countsAccessibleFromOffsetsExecSpace =
      Kokkos::SpaceAccessibility<
        offsets_memory_space, counts_memory_space>::accessible;
    if (countsAccessibleFromOffsetsExecSpace) {
      // NOTE (mfh 21 Aug 2019) Some compilers have trouble deducing
      // that operator= works if more than one template argument
      // differ.  If that should happen, introduce an intermediate
      // type here.
      counts_a = counts;
    }
    else {
      using Kokkos::view_alloc;
      using Kokkos::WithoutInitializing;
      counts_copy = counts_copy_type
        (view_alloc ("counts_copy", WithoutInitializing), numCounts);
      Kokkos::deep_copy (execSpace, counts_copy, counts);
      counts_a = counts_copy;
    }

    using functor_type =
      ComputeOffsetsFromCounts<offset_type, count_type, SizeType>;
    total = functor_type::run (execSpace, ptr_a, counts_a);
  }

  return total;
}

//! Overload that uses OffsetsViewType's execution space.
template<class OffsetsViewType,
         class CountsViewType,
         class SizeType = typename OffsetsViewType::size_type>
typename OffsetsViewType::non_const_value_type
computeOffsetsFromCounts (const OffsetsViewType& ptr,
                          const CountsViewType& counts)
{
  using execution_space = typename OffsetsViewType::execution_space;
  return computeOffsetsFromCounts (execution_space (), ptr, counts);
}

/// \brief Compute offsets from a constant count
///
/// Compute offsets from a constant count via prefix sum:
///
/// ptr[i+1] = \sum_{j=0}^{i} count
///
/// Thus, ptr[i+1] - ptr[i] = count, so that ptr[i+1] = ptr[i] +
/// count.
///
/// \return Sum of all counts; last entry of \c ptr.
///
/// \tparam OffsetsViewType Type of the Kokkos::View specialization
///   used to store the offsets; the output array of this function.
/// \tparam CountType Type of the constant count; the input argument
///   of this function.
/// \tparam SizeType The parallel loop index type; a built-in integer
///   type.  Defaults to the type of the output View's dimension.  You
///   may use a shorter type to improve performance.
///
/// The type of each entry of the \c ptr array must be able to store
/// <tt>ptr.extent (0) * count</tt>.  This functor makes no
/// attempt to check for overflow in this sum.
template<class OffsetsViewType,
         class CountType,
         class SizeType = typename OffsetsViewType::size_type>
typename OffsetsViewType::non_const_value_type
computeOffsetsFromConstantCount (const OffsetsViewType& ptr,
                                 const CountType count)
{
  static_assert (Kokkos::is_view<OffsetsViewType>::value,
                 "ptr must be a Kokkos::View.");
  static_assert (std::is_same<typename OffsetsViewType::value_type,
                   typename OffsetsViewType::non_const_value_type>::value,
                 "ptr must be a nonconst Kokkos::View.");
  static_assert (static_cast<int> (OffsetsViewType::rank) == 1,
                 "ptr must be a rank-1 Kokkos::View.");

  using offset_type = typename OffsetsViewType::non_const_value_type;
  static_assert (std::is_integral<offset_type>::value,
                 "The type of each entry of ptr must be a "
                 "built-in integer.");
  static_assert (std::is_integral<CountType>::value,
                 "CountType must be a built-in integer.");
  static_assert (std::is_integral<SizeType>::value,
                 "SizeType must be a built-in integer.");

  using device_type = typename OffsetsViewType::device_type;
  using execution_space = typename device_type::execution_space;

  offset_type total (0);
  if (ptr.extent (0) != 0) {
    using CT = CountType;
    using functor_type =
      ComputeOffsetsFromConstantCount<offset_type, CT, SizeType>;
    execution_space execSpace;
    total = functor_type::run (execSpace, ptr, count);
  }
  return total;
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_COMPUTEOFFSETS_HPP
