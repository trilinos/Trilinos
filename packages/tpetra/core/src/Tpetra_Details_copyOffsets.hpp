// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_COPYOFFSETS_HPP
#define TPETRA_DETAILS_COPYOFFSETS_HPP

/// \file Tpetra_Details_copyOffsets.hpp
/// \brief Declare and define Tpetra::Details::copyOffsets, an
///   implementation detail of Tpetra (in particular, of
///   FixedHashTable, CrsGraph, and CrsMatrix).

#include "TpetraCore_config.h"
#include "Tpetra_Details_Behavior.hpp"
#include "Kokkos_Core.hpp"
#include <limits>
#include <type_traits>

namespace Tpetra {
namespace Details {

//
// Implementation details for copyOffsets (see below).
// Users should skip over this anonymous namespace.
//
namespace { // (anonymous)

  // Implementation detail of copyOffsets (see below).  Determines
  // whether integer overflow is impossible on assignment from an
  // InputType to an OutputType.
  //
  // Implicit here is the assumption that both input and output types
  // are integers.
  template<class OutputType, class InputType>
  struct OutputCanFitInput {
  private:
    static constexpr bool output_signed = std::is_signed<OutputType>::value;
    static constexpr bool input_signed = std::is_signed<InputType>::value;

  public:
    static const bool value = sizeof (OutputType) > sizeof (InputType) ||
      (sizeof (OutputType) == sizeof (InputType) &&
       ! output_signed && input_signed);
  };

  // Avoid warnings for "unsigned integer < 0" comparisons.
  template<class InputType,
           bool input_signed = std::is_signed<InputType>::value>
  struct Negative {};

  template<class InputType>
  struct Negative<InputType, true> {
    static KOKKOS_INLINE_FUNCTION bool
    negative (const InputType src) {
      return src < InputType (0);
    }
  };

  template<class InputType>
  struct Negative<InputType, false> {
    static KOKKOS_INLINE_FUNCTION bool
    negative (const InputType /* src */) {
      return false;
    }
  };

  template<class InputType>
  KOKKOS_INLINE_FUNCTION bool negative (const InputType src) {
    return Negative<InputType>::negative (src);
  }

  template<class OutputType, class InputType>
  struct OverflowChecker {
  private:
    static constexpr bool output_signed = std::is_signed<OutputType>::value;
    static constexpr bool input_signed = std::is_signed<InputType>::value;

  public:
    // 1. Signed to unsigned could overflow due to negative numbers.
    // 2. Larger to smaller could overflow.
    // 3. Same size but unsigned to signed could overflow.
    static constexpr bool could_overflow =
      (! output_signed && input_signed) ||
      (sizeof (OutputType) < sizeof (InputType)) ||
      (sizeof (OutputType) == sizeof (InputType) &&
       output_signed && ! input_signed);

    KOKKOS_INLINE_FUNCTION bool
    overflows (const InputType src) const
    {
      if (! could_overflow) {
        return false;
      }
      else {
        // Signed to unsigned could overflow due to negative numbers.
        if (! output_signed && input_signed) {
          return negative (src);
        }
        // We're only comparing InputType with InputType here, so this
        // should not emit warnings.
        return src < minDstVal_ || src > maxDstVal_;
      }
    }

  private:
    // If InputType is unsigned and OutputType is signed, casting max
    // OutputType to InputType could overflow.  See #5548.
    InputType minDstVal_ = input_signed ?
      std::numeric_limits<OutputType>::min () : OutputType (0);
    InputType maxDstVal_ = std::numeric_limits<OutputType>::max ();
  };


  template<class OutputViewType, class InputViewType>
  void
  errorIfOverflow (const OutputViewType& dst,
                   const InputViewType& src,
                   const size_t overflowCount)
  {
    if (overflowCount == 0) {
      return;
    }

    std::ostringstream os;
    const bool plural = overflowCount != size_t (1);
    os << "copyOffsets: " << overflowCount << " value" <<
      (plural ? "s" : "") << " in src were too big (in the "
      "sense of integer overflow) to fit in dst.";

    const bool verbose = Details::Behavior::verbose ();
    if (verbose) {
      const size_t maxNumToPrint =
        Details::Behavior::verbosePrintCountThreshold();
      const size_t srcLen (src.extent (0));
      if (srcLen <= maxNumToPrint) {
        auto dst_h = Kokkos::create_mirror_view (dst);
        auto src_h = Kokkos::create_mirror_view (src);
        // DEEP_COPY REVIEW - NOT TESTED
        Kokkos::deep_copy (src_h, src);
        // DEEP_COPY REVIEW - NOT TESTED
        Kokkos::deep_copy (dst_h, dst);

        os << "  src: [";
        for (size_t k = 0; k < srcLen; ++k) {
          os << src_h[k];
          if (k + size_t (1) < srcLen) {
            os << ", ";
          }
        }
        os << "], ";

        os << " dst: [";
        for (size_t k = 0; k < srcLen; ++k) {
          os << dst_h[k];
          if (k + size_t (1) < srcLen) {
            os << ", ";
          }
        }
        os << "].";
      }
      else {
        os << "  src.extent(0) > " << maxNumToPrint << ", Tpetra's "
          "verbose print count threshold.  To increase this, set the "
          "environment variable TPETRA_VERBOSE_PRINT_COUNT_THRESHOLD "
          "to the desired threshold and rerun.  You do NOT need to "
          "rebuild Trilinos.";
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, os.str ());
  }

  // Implementation detail of copyOffsets (see below).
  //
  // Kokkos parallel_reduce functor for copying offset ("ptr") arrays.
  // Tpetra::Details::FixedHashTable uses this in its "copy"
  // constructor for converting between different Device types.  All
  // the action happens in the partial specializations for different
  // values of outputCanFitInput.  "Output can fit input" means that
  // casting the input's value type to the output's value type will
  // never result in integer overflow.
  template<class OutputViewType,
           class InputViewType,
           const bool outputCanFitInput =
             OutputCanFitInput<typename OutputViewType::non_const_value_type,
                               typename InputViewType::non_const_value_type>::value>
  class CopyOffsetsFunctor {};

  // Specialization for when overflow is possible.
  template<class OutputViewType, class InputViewType>
  class CopyOffsetsFunctor<OutputViewType, InputViewType, false> {
  public:
    using execution_space = typename OutputViewType::execution_space;
    using size_type = typename OutputViewType::size_type;
    using value_type = size_t;

    using input_value_type = typename InputViewType::non_const_value_type;
    using output_value_type = typename OutputViewType::non_const_value_type;

    CopyOffsetsFunctor (const OutputViewType& dst, const InputViewType& src) :
      dst_ (dst), src_ (src)
    {
      static_assert (Kokkos::SpaceAccessibility<
                       typename OutputViewType::memory_space,
                       typename InputViewType::memory_space>::accessible,
                     "CopyOffsetsFunctor (implements copyOffsets): Output "
                     "View's space must be able to access the input View's "
                     "memory space.");
    }

    KOKKOS_INLINE_FUNCTION void
    operator () (const size_type i, value_type& overflowCount) const {
      const input_value_type src_i = src_(i);
      if (checker_.overflows (src_i)) {
        ++overflowCount;
      }
      dst_(i) = static_cast<output_value_type> (src_i);
    }

    KOKKOS_INLINE_FUNCTION void
    operator () (const size_type i) const {
      const input_value_type src_i = src_(i);
      dst_(i) = static_cast<output_value_type> (src_i);
    }

    KOKKOS_INLINE_FUNCTION void init (value_type& overflowCount) const {
      overflowCount = 0;
    }

    KOKKOS_INLINE_FUNCTION void
    join (value_type& result,
          const value_type& current) const {
      result += current;
    }

  private:
    OutputViewType dst_;
    InputViewType src_;
    OverflowChecker<output_value_type, input_value_type> checker_;
  };

  // Specialization for when overflow is impossible.
  template<class OutputViewType, class InputViewType>
  class CopyOffsetsFunctor<OutputViewType, InputViewType, true> {
  public:
    using execution_space = typename OutputViewType::execution_space;
    using size_type = typename OutputViewType::size_type;
    using value_type = size_t;

    CopyOffsetsFunctor (const OutputViewType& dst, const InputViewType& src) :
      dst_ (dst),
      src_ (src)
    {
      static_assert (Kokkos::SpaceAccessibility<
                       typename OutputViewType::memory_space,
                       typename InputViewType::memory_space>::accessible,
                     "CopyOffsetsFunctor (implements copyOffsets): Output "
                     "View's space must be able to access the input View's "
                     "memory space.");
    }

    KOKKOS_INLINE_FUNCTION void
    operator () (const size_type i, value_type& /* overflowCount */) const {
      // Overflow is impossible in this case, so there's no need to check.
      dst_(i) = src_(i);
    }

    KOKKOS_INLINE_FUNCTION void
    operator () (const size_type i) const {
      dst_(i) = src_(i);
    }

    KOKKOS_INLINE_FUNCTION void init (value_type& overflowCount) const {
      overflowCount = 0;
    }

    KOKKOS_INLINE_FUNCTION void
    join (value_type& /* result */,
          const value_type& /* current */) const
    {}

  private:
    OutputViewType dst_;
    InputViewType src_;
  };

  // Implementation detail of copyOffsets (see below).
  //
  // We specialize copyOffsets on two different conditions:
  //
  // 1. Are the two Views' layouts the same, and do the input and
  //    output Views have the same value type?
  // 2. Can the output View's execution space access the input View's
  //    memory space?
  //
  // If (1) is true, that makes the implementation simple: just call
  // Kokkos::deep_copy (FixedHashTable always uses the same layout, no
  // matter the device type).  Otherwise, we need a custom copy
  // functor.  If (2) is true, then we can use CopyOffsetsFunctor
  // directly.  Otherwise, we have to copy the input View into the
  // output View's memory space, before we can use the functor.
  //
  template<class OutputViewType,
           class InputViewType,
           const bool sameLayoutsSameOffsetTypes =
             std::is_same<typename OutputViewType::array_layout,
                          typename InputViewType::array_layout>::value &&
             std::is_same<typename OutputViewType::non_const_value_type,
                          typename InputViewType::non_const_value_type>::value,
           const bool outputExecSpaceCanAccessInputMemSpace =
             Kokkos::SpaceAccessibility<
               typename OutputViewType::memory_space,
               typename InputViewType::memory_space>::accessible>
  struct CopyOffsetsImpl {
    static void run (const OutputViewType& dst, const InputViewType& src);
  };

  // Specialization for sameLayoutsSameOffsetTypes = true:
  //
  // If both input and output Views have the same layout, and both
  // input and output use the same type for offsets, then we don't
  // need to check for overflow, and we can use Kokkos::deep_copy
  // directly.  It doesn't matter whether the output execution space
  // can access the input memory space: Kokkos::deep_copy takes care
  // of the details.
  template<class OutputViewType,
           class InputViewType,
           const bool outputExecSpaceCanAccessInputMemSpace>
  struct CopyOffsetsImpl<OutputViewType, InputViewType,
                         true, outputExecSpaceCanAccessInputMemSpace> {
    static void run (const OutputViewType& dst, const InputViewType& src) {
      static_assert (std::is_same<typename OutputViewType::non_const_value_type,
                       typename InputViewType::non_const_value_type>::value,
                     "CopyOffsetsImpl (implementation of copyOffsets): In order"
                     " to call this specialization, the input and output must "
                     "use the same offset type.");
      static_assert (static_cast<int> (OutputViewType::rank) ==
                     static_cast<int> (InputViewType::rank),
                     "CopyOffsetsImpl (implementation of copyOffsets): In order"
                     " to call this specialization, src and dst must have the "
                     "same rank.");
      static_assert (std::is_same<typename OutputViewType::array_layout,
                       typename InputViewType::array_layout>::value,
                     "CopyOffsetsImpl (implementation of copyOffsets): In order"
                     " to call this specialization, src and dst must have the "
                     "the same array_layout.");
      // DEEP_COPY REVIEW - DEVICE-TO-DEVICE
      using execution_space = typename OutputViewType::execution_space;
      Kokkos::deep_copy (execution_space(), dst, src);
    }
  };

  // Specializations for sameLayoutsSameOffsetTypes = false:
  //
  // If input and output don't have the same layout, or use different
  // types for offsets, then we can't use Kokkos::deep_copy directly,
  // and we may have to check for overflow.

  // Specialization for sameLayoutsSameOffsetTypes = false and
  // outputExecSpaceCanAccessInputMemSpace = true:
  //
  // If the output execution space can access the input memory space,
  // then we can use CopyOffsetsFunctor directly.
  template<class OutputViewType,
           class InputViewType>
  struct CopyOffsetsImpl<OutputViewType, InputViewType,
                         false, true> {
    static void run (const OutputViewType& dst, const InputViewType& src) {
      static_assert (static_cast<int> (OutputViewType::rank) ==
                     static_cast<int> (InputViewType::rank),
                     "CopyOffsetsImpl (implementation of copyOffsets): "
                     "src and dst must have the same rank.");
      constexpr bool sameLayoutsSameOffsetTypes =
        std::is_same<typename OutputViewType::array_layout,
          typename InputViewType::array_layout>::value &&
        std::is_same<typename OutputViewType::non_const_value_type,
          typename InputViewType::non_const_value_type>::value;
      static_assert (! sameLayoutsSameOffsetTypes,
                     "CopyOffsetsImpl (implements copyOffsets): In order to "
                     "call this specialization, sameLayoutsSameOffsetTypes "
                     "must be false.  That is, either the input and output "
                     "must have different array layouts, or their value types "
                     "must differ.");
      static_assert (Kokkos::SpaceAccessibility<
                       typename OutputViewType::memory_space,
                       typename InputViewType::memory_space>::accessible,
                     "CopyOffsetsImpl (implements copyOffsets): In order to "
                     "call this specialization, the output View's space must "
                     "be able to access the input View's memory space.");
      using functor_type = CopyOffsetsFunctor<OutputViewType, InputViewType>;
      using execution_space = typename OutputViewType::execution_space;
      using size_type = typename OutputViewType::size_type;
      using range_type = Kokkos::RangePolicy<execution_space, size_type>;

      const bool debug = Details::Behavior::debug ();
      if (debug) {
        size_t overflowCount = 0; // output argument of the reduction
        Kokkos::parallel_reduce ("Tpetra::Details::copyOffsets",
                                 range_type (0, dst.extent (0)),
                                 functor_type (dst, src),
                                 overflowCount);
        errorIfOverflow (dst, src, overflowCount);
      }
      else {
        Kokkos::parallel_for ("Tpetra::Details::copyOffsets",
                              range_type (0, dst.extent (0)),
                              functor_type (dst, src));
      }
    }
  };

  // Specialization for sameLayoutsSameOffsetTypes = false and
  // outputExecSpaceCanAccessInputMemSpace = false.
  //
  // If the output execution space canNOT access the input memory
  // space, then we can't use CopyOffsetsFunctor directly.  Instead,
  // tell Kokkos to copy the input View's data into the output View's
  // memory space _first_.  Since the offset types are different for
  // this specialization, we can't just call Kokkos::deep_copy
  // directly between the input and output Views of offsets; that
  // wouldn't compile.
  //
  // This case can and does come up in practice: If the output View's
  // execution space is Cuda, it cannot currently access host memory
  // (that's the opposite direction from what UVM allows).
  // Furthermore, that case specifically requires overflow checking,
  // since (as of 28 Jan 2016 at least) Kokkos::Cuda uses a smaller
  // offset type than Kokkos' host spaces.
  template<class OutputViewType, class InputViewType>
  struct CopyOffsetsImpl<OutputViewType, InputViewType,
                         false, false> {
    static void run (const OutputViewType& dst, const InputViewType& src) {
      static_assert (static_cast<int> (OutputViewType::rank) ==
                     static_cast<int> (InputViewType::rank),
                     "CopyOffsetsImpl (implementation of copyOffsets): In order"
                     " to call this specialization, src and dst must have the "
                     "same rank.");
      constexpr bool sameLayoutsSameOffsetTypes =
        std::is_same<typename OutputViewType::array_layout,
          typename InputViewType::array_layout>::value &&
        std::is_same<typename OutputViewType::non_const_value_type,
          typename InputViewType::non_const_value_type>::value;
      static_assert (! sameLayoutsSameOffsetTypes,
                     "CopyOffsetsImpl (implements copyOffsets): In order to "
                     "call this specialization, sameLayoutsSameOffsetTypes "
                     "must be false.  That is, either the input and output "
                     "must have different array layouts, or their value types "
                     "must differ.");
      using output_space_copy_type =
        Kokkos::View<typename InputViewType::non_const_value_type*,
          Kokkos::LayoutLeft, typename OutputViewType::device_type>;
      using Kokkos::view_alloc;
      using Kokkos::WithoutInitializing;
      using execution_space = typename OutputViewType::execution_space;
      output_space_copy_type
        outputSpaceCopy (view_alloc ("outputSpace", WithoutInitializing),
                         src.extent (0));
      // DEEP_COPY REVIEW - DEVICE-TO-DEVICE
      Kokkos::deep_copy (execution_space(), outputSpaceCopy, src);

      // The output View's execution space can access
      // outputSpaceCopy's data, so we can run the functor now.
      using functor_type =
        CopyOffsetsFunctor<OutputViewType, output_space_copy_type>;
      using size_type = typename OutputViewType::size_type;
      using range_type = Kokkos::RangePolicy<execution_space, size_type>;

      const bool debug = Details::Behavior::debug ();
      if (debug) {
        size_t overflowCount = 0;
        Kokkos::parallel_reduce ("Tpetra::Details::copyOffsets",
                                 range_type (0, dst.extent (0)),
                                 functor_type (dst, outputSpaceCopy),
                                 overflowCount);
        errorIfOverflow (dst, src, overflowCount);
      }
      else {
        Kokkos::parallel_for ("Tpetra::Details::copyOffsets",
                              range_type (0, dst.extent (0)),
                              functor_type (dst, outputSpaceCopy));
      }
    }
  };
} // namespace (anonymous)

/// \brief Copy row offsets (in a sparse graph or matrix) from src
///   to dst.  The offsets may have different types.
///
/// The implementation reserves the right to do bounds checking if the
/// offsets in the two arrays have different types.
///
/// Everything above is an implementation detail of this function,
/// copyOffsets.  This function in turn is an implementation detail
/// of FixedHashTable, in particular of the "copy constructor" that
/// copies a FixedHashTable from one Kokkos device to another.
/// copyOffsets copies the array of offsets (ptr_).
template<class OutputViewType, class InputViewType>
void
copyOffsets (const OutputViewType& dst, const InputViewType& src)
{
  static_assert (Kokkos::is_view<OutputViewType>::value,
                 "OutputViewType (the type of dst) must be a Kokkos::View.");
  static_assert (Kokkos::is_view<InputViewType>::value,
                 "InputViewType (the type of src) must be a Kokkos::View.");
  static_assert (std::is_same<typename OutputViewType::value_type,
                   typename OutputViewType::non_const_value_type>::value,
                 "OutputViewType (the type of dst) must be a nonconst Kokkos::View.");
  static_assert (static_cast<int> (OutputViewType::rank) == 1,
                 "OutputViewType (the type of dst) must be a rank-1 Kokkos::View.");
  static_assert (static_cast<int> (InputViewType::rank) == 1,
                 "InputViewType (the type of src) must be a rank-1 Kokkos::View.");
  static_assert (std::is_integral<typename std::decay<decltype (dst(0)) >::type>::value,
                 "The entries of dst must be built-in integers.");
  static_assert (std::is_integral<typename std::decay<decltype (src(0)) >::type>::value,
                 "The entries of src must be built-in integers.");

  TEUCHOS_TEST_FOR_EXCEPTION
    (dst.extent (0) != src.extent (0), std::invalid_argument,
     "copyOffsets: dst.extent(0) = " << dst.extent (0)
     << " != src.extent(0) = " << src.extent (0) << ".");

  CopyOffsetsImpl<OutputViewType, InputViewType>::run (dst, src);
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_COPYOFFSETS_HPP
