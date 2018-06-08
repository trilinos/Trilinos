/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#ifndef TPETRA_DETAILS_COPYOFFSETS_HPP
#define TPETRA_DETAILS_COPYOFFSETS_HPP

/// \file Tpetra_Details_copyOffsets.hpp
/// \brief Declare and define Tpetra::Details::copyOffsets, an
///   implementation detail of Tpetra (in particular, of
///   FixedHashTable, CrsGraph, and CrsMatrix).

#include "TpetraCore_config.h"
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

  // Implementation detail of copyOffsets (see below).
  //
  // Overflow is impossible (the output can fit the input) if the
  // output type is bigger than the input type, or if the types have
  // the same size and (the output type is unsigned, or both types are
  // signed).
  //
  // Implicit here is the assumption that both input and output types
  // are integers.
  template<class T1, class T2,
           const bool T1_is_signed = std::is_signed<T1>::value,
           const bool T2_is_signed = std::is_signed<T2>::value>
  struct OutputCanFitInput {
    static const bool value = sizeof (T1) > sizeof (T2) ||
      (sizeof (T1) == sizeof (T2) &&
       (std::is_unsigned<T1>::value || (std::is_signed<T1>::value && std::is_signed<T2>::value)));
  };

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
    typedef typename OutputViewType::execution_space execution_space;
    typedef typename OutputViewType::size_type size_type;
    typedef int value_type;

    typedef typename InputViewType::non_const_value_type input_value_type;
    typedef typename OutputViewType::non_const_value_type output_value_type;

    CopyOffsetsFunctor (const OutputViewType& dst, const InputViewType& src) :
      dst_ (dst),
      src_ (src),
      // We know that output_value_type cannot fit all values of
      // input_value_type, so an input_value_type can fit all values
      // of output_value_type.  This means we can convert from
      // output_value_type to input_value_type.  This is how we test
      // whether a given input_value_type value can fit in an
      // output_value_type.
      minDstVal_ (static_cast<input_value_type> (std::numeric_limits<output_value_type>::min ())),
      maxDstVal_ (static_cast<input_value_type> (std::numeric_limits<output_value_type>::max ()))
    {
      // NOTE (mfh 29 Jan 2016): See kokkos/kokkos#178 for why we use
      // a memory space, rather than an execution space, as the first
      // argument of VerifyExecutionCanAccessMemorySpace.
      static_assert (Kokkos::Impl::VerifyExecutionCanAccessMemorySpace<
                       typename OutputViewType::memory_space,
                       typename InputViewType::memory_space>::value,
                     "CopyOffsetsFunctor (implements copyOffsets): Output "
                     "View's space must be able to access the input View's "
                     "memory space.");
    }

    KOKKOS_INLINE_FUNCTION void
    operator () (const size_type& i, value_type& noOverflow) const {
      const input_value_type src_i = src_(i);
      if (src_i < minDstVal_ || src_i > maxDstVal_) {
        noOverflow = 0;
      }
      dst_(i) = static_cast<output_value_type> (src_i);
    }

    KOKKOS_INLINE_FUNCTION void init (value_type& noOverflow) const {
      noOverflow = 1; // success (no overflow)
    }

    KOKKOS_INLINE_FUNCTION void
    join (volatile value_type& result,
          const volatile value_type& current) const {
      result = (result>0 && current>0)?1:0; // was there any overflow?
    }

  private:
    OutputViewType dst_;
    InputViewType src_;
    input_value_type minDstVal_;
    input_value_type maxDstVal_;
  };

  // Specialization for when overflow is impossible.
  template<class OutputViewType, class InputViewType>
  class CopyOffsetsFunctor<OutputViewType, InputViewType, true> {
  public:
    typedef typename OutputViewType::execution_space execution_space;
    typedef typename OutputViewType::size_type size_type;
    typedef int value_type;

    CopyOffsetsFunctor (const OutputViewType& dst, const InputViewType& src) :
      dst_ (dst),
      src_ (src)
    {
      // NOTE (mfh 29 Jan 2016): See kokkos/kokkos#178 for why we use
      // a memory space, rather than an execution space, as the first
      // argument of VerifyExecutionCanAccessMemorySpace.
      static_assert (Kokkos::Impl::VerifyExecutionCanAccessMemorySpace<
                       typename OutputViewType::memory_space,
                       typename InputViewType::memory_space>::value,
                     "CopyOffsetsFunctor (implements copyOffsets): Output "
                     "View's space must be able to access the input View's "
                     "memory space.");
    }

    KOKKOS_INLINE_FUNCTION void
    operator () (const size_type& i, value_type& /* noOverflow */) const {
      // Overflow is impossible in this case, so there's no need to check.
      dst_(i) = src_(i);
    }

    KOKKOS_INLINE_FUNCTION void init (value_type& noOverflow) const {
      noOverflow = 1; // success (no overflow)
    }

    KOKKOS_INLINE_FUNCTION void
    join (volatile value_type& result,
          const volatile value_type& current) const {
      result = (result>0 && current>0)?1:0; // was there any overflow?
    }

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
  // NOTE (mfh 29 Jan 2016): See kokkos/kokkos#178 for why we use a
  // memory space, rather than an execution space, as the first
  // argument of VerifyExecutionCanAccessMemorySpace.
  template<class OutputViewType,
           class InputViewType,
           const bool sameLayoutsSameOffsetTypes =
             std::is_same<typename OutputViewType::array_layout,
                          typename InputViewType::array_layout>::value &&
             std::is_same<typename OutputViewType::non_const_value_type,
                          typename InputViewType::non_const_value_type>::value,
           const bool outputExecSpaceCanAccessInputMemSpace =
             Kokkos::Impl::VerifyExecutionCanAccessMemorySpace<
               typename OutputViewType::memory_space,
               typename InputViewType::memory_space>::value>
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
      Kokkos::deep_copy (dst, src);
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
      // NOTE (mfh 29 Jan 2016): See kokkos/kokkos#178 for why we use
      // a memory space, rather than an execution space, as the first
      // argument of VerifyExecutionCanAccessMemorySpace.
      static_assert (Kokkos::Impl::VerifyExecutionCanAccessMemorySpace<
                       typename OutputViewType::memory_space,
                       typename InputViewType::memory_space>::value,
                     "CopyOffsetsImpl (implements copyOffsets): In order to "
                     "call this specialization, the output View's space must "
                     "be able to access the input View's memory space.");
      typedef CopyOffsetsFunctor<OutputViewType, InputViewType> functor_type;
      typedef typename OutputViewType::execution_space execution_space;
      typedef typename OutputViewType::size_type size_type;
      typedef Kokkos::RangePolicy<execution_space, size_type> range_type;

      int noOverflow = 0; // output argument of the reduction
      Kokkos::parallel_reduce (range_type (0, dst.extent (0)),
                               functor_type (dst, src),
                               noOverflow);
      TEUCHOS_TEST_FOR_EXCEPTION
        (noOverflow==0, std::runtime_error, "copyOffsets: One or more values in "
         "src were too big (in the sense of integer overflow) to fit in dst.");
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

      typedef Kokkos::View<typename InputViewType::non_const_value_type*,
                           Kokkos::LayoutLeft,
                           typename OutputViewType::device_type>
        output_space_copy_type;
      using Kokkos::ViewAllocateWithoutInitializing;
      output_space_copy_type
        outputSpaceCopy (ViewAllocateWithoutInitializing ("outputSpace"),
                         src.extent (0));
      Kokkos::deep_copy (outputSpaceCopy, src);

      // The output View's execution space can access
      // outputSpaceCopy's data, so we can run the functor now.
      typedef CopyOffsetsFunctor<OutputViewType,
                                 output_space_copy_type> functor_type;
      typedef typename OutputViewType::execution_space execution_space;
      typedef typename OutputViewType::size_type size_type;
      typedef Kokkos::RangePolicy<execution_space, size_type> range_type;

      int noOverflow = 0;
      Kokkos::parallel_reduce (range_type (0, dst.extent (0)),
                               functor_type (dst, outputSpaceCopy),
                               noOverflow);
      TEUCHOS_TEST_FOR_EXCEPTION
        (noOverflow==0, std::runtime_error, "copyOffsets: One or more values "
         "in src were too big (in the sense of integer overflow) to fit in "
         "dst.");
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
  static_assert (Kokkos::Impl::is_view<OutputViewType>::value,
                 "OutputViewType (the type of dst) must be a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<InputViewType>::value,
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
