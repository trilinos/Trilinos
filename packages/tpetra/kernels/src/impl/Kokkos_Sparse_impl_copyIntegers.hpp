/*
//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
//@HEADER
*/

#ifndef KOKKOS_SPARSE_IMPL_COPYINTEGERS_HPP_
#define KOKKOS_SPARSE_IMPL_COPYINTEGERS_HPP_

#include "TpetraKernels_config.h"
#include "Kokkos_Core.hpp"
#include <sstream>
#include <stdexcept>

namespace KokkosSparse {
namespace Impl {

namespace { // (anonymous)

// Implementation detail of copyIntegers (see below).
//
// We only need to check for overflow if overflow is possible.
// Overflow is possible under any of the following conditions:
//
//   1. Input type is bigger than output type
//   2. Output type is unsigned and input type is signed (negative
//      numbers overflow unsigned, and the max unsigned X is bigger
//      than the max signed X)
//
// Implicit here is the assumption that both input and output types
// are integers.
template<class OutputIntegerType,
         class InputIntegerType>
struct OverflowIsPossible {
  static_assert (std::is_integral<OutputIntegerType>::value,
                 "OutputIntegerType must be a built-in integer type.");
  static_assert (std::is_integral<InputIntegerType>::value,
                 "InputIntegerType must be a built-in integer type.");
private:
  static constexpr bool unsignedToSigned =
    std::is_unsigned<InputIntegerType>::value &&
    std::is_signed<OutputIntegerType>::value;
  static constexpr bool signedToUnsigned =
    std::is_signed<InputIntegerType>::value &&
    std::is_unsigned<OutputIntegerType>::value;
  static constexpr bool signedToSigned =
    std::is_signed<InputIntegerType>::value &&
    std::is_signed<OutputIntegerType>::value;

  static constexpr bool outputSameSizeAsInput =
    sizeof (OutputIntegerType) == sizeof (InputIntegerType);
  static constexpr bool outputSmallerThanInput =
    sizeof (OutputIntegerType) < sizeof (InputIntegerType);

public:
  static constexpr bool must_test_max =
    (unsignedToSigned && outputSameSizeAsInput) || outputSmallerThanInput;
  static constexpr bool must_test_min =
    signedToUnsigned || (signedToSigned && outputSmallerThanInput);
  static constexpr bool value = must_test_min || must_test_max;
};

// Implementation detail of copyIntegers (see below).
template<class OutputViewType, class InputViewType>
class CopyIntegersFunctor {
public:
  typedef typename OutputViewType::execution_space execution_space;
  typedef typename OutputViewType::size_type size_type;
  typedef int value_type;

  typedef typename InputViewType::non_const_value_type input_value_type;
  typedef typename OutputViewType::non_const_value_type output_value_type;

private:
  //! Minimum input value that we may safely cast to an output value.
  static input_value_type getMinDstVal () {
    static constexpr bool signedToUnsigned =
      std::is_signed<input_value_type>::value &&
      std::is_unsigned<output_value_type>::value;
    static constexpr bool signedToSigned =
      std::is_signed<input_value_type>::value &&
      std::is_signed<output_value_type>::value;
    static constexpr bool outputSmallerThanInput =
      sizeof (output_value_type) < sizeof (input_value_type);

    if (signedToUnsigned) {
      // Casting signed to unsigned means that we must forbid negative
      // input values.
      return static_cast<input_value_type> (0);
    }
    else if (signedToSigned && outputSmallerThanInput) {
      // If both types are signed and the output is strictly smaller,
      // than min input is min output, cast to input's type.
      //
      // We've taken care of the signed to unsigned case above.  If
      // input is unsigned, then min is zero (covered below) and no
      // need to test anyway.  Thus, only the signed to signed case
      // matters.
      return static_cast<input_value_type> (std::numeric_limits<output_value_type>::min ());
    }
    else { // no need to test; see discussion in comments above
      return std::numeric_limits<input_value_type>::min ();
    }
  }

  //! Maximum input value that we may safely cast to an output value.
  static input_value_type getMaxDstVal ()
  {
    static constexpr bool unsignedToSigned =
      std::is_unsigned<input_value_type>::value &&
      std::is_signed<output_value_type>::value;
    static constexpr bool outputSameSizeAsInput =
      sizeof (output_value_type) == sizeof (input_value_type);
    static constexpr bool outputSmallerThanInput =
      sizeof (output_value_type) < sizeof (input_value_type);

    if (unsignedToSigned && outputSameSizeAsInput) {
      // For two types with same number of bits, max unsigned > max
      // signed.
      return static_cast<input_value_type> (std::numeric_limits<output_value_type>::max ());
    }
    else if (outputSmallerThanInput) {
      return static_cast<input_value_type> (std::numeric_limits<output_value_type>::max ());
    }
    else { // no need to test
      return std::numeric_limits<input_value_type>::max ();
    }
  }

public:
  CopyIntegersFunctor (const OutputViewType& dst, const InputViewType& src) :
    dst_ (dst),
    src_ (src),
    minDstVal_ (getMinDstVal ()),
    maxDstVal_ (getMaxDstVal ()),
    testMinVal_ (OverflowIsPossible<output_value_type, input_value_type>::must_test_min),
    testMaxVal_ (OverflowIsPossible<output_value_type, input_value_type>::must_test_max)
  {
    // NOTE (mfh 29 Jan 2016): See kokkos/kokkos#178 for why we use a
    // memory space, rather than an execution space, as the first
    // argument of VerifyExecutionCanAccessMemorySpace.
    static_assert (Kokkos::Impl::VerifyExecutionCanAccessMemorySpace<
                     typename OutputViewType::memory_space,
                     typename InputViewType::memory_space>::value,
                   "CopyIntegersFunctor (implements copyIntegers): Output "
                   "View's space must be able to access the input View's "
                   "memory space.");
  }

  KOKKOS_INLINE_FUNCTION void
  operator () (const size_type& i, value_type& overflow) const {
    const input_value_type src_i = src_(i);

    if (testMinVal_ && src_i < minDstVal_) {
      overflow = 1; // underflow
    }
    else if (testMaxVal_ && src_i > maxDstVal_) {
      overflow = 2; // overflow
    }
    dst_(i) = static_cast<output_value_type> (src_i);
  }

  KOKKOS_INLINE_FUNCTION void init (value_type& overflow) const {
    overflow = 0; // success (no overflow)
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& result,
        const volatile value_type& current) const {
    result = (current > result) ? current : result;
  }

private:
  OutputViewType dst_;
  InputViewType src_;
  input_value_type minDstVal_;
  input_value_type maxDstVal_;
  bool testMinVal_; // we need this because CUDA dislikes host constexprs
  bool testMaxVal_; // we need this because CUDA dislikes host constexprs
};

// Implementation detail of copyIntegers (see below).
//
// We specialize copyIntegers on two different conditions:
//
// 1. Are the two Views' layouts the same, and do the input and
//    output Views have the same value type?
// 2. Can the output View's execution space access the input View's
//    memory space?
//
// If (1) is true, that makes the implementation simple: just call
// Kokkos::deep_copy (FixedHashTable always uses the same layout, no
// matter the device type).  Otherwise, we need a custom copy
// functor.  If (2) is true, then we can use CopyIntegersFunctor
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
struct CopyIntegersImpl {
  static void run (const OutputViewType& dst, const InputViewType& src);
};

// Specialization for sameLayoutsSameOffsetTypes = true:
//
// If both input and output Views have the same layout, and both
// input and output use the same type for integers, then we don't
// need to check for overflow, and we can use Kokkos::deep_copy
// directly.  It doesn't matter whether the output execution space
// can access the input memory space; Kokkos::deep_copy takes care
// of the details.
template<class OutputViewType,
         class InputViewType,
         const bool outputExecSpaceCanAccessInputMemSpace>
struct CopyIntegersImpl<OutputViewType, InputViewType,
                        true, outputExecSpaceCanAccessInputMemSpace> {
  static void run (const OutputViewType& dst, const InputViewType& src) {
    static_assert (std::is_same<typename OutputViewType::non_const_value_type,
                   typename InputViewType::non_const_value_type>::value,
                   "CopyIntegersImpl (implementation of copyIntegers): In order"
                   " to call this specialization, the input and output must "
                   "use the same offset type.");
    static_assert (static_cast<int> (OutputViewType::rank) ==
                   static_cast<int> (InputViewType::rank),
                   "CopyIntegersImpl (implementation of copyIntegers): In order"
                   " to call this specialization, src and dst must have the "
                   "same rank.");
    static_assert (std::is_same<typename OutputViewType::array_layout,
                   typename InputViewType::array_layout>::value,
                   "CopyIntegersImpl (implementation of copyIntegers): In order"
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
// then we can use CopyIntegersFunctor directly.
template<class OutputViewType,
         class InputViewType>
struct CopyIntegersImpl<OutputViewType, InputViewType,
                       false, true> {
  static void run (const OutputViewType& dst, const InputViewType& src) {
    static_assert (static_cast<int> (OutputViewType::rank) ==
                   static_cast<int> (InputViewType::rank),
                   "CopyIntegersImpl (implementation of copyIntegers): "
                   "src and dst must have the same rank.");
    constexpr bool sameLayoutsSameOffsetTypes =
      std::is_same<typename OutputViewType::array_layout,
      typename InputViewType::array_layout>::value &&
      std::is_same<typename OutputViewType::non_const_value_type,
      typename InputViewType::non_const_value_type>::value;
    static_assert (! sameLayoutsSameOffsetTypes,
                   "CopyIntegersImpl (implements copyIntegers): In order to "
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
                   "CopyIntegersImpl (implements copyIntegers): In order to "
                   "call this specialization, the output View's space must "
                   "be able to access the input View's memory space.");
    typedef CopyIntegersFunctor<OutputViewType, InputViewType> functor_type;
    int overflow = 0; // output argument of the reduction
    Kokkos::parallel_reduce (dst.dimension_0 (),
                             functor_type (dst, src),
                             overflow);
    if (overflow == 1) {
      throw std::runtime_error ("copyIntegers: One or more values in src underflowed on cast into dst.");
    }
    else if (overflow == 2) {
      throw std::runtime_error ("copyIntegers: One or more values in src "
                                "overflowed on cast into dst.  Some values "
                                "may also have underflowed.");
    }
  }
};

// Specialization for sameLayoutsSameOffsetTypes = false and
// outputExecSpaceCanAccessInputMemSpace = false.
//
// If the output execution space canNOT access the input memory
// space, then we can't use CopyIntegersFunctor directly.  Instead,
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
struct CopyIntegersImpl<OutputViewType, InputViewType,
                       false, false> {
  static void run (const OutputViewType& dst, const InputViewType& src) {
    static_assert (static_cast<int> (OutputViewType::rank) ==
                   static_cast<int> (InputViewType::rank),
                   "CopyIntegersImpl (implementation of copyIntegers): In order"
                   " to call this specialization, src and dst must have the "
                   "same rank.");
    constexpr bool sameLayoutsSameOffsetTypes =
      std::is_same<typename OutputViewType::array_layout,
      typename InputViewType::array_layout>::value &&
      std::is_same<typename OutputViewType::non_const_value_type,
      typename InputViewType::non_const_value_type>::value;
    static_assert (! sameLayoutsSameOffsetTypes,
                   "CopyIntegersImpl (implements copyIntegers): In order to "
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
                       src.dimension_0 ());
    Kokkos::deep_copy (outputSpaceCopy, src);

    // The output View's execution space can access
    // outputSpaceCopy's data, so we can run the functor now.
    typedef CopyIntegersFunctor<OutputViewType,
      output_space_copy_type> functor_type;
    int overflow = 0;
    Kokkos::parallel_reduce (dst.dimension_0 (),
                             functor_type (dst, outputSpaceCopy),
                             overflow);
    if (overflow == 1) {
      throw std::runtime_error ("copyIntegers: One or more values in "
                                "src underflowed on cast into dst.");
    }
    else if (overflow == 2) {
      throw std::runtime_error ("copyIntegers: One or more values in "
       "src overflowed on cast into dst.  Some values may also have "
       "underflowed.");
    }
  }
};

} // namespace (anonymous)

/// \brief Copy integers from one 1-D View to another.
///   The integers may have different types.
///
/// The implementation reserves the right to do overflow checking if
/// the integers in the two arrays have different types.
///
/// Everything above is an implementation detail of this function,
/// copyIntegers.
template<class OutputViewType, class InputViewType>
void
copyIntegers (const OutputViewType& dst, const InputViewType& src)
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

  if (dst.dimension_0 () != src.dimension_0 ()) {
    std::ostringstream os;
    os << "copyIntegers: dst.dimension_0() = " << dst.dimension_0 ()
       << " != src.dimension_0() = " << src.dimension_0 () << ".";
    throw std::invalid_argument (os.str ());
  }
  CopyIntegersImpl<OutputViewType, InputViewType>::run (dst, src);
}

} // namespace Impl
} // namespace KokkosSparse

#endif // KOKKOS_SPARSE_IMPL_COPYINTEGERS_HPP_
