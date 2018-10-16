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

#ifndef TPETRA_DETAILS_COPYCONVERT_HPP
#define TPETRA_DETAILS_COPYCONVERT_HPP

/// \file Tpetra_Details_convert.hpp
/// \brief Declare and define Tpetra::Details::copyConvert, an
///   implementation detail of Tpetra (in particular, of
///   FixedHashTable, CrsGraph, and CrsMatrix).

#include "TpetraCore_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Complex.hpp"
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace Tpetra {
namespace Details {

//
// Implementation details for copyConvert (see below).
// Users should skip over this anonymous namespace.
//
namespace { // (anonymous)

  template<class OutputValueType,
           class InputValueType>
  struct ConvertValue {
    static KOKKOS_INLINE_FUNCTION void
    convert (OutputValueType& dst, const InputValueType& src) {
      // This looks trivial, but it actually invokes OutputValueType's
      // constructor, so that needs to be marked as a __host__
      // __device__ function (e.g., via the KOKKOS_FUNCTION or
      // KOKKOS_INLINE_FUNCTION macros).
      dst = OutputValueType (src);
    }
  };

  template<class RealType>
  struct ConvertValue<RealType, Kokkos::complex<RealType> > {
    static KOKKOS_INLINE_FUNCTION void
    convert (RealType& dst, const Kokkos::complex<RealType>& src) {
      // RealType's constructor needs to be marked as a __host__
      // __device__ function (e.g., via the KOKKOS_FUNCTION or
      // KOKKOS_INLINE_FUNCTION macros).
      dst = RealType (src.real ());
    }
  };

  /// \brief Functor that helps implement copyConvert (see below).
  ///
  /// \tparam OutputViewType Type of the output Kokkos::View.
  /// \tparam InputViewType Type of the input Kokkos::View.
  template<class OutputViewType,
           class InputViewType>
  class CopyConvertFunctor {
  private:
    OutputViewType dst_;
    InputViewType src_;

  public:
    typedef typename std::decay<decltype (dst_[0])>::type output_type;
    typedef typename OutputViewType::size_type index_type;

    CopyConvertFunctor (const OutputViewType& dst, const InputViewType& src) :
      dst_ (dst),
      src_ (src)
    {
      // NOTE (mfh 29 Jan 2016): See kokkos/kokkos#178 for why we use
      // a memory space, rather than an execution space, as the first
      // argument of VerifyExecutionCanAccessMemorySpace.
      static_assert (Kokkos::Impl::VerifyExecutionCanAccessMemorySpace<
                       typename OutputViewType::memory_space,
                       typename InputViewType::memory_space>::value,
                     "CopyConvertFunctor (implements copyConvert): Output "
                     "View's space must be able to access the input View's "
                     "memory space.");
      static_assert (OutputViewType::Rank == 1 && InputViewType::Rank == 1,
                     "CopyConvertFunctor (implements copyConvert): "
                     "OutputViewType and InputViewType must be rank-1 "
                     "Kokkos::View specializations.");
    }

    KOKKOS_INLINE_FUNCTION void
    operator () (const index_type& i) const {
      using input_type = typename std::decay<decltype (src_[i])>::type;
      ConvertValue<output_type, input_type>::convert (dst_(i), src_(i));
    }
  };

  /// \brief Implementation detail of copyConvert (see below).
  ///
  /// We specialize copyConvert on two different conditions:
  ///
  /// 1. Can we just Kokkos::deep_copy from the input View to the
  ///    output View?  (That is, are the two Views' layouts the same,
  ///    and do the input and output Views have the same value type?)
  /// 2. Can the output View's execution space access the input View's
  ///    memory space?
  ///
  /// If (1) is true, that makes the implementation simple: just call
  /// Kokkos::deep_copy.  Otherwise, we need a custom copy functor
  /// (see CopyConvertFunctor above).
  ///
  /// Otherwise, if (2) is true, then we can use CopyConvertFunctor
  /// directly.  Otherwise, we have to copy the input View into the
  /// output View's memory space, before we can use the functor.
  ///
  /// NOTE (mfh 29 Jan 2016, 26 Jun 2017): See kokkos/kokkos#178 for
  /// why we use a memory space, rather than an execution space, as
  /// the first argument of VerifyExecutionCanAccessMemorySpace.
  template<class OutputViewType,
           class InputViewType,
           const bool canUseKokkosDeepCopy =
             std::is_same<typename OutputViewType::array_layout,
                          typename InputViewType::array_layout>::value &&
             std::is_same<typename OutputViewType::non_const_value_type,
                          typename InputViewType::non_const_value_type>::value,
           const bool outputExecSpaceCanAccessInputMemSpace =
             Kokkos::Impl::VerifyExecutionCanAccessMemorySpace<
               typename OutputViewType::memory_space,
               typename InputViewType::memory_space>::value>
  struct CopyConvertImpl {
    static void run (const OutputViewType& dst, const InputViewType& src);
  };

  // Specialization for canUseKokkosDeepCopy = true:
  //
  // If both input and output Views have the same layout, and both
  // input and output have the same type, then we can use
  // Kokkos::deep_copy directly.  It doesn't matter whether the output
  // execution space can access the input memory space:
  // Kokkos::deep_copy takes care of the details.
  template<class OutputViewType,
           class InputViewType,
           const bool outputExecSpaceCanAccessInputMemSpace>
  struct CopyConvertImpl<OutputViewType, InputViewType,
                         true, outputExecSpaceCanAccessInputMemSpace> {
    static void run (const OutputViewType& dst, const InputViewType& src) {
      static_assert (std::is_same<typename OutputViewType::non_const_value_type,
                       typename InputViewType::non_const_value_type>::value,
                     "CopyConvertImpl (implementation of copyConvert): In order"
                     " to call this specialization, the input and output must "
                     "use the same offset type.");
      static_assert (OutputViewType::Rank == 1 && InputViewType::Rank == 1,
                     "CopyConvertImpl (implementation of copyConvert): "
                     "OutputViewType and InputViewType must be rank-1 "
                     "Kokkos::View specializations.");
      static_assert (std::is_same<typename OutputViewType::array_layout,
                       typename InputViewType::array_layout>::value,
                     "CopyConvertImpl (implementation of copyConvert): In order"
                     " to call this specialization, src and dst must have the "
                     "the same array_layout.");
      Kokkos::deep_copy (dst, src);
    }
  };

  // Specialization for canUseKokkosDeepCopy = false and
  // outputExecSpaceCanAccessInputMemSpace = true:
  //
  // If the output execution space can access the input memory space,
  // then we can use CopyConvertFunctor directly.
  template<class OutputViewType,
           class InputViewType>
  struct CopyConvertImpl<OutputViewType,
                         InputViewType,
                         false,
                         true> {
    static void run (const OutputViewType& dst, const InputViewType& src) {
      static_assert (! std::is_same<typename OutputViewType::array_layout,
                         typename InputViewType::array_layout>::value ||
                     ! std::is_same<typename OutputViewType::non_const_value_type,
                         typename InputViewType::non_const_value_type>::value,
                     "CopyConvertImpl (implementation of copyConvert): We "
                     "should not be calling this specialization if "
                     "OutputViewType and InputViewType have the same entry "
                     "and layout types.");
      static_assert (OutputViewType::Rank == 1 && InputViewType::Rank == 1,
                     "CopyConvertImpl (implementation of copyConvert): "
                     "OutputViewType and InputViewType must both be rank-1 "
                     "Kokkos::View types.");
      // NOTE (mfh 29 Jan 2016): See kokkos/kokkos#178 for why we use
      // a memory space, rather than an execution space, as the first
      // argument of VerifyExecutionCanAccessMemorySpace.
      static_assert (Kokkos::Impl::VerifyExecutionCanAccessMemorySpace<
                       typename OutputViewType::memory_space,
                       typename InputViewType::memory_space>::value,
                     "CopyConvertImpl (implements copyConvert): In order to "
                     "call this specialization, the output View's space must "
                     "be able to access the input View's memory space.");

      typedef CopyConvertFunctor<OutputViewType, InputViewType> functor_type;
      typedef typename OutputViewType::execution_space execution_space;
      typedef typename OutputViewType::size_type index_type;
      typedef Kokkos::RangePolicy<execution_space, index_type> range_type;
      Kokkos::parallel_for (range_type (0, dst.extent (0)),
                            functor_type (dst, src));
    }
  };

  // Specialization for canUseKokkosDeepCopy = false and
  // outputExecSpaceCanAccessInputMemSpace = false.
  //
  // If the output execution space canNOT access the input memory
  // space, then we can't use CopyConvertFunctor directly.  Instead,
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
  template<class OutputViewType,
           class InputViewType>
  struct CopyConvertImpl<OutputViewType,
                         InputViewType,
                         false,
                         false> {
    static void run (const OutputViewType& dst, const InputViewType& src) {
      const bool canUseKokkosDeepCopy =
        std::is_same<typename OutputViewType::array_layout,
          typename InputViewType::array_layout>::value &&
        std::is_same<typename OutputViewType::non_const_value_type,
          typename InputViewType::non_const_value_type>::value;
      static_assert (! canUseKokkosDeepCopy,
                     "CopyConvertImpl (implementation of copyConvert): We "
                     "should not be calling this specialization if we could "
                     "have used Kokkos::deep_copy instead.");
      static_assert (OutputViewType::Rank == 1 && InputViewType::Rank == 1,
                     "CopyConvertImpl (implementation of copyConvert): "
                     "OutputViewType and InputViewType must both be rank-1 "
                     "Kokkos::View types.");

      using Kokkos::ViewAllocateWithoutInitializing;
      typedef Kokkos::View<typename InputViewType::non_const_value_type*,
        typename InputViewType::array_layout,
        typename OutputViewType::device_type> output_space_copy_type;
      output_space_copy_type
        outputSpaceCopy (ViewAllocateWithoutInitializing ("outputSpace"),
                         src.extent (0));
      Kokkos::deep_copy (outputSpaceCopy, src);

      // The output View's execution space can access
      // outputSpaceCopy's data, so we can run the functor now.
      typedef CopyConvertFunctor<OutputViewType,
        output_space_copy_type> functor_type;
      typedef typename OutputViewType::execution_space execution_space;
      typedef typename OutputViewType::size_type index_type;
      typedef Kokkos::RangePolicy<execution_space, index_type> range_type;
      Kokkos::parallel_for (range_type (0, dst.extent (0)),
                            functor_type (dst, outputSpaceCopy));
    }
  };
} // namespace (anonymous)

/// \brief Copy values from the 1-D Kokkos::View src, to the 1-D
///   Kokkos::View dst, of the same length.  The entries of src and
///   dst may have different types, but it must be possible to
///   copy-construct each entry of dst with its corresponding entry of
///   src.
///
/// Everything above is an implementation detail of this function,
/// copyConvert.
template<class OutputViewType,
         class InputViewType>
void
copyConvert (const OutputViewType& dst,
             const InputViewType& src)
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
  if (dst.extent (0) != src.extent (0)) {
    std::ostringstream os;
    os << "Tpetra::Details::copyConvert: "
       << "dst.extent(0) = " << dst.extent (0)
       << " != src.extent(0) = " << src.extent (0)
       << ".";
    throw std::invalid_argument (os.str ());
  }
  // Canonicalize the View types in order to avoid redundant instantiations.
  typedef typename OutputViewType::non_const_type output_view_type;
  typedef typename InputViewType::const_type input_view_type;
  CopyConvertImpl<output_view_type, input_view_type>::run (dst, src);
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_COPYCONVERT_HPP
