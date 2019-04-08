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
#include "Kokkos_ArithTraits.hpp"
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

  // We need separate implementations for both (T,complex) and
  // (complex,T), but we can't just overload for both cases, because
  // that would be ambiguous (e.g., (complex,complex)).
  template<class OutputValueType,
           class InputValueType,
           const bool outputIsComplex =
             Kokkos::ArithTraits<OutputValueType>::is_complex,
           const bool inputIsComplex =
             Kokkos::ArithTraits<InputValueType>::is_complex>
  struct ConvertValue
  {
    static KOKKOS_INLINE_FUNCTION void
    convert (OutputValueType& dst, const InputValueType& src)
    {
      // This looks trivial, but it actually invokes OutputValueType's
      // constructor, so that needs to be marked as a __host__
      // __device__ function (e.g., via the KOKKOS_FUNCTION or
      // KOKKOS_INLINE_FUNCTION macros).
      dst = OutputValueType (src);
    }
  };

  template<class OutputRealType, class InputComplexType>
  struct ConvertValue<OutputRealType, InputComplexType, false, true>
  {
    static KOKKOS_INLINE_FUNCTION void
    convert (OutputRealType& dst,
             const InputComplexType& src)
    {
      // OutputRealType's constructor needs to be marked with either
      // KOKKOS_FUNCTION or KOKKOS_INLINE_FUNCTION.
      using KAI = Kokkos::ArithTraits<InputComplexType>;
      dst = OutputRealType (KAI::real (src));
    }
  };

  template<class OutputComplexType, class InputRealType>
  struct ConvertValue<OutputComplexType, InputRealType, true, false>
  {
    static KOKKOS_INLINE_FUNCTION void
    convert (OutputComplexType& dst,
             const InputRealType& src)
    {
      // OutputComplexType's constructor needs to be marked with
      // either KOKKOS_FUNCTION or KOKKOS_INLINE_FUNCTION.
      using output_mag_type =
        typename Kokkos::ArithTraits<OutputComplexType>::mag_type;
      using KAM = Kokkos::ArithTraits<output_mag_type>;
      dst = OutputComplexType (src, KAM::zero ());
    }
  };
  
  template<class OutputValueType,
           class InputValueType>
  KOKKOS_INLINE_FUNCTION void
  convertValue (OutputValueType& dst, const InputValueType& src) {
    ConvertValue<OutputValueType, InputValueType>::convert (dst, src);
  }

  /// \brief Functor that helps implement copyConvert (see below).
  ///
  /// \tparam OutputViewType Type of the output Kokkos::View.
  /// \tparam InputViewType Type of the input Kokkos::View.
  template<class OutputViewType,
           class InputViewType,
           const int rank = static_cast<int> (OutputViewType::Rank)>
  class CopyConvertFunctor {};

  template<class OutputViewType,
           class InputViewType>
  class CopyConvertFunctor<OutputViewType, InputViewType, 1> {
  private:
    static_assert
    (static_cast<int> (OutputViewType::Rank) == 1 &&
     static_cast<int> (InputViewType::Rank) == 1,
     "CopyConvertFunctor (implements Tpetra::Details::copyConvert): "
     "OutputViewType and InputViewType must both have rank 1.");
    OutputViewType dst_;
    InputViewType src_;

  public:
    using index_type = typename OutputViewType::size_type;

    CopyConvertFunctor (const OutputViewType& dst,
                        const InputViewType& src) :
      dst_ (dst),
      src_ (src)
    {}

    KOKKOS_INLINE_FUNCTION void
    operator () (const index_type i) const {
      convertValue (dst_(i), src_(i));
    }
  };

  template<class OutputViewType,
           class InputViewType>
  class CopyConvertFunctor<OutputViewType, InputViewType, 2> {
  public:
    using index_type = typename OutputViewType::size_type;

  private:
    static_assert
    (static_cast<int> (OutputViewType::Rank) == 2 &&
     static_cast<int> (InputViewType::Rank) == 2,
     "CopyConvertFunctor (implements Tpetra::Details::copyConvert): "
     "OutputViewType and InputViewType must both have rank 2.");
    OutputViewType dst_;
    InputViewType src_;
    index_type numCols_;

  public:
    CopyConvertFunctor (const OutputViewType& dst,
                        const InputViewType& src) :
      dst_ (dst),
      src_ (src),
      numCols_ (dst.extent (1))
    {}

    KOKKOS_INLINE_FUNCTION void
    operator () (const index_type i) const {
      const index_type numCols = numCols_;
      for (index_type j = 0; j < numCols; ++j) {
        convertValue (dst_(i,j), src_(i,j));
      }
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
             std::is_same<typename OutputViewType::non_const_value_type,
                          typename InputViewType::non_const_value_type>::value,
           const bool outputExecSpaceCanAccessInputMemSpace =
             Kokkos::Impl::VerifyExecutionCanAccessMemorySpace<
               typename OutputViewType::memory_space,
               typename InputViewType::memory_space>::value>
  struct CopyConvertImpl {
    static void
    run (const OutputViewType& dst,
         const InputViewType& src);
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
                         true, outputExecSpaceCanAccessInputMemSpace>
  {
    static void
    run (const OutputViewType& dst,
         const InputViewType& src)
    {
      static_assert (static_cast<int> (OutputViewType::Rank) ==
                     static_cast<int> (InputViewType::Rank),
                     "CopyConvertImpl (implementation of copyConvert): "
                     "The two Views must have the same rank.");
      constexpr bool same_value_type =
        std::is_same<typename OutputViewType::non_const_value_type,
                     typename InputViewType::non_const_value_type>::value;
      static_assert (same_value_type,
                     "CopyConvertImpl (implementation of copyConvert): In "
                     "order to call this specialization, the input and output "
                     "Views must have the same value type.");
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
                         true>
  {
    static void
    run (const OutputViewType& dst,
         const InputViewType& src)
    {
      static_assert (static_cast<int> (OutputViewType::Rank) ==
                     static_cast<int> (InputViewType::Rank),
                     "CopyConvertImpl (implementation of copyConvert): "
                     "The two Views must have the same rank.");
      constexpr bool same_value_type =
        std::is_same<typename OutputViewType::non_const_value_type,
                     typename InputViewType::non_const_value_type>::value;
      static_assert (! same_value_type,
                     "CopyConvertImpl (implementation of copyConvert): We "
                     "should not be calling this specialization if the two "
                     "Views have the same value types.  It's correct to do "
                     "so, but could be slower.");
      // NOTE (mfh 29 Jan 2016): See kokkos/kokkos#178 for why we use
      // a memory space, rather than an execution space, as the first
      // argument of VerifyExecutionCanAccessMemorySpace.
      constexpr bool output_exec_space_can_access_input_mem_space =
        Kokkos::Impl::VerifyExecutionCanAccessMemorySpace<
          typename OutputViewType::memory_space,
          typename InputViewType::memory_space>::value;
      static_assert (output_exec_space_can_access_input_mem_space,
                     "CopyConvertImpl (implements copyConvert): In order to "
                     "call this specialization, the output View's execution "
                     "space must be able to access the input View's memory "
                     "space.");

      using functor_type = CopyConvertFunctor<OutputViewType, InputViewType>;
      using execution_space = typename OutputViewType::execution_space;
      using index_type = typename OutputViewType::size_type;
      using range_type = Kokkos::RangePolicy<execution_space, index_type>;
      Kokkos::parallel_for ("Tpetra::Details::copyConvert",
                            range_type (0, dst.extent (0)),
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
  struct CopyConvertImpl<OutputViewType, InputViewType, false, false>
  {
    static void
    run (const OutputViewType& dst,
         const InputViewType& src)
    {
      static_assert (static_cast<int> (OutputViewType::Rank) ==
                     static_cast<int> (InputViewType::Rank),
                     "CopyConvertImpl (implementation of copyConvert): "
                     "OutputViewType and InputViewType must have the "
                     "same rank.");
      const bool canUseKokkosDeepCopy =
        std::is_same<typename OutputViewType::non_const_value_type,
          typename InputViewType::non_const_value_type>::value;
      static_assert (! canUseKokkosDeepCopy,
                     "CopyConvertImpl (implementation of copyConvert): We "
                     "should not be calling this specialization if we could "
                     "have used Kokkos::deep_copy instead.");
      constexpr bool output_exec_space_can_access_input_mem_space =
        Kokkos::Impl::VerifyExecutionCanAccessMemorySpace<
          typename OutputViewType::memory_space,
          typename InputViewType::memory_space>::value;
      static_assert (! output_exec_space_can_access_input_mem_space,
                     "CopyConvertImpl (implements copyConvert): We should not "
                     "call this specialization if the output View's execution "
                     "space can access the input View's memory space.  It's "
                     "correct, but may be slow.");

      using output_memory_space = typename OutputViewType::memory_space;
      auto src_outputSpaceCopy =
        Kokkos::create_mirror_view (output_memory_space (), src);
      Kokkos::deep_copy (src_outputSpaceCopy, src);

      // The output View's execution space can access
      // outputSpaceCopy's data, so we can run the functor now.
      using output_space_copy_type = decltype (src_outputSpaceCopy);
      using functor_type =
        CopyConvertFunctor<OutputViewType, output_space_copy_type>;
      using execution_space = typename OutputViewType::execution_space;
      using index_type = typename OutputViewType::size_type;
      using range_type = Kokkos::RangePolicy<execution_space, index_type>;
      Kokkos::parallel_for (range_type (0, dst.extent (0)),
                            functor_type (dst, src_outputSpaceCopy));
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
                 "OutputViewType must be a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<InputViewType>::value,
                 "InputViewType must be a Kokkos::View.");
  static_assert (std::is_same<typename OutputViewType::value_type,
                   typename OutputViewType::non_const_value_type>::value,
                 "OutputViewType must be a nonconst Kokkos::View.");
  static_assert (static_cast<int> (OutputViewType::Rank) ==
                 static_cast<int> (InputViewType::Rank),
                 "src and dst must have the same rank.");

  if (dst.extent (0) != src.extent (0)) {
    std::ostringstream os;
    os << "Tpetra::Details::copyConvert: "
       << "dst.extent(0) = " << dst.extent (0)
       << " != src.extent(0) = " << src.extent (0)
       << ".";
    throw std::invalid_argument (os.str ());
  }
  if (static_cast<int> (OutputViewType::Rank) > 1 &&
      dst.extent (1) != src.extent (1)) {
    std::ostringstream os;
    os << "Tpetra::Details::copyConvert: "
       << "dst.extent(1) = " << dst.extent (1)
       << " != src.extent(1) = " << src.extent (1)
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
