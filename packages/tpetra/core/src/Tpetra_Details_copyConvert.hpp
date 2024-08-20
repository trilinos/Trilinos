// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_COPYCONVERT_HPP
#define TPETRA_DETAILS_COPYCONVERT_HPP

/// \file Tpetra_Details_copyConvert.hpp
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
           const int rank = static_cast<int> (OutputViewType::rank)>
  class CopyConvertFunctor {};

  template<class OutputViewType,
           class InputViewType>
  class CopyConvertFunctor<OutputViewType, InputViewType, 1> {
  private:
    static_assert
    (static_cast<int> (OutputViewType::rank) == 1 &&
     static_cast<int> (InputViewType::rank) == 1,
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
    (static_cast<int> (OutputViewType::rank) == 2 &&
     static_cast<int> (InputViewType::rank) == 2,
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

  //! Whether copyConvert can just use Kokkos::deep_copy.
  template<class OutputViewType, class InputViewType>
  class CanUseKokkosDeepCopy {
  private:
    static constexpr bool sameValueType =
      std::is_same<typename OutputViewType::non_const_value_type,
                   typename InputViewType::non_const_value_type>::value;
    static constexpr bool sameMemorySpace =
      std::is_same<typename OutputViewType::memory_space,
                   typename InputViewType::memory_space>::value;
    static constexpr bool sameLayout =
      std::is_same<typename OutputViewType::array_layout,
                   typename InputViewType::array_layout>::value;

  public:
    static constexpr bool value =
      sameValueType && (sameMemorySpace || sameLayout);
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
  template<class OutputViewType,
           class InputViewType,
           const bool canUseKokkosDeepCopy =
             CanUseKokkosDeepCopy<OutputViewType, InputViewType>::value,
           const bool outputExecSpaceCanAccessInputMemSpace =
             Kokkos::SpaceAccessibility<
               typename OutputViewType::memory_space,
               typename InputViewType::memory_space>::accessible>
  struct CopyConvertImpl {
    static void
    run (const OutputViewType& dst,
         const InputViewType& src);
  };

  //! Specialization for canUseKokkosDeepCopy = true.
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
      // NOTE: It's important to do the addition _inside_ the
      // reinterpret-cast.  If you reinterpret_cast the separate
      // results, you may get the wrong answer (e.g., because
      // ptrdiff_t is signed, and pointers may have arbitrary 64-bit
      // virtual addresses).  I'm speaking from experience here.
      const ptrdiff_t dst_beg =reinterpret_cast<ptrdiff_t> (dst.data ());
      const ptrdiff_t dst_end =
        reinterpret_cast<ptrdiff_t> (dst.data () + dst.span ());
      const ptrdiff_t src_beg = reinterpret_cast<ptrdiff_t> (src.data ());
      const ptrdiff_t src_end =
        reinterpret_cast<ptrdiff_t> (src.data () + src.span ());

      if (dst_end > src_beg && src_end > dst_beg) {
        // dst and src alias each other, so we can't call
        // Kokkos::deep_copy(dst,src) directly (Kokkos detects this
        // and throws, at least in debug mode).  Instead, we make
        // temporary host storage (create_mirror always makes a new
        // allocation, unlike create_mirror_view).  Use host because
        // it's cheaper to allocate.  Hopefully users aren't doing
        // aliased copies in a tight loop.
        auto src_copy = Kokkos::create_mirror (Kokkos::HostSpace (), src);
        // DEEP_COPY REVIEW - NOT TESTED
        Kokkos::deep_copy (src_copy, src);
        // DEEP_COPY REVIEW - NOT TESTED
        Kokkos::deep_copy (dst, src_copy);
      }
      else { // no aliasing
        // DEEP_COPY REVIEW - DEVICE-TO-DEVICE
        using execution_space = typename OutputViewType::execution_space;
        Kokkos::deep_copy (execution_space(), dst, src);
      }
    }
  };

  /// \brief Specialization for canUseKokkosDeepCopy = false and
  ///   outputExecSpaceCanAccessInputMemSpace = true.
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
      using functor_type = CopyConvertFunctor<OutputViewType, InputViewType>;
      using execution_space = typename OutputViewType::execution_space;
      using index_type = typename OutputViewType::size_type;
      using range_type = Kokkos::RangePolicy<execution_space, index_type>;
      Kokkos::parallel_for ("Tpetra::Details::copyConvert",
                            range_type (0, dst.extent (0)),
                            functor_type (dst, src));
    }
  };

  /// \brief Specialization for canUseKokkosDeepCopy = false and
  ///   outputExecSpaceCanAccessInputMemSpace = false.
  ///
  /// This case can and does come up in practice: If the output View's
  /// execution space is Cuda, it cannot currently access host memory
  /// (that's the opposite direction from what UVM allows).
  template<class OutputViewType,
           class InputViewType>
  struct CopyConvertImpl<OutputViewType, InputViewType, false, false>
  {
    static void
    run (const OutputViewType& dst,
         const InputViewType& src)
    {
      using output_memory_space = typename OutputViewType::memory_space;
      using output_execution_space = typename OutputViewType::execution_space;
      auto src_outputSpaceCopy =
        Kokkos::create_mirror_view (output_memory_space (), src);
      // DEEP_COPY REVIEW - DEVICE-TO-HOSTMIRROR
      Kokkos::deep_copy (output_execution_space(), src_outputSpaceCopy, src);

      // The output View's execution space can access
      // outputSpaceCopy's data, so we can run the functor now.
      using output_space_copy_type = decltype (src_outputSpaceCopy);
      using functor_type =
        CopyConvertFunctor<OutputViewType, output_space_copy_type>;
      using execution_space = typename OutputViewType::execution_space;
      using index_type = typename OutputViewType::size_type;
      using range_type = Kokkos::RangePolicy<execution_space, index_type>;
      Kokkos::parallel_for ("Tpetra::Details::copyConvert",
                            range_type (0, dst.extent (0)),
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
  static_assert (Kokkos::is_view<OutputViewType>::value,
                 "OutputViewType must be a Kokkos::View.");
  static_assert (Kokkos::is_view<InputViewType>::value,
                 "InputViewType must be a Kokkos::View.");
  static_assert (std::is_same<typename OutputViewType::value_type,
                   typename OutputViewType::non_const_value_type>::value,
                 "OutputViewType must be a nonconst Kokkos::View.");
  static_assert (static_cast<int> (OutputViewType::rank) ==
                 static_cast<int> (InputViewType::rank),
                 "src and dst must have the same rank.");

  if (dst.extent (0) != src.extent (0)) {
    std::ostringstream os;
    os << "Tpetra::Details::copyConvert: "
       << "dst.extent(0) = " << dst.extent (0)
       << " != src.extent(0) = " << src.extent (0)
       << ".";
    throw std::invalid_argument (os.str ());
  }
  if (static_cast<int> (OutputViewType::rank) > 1 &&
      dst.extent (1) != src.extent (1)) {
    std::ostringstream os;
    os << "Tpetra::Details::copyConvert: "
       << "dst.extent(1) = " << dst.extent (1)
       << " != src.extent(1) = " << src.extent (1)
       << ".";
    throw std::invalid_argument (os.str ());
  }

  // Canonicalize the View types in order to avoid redundant instantiations.
  using output_view_type =
    Kokkos::View<typename OutputViewType::non_const_data_type,
                 typename OutputViewType::array_layout,
                 typename OutputViewType::device_type>;
  using input_view_type =
    Kokkos::View<typename InputViewType::const_data_type,
                 typename InputViewType::array_layout,
                 typename InputViewType::device_type>;
  CopyConvertImpl<output_view_type, input_view_type>::run (dst, src);
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_COPYCONVERT_HPP
