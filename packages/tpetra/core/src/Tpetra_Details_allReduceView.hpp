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

#ifndef TPETRA_DETAILS_ALLREDUCEVIEW_HPP
#define TPETRA_DETAILS_ALLREDUCEVIEW_HPP

#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_isInterComm.hpp"
#include "Kokkos_Core.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <limits>
#include <type_traits>

/// \file Tpetra_Details_allReduceView.hpp
/// \brief All-reduce a 1-D or 2-D Kokkos::View

namespace { // (anonymous)

// helper for detecting views that have Cuda memory
// Technically, UVM should not be here, but it seems we
// treat UVM as CudaSpace, so I have left it. With Cuda >= 8 on Linux
// UVM is handled by page faults in the Kernel, so it should always
// be addressable.
template<class ViewType>
struct view_uses_cuda_spaces {
  static constexpr bool value =
#ifdef KOKKOS_ENABLE_CUDA
      std::is_same<typename ViewType::memory_space, Kokkos::CudaSpace>::value
   || std::is_same<typename ViewType::memory_space, Kokkos::CudaUVMSpace>::value;
#else
    false;
#endif // KOKKOS_ENABLE_CUDA
};

template<class ViewType>
struct MakeContiguousBuffer {
  static constexpr bool is_contiguous_layout =
    std::is_same<
      typename ViewType::array_layout,
      Kokkos::LayoutLeft>::value ||
    std::is_same<
      typename ViewType::array_layout,
      Kokkos::LayoutRight>::value;
  using contiguous_array_layout =
    typename std::conditional<is_contiguous_layout,
                              typename ViewType::array_layout,
                              Kokkos::LayoutLeft>::type;
  using contiguous_device_type =
    typename std::conditional<
      std::is_same<
        typename ViewType::memory_space,
        Kokkos::HostSpace>::value,
      typename ViewType::device_type,
      Kokkos::HostSpace::device_type>::type;
  using contiguous_buffer_type =
    Kokkos::View<typename ViewType::non_const_data_type,
                 contiguous_array_layout,
                 contiguous_device_type>;

  static contiguous_array_layout
  makeLayout (const ViewType& view)
  {
    // NOTE (mfh 17 Mar 2019) This would be a good chance to use if
    // constexpr, once we have C++17.
    return contiguous_array_layout (view.extent (0), view.extent (1),
                                    view.extent (2), view.extent (3),
                                    view.extent (4), view.extent (5),
                                    view.extent (6), view.extent (7));
  }

  static contiguous_buffer_type
  make (const ViewType& view)
  {
    using Kokkos::view_alloc;
    using Kokkos::WithoutInitializing;
    return contiguous_buffer_type
      (view_alloc (view.label (), WithoutInitializing),
       makeLayout (view));
  }
};

template<class ViewType>
typename MakeContiguousBuffer<ViewType>::contiguous_buffer_type
makeContiguousBuffer (const ViewType& view)
{
  return MakeContiguousBuffer<ViewType>::make (view);
}

template<class ValueType>
static void
allReduceRawContiguous (ValueType output[],
                        const ValueType input[],
                        const size_t count,
                        const Teuchos::Comm<int>& comm)
{
  using Teuchos::outArg;
  using Teuchos::REDUCE_SUM;
  using Teuchos::reduceAll;
  constexpr size_t max_int = size_t (std::numeric_limits<int>::max ());
  TEUCHOS_ASSERT( count <= size_t (max_int) );
  reduceAll<int, ValueType> (comm, REDUCE_SUM, static_cast<int> (count),
                             input, output);
}

} // namespace (anonymous)

namespace Tpetra {
namespace Details {

/// \brief All-reduce from input Kokkos::View to output Kokkos::View.
///
/// The two Views may alias one another.
template<class InputViewType, class OutputViewType>
static void
allReduceView (const OutputViewType& output,
               const InputViewType& input,
               const Teuchos::Comm<int>& comm)
{
  // If all the right conditions hold, we may all-reduce directly from
  // the input to the output.  Here are the relevant conditions:
  //
  // - assumeMpiCanAccessBuffers: May we safely assume that MPI may
  //   read from the input View and write to the output View?  (Just
  //   because MPI _can_, doesn't mean that doing so will be faster.)
  // - Do input and output Views alias each other, and is the
  //   communicator an intercommunicator?  (Intercommunicators do not
  //   permit collectives to alias input and output buffers.)
  // - Is either View noncontiguous?
  //
  // If either View is noncontiguous, we could use MPI_Type_Vector to
  // create a noncontiguous MPI_Datatype, instead of packing and
  // unpacking to resp. from a contiguous temporary buffer.  Since
  // MPI_Allreduce requires that the input and output buffers both
  // have the same MPI_Datatype, this optimization might only work if
  // the MPI communicator is an intercommunicator.  Furthermore,
  // creating an MPI_Datatype instance may require memory allocation
  // anyway.  Thus, it's probably better just to use a temporary
  // contiguous buffer.  We use a host buffer for that, since device
  // buffers are slow to allocate.

  const bool viewsAlias = output.data () == input.data ();
  if (comm.getSize () == 1) {
    if (! viewsAlias) {
      // InputViewType and OutputViewType can't be AnonymousSpace
      // Views, because deep_copy needs to know their memory spaces.
      Kokkos::deep_copy (output, input);
    }
    return;
  }

  // we must esnure MPI can handle the pointers we pass it
  // if CudaAware, we are done
  // otherwise, if the views use Cuda, then we should copy them
  const bool mpiCannotAccessBuffers =
    // if assumeMpiIsCudaAware, then we can access cuda buffers
    ! ::Tpetra::Details::Behavior::assumeMpiIsCudaAware ()
    && (
         view_uses_cuda_spaces<OutputViewType>::value
         ||
         view_uses_cuda_spaces<InputViewType>::value 
        );

  const bool needContiguousTemporaryBuffers =
    // we must alloc/copy if MPI cannot access the buffers
    mpiCannotAccessBuffers                ||
    // If the comm is Inter and the views alias we must alloc/copy
    (::Tpetra::Details::isInterComm (comm)
     &&
     viewsAlias)                          ||
    // if either view is not contiguous then we must alloc/copy
    ! output.span_is_contiguous ()        ||
    ! input.span_is_contiguous ();

  if (needContiguousTemporaryBuffers) {
    auto output_tmp = makeContiguousBuffer (output);
    auto input_tmp = makeContiguousBuffer (input);
    Kokkos::deep_copy (input_tmp, input);
    // It's OK if LayoutLeft allocations have padding at the end of
    // each row.  MPI might write to those padding bytes, but it's
    // undefined behavior for users to use Kokkos to access whatever
    // bytes are there, and the bytes there don't need to define valid
    // ValueType instances.
    allReduceRawContiguous (output_tmp.data (), input_tmp.data (),
                            output_tmp.span (), comm);
    Kokkos::deep_copy (output, output_tmp);
  }
  else {
    allReduceRawContiguous (output.data (), input.data (),
                            output.span (), comm);
  }
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_ALLREDUCEVIEW_HPP
