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
#include "Tpetra_Details_temporaryViewUtils.hpp"
#include <limits>
#include <type_traits>

/// \file Tpetra_Details_allReduceView.hpp
/// \brief All-reduce a 1-D or 2-D Kokkos::View

namespace Tpetra {
namespace Details {

template<typename InputViewType, typename OutputViewType>
static void
allReduceRawContiguous (const OutputViewType& output,
                        const InputViewType& input,
                        const Teuchos::Comm<int>& comm)
{
  using Teuchos::outArg;
  using Teuchos::REDUCE_SUM;
  using Teuchos::reduceAll;
  using ValueType = typename InputViewType::non_const_value_type;
  size_t count = input.span();
  TEUCHOS_ASSERT( count <= size_t (INT_MAX) );
  if(isInterComm(comm) && input.data() == output.data())
  {
    //Can't do in-place collective on an intercomm,
    //so use a separate copy as the input.
    typename InputViewType::array_layout layout(input.extent(0), input.extent(1), input.extent(2), input.extent(3), input.extent(4), input.extent(5), input.extent(6), input.extent(7)); 
    Kokkos::View<typename InputViewType::non_const_data_type, typename InputViewType::array_layout, typename InputViewType::device_type>
      tempInput(Kokkos::ViewAllocateWithoutInitializing("tempInput"), layout);
    Kokkos::deep_copy(tempInput, input);
    reduceAll<int, ValueType> (comm, REDUCE_SUM, static_cast<int> (count),
        tempInput.data(), output.data());
  }
  else
    reduceAll<int, ValueType> (comm, REDUCE_SUM, static_cast<int> (count),
        input.data(), output.data());
}

/// \brief All-reduce from input Kokkos::View to output Kokkos::View.
///
/// The two Views may alias one another.
template<class InputViewType, class OutputViewType>
static void
allReduceView (const OutputViewType& output,
               const InputViewType& input,
               const Teuchos::Comm<int>& comm)
{
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
  using Layout = typename TempView::UnifiedContiguousLayout<InputViewType, OutputViewType>::type;
  //if one or both is already in the correct layout, toLayout returns the same view
  auto inputContig = TempView::toLayout<InputViewType, Layout>(input);
  auto outputContig = TempView::toLayout<InputViewType, Layout>(output);
  if(Tpetra::Details::Behavior::assumeMpiIsCudaAware())
  {
    allReduceRawContiguous(outputContig, inputContig, comm);
  }
  else
  {
    auto inputMPI = TempView::toMPISafe<decltype(inputContig), false>(inputContig);
    auto outputMPI = TempView::toMPISafe<decltype(outputContig), false>(outputContig);
    allReduceRawContiguous(outputMPI, inputMPI, comm);
    Kokkos::deep_copy(outputContig, outputMPI);
  }
  Kokkos::deep_copy(output, outputContig);
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_ALLREDUCEVIEW_HPP
