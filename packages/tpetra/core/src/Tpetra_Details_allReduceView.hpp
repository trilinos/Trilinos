// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    // DEEP_COPY REVIEW - This could be either DEVICE-TO-DEVICE or HOST-TO-HOST
    // Either way, MPI is called right afterwards, meaning we'd need a sync on device
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

  // using execution_space = typename OutputViewType::execution_space;
  const bool viewsAlias = output.data () == input.data ();
  if (comm.getSize () == 1) {
    if (! viewsAlias) {
      // InputViewType and OutputViewType can't be AnonymousSpace
      // Views, because deep_copy needs to know their memory spaces.
      // DEEP_COPY REVIEW - NOT TESTED
      Kokkos::deep_copy (output, input);
    }
    return;
  }

  // we must ensure MPI can handle the pointers we pass it
  // if GPUAware, we are done
  // otherwise, if the views use GPUs, then we should copy them
  using Layout = typename TempView::UnifiedContiguousLayout<InputViewType, OutputViewType>::type;
  //if one or both is already in the correct layout, toLayout returns the same view
  auto inputContig = TempView::toLayout<InputViewType, Layout>(input);
  auto outputContig = TempView::toLayout<InputViewType, Layout>(output);
  if(Tpetra::Details::Behavior::assumeMpiIsGPUAware())
  {
    allReduceRawContiguous(outputContig, inputContig, comm);

  }
  else
  {
    auto inputMPI = TempView::toMPISafe<decltype(inputContig), false>(inputContig);
    auto outputMPI = TempView::toMPISafe<decltype(outputContig), false>(outputContig);
    allReduceRawContiguous(outputMPI, inputMPI, comm);
    // DEEP_COPY REVIEW - Could be either
    Kokkos::deep_copy(outputContig, outputMPI);
  }
    // DEEP_COPY REVIEW - Could be either
  Kokkos::deep_copy(output, outputContig);
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_ALLREDUCEVIEW_HPP
