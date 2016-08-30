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

#ifndef KOKKOS_SPARSE_IMPL_ASSIGNORCOPYINTEGERS_HPP_
#define KOKKOS_SPARSE_IMPL_ASSIGNORCOPYINTEGERS_HPP_

#include "TpetraKernels_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Sparse_impl_copyIntegers.hpp"

namespace KokkosSparse {
namespace Impl {

/// \brief Either copy or assign input to output, where both input and
///   output are 1-D Kokkos Views of integers
///
/// \param out [in/out] Output: 1-D View of integers
/// \param in [in] Input: 1-D View of integers
/// \param reuseGraph [in] Whether to reuse data in the output array
///
/// If it is legal to assign the input View to the output View (e.g.,
/// if they are Views of the same type, with \c out possibly being a
/// const View), then do so.  Otherwise, allocate the output View if
/// necessary, and copy the contents of the input View into the output
/// View.
///
/// If the output View is a const View and we can't assign directly,
/// allocate a new nonconst View, copy the contents of the input into
/// the new nonconst View, then assign the result to the output View.
///
/// This function will attempt to reuse output storage if it has
/// already been allocated.  However, it will always copy the data,
/// unless you tell it to "reuse the graph."  If you tell it to reuse
/// the graph, then it will not copy if the output allocation already
/// exists.
///
/// This function reserves the right to check for integer overflow if
/// necessary, and to throw on overflow.  This is useful if the output
/// View stores smaller integers than the input View, for example, int
/// vs. size_t.
template<class OutputType, class InputType>
void
assignOrCopyIntegers (OutputType& out,
                      const InputType& in,
                      const bool reuseGraph)
{
  constexpr bool outputViewIsNonConst =
    std::is_same<typename OutputType::non_const_value_type,
      typename OutputType::value_type>::value;
  constexpr bool inputViewIsConst =
    std::is_same<typename InputType::const_value_type,
      typename InputType::value_type>::value;
  // FIXME (mfh 27 Aug 2016) Kokkos lacks an "assignable" predicate.
  // Also, this does not account for memory traits.
  constexpr bool assignable =
    std::is_same<typename OutputType::non_const_value_type,
      typename InputType::non_const_value_type>::value &&
    std::is_same<typename OutputType::array_layout,
      typename InputType::array_layout>::value &&
    std::is_same<typename OutputType::device_type,
      typename InputType::device_type>::value &&
    ! (outputViewIsNonConst && inputViewIsConst);

  // This achieves the desired goal if 'in' is assignable to 'out'.
  // Otherwise, it does nothing (just assigns 'out' to 'out').
  out = Kokkos::Impl::if_c<assignable,
    InputType, OutputType>::select (in, out);
  if (assignable) {
    return; // of course no overflow; we just assigned 'in' to 'out'
  }
  else {
    if (reuseGraph && out.dimension_0 () == in.dimension_0 ()) {
      return; // just reuse the existing data in 'out'
    }
    else {
      // 'out' may be const.  In that case, we need to make a temp
      // output array, copy into that, then assign to 'out'.  If 'out'
      // is NOT const, we can use it directly, as long as it is long
      // enough.
      typedef typename OutputType::non_const_type nc_output_type;
      nc_output_type out_tmp =
        Kokkos::Impl::if_c<outputViewIsNonConst,
          OutputType,
          nc_output_type>::select (out, nc_output_type (out.label (), in.dimension_0 ()));
      if (out_tmp.dimension_0 () != in.dimension_0 ()) {
        // If we get here, that means that out is NOT const, and is
        // not long enough.  Repeat the above allocation step.  In
        // this case, assume that 'out' may not have a label yet.
        const std::string label = (out.label () == "") ?
          in.label () : out.label ();
        out_tmp = nc_output_type (label, in.dimension_0 ());
      }
      // At this point, 'out' had better have the right length.
      if (out.dimension_0 () != in.dimension_0 ()) {
        throw std::logic_error ("assignOrCopyIntegers: Output array still does "
                                "not have the right length.  Please report this"
                                " bug to the KokkosKernels developers.");
      }
      copyIntegers (out_tmp, in);
      out = out_tmp;
    }
  }
}

} // namespace Impl
} // namespace KokkosSparse

#endif // KOKKOS_SPARSE_IMPL_ASSIGNORCOPYINTEGERS_HPP_
