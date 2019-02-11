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

#ifndef TPETRA_DETAILS_DUALVIEWUTIL_HPP
#define TPETRA_DETAILS_DUALVIEWUTIL_HPP

#include "TpetraCore_config.h"
#include "Kokkos_DualView.hpp"
#include "Teuchos_ArrayView.hpp"
#include <ostream>
#include <string>

//! Namespace for Tpetra classes and methods
namespace Tpetra {

/// \brief Namespace for Tpetra implementation details.
/// \warning Do NOT rely on the contents of this namespace.
namespace Details {

/// \brief Namespace for implementation details of Tpetra
///   implementation details.
/// \warning Do NOT rely on the contents of this namespace.
namespace Impl {

void throw_if_false (const bool cond, const char msg[]);

} // namespace Impl


auto view_alloc_no_init (const std::string& label) ->
  decltype (Kokkos::view_alloc (label, Kokkos::WithoutInitializing));

template<class ElementType, class DeviceType>
void
makeDualViewFromOwningHostView
  (Kokkos::DualView<ElementType*, DeviceType>& dv,
   const typename Kokkos::DualView<ElementType*, DeviceType>::t_host& hostView)
{
  using dual_view_type = Kokkos::DualView<ElementType*, DeviceType>;
  using dev_view_type = typename dual_view_type::t_dev;
  using host_view_type = typename dual_view_type::t_host;

  const auto size = hostView.size ();
  auto devView = Kokkos::create_mirror_view (DeviceType (), hostView);

#if defined(KOKKOS_ENABLE_CUDA)
  constexpr bool is_cuda = std::is_same<typename DeviceType::execution_space, Kokkos::Cuda>::value;
#else
  constexpr bool is_cuda = false;
#endif
  using Impl::throw_if_false;
  throw_if_false (! is_cuda || devView.data () != hostView.data (),
                  "If running with CUDA, then create_mirror_view needs to "
                  "return a View with a different pointer.  Please report "
                  "this bug to the Tpetra developers." );
  Kokkos::deep_copy (devView, hostView);
  dv = dual_view_type (devView, hostView);
}

template<class ElementType, class DeviceType>
void
makeDualViewFromArrayView (Kokkos::DualView<ElementType*, DeviceType>& dv,
                           const Teuchos::ArrayView<const ElementType>& av,
                           const std::string& label)
{
  using dual_view_type = Kokkos::DualView<ElementType*, DeviceType>;
  using dev_view_type = typename dual_view_type::t_dev;
  using host_view_type = typename dual_view_type::t_host;
  using const_host_view_type = typename host_view_type::const_type;

  const auto size = av.size ();
  const ElementType* ptr = (size == 0) ? nullptr : av.getRawPtr ();
  const_host_view_type inView (ptr, size);
  host_view_type hostView (view_alloc_no_init (label), size);
  Kokkos::deep_copy (hostView, inView);

  makeDualViewFromOwningHostView (dv, hostView);
}

template<class ElementType, class DeviceType>
void
makeDualViewFromVector (Kokkos::DualView<ElementType*, DeviceType>& dv,
                        const std::vector<ElementType>& vec,
                        const std::string& label)
{
  using dual_view_type = Kokkos::DualView<ElementType*, DeviceType>;
  using dev_view_type = typename dual_view_type::t_dev;
  using host_view_type = typename dual_view_type::t_host;
  using const_host_view_type = typename host_view_type::const_type;

  const auto size = vec.size ();
  const ElementType* ptr = (size == 0) ? nullptr : vec.data ();
  const_host_view_type inView (ptr, size);
  host_view_type hostView (view_alloc_no_init (label), size);
  Kokkos::deep_copy (hostView, inView);

  makeDualViewFromOwningHostView (dv, hostView);
}

template<class ElementType, class DeviceType>
void
printDualView (std::ostream& out,
               const Kokkos::DualView<ElementType*, DeviceType>& dv,
               const std::string& name)
{
  out << name << ": ";
  const size_t size = size_t (dv.extent (0));
  const auto hostView = dv.view_host ();

  out << "[";
  for (size_t k = 0; k < size; ++k) {
    out << hostView[k];
    if (k + size_t (1) < size) {
      out << ",";
    }
  }
  out << "]";
}

} // namespace Details

} // namespace Tpetra

#endif // TPETRA_DETAILS_DUALVIEWUTIL_HPP
