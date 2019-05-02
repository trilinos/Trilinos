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

#ifndef TPETRA_WITHLOCALACCESS_MULTIVECTOR_HPP
#define TPETRA_WITHLOCALACCESS_MULTIVECTOR_HPP

#include "Tpetra_withLocalAccess.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include <memory>

/// \file Tpetra_withLocalAccess_MultiVector.hpp
/// \brief Include this file to make Tpetra::MultiVector and
///   Tpetra::Vector work with Tpetra::withLocalAccess.

namespace Tpetra {
  namespace Details {

    //! Specialization of GetMasterLocalObject for Tpetra::MultiVector.
    template<class SC, class LO, class GO, class NT,
             class MemorySpace,
             const AccessMode am>
    struct GetMasterLocalObject<
      LocalAccess<
        Tpetra::MultiVector<SC, LO, GO, NT>, MemorySpace, am> > {
    public:
      using local_access_type =
        LocalAccess<Tpetra::MultiVector<SC, LO, GO, NT>, MemorySpace, am>;
    private:
      using global_object_type =
        typename local_access_type::global_object_type;
      using memory_space = typename local_access_type::memory_space;
      static constexpr AccessMode access_mode =
        local_access_type::access_mode;
      using non_const_value_type =
        typename global_object_type::impl_scalar_type;
      using value_type = typename std::conditional<
          access_mode == AccessMode::ReadOnly,
          const non_const_value_type,
          non_const_value_type
        >::type;

      // FIXME (mfh 22 Oct 2018, 25 Apr 2019) Need to make sure that
      // the execution space matches.  If not, we would need to
      // allocate a new View, and then we should actually make the
      // std::unique_ptr's destructor "copy back."  This is why
      // master_local_object_type is a std::unique_ptr<view_type>, not
      // just a view_type.
      //
      // mfh 01 May 2019: For now, we avoid allocation and copy back,
      // by using only the Views available in the MV's DualView.
      using dual_view_type = typename global_object_type::dual_view_type;
      static constexpr bool is_host =
        std::is_same<memory_space, Kokkos::HostSpace>::value;
      using result_device_type = typename std::conditional<
        is_host,
        typename dual_view_type::t_host::device_type,
        typename dual_view_type::t_dev::device_type>::type;
      using view_type = Kokkos::View<
        value_type**,
        typename dual_view_type::t_dev::array_layout,
        result_device_type>;

    public:
      using master_local_object_type = std::unique_ptr<view_type>;

      static master_local_object_type
      get (local_access_type LA)
      {
        if (LA.isValid ()) {
          if (access_mode == Details::AccessMode::WriteOnly) {
            LA.G_.clear_sync_state ();
          }

          // The various templated methods want an execution space
          // rather than a memory space.  Otherwise, DualView of
          // CudaUVMSpace complains that HostSpace is not one of its
          // two memory spaces.  (Both the device and the host Views
          // of a DualView of CudaUVMSpace have memory_space =
          // CudaUVMSpace.)
          using execution_space = typename memory_space::execution_space;

          if (LA.G_.template need_sync<execution_space> ()) {
            LA.G_.template sync<execution_space> ();
          }
          if (access_mode != Details::AccessMode::ReadOnly) {
            LA.G_.template modify<execution_space> ();
          }

          // See note about "copy back" above.
          auto G_lcl_2d = LA.G_.template getLocalView<execution_space> ();
          // This converts the View to const if applicable.
          return std::unique_ptr<view_type> (new view_type (G_lcl_2d));
        }
        else { // invalid; return "null" Kokkos::View
          return std::unique_ptr<view_type> (new view_type ());
        }
      }
    };

    //! Specialization of GetMasterLocalObject for Tpetra::Vector.
    template<class SC, class LO, class GO, class NT,
             class MemorySpace,
             const Details::AccessMode am>
    struct GetMasterLocalObject<
      LocalAccess<
        Tpetra::Vector<SC, LO, GO, NT>, MemorySpace, am> > {
    public:
      using local_access_type =
        LocalAccess<Tpetra::Vector<SC, LO, GO, NT>, MemorySpace, am>;
    private:
      using global_object_type =
        typename local_access_type::global_object_type;
      using memory_space = typename local_access_type::memory_space;
      static constexpr AccessMode access_mode =
        local_access_type::access_mode;
      using non_const_value_type =
        typename global_object_type::impl_scalar_type;
      using value_type = typename std::conditional<
          access_mode == AccessMode::ReadOnly,
          const non_const_value_type,
          non_const_value_type
        >::type;

      // FIXME (mfh 22 Oct 2018, 25 Apr 2019) Need to make sure that
      // the execution space matches.  If not, we would need to
      // allocate a new View, and then we should actually make the
      // std::unique_ptr's destructor "copy back."  This is why
      // master_local_object_type is a std::unique_ptr<view_type>, not
      // just a view_type.
      //
      // mfh 01 May 2019: For now, we avoid allocation and copy back,
      // by using only the Views available in the MV's DualView.
      using dual_view_type = typename global_object_type::dual_view_type;
      static constexpr bool is_host =
        std::is_same<memory_space, Kokkos::HostSpace>::value;
      using result_device_type = typename std::conditional<
        is_host,
        typename dual_view_type::t_host::device_type,
        typename dual_view_type::t_dev::device_type>::type;
      using view_type = Kokkos::View<
        value_type*,
        typename dual_view_type::t_dev::array_layout,
        result_device_type>;

    public:
      using master_local_object_type = std::unique_ptr<view_type>;

      static master_local_object_type
      get (local_access_type LA)
      {
        if (LA.isValid ()) {
          if (access_mode == Details::AccessMode::WriteOnly) {
            LA.G_.clear_sync_state ();
          }

          // The various templated methods want an execution space
          // rather than a memory space.  Otherwise, DualView of
          // CudaUVMSpace complains that HostSpace is not one of its
          // two memory spaces.  (Both the device and the host Views
          // of a DualView of CudaUVMSpace have memory_space =
          // CudaUVMSpace.)
          using execution_space = typename memory_space::execution_space;

          if (LA.G_.template need_sync<execution_space> ()) {
            LA.G_.template sync<execution_space> ();
          }
          if (access_mode != Details::AccessMode::ReadOnly) {
            LA.G_.template modify<execution_space> ();
          }

          // See note about "copy back" above.
          auto G_lcl_2d = LA.G_.template getLocalView<execution_space> ();
          auto G_lcl_1d = Kokkos::subview (G_lcl_2d, Kokkos::ALL (), 0);
          // This converts the View to const if applicable.
          return std::unique_ptr<view_type> (new view_type (G_lcl_1d));
        }
        else { // invalid; return "null" Kokkos::View
          return std::unique_ptr<view_type> (new view_type ());
        }
      }
    };

    /// \brief Specialization of GetNonowningLocalObject for Kokkos::View.
    ///
    /// This is meant for the result of GetMasterLocalObject::get for
    /// a Tpetra::MultiVector or Tpetra::Vector.  DataType is
    /// <tt>impl_scalar_type**</tt> for a MultiVector, and
    /// <tt>impl_scalar_type*</tt> for a Vector.
    template<class DataType,
             class LayoutType,
             class MemorySpace>
    struct GetNonowningLocalObject<
      std::unique_ptr<
        Kokkos::View<DataType, LayoutType, MemorySpace>>>
    {
    private:
      using input_view_type =
        Kokkos::View<DataType, LayoutType, MemorySpace>;
      using output_view_type =
        Kokkos::View<DataType,
                     LayoutType,
                     MemorySpace,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    public:
      using master_local_object_type = std::unique_ptr<input_view_type>;
      using nonowning_local_object_type = output_view_type;

      static nonowning_local_object_type
      get (const master_local_object_type& M)
      {
        input_view_type* viewPtr = M.get ();
        return viewPtr == nullptr ?
          nonowning_local_object_type () :
          nonowning_local_object_type (*viewPtr);
      }
    };
  } // namespace Details
} // namespace Tpetra

#endif // TPETRA_WITHLOCALACCESS_MULTIVECTOR_HPP

