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
        Tpetra::MultiVector<SC, LO, GO, NT>, MemorySpace, am> >
    {
    public:
      using local_access_type =
        LocalAccess<Tpetra::MultiVector<SC, LO, GO, NT>, MemorySpace, am>;
    private:
      using global_object_type = Tpetra::MultiVector<SC, LO, GO, NT>;

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

      // MemorySpace=CudaSpace: false.
      // MemorySpace=CudaUVMSpace: false.
      // MemorySpace=CudaHostPinnedSpace: true.
      // MemorySpace=HostSpace: true.
      static constexpr bool is_host =
        std::is_same<
          typename MemorySpace::execution_space::memory_space,
          Kokkos::HostSpace>::value;

    public:
      // This alias is for the MultiVector specialization of
      // GetNonowningLocalObject.  withLocalAccess itself does not
      // need this to be public.
      //
      // Owning View type is always a View of nonconst.  If you own
      // the data, you need to be able to modify them.
      using master_local_view_type = typename std::conditional<
        is_host,
        typename dual_view_type::t_host,
        typename dual_view_type::t_dev>::type;

      static_assert
      (static_cast<int> (master_local_view_type::Rank) == 2,
       "Rank of master_local_view_type must be 2.  "
       "Please report this bug to the Tpetra developers.");

      // This alias is required by withLocalAccess.
      using master_local_object_type =
        std::unique_ptr<master_local_view_type>;

      // This method is required by withLocalAccess.
      static master_local_object_type
      get (local_access_type LA)
      {
        if (LA.isValid ()) {
          // Intel 17.0.1 requires the static_cast.  Otherwise, you'll
          // get build errors of the form "error: a built-in binary
          // operator applied to a scoped enumeration requires two
          // operands of the same type."
          if (static_cast<Details::AccessMode> (am) ==
              Details::AccessMode::WriteOnly) {
            LA.G_.clear_sync_state ();
          }

          // The various templated methods want an execution space
          // rather than a memory space.  Otherwise, DualView of
          // CudaUVMSpace complains that HostSpace is not one of its
          // two memory spaces.  (Both the device and the host Views
          // of a DualView of CudaUVMSpace have memory_space =
          // CudaUVMSpace.)
          using execution_space = typename MemorySpace::execution_space;

          if (LA.G_.template need_sync<execution_space> ()) {
            LA.G_.template sync<execution_space> ();
          }
          // Intel 17.0.1 requires the static_cast.  Otherwise, you'll
          // get build errors of the form "error: a built-in binary
          // operator applied to a scoped enumeration requires two
          // operands of the same type."
          if (static_cast<Details::AccessMode> (am) !=
              Details::AccessMode::ReadOnly) {
            LA.G_.template modify<execution_space> ();
          }

          // See note about "copy back" above.
          auto G_lcl_2d = LA.G_.template getLocalView<execution_space> ();
          // This converts the View to const if applicable.
          // Once we can use C++14, switch to std::make_unique.
          return std::unique_ptr<master_local_view_type>
            (new master_local_view_type (G_lcl_2d));
        }
        else { // invalid; return "null" Kokkos::View
          return std::unique_ptr<master_local_view_type>
            (new master_local_view_type ());
        }
      }
    };

    //! Specialization of GetMasterLocalObject for Tpetra::Vector.
    template<class SC, class LO, class GO, class NT,
             class MemorySpace,
             const Details::AccessMode am>
    struct GetMasterLocalObject<
      LocalAccess<
        Tpetra::Vector<SC, LO, GO, NT>, MemorySpace, am> >
    {
    private:
      using global_object_type = Tpetra::Vector<SC, LO, GO, NT>;
      using parent_global_object_type = Tpetra::MultiVector<SC, LO, GO, NT>;
      using parent_local_access_type =
        LocalAccess<parent_global_object_type, MemorySpace, am>;
      using mv_gmlo = GetMasterLocalObject<parent_local_access_type>;
      using parent_master_local_view_type =
        typename mv_gmlo::master_local_view_type;

    public:
      using local_access_type =
        LocalAccess<global_object_type, MemorySpace, am>;

    public:
      // This alias is for the Vector specialization of
      // GetNonowningLocalObject.  withLocalAccess itself does not
      // need this to be public.
      using master_local_view_type = decltype (Kokkos::subview
        (parent_master_local_view_type (), Kokkos::ALL (), 0));

      static_assert
      (static_cast<int> (master_local_view_type::Rank) == 1,
       "Rank of master_local_view_type must be 1.  "
       "Please report this bug to the Tpetra developers.");

      // This alias is required by withLocalAccess.
      using master_local_object_type =
        std::unique_ptr<master_local_view_type>;

      // This method is required by withLocalAccess.
      static master_local_object_type
      get (local_access_type LA)
      {
        if (LA.isValid ()) {
          // Intel 17.0.1 requires the static_cast.  Otherwise, you'll
          // get build errors of the form "error: a built-in binary
          // operator applied to a scoped enumeration requires two
          // operands of the same type."
          if (static_cast<Details::AccessMode> (am) ==
              Details::AccessMode::WriteOnly) {
            LA.G_.clear_sync_state ();
          }

          // The various templated methods want an execution space
          // rather than a memory space.  Otherwise, DualView of
          // CudaUVMSpace complains that HostSpace is not one of its
          // two memory spaces.  (Both the device and the host Views
          // of a DualView of CudaUVMSpace have memory_space =
          // CudaUVMSpace.)
          using execution_space = typename MemorySpace::execution_space;

          if (LA.G_.template need_sync<execution_space> ()) {
            LA.G_.template sync<execution_space> ();
          }
          // Intel 17.0.1 requires the static_cast.  Otherwise, you'll
          // get build errors of the form "error: a built-in binary
          // operator applied to a scoped enumeration requires two
          // operands of the same type."
          if (static_cast<Details::AccessMode> (am) !=
              Details::AccessMode::ReadOnly) {
            LA.G_.template modify<execution_space> ();
          }

          // See note about "copy back" above.
          auto G_lcl_2d = LA.G_.template getLocalView<execution_space> ();
          auto G_lcl_1d = Kokkos::subview (G_lcl_2d, Kokkos::ALL (), 0);

          // This converts the View to const if applicable.
          // Once we can use C++14, switch to std::make_unique.
          return std::unique_ptr<master_local_view_type>
            (new master_local_view_type (G_lcl_1d));
        }
        else { // invalid; return "null" Kokkos::View
          return std::unique_ptr<master_local_view_type>
            (new master_local_view_type ());
        }
      }
    };

    /// \brief Specialization of GetNonowningLocalObject for
    ///   Tpetra::MultiVector.
    template<class SC, class LO, class GO, class NT,
             class MemorySpace,
             const AccessMode am>
    struct GetNonowningLocalObject<
      LocalAccess<
        Tpetra::MultiVector<SC, LO, GO, NT>, MemorySpace, am> >
    {
    public:
      using local_access_type = LocalAccess<
        Tpetra::MultiVector<SC, LO, GO, NT>, MemorySpace, am>;

    private:
      using input_view_type =
        typename GetMasterLocalObject<local_access_type>::master_local_view_type;
      // input_view_type::non_const_data_type is
      // MV::impl_scalar_type**, where
      // MV = Tpetra::MultiVector<SC, LO, GO, NT>.
      //
      // Intel 17.0.1 requires the static_cast.  Otherwise, you'll get
      // build errors of the form "error: a built-in binary operator
      // applied to a scoped enumeration requires two operands of the
      // same type."
      using output_data_type = typename std::conditional<
        static_cast<AccessMode> (am) == AccessMode::ReadOnly,
        typename input_view_type::const_data_type,
        typename input_view_type::non_const_data_type>::type;
      using output_view_type =
        Kokkos::View<output_data_type,
                     typename input_view_type::array_layout,
                     typename input_view_type::device_type,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
    public:
      using master_local_object_type =
        typename GetMasterLocalObject<local_access_type>::master_local_object_type;
      using nonowning_local_object_type = output_view_type;

      static nonowning_local_object_type
      get (local_access_type /* LA */,
           const master_local_object_type& M)
      {
        input_view_type* viewPtr = M.get ();
        return viewPtr == nullptr ?
          nonowning_local_object_type () :
          nonowning_local_object_type (*viewPtr);
      }
    };

    /// \brief Specialization of GetNonowningLocalObject for
    ///   Tpetra::Vector.
    template<class SC, class LO, class GO, class NT,
             class MemorySpace,
             const AccessMode am>
    struct GetNonowningLocalObject<
      LocalAccess<
        ::Tpetra::Vector<SC, LO, GO, NT>, MemorySpace, am> >
    {
    public:
      using local_access_type = LocalAccess<
        Tpetra::Vector<SC, LO, GO, NT>, MemorySpace, am>;

    private:
      using input_view_type =
        typename GetMasterLocalObject<local_access_type>::master_local_view_type;
      // input_view_type::non_const_data_type is V::impl_scalar_type*,
      // where V = Tpetra::Vector<SC, LO, GO, NT>.
      //
      // Intel 17.0.1 requires the static_cast.  Otherwise, you'll get
      // build errors of the form "error: a built-in binary operator
      // applied to a scoped enumeration requires two operands of the
      // same type."
      using output_data_type = typename std::conditional<
        static_cast<AccessMode> (am) == AccessMode::ReadOnly,
        typename input_view_type::const_data_type,
        typename input_view_type::non_const_data_type>::type;
      using output_view_type =
        Kokkos::View<output_data_type,
                     typename input_view_type::array_layout,
                     typename input_view_type::device_type,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
    public:
      using master_local_object_type =
        typename GetMasterLocalObject<local_access_type>::master_local_object_type;
      using nonowning_local_object_type = output_view_type;

      static nonowning_local_object_type
      get (local_access_type /* LA */,
           const master_local_object_type& M)
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

