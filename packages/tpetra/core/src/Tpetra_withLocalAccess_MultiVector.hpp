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

    // We need these forward declarations so that LocalAccess knows
    // these partial specializations exist, and doesn't just defer to
    // the generic empty version.

    template<class SC, class LO, class GO, class NT, class ... Args>
    struct GetMasterLocalObject<
      LocalAccess<Tpetra::MultiVector<SC, LO, GO, NT>, Args...>
      >;
    template<class SC, class LO, class GO, class NT, class ... Args>
    struct GetMasterLocalObject<
      LocalAccess<Tpetra::Vector<SC, LO, GO, NT>, Args...>
      >;

    template<class SC, class LO, class GO, class NT, class ... Args>
    struct GetNonowningLocalObject<
      LocalAccess<Tpetra::MultiVector<SC, LO, GO, NT>, Args...>
      >;
    template<class SC, class LO, class GO, class NT, class ... Args>
    struct GetNonowningLocalObject<
      LocalAccess<Tpetra::Vector<SC, LO, GO, NT>, Args...>
      >;

    //////////////////////////////////////////////////////////////////

    template<class Space>
    struct is_host_space {
      // Space=CudaSpace: false.
      // Space=CudaUVMSpace: false.
      // Space=CudaHostPinnedSpace: true.
      // Space=HostSpace: true.
      //
      // Space=Cuda: true.
      // Space=OpenMP: false.
      // Space=Serial: false.
      // space=Threads: false.
      static constexpr bool value =
        std::is_same<typename Space::execution_space::memory_space,
                     Kokkos::HostSpace>::value;
    };

    //////////////////////////////////////////////////////////////////

    //! Specialization of GetMasterLocalObject for Tpetra::MultiVector.
    template<class SC, class LO, class GO, class NT, class ... Args>
    struct GetMasterLocalObject<
      LocalAccess<Tpetra::MultiVector<SC, LO, GO, NT>, Args...>
      >
    {
    private:
      using global_object_type = Tpetra::MultiVector<SC, LO, GO, NT>;

    public:
      using local_access_type =
        LocalAccess<global_object_type, Args...>;

    private:
      using execution_space =
        typename local_access_type::execution_space;
      static_assert(
        Kokkos::Impl::is_execution_space<execution_space>::value,
        "LocalAccess<Args...>::execution_space is not a valid "
        "Kokkos execution space.");

      using memory_space = typename local_access_type::memory_space;
      static_assert(
        Kokkos::Impl::is_memory_space<memory_space>::value,
        "LocalAccess<Args...>::memory_space is not a valid "
        "Kokkos memory space.");

      using access_mode = typename local_access_type::access_mode;
      static_assert(
        is_access_mode<access_mode>::value,
        "LocalAccess<Args...>::access_mode is not an Access "
        "type.");

      // FIXME (mfh 22 Oct 2018, 25 Apr 2019) Need to make sure that
      // the execution space matches.  If not, we would need to
      // allocate a new View, and then we should actually make the
      // std::unique_ptr's destructor "copy back."  This is why
      // master_local_object_type (see below) is a
      // std::unique_ptr<Kokkos::View<...>>, not just a
      // Kokkos::View<...>.
      //
      // mfh 01 May 2019: For now, we avoid allocation and copy back,
      // by using only the Views available in the MV's DualView.
      using dual_view_type =
        typename global_object_type::dual_view_type;

    public:
      // This alias is for the Tpetra::MultiVector specialization of
      // GetNonowningLocalObject.  withLocalAccess itself does not
      // need this to be public.
      //
      // Owning View type is always a View of nonconst.  If you own
      // the data, you need to be able to modify them.
      using master_local_view_type = typename std::conditional<
        is_host_space<execution_space>::value,
        typename dual_view_type::t_host,
        typename dual_view_type::t_dev>::type;

      static_assert(
        static_cast<int>(master_local_view_type::Rank) == 2,
        "Rank of master_local_view_type must be 2.  "
        "Please report this bug to the Tpetra developers.");

    private:
      static master_local_view_type
      getOwningPreActions(local_access_type& LA)
      {
        if (! LA.isValid()) { // return "null" Kokkos::View
          return master_local_view_type();
        }
        else {
          using access_mode = typename local_access_type::access_mode;
          constexpr bool is_write_only =
            std::is_same<access_mode, write_only>::value;
          if (is_write_only) {
            LA.G_.clear_sync_state();
          }

          // Given that Tpetra::(Multi)Vector currently uses
          // Kokkos::DualView, here is how get() must behave:
          //
          // - LA's memory space tells us which allocation to view.
          //
          // - LA's execution space tells us whether we need to fence.
          //
          //   - If LA's execution space equals the MultiVector's
          //     execution space, then there is no need to fence, no
          //     matter what the requested memory space is.
          //
          //   - Else, if LA's execution space is a host space, but
          //     the MultiVector needs sync to host, then we must
          //     fence the MultiVector's (device) execution space, to
          //     ensure that device kernels aren't concurrently
          //     modifying the MultiVector's local data.
          //     Tpetra::MultiVector::fence (should) do that for us,
          //     via Kokkos::DualView::fence.

          // It's easier to use an execution space than a memory space
          // in sync<Space>.  Otherwise, DualView of CudaUVMSpace
          // complains that HostSpace is not one of its two memory
          // spaces.  (Both the device and the host Views of a
          // DualView of CudaUVMSpace have memory_space =
          // CudaUVMSpace.)  Furthermore, this handles the case where
          // Kokkos::DefaultHostExecutionSpace != the MultiVector's
          // execution space, but the latter is still a host space
          // (e.g., Kokkos::OpenMP vs. Kokkos::Serial).

          if (LA.G_.template need_sync<execution_space>()) {
            LA.G_.template sync<execution_space>();
          }

          constexpr bool is_read_only =
            std::is_same<access_mode, read_only>::value;
          if (! is_read_only) {
            LA.G_.template modify<execution_space>();
          }

          // See note about "copy back" above.
          auto G_lcl_2d =
            LA.G_.template getLocalView<execution_space>();
          // This converts the View to const if applicable.
          return master_local_view_type(G_lcl_2d);
        }
      }

    public:
      // This alias is required by withLocalAccess.
      // using master_local_object_type = std::unique_ptr<
      //   master_local_view_type,
      //   typename impl_type::deleter_type<master_local_view_type>>;
      using master_local_object_type = std::unique_ptr<
        master_local_view_type>;

      // This method is required by withLocalAccess.
      static master_local_object_type
      get(local_access_type LA)
      {
        auto G_lcl_2d = getOwningPreActions(LA);
        // Once we can use C++14, switch to std::make_unique.
        // return master_local_object_type(
        //   new master_local_view_type(G_lcl_2d),
        //   impl_type::getOwningPostActions(LA, G_lcl_2d));
        return master_local_object_type(
          new master_local_view_type(G_lcl_2d));
      }
    };

    static_assert(
      Kokkos::is_view<
        GetMasterLocalObject<
          LocalAccess<Tpetra::MultiVector<>, read_only>
        >::master_local_view_type
      >::value, "Missing GetMasterLocalObject specialization");

    static_assert(
      Kokkos::is_view<
        GetMasterLocalObject<
          LocalAccess<Tpetra::MultiVector<>, Kokkos::HostSpace, read_only>
        >::master_local_view_type
      >::value, "Missing GetMasterLocalObject specialization");

    //////////////////////////////////////////////////////////////////

    //! Specialization of GetMasterLocalObject for Tpetra::Vector.
    template<class SC, class LO, class GO, class NT, class ... Args>
    struct GetMasterLocalObject<
      LocalAccess<Tpetra::Vector<SC, LO, GO, NT>, Args...>
      >
    {
    private:
      using global_object_type = Tpetra::Vector<SC, LO, GO, NT>;

    public:
      using local_access_type =
        LocalAccess<global_object_type, Args...>;

    private:
      using execution_space =
        typename local_access_type::execution_space;
      static_assert(
        Kokkos::Impl::is_execution_space<execution_space>::value,
        "LocalAccess<Args...>::execution_space is not a valid "
        "Kokkos execution space.");

      using memory_space = typename local_access_type::memory_space;
      static_assert(
        Kokkos::Impl::is_memory_space<memory_space>::value,
        "LocalAccess<Args...>::memory_space is not a valid "
        "Kokkos memory space.");

      using access_mode = typename local_access_type::access_mode;
      static_assert(
        is_access_mode<access_mode>::value,
        "LocalAccess<Args...>::access_mode is not an Access "
        "type.");

      // FIXME (mfh 22 Oct 2018, 25 Apr 2019) Need to make sure that
      // the execution space matches.  If not, we would need to
      // allocate a new View, and then we should actually make the
      // std::unique_ptr's destructor "copy back."  This is why
      // master_local_object_type (see below) is a
      // std::unique_ptr<Kokkos::View<...>>, not just a
      // Kokkos::View<...>.
      //
      // mfh 01 May 2019: For now, we avoid allocation and copy back,
      // by using only the Views available in the MV's DualView.
      using dual_view_type =
        typename global_object_type::dual_view_type;

      // Owning View type is always a View of nonconst.  If you own
      // the data, you need to be able to modify them.
      using master_local_view_2d_type = typename std::conditional<
        is_host_space<execution_space>::value,
        typename dual_view_type::t_host,
        typename dual_view_type::t_dev>::type;

    public:
      // This alias is for the Tpetra::Vector specialization of
      // GetNonowningLocalObject.  withLocalAccess itself does not
      // need this to be public.
      using master_local_view_type = decltype(
        Kokkos::subview(master_local_view_2d_type(),
          Kokkos::ALL(), 0));

      static_assert(
        static_cast<int>(master_local_view_type::Rank) == 1,
        "Rank of master_local_view_type must be 1.  "
        "Please report this bug to the Tpetra developers.");

    private:
      static master_local_view_2d_type
      getOwningPreActions(local_access_type& LA)
      {
        if (! LA.isValid()) { // return "null" Kokkos::View
          return master_local_view_2d_type();
        }
        else {
          using access_mode = typename local_access_type::access_mode;
          constexpr bool is_write_only =
            std::is_same<access_mode, write_only>::value;
          if (is_write_only) {
            LA.G_.clear_sync_state();
          }

          if (LA.G_.template need_sync<execution_space>()) {
            LA.G_.template sync<execution_space>();
          }

          constexpr bool is_read_only =
            std::is_same<access_mode, read_only>::value;
          if (! is_read_only) {
            LA.G_.template modify<execution_space>();
          }

          // See note about "copy back" above.
          auto G_lcl_2d =
            LA.G_.template getLocalView<execution_space>();
          // This converts the View to const if applicable.
          return master_local_view_2d_type(G_lcl_2d);
        }
      }

    public:
      // This alias is required by withLocalAccess.
      // using master_local_object_type = std::unique_ptr<
      //   master_local_view_type,
      //   typename impl_type::deleter_type<master_local_view_type>>;
      using master_local_object_type = std::unique_ptr<
        master_local_view_type>;

      // This method is required by withLocalAccess.
      static master_local_object_type
      get(local_access_type LA)
      {
        master_local_view_2d_type G_lcl_2d = getOwningPreActions(LA);
        master_local_view_type G_lcl_1d =
          Kokkos::subview(G_lcl_2d, Kokkos::ALL(), 0);
        // Once we can use C++14, switch to std::make_unique.
        // return master_local_object_type(
        //   new master_local_view_type(G_lcl_1d),
        //   impl_type::getOwningPostActions(LA, G_lcl_2d));

        return master_local_object_type(
          new master_local_view_type(G_lcl_1d));
      }
    };

    static_assert(
      Kokkos::is_view<
        GetMasterLocalObject<
          LocalAccess<Tpetra::Vector<>, read_only>
        >::master_local_view_type
      >::value, "Missing GetMasterLocalObject specialization");

    static_assert(
      Kokkos::is_view<
        GetMasterLocalObject<
          LocalAccess<Tpetra::Vector<>, Kokkos::HostSpace, read_only>
        >::master_local_view_type
      >::value, "Missing GetMasterLocalObject specialization");


    //////////////////////////////////////////////////////////////////

    /// \brief Specialization of GetNonowningLocalObject for
    ///   Tpetra::MultiVector.
    template<class SC, class LO, class GO, class NT, class ... Args>
    struct GetNonowningLocalObject<
      LocalAccess<Tpetra::MultiVector<SC, LO, GO, NT>, Args...>
      >
    {
    private:
      using global_object_type = Tpetra::MultiVector<SC, LO, GO, NT>;
    public:
      using local_access_type =
        LocalAccess<global_object_type, Args...>;

    private:
      using access_mode = typename local_access_type::access_mode;
      static_assert(is_access_mode<access_mode>::value,
        "Please report this bug to the Tpetra developers.");

      using master_local_view_type =
        typename GetMasterLocalObject<local_access_type>::
          master_local_view_type;
      static_assert(
        static_cast<int>(master_local_view_type::Rank) == 2,
        "Rank of master_local_view_type must be 2.");

      // master_local_view_type::non_const_data_type is
      // MV::impl_scalar_type**, where
      // MV = Tpetra::MultiVector<SC, LO, GO, NT>.
      using output_data_type = typename std::conditional<
        std::is_same<access_mode, read_only>::value,
        typename master_local_view_type::const_data_type,
        typename master_local_view_type::non_const_data_type>::type;

    public:
      using master_local_object_type =
        typename GetMasterLocalObject<local_access_type>::
          master_local_object_type;
      using nonowning_local_object_type =
        Kokkos::View<output_data_type,
                     typename master_local_view_type::array_layout,
                     typename master_local_view_type::device_type,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

      static nonowning_local_object_type
      get(local_access_type /* LA */,
          const master_local_object_type& M)
      {
        master_local_view_type* viewPtr = M.get();
        return viewPtr == nullptr ?
          nonowning_local_object_type() :
          nonowning_local_object_type(*viewPtr);
      }
    };

    //////////////////////////////////////////////////////////////////

    /// \brief Specialization of GetNonowningLocalObject for
    ///   Tpetra::Vector.
    template<class SC, class LO, class GO, class NT, class ... Args>
    struct GetNonowningLocalObject<
      LocalAccess<Tpetra::Vector<SC, LO, GO, NT>, Args...>
      >
    {
    private:
      using global_object_type = Tpetra::Vector<SC, LO, GO, NT>;
    public:
      using local_access_type =
        LocalAccess<global_object_type, Args...>;

    private:
      using access_mode = typename local_access_type::access_mode;
      static_assert(is_access_mode<access_mode>::value,
        "Please report this bug to the Tpetra developers.");

      using master_local_view_type =
        typename GetMasterLocalObject<local_access_type>::
          master_local_view_type;
      static_assert(
        static_cast<int>(master_local_view_type::Rank) == 1,
        "Rank of master_local_view_type must be 1.");

      // input_view_type::non_const_data_type is V::impl_scalar_type*,
      // where V = Tpetra::Vector<SC, LO, GO, NT>.
      using output_data_type = typename std::conditional<
        std::is_same<access_mode, read_only>::value,
        typename master_local_view_type::const_data_type,
        typename master_local_view_type::non_const_data_type>::type;

    public:
      using master_local_object_type =
        typename GetMasterLocalObject<local_access_type>::
          master_local_object_type;
      using nonowning_local_object_type =
        Kokkos::View<output_data_type,
                     typename master_local_view_type::array_layout,
                     typename master_local_view_type::device_type,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

      static nonowning_local_object_type
      get (local_access_type /* LA */,
           const master_local_object_type& M)
      {
        master_local_view_type* viewPtr = M.get();
        return viewPtr == nullptr ?
          nonowning_local_object_type() :
          nonowning_local_object_type(*viewPtr);
      }
    };

  } // namespace Details
} // namespace Tpetra

#endif // TPETRA_WITHLOCALACCESS_MULTIVECTOR_HPP
