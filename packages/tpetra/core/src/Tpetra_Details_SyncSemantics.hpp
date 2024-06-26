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
#ifndef TPETRA_DETAILS_SYNC_SEMANTICS_HPP
#define TPETRA_DETAILS_SYNC_SEMANTICS_HPP
/// \file Tpetra_Details_SyncSemantics.hpp
/// \brief Declaration of various tools controlling the device sync semantics
/// of Tpetra

#include <string>
#include "Kokkos_Core.hpp"

namespace Tpetra {

  namespace Details {
    /// \brief  Relax the Tpetra syncronization semantics.
    ///
    /// This will trigger a Kokkos::fence() immediately on call
    /// WARNING: This is not intended for use outside of Trilinos.
    void enableRelaxedSyncs();

    /// \brief  Use the default Tpetra syncronization semantics.
    ///
    /// This will trigger a Kokkos::fence() immediately on call
    /// WARNING: This is not intended for use outside of Trilinos.
    void disableRelaxedSyncs();

    /// \brief  Query the state of Tpetra synchronization semantics
    ///
    /// WARNING: This is not intended for use outside of Trilinos.
    bool areRelaxedSyncsEnabled();

  };


  // Tpetra wrapper of Kokkos::fence();
  void fence();
  void fence(const std::string& label);

  // Tpetra wrappers ofKokkos::deep_copy();
  template<class ExecSpace, class ViewDest, class ViewSrc>
  void deep_copy(const ExecSpace &exec_space, const ViewDest &dest, const ViewSrc &src) {
    Kokkos::deep_copy(exec_space,dest,src);
  }

  template<class ExecSpace, class ViewDest>
  void deep_copy(const ExecSpace &exec_space, const ViewDest &dest, const typename ViewDest::value_type &src) {
    Kokkos::deep_copy(exec_space,dest,src);
  }

  template<class ExecSpace, class ViewSrc>
  void deep_copy(const ExecSpace &exec_space, typename ViewSrc::value_type &dest, const ViewSrc &src) {
    Kokkos::deep_copy(exec_space,dest,src);
  }

  template<class ViewDest, class ViewSrc>
  void deep_copy(const ViewDest &dest, const ViewSrc &src) {
    if(Details::areRelaxedSyncsEnabled()) {
      auto exec_space = Kokkos::DefaultExecutionSpace();
      Kokkos::deep_copy(exec_space,dest,src);
      exec_space.fence();
    }
    else {
      Kokkos::deep_copy(dest,src);
    }
  }

  template<class ViewDest>
  void deep_copy(const ViewDest &dest, const typename ViewDest::value_type &src) {
    if(Details::areRelaxedSyncsEnabled()) {
      auto exec_space = Kokkos::DefaultExecutionSpace();
      Kokkos::deep_copy(exec_space,dest,src);
      exec_space.fence();
    }
    else {
      Kokkos::deep_copy(dest,src);
    }
  }

  template<class ViewSrc>
  void deep_copy(typename ViewSrc::value_type &dest, const ViewSrc &src) {
    if(Details::areRelaxedSyncsEnabled()) {
      auto exec_space = Kokkos::DefaultExecutionSpace();
      Kokkos::deep_copy(exec_space,dest,src);
      exec_space.fence();
    }
    else {
      Kokkos::deep_copy(dest,src);
    }
  }



} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_SYNC_SEMANTICS_HPP
