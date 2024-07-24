// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_TESTINGUTILITIES_HPP_
#define TPETRA_TESTINGUTILITIES_HPP_

/// \file Tpetra_TestingUtilities.hpp
/// \brief Internal utilities for testing Tpetra.
///
/// \warning This header file and its contents are implementation
///   details of Tpetra.  Users must not rely on this file existing,
///   or on any contents of this file.

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "KokkosCompat_View.hpp"
#include "Tpetra_Core.hpp"
#include "TpetraCore_ETIHelperMacros.h"
#include "Teuchos_DefaultSerialComm.hpp"
#include "Teuchos_CommHelpers.hpp"

#define TPETRA_GLOBAL_SUCCESS_CHECK(out,comm,success)  \
  { \
    int tgscLclSuccess = success ? 1 : 0; \
    int tgscGblSuccess = 1; \
    Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_MIN, tgscLclSuccess, Teuchos::outArg (tgscGblSuccess)); \
    if (tgscGblSuccess == 1) { \
      out << "Succeeded on all processes!" << std::endl; \
    } else { \
      out << "FAILED on at least one process!" << std::endl; \
    } \
    TEST_EQUALITY_CONST(tgscGblSuccess, 1);  \
    success = (bool) tgscGblSuccess; \
  }


namespace Tpetra {
  namespace TestingUtilities {

    /// \brief Whether to test MPI functionality.
    ///
    /// \warning This variable is an implementation detail of Tpetra.
    ///   Users must not rely on this variable's name or contents, set
    ///   this variable, or rely on behavior resulting from setting
    ///   this variable.
    ///
    /// This variable is true by default.  Its value affects the
    /// behavior of getDefaultComm() (see below in this header file).
    bool testMpi = true;

    /// \brief Get the default communicator for Tpetra tests.
    ///
    /// \warning This function is an implementation detail of Tpetra.
    ///   Users must not call this function or rely on its behavior.
    ///
    /// If testMpi (see above in this header file) false, this
    /// function will return a Teuchos::SerialComm.  Otherwise, this
    /// fucntion will return the default communicator from
    /// Tpetra::getDefaultComm.  If Trilinos was built with MPI
    /// support, the resulting communicator will be a Teuchos::MpiComm
    /// that wraps <tt>MPI_COMM_WORLD</tt>.  Otherwise, it will be
    /// either a Teuchos::SerialComm, or a Teuchos::MpiComm that wraps
    /// <tt>MPI_COMM_SELF</tt>.
    Teuchos::RCP<const Teuchos::Comm<int> >
    getDefaultComm ()
    {
      if (testMpi) {
        return Tpetra::getDefaultComm ();
      }
      else {
        // Always return the same Comm instance.  Create it if it
        // hasn't already been created.  A function-static RCP has an
        // initial value of Teuchos::null.
        static Teuchos::RCP<const Teuchos::SerialComm<int> > serialComm_;
        if (serialComm_.is_null ()) {
          serialComm_ = Teuchos::rcp (new Teuchos::SerialComm<int> ());
        }
        return serialComm_;
      }
    }
    
    // Wrap a subview of Kokkos::View into an ArrayRCP
    template<class ViewType>
    Teuchos::ArrayRCP<typename ViewType::value_type> arcp_from_view(ViewType &view, int size=-1) {
      if(size == -1) size=view.extent(0);
      return Kokkos::Compat::persistingView<ViewType>(view, 0, size);
    }

     // Copy a subview of Kokkos::View into an Teuchos::Array
    template<class T>
    Kokkos::View<T*,Kokkos::LayoutLeft,Kokkos::HostSpace> copy_view_from_array(const Teuchos::Array<T> &v_in, int size=-1) {
      if(size == -1) size=v_in.size();
      Kokkos::View<const T*,Kokkos::LayoutLeft,Kokkos::HostSpace> v_in_wrap(v_in.data(),size);
      Kokkos::View<T*,Kokkos::LayoutLeft,Kokkos::HostSpace> v_out("array_copy",size);
      Kokkos::deep_copy(v_out,v_in_wrap);
      return v_out;
    }

  } // namespace TestingUtilities
} // namespace Tpetra

#endif // TPETRA_TESTINGUTILITIES_HPP_
