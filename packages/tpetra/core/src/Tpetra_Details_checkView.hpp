// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_CHECKVIEW_HPP
#define TPETRA_DETAILS_CHECKVIEW_HPP

/// \file Tpetra_Details_checkView.hpp
/// \brief Declaration of functions for checking whether a given
///   pointer is accessible from a given Kokkos execution space.
///
/// \warning This header file and its contents are implementation
///   details of Tpetra.

#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_Details_WrappedDualView.hpp"
#include "Kokkos_DualView.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <sstream>

namespace Tpetra {
namespace Details {

std::string memorySpaceName (const void* ptr);

/// \brief Is the given View valid?
///
/// "Valid" means one of the following
///
/// <ol>
/// <li> The View has zero size; or </li>
/// <li> the View has nonzero size, a nonnull pointer, and the pointer
///      is accessible from ViewType::execution_space (e.g., not a
///      host pointer with execution_space Kokkos::Cuda). </li>
/// </ol>
///
/// \param lclErrStrm [out] If the View is invalid, and this pointer
///   is nonnull, then write a human-readable explanation of what's
///   wrong with the View to the stream.
///
/// \param myMpiProcessRank [in] The rank of the calling MPI process,
///   in whatever communicator is relevant to the caller.  Only used
///   as part of human-readable error output to <tt>*lclErrStrm</tt>.
///
/// \param view [in] The Kokkos::View to investigate.
template<class DataType, class ... Properties>
bool
checkLocalViewValidity
  (std::ostream* lclErrStrm,
   const int myMpiProcessRank,
   const Kokkos::View<DataType, Properties...>& view)
{
  using Teuchos::TypeNameTraits;
  using std::endl;
  using view_type = Kokkos::View<DataType, Properties...>;

  if (view.size () == 0) {
    // Kokkos::View can be zero size with a nonnull pointer.
    // Even std::vector can have this behavior.
    return true;
  }
  else { // nonzero size
    auto ptr = view.data ();

    if (ptr == nullptr) {
      if (lclErrStrm != nullptr) {
        const std::string viewName = TypeNameTraits<view_type>::name ();
        *lclErrStrm << "Proc " << myMpiProcessRank << ": Kokkos::View "
          "of type " << viewName << " has nonzero size " << view.size ()
          << " but a null pointer." << endl;
      }
      return false;
    }
    else 
      return true;
  }
}

/// \brief Is the given Kokkos::DualView valid?
///
/// A DualView is valid if both of its constituent Views are valid.
template<class DataType ,
         class... Args>
bool
checkLocalDualViewValidity
  (std::ostream* const lclErrStrm,
   const int myMpiProcessRank,
   const Kokkos::DualView<DataType, Args...>& dv)
{
  const bool dev_good =
    checkLocalViewValidity (lclErrStrm, myMpiProcessRank,
                            dv.view_device ());
  const bool host_good =
    checkLocalViewValidity (lclErrStrm, myMpiProcessRank,
                            dv.view_host ());
  const bool good = dev_good && host_good;
  if (! good && lclErrStrm != nullptr) {
    using Teuchos::TypeNameTraits;
    using std::endl;
    using dv_type =
      Kokkos::DualView<DataType, Args...>;

    const std::string dvName = TypeNameTraits<dv_type>::name ();
    *lclErrStrm << "Proc " << myMpiProcessRank << ": Kokkos::DualView "
      "of type " << dvName << " has one or more invalid Views.  See "
      "above error messages from this MPI process for details." << endl;
  }
  return good;
}

template<class DataType ,
         class... Args>
bool
checkGlobalDualViewValidity
(std::ostream* const gblErrStrm,
 const Kokkos::DualView<DataType, Args...>& dv,
 const bool verbose,
 const Teuchos::Comm<int>* const comm)
{
  using std::endl;
  const int myRank = comm == nullptr ? 0 : comm->getRank ();
  std::ostringstream lclErrStrm;
  int lclSuccess = 1;

  try {
    const bool lclValid =
      checkLocalDualViewValidity (&lclErrStrm, myRank, dv);
    lclSuccess = lclValid ? 1 : 0;
  }
  catch (std::exception& e) {
    lclErrStrm << "Proc " << myRank << ": checkLocalDualViewValidity "
      "threw an exception: " << e.what () << endl;
    lclSuccess = 0;
  }
  catch (...) {
    lclErrStrm << "Proc " << myRank << ": checkLocalDualViewValidity "
      "threw an exception not a subclass of std::exception." << endl;
    lclSuccess = 0;
  }

  int gblSuccess = 0; // output argument
  if (comm == nullptr) {
    gblSuccess = lclSuccess;
  }
  else {
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  }

  if (gblSuccess != 1 && gblErrStrm != nullptr) {
    *gblErrStrm << "On at least one (MPI) process, the "
      "Kokkos::DualView has "
      "either the device or host pointer in the "
      "DualView equal to null, but the DualView has a nonzero number of "
      "rows.  For more detailed information, please rerun with the "
      "TPETRA_VERBOSE environment variable set to 1. ";
    if (verbose) {
      *gblErrStrm << "  Here are error messages from all "
        "processes:" << endl;
      if (comm == nullptr) {
        *gblErrStrm << lclErrStrm.str ();
      }
      else {
        using Tpetra::Details::gathervPrint;
        gathervPrint (*gblErrStrm, lclErrStrm.str (), *comm);
      }
    }
   *gblErrStrm << endl;
  }
  return gblSuccess == 1;
}


/// \brief Is the given Tpetra::WrappedDualView valid?
///
/// A WrappedDualView is valid if both of its constituent Views are valid.
template<class DataType ,
         class... Args>
bool
checkLocalWrappedDualViewValidity
  (std::ostream* const lclErrStrm,
   const int myMpiProcessRank,
   const Tpetra::Details::WrappedDualView<Kokkos::DualView<DataType, Args...> >& dv)
{
  const bool dev_good  = dv.is_valid_device();
  const bool host_good = dv. is_valid_host();
  const bool good = dev_good && host_good;
  if (! good && lclErrStrm != nullptr) {
    using Teuchos::TypeNameTraits;
    using std::endl;
    using dv_type =
      Tpetra::Details::WrappedDualView<Kokkos::DualView<DataType, Args...> >;

    const std::string dvName = TypeNameTraits<dv_type>::name ();
    *lclErrStrm << "Proc " << myMpiProcessRank << ": Tpetra::WrappedDualView "
      "of type " << dvName << " has one or more invalid Views.  See "
      "above error messages from this MPI process for details." << endl;
  }
  return good;
}

template<class DataType ,
         class... Args>
bool
checkGlobalWrappedDualViewValidity
(std::ostream* const gblErrStrm,
 const Tpetra::Details::WrappedDualView<Kokkos::DualView<DataType, Args...> >& dv,
 const bool verbose,
 const Teuchos::Comm<int>* const comm)
{
  using std::endl;
  const int myRank = comm == nullptr ? 0 : comm->getRank ();
  std::ostringstream lclErrStrm;
  int lclSuccess = 1;

  try {
    const bool lclValid =
      checkLocalWrappedDualViewValidity (&lclErrStrm, myRank, dv);
    lclSuccess = lclValid ? 1 : 0;
  }
  catch (std::exception& e) {
    lclErrStrm << "Proc " << myRank << ": checkLocalDualViewValidity "
      "threw an exception: " << e.what () << endl;
    lclSuccess = 0;
  }
  catch (...) {
    lclErrStrm << "Proc " << myRank << ": checkLocalDualViewValidity "
      "threw an exception not a subclass of std::exception." << endl;
    lclSuccess = 0;
  }

  int gblSuccess = 0; // output argument
  if (comm == nullptr) {
    gblSuccess = lclSuccess;
  }
  else {
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  }

  if (gblSuccess != 1 && gblErrStrm != nullptr) {
    *gblErrStrm << "On at least one (MPI) process, the "
      "Kokkos::DualView has "
      "either the device or host pointer in the "
      "DualView equal to null, but the DualView has a nonzero number of "
      "rows.  For more detailed information, please rerun with the "
      "TPETRA_VERBOSE environment variable set to 1. ";
    if (verbose) {
      *gblErrStrm << "  Here are error messages from all "
        "processes:" << endl;
      if (comm == nullptr) {
        *gblErrStrm << lclErrStrm.str ();
      }
      else {
        using Tpetra::Details::gathervPrint;
        gathervPrint (*gblErrStrm, lclErrStrm.str (), *comm);
      }
    }
   *gblErrStrm << endl;
  }
  return gblSuccess == 1;
}


} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_CHECKVIEW_HPP
