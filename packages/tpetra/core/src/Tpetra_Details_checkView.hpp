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

#ifndef TPETRA_DETAILS_CHECKVIEW_HPP
#define TPETRA_DETAILS_CHECKVIEW_HPP

/// \file Tpetra_Details_checkPointer.hpp
/// \brief Declaration of functions for checking whether a given
///   pointer is accessible from a given Kokkos execution space.
///
/// \warning This header file and its contents are implementation
///   details of Tpetra.

#include "Tpetra_Details_checkPointer.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Kokkos_DualView.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <sstream>

namespace Tpetra {
namespace Details {

template<class ExecutionSpace>
bool
pointerAccessibleFromExecutionSpace (const void* ptr,
                                     const ExecutionSpace& space);

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
/// This function doesn't promise exact results for anything other
/// than CudaSpace, CudaUVMSpace, or CudaHostPinnedSpace.  The point
/// of this function is for Tpetra classes to debug user error in
/// which users create a Kokkos::View with a raw pointer in the wrong
/// memory space (e.g., a host pointer, when they should have used a
/// UVM pointer).
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
  using ES = typename view_type::execution_space;

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
    else { // nonnull pointer, nonzero size
      const bool canAcc = pointerAccessibleFromExecutionSpace (ptr, ES ());
      if (! canAcc && lclErrStrm != nullptr) {
        const std::string viewName = TypeNameTraits<view_type>::name ();
        const std::string esName = TypeNameTraits<ES>::name ();
        *lclErrStrm << "Proc " << myMpiProcessRank << ": Kokkos::View "
          "of type " << viewName << " and nonzero size " << view.size ()
          << " has a pointer " << ptr << " which is nonnull, but not "
          "accessible from the View's claimed execution space "
          << esName << ".  As far as I can tell, the View's pointer "
          "(view.data()) lives in memory space " << memorySpaceName (ptr)
          << "." << endl;
      }
      return canAcc;
    }
  }
}

/// \brief Is the given Kokkos::DualView valid?
///
/// A DualView is valid if both of its constituent Views are valid.
template<class DataType ,
         class Arg1Type = void ,
         class Arg2Type = void ,
         class Arg3Type = void>
bool
checkLocalDualViewValidity
  (std::ostream* const lclErrStrm,
   const int myMpiProcessRank,
   const Kokkos::DualView<DataType, Arg1Type, Arg2Type, Arg3Type>& dv)
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
      Kokkos::DualView<DataType, Arg1Type, Arg2Type, Arg3Type>;

    const std::string dvName = TypeNameTraits<dv_type>::name ();
    *lclErrStrm << "Proc " << myMpiProcessRank << ": Kokkos::DualView "
      "of type " << dvName << " has one or more invalid Views.  See "
      "above error messages from this MPI process for details." << endl;
  }
  return good;
}

template<class DataType ,
         class Arg1Type = void ,
         class Arg2Type = void ,
         class Arg3Type = void>
bool
checkGlobalDualViewValidity
(std::ostream* const gblErrStrm,
 const Kokkos::DualView<DataType, Arg1Type, Arg2Type, Arg3Type>& dv,
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
      "Kokkos::DualView has either device or host pointer that "
      "is not accessible from the corresponding execution space.  "
      "This can happen if you, the user, created the DualView "
      "with raw pointer(s) that is/are not in the correct memory "
      "space.  For example, you may have a Kokkos::DualView whose "
      "device memory_space is Kokkos::CudaUVMSpace, but you gave "
      "the device View's constructor a raw host pointer.  It may "
      "also happen if either the device or host pointer in the "
      "DualView is null, but the DualView has a nonzero number of "
      "rows.  For more detailed information, please rerun with the "
      "TPETRA_VERBOSE environment variable set to 1.  (You do not "
      "need to recompile.)";
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
