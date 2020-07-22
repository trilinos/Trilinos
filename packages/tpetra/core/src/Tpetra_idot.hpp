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

#ifndef TPETRA_DETAILS_IDOT_HPP
#define TPETRA_DETAILS_IDOT_HPP

/// \file Tpetra_idot.hpp
/// \brief Declaration of Tpetra::idot, a nonblocking dot product.
///
/// Tpetra::idot implements a nonblocking dot product.  It takes two
/// Tpetra::Vector or Tpetra::MultiVector inputs, and computes their
/// vector-wise dot product(s).
///
/// The overloads of Tpetra::idot in this file are the only contents
/// of this file upon which users should rely.  Users should not rely
/// on anything else in this file.  In particular, anything in the
/// Tpetra::Details namespace is an implementation detail and is NOT
/// for users.
///
/// If you want to find the only functions in this file that you are
/// supposed to use, search for "SKIP DOWN TO HERE" (omit the quotes).

#include "Tpetra_Details_iallreduce.hpp"
#include "Tpetra_Details_isInterComm.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "KokkosBlas1_dot.hpp"
#include <stdexcept>
#include <sstream>

namespace Tpetra {
namespace Details {

/// \brief Compute dot product locally, with the kernel operating on Views in memory_space.
///   Preconditions:
///    - localResult lives in Device::memory_space.
///    - X's underlying view is up-to-date in Device::memory_space.
template<class SC, class LO, class GO, class NT, class ResultView, class memory_space>
void idotLocal(const ResultView& localResult,
               const ::Tpetra::MultiVector<SC, LO, GO, NT>& X,
               const ::Tpetra::MultiVector<SC, LO, GO, NT>& Y)
{
  using MV = ::Tpetra::MultiVector<SC, LO, GO, NT>;
  using dot_type = typename MV::dot_type;
  using dual_view_type = typename MV::dual_view_type;
  using dev_mem_space = typename dual_view_type::t_dev::memory_space;
  using host_mem_space = typename dual_view_type::t_host::memory_space;
  //Get the other memory space for MV's DualView, whatever that is (it may be the same in case of UVM)
  using other_memory_space = typename std::conditional<std::is_same<memory_space, dev_mem_space>::value,
        host_mem_space, dev_mem_space>::type;
  using pair_type = Kokkos::pair<size_t, size_t>;
  const size_t numRows = X.getLocalLength ();
  const pair_type rowRange (0, numRows);
  const size_t X_numVecs = X.getNumVectors ();
  const size_t Y_numVecs = Y.getNumVectors ();
  const size_t numVecs = X_numVecs > Y_numVecs ? X_numVecs : Y_numVecs;
  // Check compatibility of number of columns; allow special cases of
  // a multivector dot a single vector, or vice versa.
  if (X_numVecs != Y_numVecs &&
      X_numVecs != size_t (1) &&
      Y_numVecs != size_t (1)) {
    std::ostringstream os;
    os << "Tpetra::idot: X.getNumVectors() = " << X_numVecs
       << " != Y.getNumVectors() = " << Y_numVecs
       << ", but neither is 1.";
    throw std::invalid_argument (os.str ());
  }
  //Can the multivector dot kernel be used?
  bool useMVDot = X.isConstantStride() && Y.isConstantStride() && X_numVecs == Y_numVecs;
  //Does Y need to be copied to the same memory_space as X and localResult?
  bool copyY = Y.template need_sync<memory_space>();
  if(useMVDot)
  {
    //Get local X (we know it doesn't need a sync to memory_space)
    auto X_lcl = X.template getLocalView<memory_space>();
    //Get local Y, and make a copy if it needs sync. Y is const so we can't just sync to memory_space.
    decltype (X_lcl) Y_lcl;
    if (copyY) {
      auto Y_lcl_other = Y.template getLocalView<other_memory_space>();
      typename decltype (Y_lcl)::non_const_type
        Y_lcl_nc (Kokkos::ViewAllocateWithoutInitializing("Y_lcl"), Y_lcl_other.extent(0), Y_numVecs);
      Kokkos::deep_copy (Y_lcl_nc, Y_lcl_other);
      Y_lcl = Y_lcl_nc;
    }
    else { // Y already sync'd to dev; just use its dev view
      Y_lcl = Y.template getLocalView<memory_space>();
    }
    if (numVecs == size_t (1)) {
      auto X_lcl_1d = Kokkos::subview (X_lcl, rowRange, 0);
      auto Y_lcl_1d = Kokkos::subview (Y_lcl, rowRange, 0);
      auto result_0d = Kokkos::subview (localResult, 0);
      KokkosBlas::dot (result_0d, X_lcl_1d, Y_lcl_1d);
    }
    else {
      auto X_lcl_2d = Kokkos::subview (X_lcl, rowRange, pair_type (0, X_numVecs));
      auto Y_lcl_2d = Kokkos::subview (Y_lcl, rowRange, pair_type (0, Y_numVecs));
      KokkosBlas::dot (localResult, X_lcl_2d, Y_lcl_2d);
    }
  }
  else
  {
    //Need to compute each dot individually, using 1D subviews of X_lcl and Y_lcl
    Kokkos::View<dot_type*, memory_space> YLocalCol_nc;
    if(copyY) {
      //Need to allocate Ycol, as a place to store each column of Y temporarily
      YLocalCol_nc = Kokkos::View<dot_type*, memory_space>(
          Kokkos::ViewAllocateWithoutInitializing("YLocalCol_nc"), numRows);
    }
    for(size_t vec = 0; vec < numVecs; vec++) {
      //Get the Vectors for each column of X and Y that will be dotted together.
      size_t Xvec = (numVecs == X_numVecs) ? vec : size_t(0);
      size_t Yvec = (numVecs == Y_numVecs) ? vec : size_t(0);
      auto Xcol = X.getVector(Xvec);
      auto Ycol = Y.getVector(Yvec);
      auto XLocalCol = Kokkos::subview(Xcol->template getLocalView<memory_space>(), rowRange, 0);
      decltype(XLocalCol) YLocalCol;
      if(copyY) {
        //Get a copy of Y's column if needed. If X has k columns but Y has 1, only copy to YLocalCol_nc once.
        if(vec == 0 || Y_numVecs > 1) {
          auto YLocalCol_other_space = Kokkos::subview(Ycol->template getLocalView<other_memory_space>(), rowRange, 0);
          Kokkos::deep_copy(YLocalCol_nc, YLocalCol_other_space);
        }
        YLocalCol = YLocalCol_nc;
      }
      else {
        YLocalCol = Kokkos::subview(Ycol->template getLocalView<memory_space>(), rowRange, 0);
      }
      //Compute the rank-1 dot product, and place the result in an element of localResult
      KokkosBlas::dot(Kokkos::subview(localResult, vec), XLocalCol, YLocalCol);
    }
  }
}

//Fallback version: produces same result, but is synchronous.
//Used only when MPI is not CUDA-aware, and the global result lives on device.
//It's not allowed to use a CudaSpace or CudaUVMSpace view as a target buffer for MPI
//unless CUDA-aware.
template<class SC, class LO, class GO, class NT, class ResultView>
void blockingDotImpl(
    const ResultView& globalResult,
    const ::Tpetra::MultiVector<SC, LO, GO, NT>& X,
    const ::Tpetra::MultiVector<SC, LO, GO, NT>& Y,
    bool mpiReduceInPlace)
{
  using MV = ::Tpetra::MultiVector<SC, LO, GO, NT>;
  using dot_type = typename MV::dot_type;
  using result_dev_view_type = Kokkos::View<dot_type*, typename NT::device_type>;
  using result_mirror_view_type = typename result_dev_view_type::HostMirror;
  using result_host_view_type = Kokkos::View<dot_type*, Kokkos::HostSpace>;
  using dev_mem_space = typename result_dev_view_type::memory_space;
  using mirror_mem_space = typename result_mirror_view_type::memory_space;
  using unmanaged_result_dev_view_type = Kokkos::View<dot_type*, dev_mem_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using unmanaged_result_host_view_type = Kokkos::View<dot_type*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  const size_t numVecs = globalResult.extent(0);
  //Logic to compute the local dot is no different than when CUDA-aware MPI.
  result_host_view_type localHostResult(Kokkos::ViewAllocateWithoutInitializing("HostLocalDotResult"), numVecs);
  if(X.template need_sync<dev_mem_space>())
  {
    //compute local result on host.
    //NOTE: UVM never takes this branch because then dev_mem_space == mirror_mem_space.
    //This means that the dot is actually computed on host.
    idotLocal<SC, LO, GO, NT, result_host_view_type, mirror_mem_space>(localHostResult, X, Y);
  }
  else
  {
    //compute local result on temporary device view, then copy that to host.
    result_dev_view_type localDeviceResult(Kokkos::ViewAllocateWithoutInitializing("DeviceLocalDotResult"), numVecs);
    idotLocal<SC, LO, GO, NT, result_dev_view_type, mirror_mem_space>(localDeviceResult, X, Y);
    //NOTE: no fence is required: deep_copy will fence.
    Kokkos::deep_copy(localHostResult, localDeviceResult);
  }
  result_host_view_type globalHostResultOwning;
  unmanaged_result_host_view_type globalHostResultNonowning;
  if(mpiReduceInPlace)
  {
    globalHostResultNonowning = unmanaged_result_host_view_type(localHostResult.data(), numVecs);
  }
  else
  {
    globalHostResultOwning = result_host_view_type(Kokkos::ViewAllocateWithoutInitializing("HostGlobalDotResult"), numVecs);
    globalHostResultNonowning = unmanaged_result_host_view_type(globalHostResultOwning.data(), numVecs);
  }
  Teuchos::reduceAll<int, dot_type> (*X.getMap()->getComm(), Teuchos::REDUCE_SUM, numVecs, localHostResult.data(), globalHostResultNonowning.data());
  Kokkos::deep_copy(globalResult, globalHostResultNonowning);
}

/// \brief Internal (common) version of idot, a global dot product
/// that uses a non-blocking MPI reduction.
///
/// Will always try to minimize the amount of host-device copies.
/// However, if the globalResult is in a CUDA device space (CudaSpace
/// or CudaUVMSpace) and MPI is not CUDA-aware, a synchronous reduction
/// will be used instead to ensure correct execution.
template<class SC, class LO, class GO, class NT, class ResultView>
std::shared_ptr< ::Tpetra::Details::CommRequest>
idotImpl(const ResultView& globalResult,
         const ::Tpetra::MultiVector<SC, LO, GO, NT>& X,
         const ::Tpetra::MultiVector<SC, LO, GO, NT>& Y)
{
  using pair_type = ::Kokkos::pair<size_t, size_t>;
  using MV = ::Tpetra::MultiVector<SC, LO, GO, NT>;
  using dot_type = typename MV::dot_type;
  using global_result_memspace = typename ResultView::memory_space;
  using result_dev_view_type = Kokkos::View<dot_type*, typename NT::device_type>;
  using result_mirror_view_type = typename result_dev_view_type::HostMirror;
  using result_host_view_type = Kokkos::View<dot_type*, Kokkos::HostSpace>;
  using dev_mem_space = typename result_dev_view_type::memory_space;
  using mirror_mem_space = typename result_mirror_view_type::memory_space;
  using unmanaged_result_dev_view_type = Kokkos::View<dot_type*, dev_mem_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using unmanaged_result_mirror_view_type = Kokkos::View<dot_type*, mirror_mem_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using unmanaged_result_host_view_type = Kokkos::View<dot_type*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  auto comm = X.getMap()->getComm();
  bool mpiReduceInPlace = !Details::isInterComm(*comm);
  //note: numVecs is verified in idotLocal
  const size_t numVecs = globalResult.extent(0);
  //BMK: If either local or global results live in a CUDA device space but
  //MPI is not CUDA-aware, we are forced to perform (synchronous) communication
  //using temporary HostSpace view(s).
#ifdef KOKKOS_ENABLE_CUDA
  //If globalResult's memory space is a Cuda space, then the user's global result lives on device.
  //For this case, there is currently no way to do the asynchronous MPI reduction to a temporary HostSpace
  //view, and then deep copy to the original globalResult. 
  //However, this could also be solved
  //(and still be globally asynchronous) by adding a
  //deep_copy callback to CommRequest.
  if(!Tpetra::Details::Behavior::assumeMpiIsCudaAware() &&
     (std::is_same<global_result_memspace, Kokkos::CudaSpace>::value || 
     std::is_same<global_result_memspace, Kokkos::CudaUVMSpace>::value))
  {
    blockingDotImpl<SC, LO, GO, NT, ResultView>(globalResult, X, Y, mpiReduceInPlace);
    return Tpetra::Details::Impl::emptyCommRequest();
  }
#endif
  bool X_latestOnHost = X.template need_sync<dev_mem_space>();
  bool runOnHost = !std::is_same<dev_mem_space, mirror_mem_space>::value && X_latestOnHost;
  //Wherever the kernel runs, allocate a separate local result if either:
  //  * comm is an "InterComm" (can't do in-place collectives)
  //  * globalResult is not accessible from the space of X's latest data
  const bool allocateLocalResult = !mpiReduceInPlace ||
    (X_latestOnHost ?
      !Kokkos::Impl::MemorySpaceAccess<mirror_mem_space, global_result_memspace>::accessible :
      !Kokkos::Impl::MemorySpaceAccess<dev_mem_space, global_result_memspace>::accessible);
  if(runOnHost) {
    //Note BMK: cannot get here if the device memory space is CudaUVMSpace,
    //because in that case the mirror space of MV's DualView is also CudaUVMSpace. Verified this by printing typeid.
    unmanaged_result_host_view_type nonowningLocalResult;
    result_host_view_type localResult;
    if(allocateLocalResult) {
      localResult = result_host_view_type(Kokkos::ViewAllocateWithoutInitializing("localResult"), numVecs);
      nonowningLocalResult = unmanaged_result_host_view_type(localResult.data(), numVecs);
    }
    else
      nonowningLocalResult = unmanaged_result_host_view_type(globalResult.data(), numVecs);
    idotLocal<SC, LO, GO, NT, unmanaged_result_host_view_type, mirror_mem_space>(nonowningLocalResult, X, Y);
    return iallreduce(nonowningLocalResult, globalResult, ::Teuchos::REDUCE_SUM, *comm);
  }
  else {
    //running on device
    unmanaged_result_dev_view_type nonowningLocalResult;
    result_dev_view_type localResult;
    if(allocateLocalResult) {
      localResult = result_dev_view_type(Kokkos::ViewAllocateWithoutInitializing("localResult"), numVecs);
      nonowningLocalResult = unmanaged_result_dev_view_type(localResult.data(), numVecs);
    }
    else
      nonowningLocalResult = unmanaged_result_dev_view_type(globalResult.data(), numVecs);
    idotLocal<SC, LO, GO, NT, unmanaged_result_dev_view_type, dev_mem_space>(nonowningLocalResult, X, Y);
    //If MPI not CUDA-aware, must copy local result to HostSpace before iallreduce.
    //That view will persistent in the CommRequest, and localResult won't.
    bool communicateFromHost = false;
#ifdef KOKKOS_ENABLE_CUDA
    communicateFromHost = !Tpetra::Details::Behavior::assumeMpiIsCudaAware();
#endif
    if(communicateFromHost)
    {
      result_host_view_type hostLocalResult(Kokkos::ViewAllocateWithoutInitializing("hostLocalResult"), numVecs);
      //The deep copy fences.
      Kokkos::deep_copy(hostLocalResult, nonowningLocalResult);
      return iallreduce(hostLocalResult, globalResult, ::Teuchos::REDUCE_SUM, *comm);
    }
    else
    {
      //(Only fence in idot) required because the device-space result of idotLocal will be accessed directly by MPI.
      typename dev_mem_space::execution_space().fence();
      return iallreduce(nonowningLocalResult, globalResult, ::Teuchos::REDUCE_SUM, *comm);
    }
  }
}
} // namespace Details

//
// SKIP DOWN TO HERE
//

/// \brief Nonblocking dot product, with either Tpetra::MultiVector or
///   Tpetra::Vector inputs, and raw pointer or raw array output.
///
/// This implements the following cases:
/// <ul>
/// <li> If X and Y have the same number of columns (zero or more),
///      compute the dot products of corresponding columns of X and
///      Y.  That is, resultRaw[j] = dot(X(:,j), Y(:,j)). </li>
/// <li> If X has one column and Y has more than one column, compute
///      the dot products of X and each column of Y in turn.
///      That is, resultRaw[j] = dot(X(:,0), Y(:,j)). </li>
/// <li> If X has more than one column and Y has one column, compute
///      the dot products of each column of X in turn with X.
///      That is, resultRaw[j] = dot(X(:,j), Y(:,0)). </li>
/// </ul>
///
/// \param resultRaw [out] Output; raw pointer or raw array to the
///   return value(s).  It is only valid to read this after calling
///   wait() on the return value.  It must be legal to write to the
///   first <tt>std::max(X.getNumVectors(), Y.getNumVectors())</tt>
///   entries of this array.  This memory must be accessible from the
///   host CPU.
///
/// \param X [in] First input Tpetra::MultiVector or Tpetra::Vector.
///   This must have same number of rows (globally, and on each (MPI)
///   process) as Y.  If this is a Tpetra::MultiVector, then this must
///   either have the same number of columns as Y, or have one column.
///
/// \param Y [in] Second input Tpetra::MultiVector or Tpetra::Vector.
///   This must have same number of rows (globally, and on each (MPI)
///   process) as X.  If this is a Tpetra::MultiVector, then this must
///   either have the same number of columns as X, or have one column.
///
/// \return Pointer to an object representing the nonblocking
///   collective (communication operation).  Call wait() on this
///   object to complete the collective.  After calling wait(), you
///   may read the result.
///
/// \tparam SC Same as the first template parameter of Tpetra::Vector.
/// \tparam LO Same as the second template parameter of Tpetra::Vector.
/// \tparam GO Same as the third template parameter of Tpetra::Vector.
/// \tparam NT Same as the fourth template parameter of Tpetra::Vector.
///
/// In this version of the function, the dot product result goes into
/// an array (just a raw pointer).  The \c dot_type typedef is the
/// type of a dot product result, for a Tpetra::Vector whose entries
/// have type \c SC (the "Scalar" type).  For most \c SC types,
/// <tt>dot_type == SC</tt>.  However, once you start dipping into
/// more interesting Scalar types, such as those found in the Sacado
/// or Stokhos packages, \c dot_type may be a different type.  Most
/// users should not have to worry about this, but Tpetra developers
/// may need to worry about this.
template<class SC, class LO, class GO, class NT>
std::shared_ptr< ::Tpetra::Details::CommRequest>
idot (typename ::Tpetra::MultiVector<SC, LO, GO, NT>::dot_type* resultRaw,
      const ::Tpetra::MultiVector<SC, LO, GO, NT>& X,
      const ::Tpetra::MultiVector<SC, LO, GO, NT>& Y)
{
  using dot_type = typename ::Tpetra::Vector<SC, LO, GO, NT>::dot_type;
  const size_t X_numVecs = X.getNumVectors ();
  const size_t Y_numVecs = Y.getNumVectors ();
  const size_t numVecs = (X_numVecs > Y_numVecs) ? X_numVecs : Y_numVecs;
  Kokkos::View<dot_type*, Kokkos::HostSpace> resultView(resultRaw, numVecs);
  return Details::idotImpl<SC,LO,GO,NT>(resultView, X, Y);
}

/// \brief Nonblocking dot product, with Tpetra::MultiVector inputs,
///   and rank-1 (one-dimensional array) Kokkos::View output.
///
/// This implements the following cases:
/// <ul>
/// <li> If X and Y have the same number of columns (zero or more),
///      compute the dot products of corresponding columns of X and
///      Y.  That is, result[j] = dot(X(:,j), Y(:,j)). </li>
/// <li> If X has one column and Y has more than one column, compute
///      the dot products of X and each column of Y in turn.
///      That is, result[j] = dot(X(:,0), Y(:,j)). </li>
/// <li> If X has more than one column and Y has one column, compute
///      the dot products of each column of X in turn with X.
///      That is, result[j] = dot(X(:,j), Y(:,0)). </li>
/// </ul>
///
/// \param result [out] Output; rank-1 Kokkos::View.  It is only valid
///   to read the entries of this after calling wait() on the return
///   value.  This must have length
///   <tt>std::max(X.getNumVectors(), Y.getNumVectors())</tt>.
///
/// \param X [in] First input Tpetra::MultiVector or Tpetra::Vector.
///   This must have same number of rows (globally, and on each (MPI)
///   process) as Y.  If this is a Tpetra::MultiVector, then this must
///   either have the same number of columns as Y, or have one column.
///
/// \param Y [in] Second input Tpetra::MultiVector or Tpetra::Vector.
///   This must have same number of rows (globally, and on each (MPI)
///   process) as X.  If this is a Tpetra::MultiVector, then this must
///   either have the same number of columns as X, or have one column.
///
/// \return Pointer to an object representing the nonblocking
///   collective (communication operation).  Call wait() on this
///   object to complete the collective.  After calling wait(), you
///   may read the results.
///
/// \tparam SC Same as the first template parameter of Tpetra::MultiVector.
/// \tparam LO Same as the second template parameter of Tpetra::MultiVector.
/// \tparam GO Same as the third template parameter of Tpetra::MultiVector.
/// \tparam NT Same as the fourth template parameter of Tpetra::MultiVector.
///
/// Compute the dot product of each column of X, with each
/// corresponding column of Y.  If X has a single column, then compute
/// the dot product of X with each column of Y in turn.  If Y has a
/// single column, then compute the dot product of each column of X in
/// turn with Y.
///
/// In this version of the function, the dot product results go into a
/// rank-1 (one-dimensional) Kokkos::View.  We prefer that you use the
/// versions of idot() that take a Kokkos::View as their output
/// argument, because this ensures that the output will still exist
/// (not be deallocated or fall out of scope).  The versions of idot()
/// that take a raw pointer cannot promise that the memory will
/// continue to exist until the dot product is done.
///
/// The \c dot_type typedef is the type of a dot product result, for a
/// Tpetra::Vector whose entries have type \c SC (the "Scalar" type).
/// For most \c SC types, <tt>dot_type == SC</tt>.  However, once you
/// start dipping into more interesting Scalar types, such as those
/// found in the Sacado or Stokhos packages, \c dot_type may be a
/// different type.  Most users should not have to worry about this,
/// but Tpetra developers may need to worry about this.
template<class SC, class LO, class GO, class NT>
std::shared_ptr< ::Tpetra::Details::CommRequest>
idot (const Kokkos::View<typename ::Tpetra::MultiVector<SC, LO, GO, NT>::dot_type*,
        typename ::Tpetra::MultiVector<SC, LO, GO, NT>::device_type>& result,
      const ::Tpetra::MultiVector<SC, LO, GO, NT>& X,
      const ::Tpetra::MultiVector<SC, LO, GO, NT>& Y)
{
  return Details::idotImpl<SC,LO,GO,NT>(result, X, Y);
}

/// \brief Nonblocking dot product, with Tpetra::Vector inputs, and
///   rank-0 (single value) Kokkos::View output.
///
/// This function computes result() = dot(X,Y).
///
/// \param result [out] Output; rank-0 Kokkos::View of the return
///   value.  It is only valid to read this after calling wait() on
///   the return value.
///
/// \param X [in] First input Tpetra::Vector.  This must have same
///   number of rows (globally, and on each (MPI) process) as Y.
///
/// \param Y [in] Second input Tpetra::Vector.  This must have same
///   number of rows (globally, and on each (MPI) process) as X.
///
/// \return Pointer to an object representing the nonblocking
///   collective (communication operation).  Call wait() on this
///   object to complete the collective.  After calling wait(), you
///   may read the result.
///
/// \tparam SC Same as the first template parameter of Tpetra::Vector.
/// \tparam LO Same as the second template parameter of Tpetra::Vector.
/// \tparam GO Same as the third template parameter of Tpetra::Vector.
/// \tparam NT Same as the fourth template parameter of Tpetra::Vector.
///
/// In this version of the function, the dot product result goes into
/// a rank-0 ("zero-dimensional") Kokkos::View.  Rank-0 Views just
/// view a single value.  We prefer that you use the versions of
/// idot() that take a Kokkos::View as their output argument, because
/// this ensures that the output will still exist (not be deallocated
/// or fall out of scope).  The versions of idot() that take a raw
/// pointer cannot promise that the memory will continue to exist
/// until the dot product is done.
///
/// The \c dot_type typedef is the type of a dot product result, for a
/// Tpetra::Vector whose entries have type \c SC (the "Scalar" type).
/// For most \c SC types, <tt>dot_type == SC</tt>.  However, once you
/// start dipping into more interesting Scalar types, such as those
/// found in the Sacado or Stokhos packages, \c dot_type may be a
/// different type.  Most users should not have to worry about this,
/// but Tpetra developers may need to worry about this.
template<class SC, class LO, class GO, class NT>
std::shared_ptr< ::Tpetra::Details::CommRequest>
idot (const Kokkos::View<typename ::Tpetra::Vector<SC, LO, GO, NT>::dot_type,
        typename ::Tpetra::Vector<SC, LO, GO, NT>::device_type>& result,
      const ::Tpetra::Vector<SC, LO, GO, NT>& X,
      const ::Tpetra::Vector<SC, LO, GO, NT>& Y)
{
  using dot_type = typename ::Tpetra::Vector<SC, LO, GO, NT>::dot_type;
  using result_device_t = typename ::Tpetra::Vector<SC, LO, GO, NT>::device_type;
  Kokkos::View<dot_type*, result_device_t, Kokkos::MemoryTraits<Kokkos::Unmanaged>> result1D(result.data(), 1);
  return Details::idotImpl<SC,LO,GO,NT>(result1D, X, Y);
}

} // namespace Tpetra

#endif // TPETRA_DETAILS_IDOT_HPP
