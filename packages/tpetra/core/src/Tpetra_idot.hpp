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
#include "KokkosBlas1_dot.hpp"
#include <stdexcept>
#include <sstream>

namespace Tpetra {
namespace Details {

/// \brief Compute the local dot product(s), columnwise, of X and Y.
///
/// \warning This is an implementation detail of Tpetra::idot.
///   Users should never call this function.
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
/// \tparam SC Same as the first template parameter of Tpetra::MultiVector.
/// \tparam LO Same as the second template parameter of Tpetra::MultiVector.
/// \tparam GO Same as the third template parameter of Tpetra::MultiVector.
/// \tparam NT Same as the fourth template parameter of Tpetra::MultiVector.
///
/// \param resultRaw [out] Raw pointer to output array of dot products.
/// \param X [in] First input MultiVector.
/// \param Y [in] Second input MultiVector.
/// \param resultOnDevice [in] Whether \c resultRaw points to memory
///   accessible from device.  If not, it may only be accessed from a
///   host execution space.
template<class SC, class LO, class GO, class NT>
void
lclDotRaw (typename ::Tpetra::MultiVector<SC, LO, GO, NT>::dot_type* const resultRaw,
           const ::Tpetra::MultiVector<SC, LO, GO, NT>& X,
           const ::Tpetra::MultiVector<SC, LO, GO, NT>& Y,
           const bool resultOnDevice)
{
  using ::Kokkos::Impl::MemorySpaceAccess;
  using ::Kokkos::deep_copy;
  using ::Kokkos::HostSpace;
  using ::Kokkos::subview;
  using ::Kokkos::View;
  typedef ::Kokkos::pair<size_t, size_t> pair_type;
  typedef ::Tpetra::MultiVector<SC, LO, GO, NT> MV;
  typedef typename MV::dot_type dot_type;
  typedef typename MV::device_type device_type;
  typedef typename MV::dual_view_type::t_dev::memory_space dev_memory_space;
  typedef typename MV::dual_view_type::t_host::memory_space host_memory_space;
  typedef View<dot_type*, device_type> result_dev_view_type;
  typedef typename result_dev_view_type::HostMirror result_host_view_type;

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

  const bool X_latestOnHost = X.template need_sync<dev_memory_space> () &&
    ! X.template need_sync<host_memory_space> ();
  const bool Y_latestOnHost = Y.template need_sync<dev_memory_space> () &&
    ! Y.template need_sync<host_memory_space> ();
  // Let X guide whether to execute on device or host.
  const bool runOnHost = X_latestOnHost;
  // Assume that we have to copy Y to where we need to run.  idot()
  // takes the input MultiVectors as const, so it can't sync them.  We
  // just have to copy their data, if needed, to the memory space
  // where we want to run.
  const bool Y_copyToHost = (X_latestOnHost && ! Y_latestOnHost);
  const bool Y_copyToDev = (! X_latestOnHost && Y_latestOnHost);

  result_dev_view_type result;
  result_host_view_type result_h;
  if (resultOnDevice) {
    result = result_dev_view_type (resultRaw, numVecs);
    if (runOnHost) {
      // create_mirror_view doesn't promise to create a deep copy, so we
      // can't rely on this case to avoid an extra buffer if the input
      // communicator is an intercomm.
      result_h = ::Kokkos::create_mirror_view (result);
    }
    // If running on device, not host, then we don't use result_h.
  }
  else {
    result_h = result_host_view_type (resultRaw, numVecs);
    if (! runOnHost) { // running on device, not host
      // We have to do a separate allocation here, since resultRaw is
      // not accessible from device.
      result = result_dev_view_type ("result", numVecs);
    }
    // If running on host, not device, then we don't use result.
  }

  const bool X_or_Y_have_nonconstant_stride =
    ! X.isConstantStride () || ! Y.isConstantStride ();
  if (X_or_Y_have_nonconstant_stride && numVecs > size_t (1)) {
    // The "nonconstant stride" case relates to the ability of a
    // Tpetra::MultiVector to view multiple, noncontiguous columns of
    // another Tpetra::MultiVector.  In this case, we can't just get
    // the Kokkos::View out of X and Y, since the Kokkos::View has all
    // the columns, not just the columns that X or Y view.  The
    // Kokkos::View of X or Y would view a contiguous subset of
    // columns of the latter, "parent" MultiVector, not of the former,
    // "child" MultiVector.  (Kokkos::View doesn't have a layout to
    // cover this case where consecutive entries within a column are
    // contiguous, but consecutive entries within a row may have
    // varying strides.)  Instead, we work one column at a time.
    for (size_t j = 0; j < numVecs; ++j) {
      // numVecs = max(X_numVecs, Y_numVecs).  Allow special case of
      // X_numVecs != Y_numVecs && (X_numVecs == 1 || Y_numVecs == 1).
      const size_t X_col = (X_numVecs == size_t (1)) ? size_t (0) : j;
      const size_t Y_col = (Y_numVecs == size_t (1)) ? size_t (0) : j;
      auto X_j = X.getVector (X_col);
      auto Y_j = Y.getVector (Y_col);

      if (runOnHost) {
        auto X_j_2d_h = X_j->template getLocalView<host_memory_space> ();
        auto X_j_1d_h = subview (X_j_2d_h, rowRange, 0);
        decltype (X_j_1d_h) Y_j_1d_h;

        if (Y_copyToHost) {
          auto Y_j_2d = Y_j->template getLocalView<dev_memory_space> ();
          auto Y_j_1d = subview (Y_j_2d, rowRange, 0);
          typename decltype (Y_j_1d_h)::non_const_type
            Y_j_1d_h_nc ("Y_j_1d_h", Y_j_1d.extent (0));
          deep_copy (Y_j_1d_h_nc, Y_j_1d);
          Y_j_1d_h = Y_j_1d_h_nc;
        }
        else { // Y already sync'd to host; just use its host view
          auto Y_j_2d_h = Y_j->template getLocalView<host_memory_space> ();
          Y_j_1d_h = subview (Y_j_2d_h, rowRange, 0);
        }
        auto result_h_j = subview (result_h, j);
        KokkosBlas::dot (result_h_j, X_j_1d_h, Y_j_1d_h);
      }
      else { // run on device
        auto X_j_2d = X_j->template getLocalView<dev_memory_space> ();
        auto X_j_1d = subview (X_j_2d, rowRange, 0);
        decltype (X_j_1d) Y_j_1d;

        if (Y_copyToDev) {
          auto Y_j_2d_h = Y_j->template getLocalView<host_memory_space> ();
          auto Y_j_1d_h = subview (Y_j_2d_h, rowRange, 0);
          typename decltype (Y_j_1d)::non_const_type
            Y_j_1d_nc ("Y_j_1d", Y_j_1d_h.extent (0));
          deep_copy (Y_j_1d_nc, Y_j_1d_h);
          Y_j_1d = Y_j_1d_nc;
        }
        else { // Y already sync'd to dev; just use its dev view
          auto Y_j_2d = Y_j->template getLocalView<dev_memory_space> ();
          Y_j_1d = subview (Y_j_2d, rowRange, 0);
        }
        auto result_j = subview (result, j);
        KokkosBlas::dot (result_j, X_j_1d, Y_j_1d);
      }
    } // for each column j of X and Y
  }
  else { // numVecs == 1 || ! X_or_Y_have_nonconstant_stride
    if (runOnHost) {
      auto X_lcl_h = X.template getLocalView<host_memory_space> ();
      decltype (X_lcl_h) Y_lcl_h;

      if (Y_copyToHost) {
        auto Y_lcl = Y.template getLocalView<dev_memory_space> ();
        typename decltype (Y_lcl_h)::non_const_type
          Y_lcl_h_nc ("Y_lcl_h", Y_lcl.extent (0), Y_numVecs);
        deep_copy (Y_lcl_h_nc, Y_lcl);
        Y_lcl_h = Y_lcl_h_nc;
      }
      else { // Y already sync'd to host; just use its host view
        Y_lcl_h = Y.template getLocalView<host_memory_space> ();
      }

      if (numVecs == size_t (1)) {
        auto X_lcl_h_1d = subview (X_lcl_h, rowRange, 0);
        auto Y_lcl_h_1d = subview (Y_lcl_h, rowRange, 0);
        auto result_h_0d = subview (result_h, 0);
        KokkosBlas::dot (result_h_0d, X_lcl_h_1d, Y_lcl_h_1d);
      }
      else {
        auto X_lcl_h_2d = subview (X_lcl_h, rowRange,
                                   pair_type (0, X_numVecs));
        auto Y_lcl_h_2d = subview (Y_lcl_h, rowRange,
                                   pair_type (0, Y_numVecs));
        KokkosBlas::dot (result_h, X_lcl_h_2d, Y_lcl_h_2d);
      }
    }
    else { // run on device
      auto X_lcl = X.template getLocalView<dev_memory_space> ();
      decltype (X_lcl) Y_lcl;

      if (Y_copyToDev) {
        auto Y_lcl_h = Y.template getLocalView<host_memory_space> ();
        typename decltype (Y_lcl)::non_const_type
          Y_lcl_nc ("Y_lcl", Y_lcl_h.extent (0), Y_numVecs);
        deep_copy (Y_lcl_nc, Y_lcl_h);
        Y_lcl = Y_lcl_nc;
      }
      else { // Y already sync'd to dev; just use its dev view
        Y_lcl = Y.template getLocalView<dev_memory_space> ();
      }

      if (numVecs == size_t (1)) {
        auto X_lcl_1d = subview (X_lcl, rowRange, 0);
        auto Y_lcl_1d = subview (Y_lcl, rowRange, 0);
        auto result_0d = subview (result, 0);
        KokkosBlas::dot (result_0d, X_lcl_1d, Y_lcl_1d);
      }
      else {
        auto X_lcl_2d = subview (X_lcl, rowRange,
                                 pair_type (0, X_numVecs));
        auto Y_lcl_2d = subview (Y_lcl, rowRange,
                                 pair_type (0, Y_numVecs));
        KokkosBlas::dot (result, X_lcl_2d, Y_lcl_2d);
      }
    } // host or device?
  } // multivector with nonconstant stride?

  if (runOnHost && resultOnDevice) {
    // Copy result from host to device, where the user wanted it.
    deep_copy (result, result_h);
  }
  else if (! runOnHost && ! resultOnDevice) {
    // Copy result from device to host, where the user wanted it.
    deep_copy (result_h, result);
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
  using ::Kokkos::Impl::MemorySpaceAccess;
  using ::Kokkos::HostSpace;
  using ::Kokkos::View;
  using ::Tpetra::Details::iallreduce;
  typedef ::Tpetra::MultiVector<SC, LO, GO, NT> MV;
  typedef typename MV::dot_type dot_type;
  typedef typename MV::device_type device_type;
  typedef View<dot_type*, device_type> result_dev_view_type;
  typedef typename result_dev_view_type::HostMirror result_host_view_type;
  typedef typename device_type::memory_space dev_memory_space;

  auto map = X.getMap ();
  auto comm = map.is_null () ? Teuchos::null : map->getComm ();
  if (comm.is_null ()) { // calling process does not participate
    return std::shared_ptr< ::Tpetra::Details::CommRequest> (NULL);
  }

  const size_t X_numVecs = X.getNumVectors ();
  const size_t Y_numVecs = Y.getNumVectors ();
  const size_t numVecs = (X_numVecs > Y_numVecs) ? X_numVecs : Y_numVecs;

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

  // If the input communicator is an intercomm, then the input and
  // output buffers of iallreduce may not alias one another, so we
  // must allocate a temporary buffer for the local dot product(s).
  const bool needCopy = Details::isInterComm (*comm);

  // We know that resultRaw is host memory.  Can the device access it?
  // If so, lclDotRaw may be able to avoid a copy.
  const bool resultOnDevice =
    MemorySpaceAccess<dev_memory_space, HostSpace>::accessible;
  if (resultOnDevice) {
    result_dev_view_type gblResult (resultRaw, numVecs);
    result_dev_view_type lclResult = needCopy ?
      result_dev_view_type ("lclResult", numVecs) :
      gblResult;
    Details::lclDotRaw (lclResult.data (), X, Y, resultOnDevice);
    return iallreduce (lclResult, gblResult, ::Teuchos::REDUCE_SUM, *comm);
  }
  else {
    result_host_view_type gblResult (resultRaw, numVecs);
    result_host_view_type lclResult = needCopy ?
      result_host_view_type ("lclResult", numVecs) :
      gblResult;
    Details::lclDotRaw (lclResult.data (), X, Y, resultOnDevice);
    return iallreduce (lclResult, gblResult, ::Teuchos::REDUCE_SUM, *comm);
  }
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
  using ::Kokkos::Impl::MemorySpaceAccess;
  using ::Kokkos::HostSpace;
  using ::Kokkos::View;
  using ::Tpetra::Details::iallreduce;
  typedef ::Tpetra::MultiVector<SC, LO, GO, NT> MV;
  typedef typename MV::dot_type dot_type;
  typedef typename MV::device_type device_type;
  typedef View<dot_type, device_type> result_view_type;

  auto map = X.getMap ();
  auto comm = map.is_null () ? Teuchos::null : map->getComm ();
  if (comm.is_null ()) { // calling process does not participate
    return std::shared_ptr< ::Tpetra::Details::CommRequest> (NULL);
  }

  // If the input communicator is an intercomm, then the input and
  // output buffers of iallreduce may not alias one another, so we
  // must allocate a temporary buffer for the local dot product.
  const bool needCopy = Details::isInterComm (*comm);

  constexpr bool resultOnDevice = true; // 'result' is a device View
  result_view_type gblResult = result;
  result_view_type lclResult = needCopy ?
    result_view_type ("lclResult") :
    gblResult;
  Details::lclDotRaw (lclResult.data (), X, Y, resultOnDevice);
  return iallreduce (lclResult, gblResult, ::Teuchos::REDUCE_SUM, *comm);
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
  using ::Kokkos::Impl::MemorySpaceAccess;
  using ::Kokkos::HostSpace;
  using ::Kokkos::View;
  using ::Tpetra::Details::iallreduce;
  typedef ::Tpetra::MultiVector<SC, LO, GO, NT> MV;
  typedef typename MV::dot_type dot_type;
  typedef typename MV::device_type device_type;
  typedef View<dot_type*, device_type> result_view_type;

  auto map = X.getMap ();
  auto comm = map.is_null () ? ::Teuchos::null : map->getComm ();
  if (comm.is_null ()) { // calling process does not participate
    return std::shared_ptr< ::Tpetra::Details::CommRequest> (NULL);
  }

  // Check compatibility of number of columns; allow special cases of
  // a multivector dot a single vector, or vice versa.  It's OK to
  // throw without MPI communication, since the number of columns of a
  // Tpetra::MultiVector must be the same on all processes of its
  // communicator.
  const size_t X_numVecs = X.getNumVectors ();
  const size_t Y_numVecs = Y.getNumVectors ();
  if (X_numVecs != Y_numVecs &&
      X_numVecs != size_t (1) &&
      Y_numVecs != size_t (1)) {
    std::ostringstream os;
    os << "Tpetra::idot: X.getNumVectors() = " << X_numVecs
       << " != Y.getNumVectors() = " << Y_numVecs
       << ", but neither is 1.";
    throw std::invalid_argument (os.str ());
  }

  // If the input communicator is an intercomm, then the input and
  // output buffers of iallreduce may not alias one another, so we
  // must allocate a temporary buffer for the local dot product(s).
  const bool needCopy = Details::isInterComm (*comm);

  constexpr bool resultOnDevice = true; // 'result' is a device View
  result_view_type gblResult = result;
  result_view_type lclResult = needCopy ?
    result_view_type ("lclResult", result.extent (0)) :
    gblResult;
  Details::lclDotRaw (lclResult.data (), X, Y, resultOnDevice);
  return iallreduce (lclResult, gblResult, ::Teuchos::REDUCE_SUM, *comm);
}

} // namespace Tpetra

#endif // TPETRA_DETAILS_IDOT_HPP
