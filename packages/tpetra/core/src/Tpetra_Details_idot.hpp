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

/// \file Tpetra_Details_idot.hpp
/// \brief Declaration of Tpetra::Details::idot, a nonblocking dot product.
///
/// \warning This file and its contents are implementation details of
///   Tpetra.  Users must not rely on them.
///
/// Tpetra::Details::idot implements a nonblocking dot product.  That
/// is the only thing in this file upon which Tpetra developers should
/// rely.  Tpetra developers should not rely on anything else in this
/// file.  <i>Users</i> may not rely on <i>anything</i> in this file!
///
/// If you want to find the only thing in this file that you are
/// supposed to use, search for "SKIP DOWN TO HERE" (no quotes).
/// "You" only refers to Tpetra developers.  Users, this file is not
/// for you!

#include "Tpetra_Details_iallreduce.hpp"
#include "Tpetra_Details_isInterComm.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Kokkos_Blas1_MV.hpp"

namespace Tpetra {
namespace Details {
namespace Impl {

template<class ViewType,
         const bool isView = Kokkos::Impl::is_view<ViewType>::value>
struct RankOfView {};

template<class T>
struct RankOfView<T*, false> {
  static const int rank = 1;
};

template<class T>
struct RankOfView<const T*, false> {
  static const int rank = 1;
};

template<class T>
struct RankOfView<T, true> {
  static const int rank = static_cast<int> (T::rank);
};

/// \brief Struct implementing Tpetra::Details::idot.
///
/// We use this struct to do partial specialization on the different
/// kinds of arguments that the overloads of idot can accept.
template<class ResultViewType,
         class VecViewType,
         const int resultViewRank = RankOfView<ResultViewType>::rank,
         const int vecViewRank = static_cast<int> (VecViewType::rank),
         const bool resultViewTypeIsView = Kokkos::Impl::is_view<ResultViewType>::value>
struct Idot {
  static std::shared_ptr<CommRequest>
  idot (const ResultViewType& result,
        const VecViewType& X,
        const VecViewType& Y,
        const ::Teuchos::Comm<int>& comm);
};

/// \brief Specialization for rank-1 result and rank-2 (multi)vectors.
///
/// Any LayoutLeft or LayoutRight rank-1 View can be assigned to a
/// View of the same rank and either of these layouts, so we don't
/// need two specializations.
template<class ResultViewType, class VecViewType>
struct Idot<ResultViewType, VecViewType, 1, 2, true> {
  typedef ResultViewType dot_view_type;
  typedef VecViewType multivec_view_type;

  static std::shared_ptr<CommRequest>
  idot (const dot_view_type& result,
        const multivec_view_type& X,
        const multivec_view_type& Y,
        const ::Teuchos::Comm<int>& comm)
  {
    static_assert (Kokkos::Impl::is_view<ResultViewType>::value,
                   "ResultViewType must be a Kokkos::View specialization.");
    static_assert (Kokkos::Impl::is_view<VecViewType>::value,
                   "VecViewType must be a Kokkos::View specialization.");
    static_assert (static_cast<int> (ResultViewType::rank) == 1,
                   "ResultViewType must be a rank-1 Kokkos::View.");
    static_assert (static_cast<int> (VecViewType::rank) == 2,
                   "VecViewType must be a rank-2 Kokkos::View.");
    using ::Tpetra::Details::iallreduce;
    using ::Teuchos::REDUCE_SUM;
    typedef typename dot_view_type::device_type device_type;

    // TODO (mfh 19 Nov 2016) Once we get asynchronous kernel launch, it
    // would make sense to attach launching the iallreduce to the dot
    // product as a continuation task.  This would make the returned
    // "request" not actually a "CommRequest"; rather, it would be the
    // future of the asynchronous task.  The only issue is that this
    // approach would require MPI_THREAD_MULTIPLE (or issuing MPI calls
    // through an intermediary service that funnels or serializes them).
    ::KokkosBlas::dot (result, X, Y);

    if (isInterComm (comm)) {
      auto resultCopy = ::Kokkos::create_mirror (device_type (), result);
      ::Kokkos::deep_copy (resultCopy, result);
      return iallreduce (resultCopy, result, REDUCE_SUM, comm);
    }
    else {
      return iallreduce (result, result, REDUCE_SUM, comm);
    }
  }
};

//! Specialization for raw pointer result and rank-2 (multi)vectors.
template<class ResultViewType, class VecViewType>
struct Idot<ResultViewType, VecViewType, 1, 2, false> {
  typedef ResultViewType dot_view_type;
  typedef VecViewType multivec_view_type;

  static std::shared_ptr<CommRequest>
  idot (dot_view_type result,
        const multivec_view_type& X,
        const multivec_view_type& Y,
        const ::Teuchos::Comm<int>& comm)
  {
    static_assert (! Kokkos::Impl::is_view<ResultViewType>::value,
                   "ResultViewType must NOT be a Kokkos::View specialization.");
    static_assert (Kokkos::Impl::is_view<VecViewType>::value,
                   "VecViewType must be a Kokkos::View specialization.");
    static_assert (static_cast<int> (VecViewType::rank) == 2,
                   "VecViewType must be a rank-2 Kokkos::View.");
    static_assert (std::is_pointer<ResultViewType>::value ||
                   std::is_array<ResultViewType>::value,
                   "ResultViewType must either be a pointer or an array.");

    // Assume that we can access dot_view_type on host, but use
    // multivec_view_type's device_type for the Kokkos kernel.
    typedef typename multivec_view_type::HostMirror::execution_space host_execution_space;
    typedef Kokkos::Device<host_execution_space, Kokkos::HostSpace> host_device_type;
    typedef Kokkos::View<ResultViewType, host_device_type> result_view_type;

    const size_t numVecs = X.dimension_1 ();
    result_view_type resultView (result, numVecs);
    return Idot<result_view_type, VecViewType, 1, 2, true>::idot (resultView, X, Y, comm);
  }
};

//! Specialization for rank-0 result and rank-1 vectors.
template<class ResultViewType, class VecViewType>
struct Idot<ResultViewType, VecViewType, 0, 1, true> {
  typedef ResultViewType dot_view_type;
  typedef VecViewType vec_view_type;

  static std::shared_ptr<CommRequest>
  idot (const dot_view_type& result,
        const vec_view_type& X,
        const vec_view_type& Y,
        const ::Teuchos::Comm<int>& comm)
  {
    static_assert (Kokkos::Impl::is_view<ResultViewType>::value,
                   "ResultViewType must be a Kokkos::View specialization.");
    static_assert (Kokkos::Impl::is_view<VecViewType>::value,
                   "VecViewType must be a Kokkos::View specialization.");
    static_assert (static_cast<int> (ResultViewType::rank) == 0,
                   "ResultViewType must be a rank-0 Kokkos::View.");
    static_assert (static_cast<int> (VecViewType::rank) == 1,
                   "VecViewType must be a rank-1 Kokkos::View.");
    using ::Tpetra::Details::iallreduce;
    using ::Teuchos::REDUCE_SUM;
    typedef typename dot_view_type::device_type device_type;

    // TODO (mfh 19 Nov 2016) Once we get asynchronous kernel launch, it
    // would make sense to attach launching the iallreduce to the dot
    // product as a continuation task.  This would make the returned
    // "request" not actually a "CommRequest"; rather, it would be the
    // future of the asynchronous task.  The only issue is that this
    // approach would require MPI_THREAD_MULTIPLE (or issuing MPI calls
    // through an intermediary service that funnels or serializes them).
    ::KokkosBlas::dot (result, X, Y);

    if (isInterComm (comm)) {
      auto resultCopy = ::Kokkos::create_mirror (device_type (), result);
      ::Kokkos::deep_copy (resultCopy, result);
      return iallreduce (resultCopy, result, REDUCE_SUM, comm);
    }
    else {
      return iallreduce (result, result, REDUCE_SUM, comm);
    }
  }
};

/// \brief Specialization for raw pointer result and rank-1 vectors.
///
/// Raw pointer result counts as a "rank-1 View."
template<class ResultViewType, class VecViewType>
struct Idot<ResultViewType, VecViewType, 1, 1, false> {
  typedef ResultViewType dot_view_type;
  typedef VecViewType vec_view_type;

  static std::shared_ptr<CommRequest>
  idot (dot_view_type result,
        const vec_view_type& X,
        const vec_view_type& Y,
        const ::Teuchos::Comm<int>& comm)
  {
    static_assert (! Kokkos::Impl::is_view<ResultViewType>::value,
                   "ResultViewType must NOT be a Kokkos::View specialization.");
    static_assert (Kokkos::Impl::is_view<VecViewType>::value,
                   "VecViewType must be a Kokkos::View specialization.");
    static_assert (static_cast<int> (VecViewType::rank) == 1,
                   "VecViewType must be a rank-1 Kokkos::View.");
    static_assert (std::is_pointer<ResultViewType>::value ||
                   std::is_array<ResultViewType>::value,
                   "ResultViewType must either be a pointer or an array.");

    // Assume that we can access dot_view_type on host, but use
    // vec_view_type's device_type for the Kokkos kernel.
    typedef typename vec_view_type::HostMirror::execution_space host_execution_space;
    typedef Kokkos::Device<host_execution_space, Kokkos::HostSpace> host_device_type;
    // We must ensure that the result_view_type is a rank-0 View.
    // Furthermore, ResultViewType may be either a pointer (dot_type*)
    // or an array (dot_type[], dot_type[1], etc.).  This is why we go
    // through a bit of trouble to extract the dot product value type.
    typedef typename Kokkos::View<ResultViewType, host_device_type>::value_type dot_type;
    typedef Kokkos::View<dot_type, host_device_type> result_view_type;

    result_view_type resultView (result);
    return Idot<result_view_type, VecViewType, 0, 1, true>::idot (resultView, X, Y, comm);
  }
};

} // namespace Impl

//
// SKIP DOWN TO HERE
//

/// \brief Nonblocking dot product, with either Tpetra::MultiVector or
///   Tpetra::Vector inputs, and raw pointer output.
///
/// \param result [out] Output; raw array (/ pointer) to the return
///   value(s).  It is only valid to read this after calling wait() on
///   the return value.  It must be legal to write to the first
///   <tt>X.getNumVectors()</tt> entries of this array.
///
/// \param X [in] First input Tpetra::MultiVector or Tpetra::Vector.
///   This must have same number of rows (globally, and on each (MPI)
///   process) as Y.  If this is a Tpetra::MultiVector, then this must
///   have the same number of columns as Y.
///
/// \param Y [in] Second input Tpetra::MultiVector or Tpetra::Vector.
///   This must have same number of rows (globally, and on each (MPI)
///   process) as X.  If this is a Tpetra::MultiVector, then this must
///   have the same number of columns as X.
///
/// \return Pointer to an object representing the nonblocking
///   collective (communication operation).  Call wait() on this
///   object to complete the collective.  After calling wait(), you
///   may read the result.
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
std::shared_ptr<CommRequest>
idot (typename ::Tpetra::Vector<SC, LO, GO, NT>::dot_type* result,
      const ::Tpetra::MultiVector<SC, LO, GO, NT>& X,
      const ::Tpetra::MultiVector<SC, LO, GO, NT>& Y)
{
  using ::Teuchos::Comm;
  using ::Teuchos::RCP;
  typedef ::Tpetra::MultiVector<SC, LO, GO, NT> mv_type;
  typedef typename mv_type::device_type device_type;
  typedef typename device_type::memory_space dev_memory_space;
  typedef typename ::Kokkos::View<SC*, device_type>::host_mirror_space::memory_space host_memory_space;
  typedef typename ::Tpetra::Vector<SC, LO, GO, NT>::dot_type* result_view_type;

  auto map = X.getMap ();
  RCP<const Comm<int> > comm = map.is_null () ? Teuchos::null : map->getComm ();
  if (! comm.is_null ()) {
    if (X.template need_sync<dev_memory_space> () &&
        ! X.template need_sync<host_memory_space> ()) { // use host version
      auto X_lcl = X.template getLocalView<host_memory_space> ();
      auto Y_lcl = Y.template getLocalView<host_memory_space> ();

      if (X.getNumVectors () == 1) {
        auto X_lcl_1d = Kokkos::subview (X_lcl, Kokkos::pair<size_t, size_t> (0, X.getLocalLength ()), 0);
        auto Y_lcl_1d = Kokkos::subview (Y_lcl, Kokkos::pair<size_t, size_t> (0, Y.getLocalLength ()), 0);
        typedef typename decltype (X_lcl_1d)::const_type vec_view_type;
        return Impl::Idot<result_view_type, vec_view_type>::idot (result, X_lcl_1d, Y_lcl_1d, *comm);
      }
      else {
        auto X_lcl_2d = Kokkos::subview (X_lcl, Kokkos::pair<size_t, size_t> (0, X.getLocalLength ()),
                                         Kokkos::pair<size_t, size_t> (0, X.getNumVectors ()));
        auto Y_lcl_2d = Kokkos::subview (Y_lcl, Kokkos::pair<size_t, size_t> (0, Y.getLocalLength ()),
                                         Kokkos::pair<size_t, size_t> (0, X.getNumVectors ()));
        typedef typename decltype (X_lcl_2d)::const_type vec_view_type;
        return Impl::Idot<result_view_type, vec_view_type>::idot (result, X_lcl_2d, Y_lcl_2d, *comm);
      }
    }
    else { // use device version
      auto X_lcl = X.template getLocalView<dev_memory_space> ();
      auto Y_lcl = Y.template getLocalView<dev_memory_space> ();

      if (X.getNumVectors () == 1) {
        auto X_lcl_1d = Kokkos::subview (X_lcl, Kokkos::pair<size_t, size_t> (0, X.getLocalLength ()), 0);
        auto Y_lcl_1d = Kokkos::subview (Y_lcl, Kokkos::pair<size_t, size_t> (0, Y.getLocalLength ()), 0);
        typedef typename decltype (X_lcl_1d)::const_type vec_view_type;
        return Impl::Idot<result_view_type, vec_view_type>::idot (result, X_lcl_1d, Y_lcl_1d, *comm);
      }
      else {
        auto X_lcl_2d = Kokkos::subview (X_lcl, Kokkos::pair<size_t, size_t> (0, X.getLocalLength ()),
                                         Kokkos::pair<size_t, size_t> (0, X.getNumVectors ()));
        auto Y_lcl_2d = Kokkos::subview (Y_lcl, Kokkos::pair<size_t, size_t> (0, Y.getLocalLength ()),
                                         Kokkos::pair<size_t, size_t> (0, X.getNumVectors ()));
        typedef typename decltype (X_lcl_2d)::const_type vec_view_type;
        return Impl::Idot<result_view_type, vec_view_type>::idot (result, X_lcl_2d, Y_lcl_2d, *comm);
      }
    }
  }
  else { // calling process does not participate
    return std::shared_ptr<CommRequest> (NULL);
  }
}

/// \brief Nonblocking dot product, with Tpetra::Vector inputs, and
///   rank-0 (single value) Kokkos::View output.
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
std::shared_ptr<CommRequest>
idot (const Kokkos::View<typename ::Tpetra::Vector<SC, LO, GO, NT>::dot_type,
        typename ::Tpetra::Vector<SC, LO, GO, NT>::device_type>& result,
      const ::Tpetra::Vector<SC, LO, GO, NT>& X,
      const ::Tpetra::Vector<SC, LO, GO, NT>& Y)
{
  using ::Teuchos::Comm;
  using ::Teuchos::RCP;
  typedef ::Tpetra::Vector<SC, LO, GO, NT> vec_type;
  typedef typename vec_type::device_type device_type;
  typedef typename device_type::memory_space dev_memory_space;
  typedef typename ::Kokkos::View<SC*, device_type>::host_mirror_space::memory_space host_memory_space;
  typedef typename ::Tpetra::Vector<SC, LO, GO, NT>::dot_type dot_type;
  typedef Kokkos::View<dot_type, device_type> result_view_type;

  auto map = X.getMap ();
  RCP<const Comm<int> > comm = map.is_null () ? Teuchos::null : map->getComm ();
  if (! comm.is_null ()) {
    if (X.template need_sync<dev_memory_space> () &&
        ! X.template need_sync<host_memory_space> ()) { // use host version
      auto X_lcl = X.template getLocalView<host_memory_space> ();
      auto Y_lcl = Y.template getLocalView<host_memory_space> ();
      auto X_lcl_1d = Kokkos::subview (X_lcl, Kokkos::pair<size_t, size_t> (0, X.getLocalLength ()), 0);
      auto Y_lcl_1d = Kokkos::subview (Y_lcl, Kokkos::pair<size_t, size_t> (0, Y.getLocalLength ()), 0);
      typedef typename decltype (X_lcl_1d)::const_type vec_view_type;
      return Impl::Idot<result_view_type, vec_view_type>::idot (result, X_lcl_1d, Y_lcl_1d, *comm);
    }
    else { // use device version
      auto X_lcl = X.template getLocalView<dev_memory_space> ();
      auto Y_lcl = Y.template getLocalView<dev_memory_space> ();
      auto X_lcl_1d = Kokkos::subview (X_lcl, Kokkos::pair<size_t, size_t> (0, X.getLocalLength ()), 0);
      auto Y_lcl_1d = Kokkos::subview (Y_lcl, Kokkos::pair<size_t, size_t> (0, Y.getLocalLength ()), 0);
      typedef typename decltype (X_lcl_1d)::const_type vec_view_type;
      return Impl::Idot<result_view_type, vec_view_type>::idot (result, X_lcl_1d, Y_lcl_1d, *comm);
    }
  }
  else { // calling process does not participate
    return std::shared_ptr<CommRequest> (NULL);
  }
}

/// \brief Nonblocking dot product, with Tpetra::MultiVector inputs,
///   and rank-1 (one-dimensional array) Kokkos::View output.
///
/// \param result [out] Output; rank-1 Kokkos::View.  It is only valid
///   to read the entries of this after calling wait() on the return
///   value.
///
/// \param X [in] First input Tpetra::MultiVector or Tpetra::Vector.
///   This must have same number of rows (globally, and on each (MPI)
///   process) as Y.  If this is a Tpetra::MultiVector, then this must
///   have the same number of columns as Y.
///
/// \param Y [in] Second input Tpetra::MultiVector or Tpetra::Vector.
///   This must have same number of rows (globally, and on each (MPI)
///   process) as X.  If this is a Tpetra::MultiVector, then this must
///   have the same number of columns as X.
///
/// \return Pointer to an object representing the nonblocking
///   collective (communication operation).  Call wait() on this
///   object to complete the collective.  After calling wait(), you
///   may read the results.
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
std::shared_ptr<CommRequest>
idot (const Kokkos::View<typename ::Tpetra::MultiVector<SC, LO, GO, NT>::dot_type*,
        typename ::Tpetra::MultiVector<SC, LO, GO, NT>::device_type>& result,
      const ::Tpetra::MultiVector<SC, LO, GO, NT>& X,
      const ::Tpetra::MultiVector<SC, LO, GO, NT>& Y)
{
  using ::Teuchos::Comm;
  using ::Teuchos::RCP;
  typedef ::Tpetra::MultiVector<SC, LO, GO, NT> mv_type;
  typedef typename mv_type::device_type device_type;
  typedef typename device_type::memory_space dev_memory_space;
  typedef typename ::Kokkos::View<SC*, device_type>::host_mirror_space::memory_space host_memory_space;
  typedef Kokkos::View<typename ::Tpetra::MultiVector<SC, LO, GO, NT>::dot_type*,
    typename ::Tpetra::MultiVector<SC, LO, GO, NT>::device_type> result_view_type;

  auto map = X.getMap ();
  RCP<const Comm<int> > comm = map.is_null () ? Teuchos::null : map->getComm ();
  if (! comm.is_null ()) {
    if (X.template need_sync<dev_memory_space> () &&
        ! X.template need_sync<host_memory_space> ()) { // use host version
      auto X_lcl = X.template getLocalView<host_memory_space> ();
      auto Y_lcl = Y.template getLocalView<host_memory_space> ();

      if (X.getNumVectors () == 1) {
        auto X_lcl_1d = Kokkos::subview (X_lcl, Kokkos::pair<size_t, size_t> (0, X.getLocalLength ()), 0);
        auto Y_lcl_1d = Kokkos::subview (Y_lcl, Kokkos::pair<size_t, size_t> (0, Y.getLocalLength ()), 0);
        typedef typename decltype (X_lcl_1d)::const_type vec_view_type;
        return Impl::Idot<result_view_type, vec_view_type>::idot (result, X_lcl_1d, Y_lcl_1d, *comm);
      }
      else {
        auto X_lcl_2d = Kokkos::subview (X_lcl, Kokkos::pair<size_t, size_t> (0, X.getLocalLength ()),
                                         Kokkos::pair<size_t, size_t> (0, X.getNumVectors ()));
        auto Y_lcl_2d = Kokkos::subview (Y_lcl, Kokkos::pair<size_t, size_t> (0, Y.getLocalLength ()),
                                         Kokkos::pair<size_t, size_t> (0, X.getNumVectors ()));
        typedef typename decltype (X_lcl_2d)::const_type vec_view_type;
        return Impl::Idot<result_view_type, vec_view_type>::idot (result, X_lcl_2d, Y_lcl_2d, *comm);
      }
    }
    else { // use device version
      auto X_lcl = X.template getLocalView<dev_memory_space> ();
      auto Y_lcl = Y.template getLocalView<dev_memory_space> ();

      if (X.getNumVectors () == 1) {
        auto X_lcl_1d = Kokkos::subview (X_lcl, Kokkos::pair<size_t, size_t> (0, X.getLocalLength ()), 0);
        auto Y_lcl_1d = Kokkos::subview (Y_lcl, Kokkos::pair<size_t, size_t> (0, Y.getLocalLength ()), 0);
        typedef typename decltype (X_lcl_1d)::const_type vec_view_type;
        return Impl::Idot<result_view_type, vec_view_type>::idot (result, X_lcl_1d, Y_lcl_1d, *comm);
      }
      else {
        auto X_lcl_2d = Kokkos::subview (X_lcl, Kokkos::pair<size_t, size_t> (0, X.getLocalLength ()),
                                         Kokkos::pair<size_t, size_t> (0, X.getNumVectors ()));
        auto Y_lcl_2d = Kokkos::subview (Y_lcl, Kokkos::pair<size_t, size_t> (0, Y.getLocalLength ()),
                                         Kokkos::pair<size_t, size_t> (0, X.getNumVectors ()));
        typedef typename decltype (X_lcl_2d)::const_type vec_view_type;
        return Impl::Idot<result_view_type, vec_view_type>::idot (result, X_lcl_2d, Y_lcl_2d, *comm);
      }
    }
  }
  else { // calling process does not participate
    return std::shared_ptr<CommRequest> (NULL);
  }
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_IDOT_HPP
