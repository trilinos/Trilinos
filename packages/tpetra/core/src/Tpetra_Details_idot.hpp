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
#include <stdexcept>
#include <sstream>

namespace Tpetra {
namespace Details {
namespace Impl {

/// \brief Implementation detail of Idot specializations (see below in
///   this header file) with raw pointer or array DotViewType.
///
/// \tparam DotViewType Type of the idot() result argument.  Must be
///   either a raw pointer to nonconst (T* for some T) or a raw array
///   of nonconst (T[] for some T).
/// \tparam VecViewType Type of the vector or multivector View
///   arguments of idot().
/// \tparam dotViewRank Rank of view_type, the Kokkos::View that wraps
///   the input raw pointer or raw array given to getView().  Please
///   only use the default value here.
template<class DotViewType,
         class VecViewType,
         const int dotViewRank = static_cast<int> (VecViewType::rank) - 1>
struct GetDotView {
  /// \brief Kokkos::Device type of the Kokkos::View returned by getView().
  ///
  /// Assume that we can access DotViewType on host, but use
  /// VecViewType's device_type for the Kokkos kernel.
  typedef ::Kokkos::Device<typename VecViewType::HostMirror::execution_space,
                           ::Kokkos::HostSpace> host_device_type;

  /// \brief Array layout of the Kokkos::View returned by getView().
  ///
  /// We want view_type to have the same layout as VecViewType, since
  /// that is what both KokkosBlas::dot and iallreduce expect.
  /// However, only Kokkos::LayoutLeft and Kokkos::LayoutRight work
  /// for rank-0 Views, so we have to make sure that the layout is one
  /// of these.
  typedef typename std::conditional<std::is_same<typename VecViewType::array_layout,
                                                 ::Kokkos::LayoutLeft>::value,
                                    ::Kokkos::LayoutLeft,
                                    ::Kokkos::LayoutRight>::type array_layout;

  /// \brief Type of the Kokkos::View that wraps the raw pointer or
  ///   raw array input to getView().
  typedef ::Kokkos::View<DotViewType, array_layout, host_device_type> view_type;

  /// \brief Return a Kokkos::View that wraps the raw pointer or raw
  ///   array input argument \c raw.
  ///
  /// \param raw [in] Raw pointer or raw array to wrap.
  /// \param numVecs [in] Number of entries in the input raw array.
  ///   If the input is a raw pointer to a single entry, this may only
  ///   be one.
  static view_type getView (DotViewType raw, const size_t numVecs);
};

//! Specialization for returning a rank-0 Kokkos::View.
template<class DotViewType,
         class VecViewType>
struct GetDotView<DotViewType, VecViewType, 0> {
  /// \brief Kokkos::Device type of the Kokkos::View returned by getView().
  ///
  /// Assume that we can access DotViewType on host, but use
  /// VecViewType's device_type for the Kokkos kernel.
  typedef ::Kokkos::Device<typename VecViewType::HostMirror::execution_space,
                           ::Kokkos::HostSpace> host_device_type;

  /// \brief Type of the one entry of the Kokkos::View returned by getView().
  ///
  /// view_type must be a rank-0 Kokkos::View.  This means we can't
  /// just use DotViewType as the first template parameter of
  /// Kokkos::View.  However, we can use Kokkos::View<DotViewType> to
  /// get the type of the (only) entry of DotViewType.
  typedef typename ::Kokkos::View<DotViewType,
    host_device_type>::non_const_value_type dot_type;

  /// \brief Array layout of the Kokkos::View returned by getView().
  ///
  /// We want view_type to have the same layout as VecViewType, since
  /// that is what both KokkosBlas::dot and iallreduce expect.
  /// However, only Kokkos::LayoutLeft and Kokkos::LayoutRight work
  /// for rank-0 Views, so we have to make sure that the layout is one
  /// of these.
  typedef typename std::conditional<std::is_same<typename VecViewType::array_layout,
                                                 ::Kokkos::LayoutLeft>::value,
                                    ::Kokkos::LayoutLeft,
                                    ::Kokkos::LayoutRight>::type array_layout;

  /// \brief Type of the Kokkos::View that wraps the raw pointer or
  ///   raw array input to getView().
  typedef ::Kokkos::View<dot_type, array_layout, host_device_type> view_type;

  /// \brief Return a Kokkos::View that wraps the raw pointer or raw
  ///   array input argument \c raw.
  ///
  /// \param raw [in] Raw pointer or raw array to wrap.
  /// \param numVecs [in] Number of entries; must be exactly 1 in this case.
  static view_type getView (DotViewType raw, const size_t numVecs) {
    if (numVecs != 1) {
      std::ostringstream os;
      os << "You cannot create a Kokkos::View of rank 0 from a raw pointer or "
        "array, if the number of entries is not exactly 1.  You specified " <<
        numVecs << " entries.";
      throw std::invalid_argument (os.str ());
    }
    return view_type (raw);
  }
};

//! Specialization for returning a rank-1 Kokkos::View.
template<class DotViewType,
         class VecViewType>
struct GetDotView<DotViewType, VecViewType, 1> {
  /// \brief Kokkos::Device type of the Kokkos::View returned by getView().
  ///
  /// Assume that we can access DotViewType on host, but use
  /// VecViewType's device_type for the Kokkos kernel.
  typedef ::Kokkos::Device<typename VecViewType::HostMirror::execution_space,
                           ::Kokkos::HostSpace> host_device_type;

  /// \brief Array layout of the Kokkos::View returned by getView().
  ///
  /// We want view_type to have the same layout as VecViewType, since
  /// that is what both KokkosBlas::dot and iallreduce expect.
  /// However, if we're wrapping a raw pointer or raw array, we need
  /// to make sure that we only use either Kokkos::LayoutLeft or
  /// Kokkos::LayoutRight.
  typedef typename std::conditional<std::is_same<typename VecViewType::array_layout,
                                                 ::Kokkos::LayoutLeft>::value,
                                    ::Kokkos::LayoutLeft,
                                    ::Kokkos::LayoutRight>::type array_layout;

  /// \brief Type of the Kokkos::View that wraps the raw pointer or
  ///   raw array input to getView().
  typedef ::Kokkos::View<DotViewType, array_layout, host_device_type> view_type;

  /// \brief Return a Kokkos::View that wraps the raw pointer or raw
  ///   array input argument \c raw.
  ///
  /// \param raw [in] Raw pointer or raw array to wrap.
  /// \param numVecs [in] Number of entries in the input raw array.
  ///   If the input is a raw pointer to a single entry, this may only
  ///   be one.
  static view_type getView (DotViewType raw, const size_t numVecs) {
    static_assert (std::is_same<array_layout, ::Kokkos::LayoutLeft>::value ||
                   std::is_same<array_layout, ::Kokkos::LayoutRight>::value,
                   "This function requires that VecViewType be either "
                   "LayoutLeft or LayoutRight.");
    return view_type (raw, numVecs);
  }
};

/// \brief Struct implementing Tpetra::Details::idot.
///
/// We use this struct to do partial specialization on the different
/// kinds of arguments that the overloads of idot can accept.
///
/// \tparam DotViewType Type of the result of the dot product.  This
///   may be any of the following: a rank-0 or rank-1 Kokkos::View, a
///   raw pointer, or a raw array.  If a Kokkos::View, the rank must
///   be exactly one less than the rank of VecViewType.  If a raw
///   pointer or raw array, we assume that this is a pointer to
///   host-accessible memory.
/// \tparam VecViewType Type of the vector / multivector input
///   arguments.  Must be a rank-1 or rank-2 Kokkos::View.
/// \tparam dotViewTypeIsView Whether DotViewType is a Kokkos::View.
template<class DotViewType,
         class VecViewType,
         const bool dotViewTypeIsView = Kokkos::Impl::is_view<DotViewType>::value>
struct Idot {
  static std::shared_ptr<CommRequest>
  idot (const DotViewType& result,
        const VecViewType& X,
        const VecViewType& Y,
        const ::Teuchos::Comm<int>& comm);
};

//! Specialization for DotViewType and VecViewType both Kokkos::View.
template<class DotViewType,
         class VecViewType>
struct Idot<DotViewType, VecViewType, true> {
  typedef DotViewType dot_view_type;
  typedef VecViewType vec_view_type;

  static std::shared_ptr<CommRequest>
  idot (const dot_view_type& result,
        const vec_view_type& X,
        const vec_view_type& Y,
        const ::Teuchos::Comm<int>& comm)
  {
    static_assert (Kokkos::Impl::is_view<DotViewType>::value,
                   "DotViewType must be a Kokkos::View specialization.");
    static_assert (Kokkos::Impl::is_view<VecViewType>::value,
                   "VecViewType must be a Kokkos::View specialization.");
    static_assert (static_cast<int> (VecViewType::rank) != 1 ||
                   static_cast<int> (DotViewType::rank) == 0,
                   "If VecViewType has rank 1, then DotViewType must have "
                   "rank 0.");
    static_assert (static_cast<int> (VecViewType::rank) != 2 ||
                   static_cast<int> (DotViewType::rank) == 1,
                   "If VecViewType has rank 2, then DotViewType must have "
                   "rank 1.");
    using ::Tpetra::Details::iallreduce;
    using ::Teuchos::REDUCE_SUM;
    using ::Kokkos::Impl::MemorySpaceAccess;
    typedef typename dot_view_type::memory_space dot_mem_space;
    typedef typename vec_view_type::memory_space vec_mem_space;

    // TODO (mfh 19 Nov 2016) Once we get asynchronous kernel launch,
    // it would make sense to attach launching the iallreduce to the
    // dot product as a continuation task.  This would make the
    // returned "request" not actually a "CommRequest"; rather, it
    // would be the future of the asynchronous task.  The only issue
    // is that this approach would require MPI_THREAD_MULTIPLE (or
    // issuing MPI calls through an intermediary service that funnels
    // or serializes them).
    //
    // KokkosBlas::dot uses the execution space of its vector inputs
    // (in particular, the first vector input) to run the kernel.  If
    // the vectors' execution space cannot access the result's memory
    // space, then we need to make a deep copy of the result buffer.
    // We also need to do this if the input comm is an intercomm, for
    // a different reason.  The reasons differ for each of these
    // cases, but the code is the same.  (MemorySpaceAccess wants
    // memory spaces for both of its arguments.)
    const bool needCopy =
      ! MemorySpaceAccess<vec_mem_space, dot_mem_space>::accessible ||
      isInterComm (comm);

    if (needCopy) {
      typedef typename vec_view_type::device_type vec_device_type;
      auto resultCopy = ::Kokkos::create_mirror (vec_device_type (), result);
      static_assert (std::is_same<typename dot_view_type::array_layout,
                     typename decltype (resultCopy)::array_layout>::value,
                     "result and resultCopy must have the same array_layout.");
      try {
        ::KokkosBlas::dot (resultCopy, X, Y);
      }
      catch (std::exception& e) {
        std::ostringstream os;
        os << "Tpetra::Details::idot: KokkosBlas::dot threw an exception: "
           << e.what ();
        throw std::runtime_error (os.str ());
      }
      return iallreduce (resultCopy, result, REDUCE_SUM, comm);
    }
    else {
      try {
        ::KokkosBlas::dot (result, X, Y);
      }
      catch (std::exception& e) {
        std::ostringstream os;
        os << "Tpetra::Details::idot: KokkosBlas::dot threw an exception: "
           << e.what ();
        throw std::runtime_error (os.str ());
      }
      return iallreduce (result, result, REDUCE_SUM, comm);
    }
  }
};

//! Specialization for raw pointer dot product result type.
template<class DotViewType, class VecViewType>
struct Idot<DotViewType, VecViewType, false> {
  typedef DotViewType dot_view_type;
  typedef VecViewType vec_view_type;

  static std::shared_ptr<CommRequest>
  idot (dot_view_type result,
        const vec_view_type& X,
        const vec_view_type& Y,
        const ::Teuchos::Comm<int>& comm)
  {
    static_assert (! Kokkos::Impl::is_view<DotViewType>::value,
                   "DotViewType must NOT be a Kokkos::View specialization.");
    static_assert (Kokkos::Impl::is_view<VecViewType>::value,
                   "VecViewType must be a Kokkos::View specialization.");
    static_assert (static_cast<int> (VecViewType::rank) == 1 ||
                   static_cast<int> (VecViewType::rank) == 2,
                   "VecViewType must have rank 1 or rank 2.");
    static_assert (std::is_pointer<DotViewType>::value ||
                   std::is_array<DotViewType>::value,
                   "DotViewType must either be a pointer or an array.");

    // Even if VecViewType is rank 1, it's still syntactically correct
    // to call X.dimension_1().
    const size_t numVecs =
      (static_cast<int> (VecViewType::rank) == 1) ?
      static_cast<size_t> (1) :
      static_cast<size_t> (X.dimension_1 ());
    // tmp_dot_view_type must have rank exactly 1 less than vec_view_type.
    typedef typename GetDotView<DotViewType, VecViewType>::view_type
      tmp_dot_view_type;
    tmp_dot_view_type resultView =
      GetDotView<DotViewType, VecViewType>::getView (result, numVecs);
    typedef Idot<tmp_dot_view_type, vec_view_type, true> impl_type;
    return impl_type::idot (resultView, X, Y, comm);
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
  using ::Kokkos::subview;
  using ::Teuchos::Comm;
  using ::Teuchos::RCP;
  typedef ::Tpetra::MultiVector<SC, LO, GO, NT> mv_type;
  typedef typename mv_type::device_type DT;
  typedef typename DT::memory_space dev_memory_space;
  typedef typename ::Kokkos::View<SC*, DT>::host_mirror_space::memory_space
    host_memory_space;
  typedef typename ::Tpetra::Vector<SC, LO, GO, NT>::dot_type* result_view_type;
  typedef typename ::Kokkos::pair<size_t, size_t> pair_type;

  auto map = X.getMap ();
  RCP<const Comm<int> > comm = map.is_null () ? Teuchos::null : map->getComm ();
  if (! comm.is_null ()) {
    if (X.template need_sync<dev_memory_space> () &&
        ! X.template need_sync<host_memory_space> ()) { // use host version
      auto X_lcl = X.template getLocalView<host_memory_space> ();
      auto Y_lcl = Y.template getLocalView<host_memory_space> ();

      if (X.getNumVectors () == 1) {
        auto X_lcl_1d = subview (X_lcl, pair_type (0, X.getLocalLength ()), 0);
        auto Y_lcl_1d = subview (Y_lcl, pair_type (0, Y.getLocalLength ()), 0);
        typedef typename decltype (X_lcl_1d)::const_type vec_view_type;
        typedef Impl::Idot<result_view_type, vec_view_type> impl_type;
        return impl_type::idot (result, X_lcl_1d, Y_lcl_1d, *comm);
      }
      else {
        auto X_lcl_2d = subview (X_lcl, pair_type (0, X.getLocalLength ()),
                                 pair_type (0, X.getNumVectors ()));
        auto Y_lcl_2d = subview (Y_lcl, pair_type (0, Y.getLocalLength ()),
                                 pair_type (0, X.getNumVectors ()));
        typedef typename decltype (X_lcl_2d)::const_type vec_view_type;
        typedef Impl::Idot<result_view_type, vec_view_type> impl_type;
        return impl_type::idot (result, X_lcl_2d, Y_lcl_2d, *comm);
      }
    }
    else { // use device version
      auto X_lcl = X.template getLocalView<dev_memory_space> ();
      auto Y_lcl = Y.template getLocalView<dev_memory_space> ();

      if (X.getNumVectors () == 1) {
        auto X_lcl_1d = subview (X_lcl, pair_type (0, X.getLocalLength ()), 0);
        auto Y_lcl_1d = subview (Y_lcl, pair_type (0, Y.getLocalLength ()), 0);
        typedef typename decltype (X_lcl_1d)::const_type vec_view_type;
        typedef Impl::Idot<result_view_type, vec_view_type> impl_type;
        return impl_type::idot (result, X_lcl_1d, Y_lcl_1d, *comm);
      }
      else {
        auto X_lcl_2d = subview (X_lcl, pair_type (0, X.getLocalLength ()),
                                 pair_type (0, X.getNumVectors ()));
        auto Y_lcl_2d = subview (Y_lcl, pair_type (0, Y.getLocalLength ()),
                                 pair_type (0, X.getNumVectors ()));
        typedef typename decltype (X_lcl_2d)::const_type vec_view_type;
        typedef Impl::Idot<result_view_type, vec_view_type> impl_type;
        return impl_type::idot (result, X_lcl_2d, Y_lcl_2d, *comm);
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
  using ::Kokkos::subview;
  using ::Teuchos::Comm;
  using ::Teuchos::RCP;
  typedef ::Tpetra::Vector<SC, LO, GO, NT> vec_type;
  typedef typename vec_type::device_type DT;
  typedef typename DT::memory_space dev_memory_space;
  typedef typename ::Kokkos::View<SC*, DT>::host_mirror_space::memory_space
    host_memory_space;
  typedef typename ::Tpetra::Vector<SC, LO, GO, NT>::dot_type dot_type;
  typedef ::Kokkos::View<dot_type, DT> result_view_type;
  typedef ::Kokkos::pair<size_t, size_t> pair_type;

  auto map = X.getMap ();
  RCP<const Comm<int> > comm = map.is_null () ? Teuchos::null : map->getComm ();
  if (! comm.is_null ()) {
    if (X.template need_sync<dev_memory_space> () &&
        ! X.template need_sync<host_memory_space> ()) { // use host version
      auto X_lcl = X.template getLocalView<host_memory_space> ();
      auto Y_lcl = Y.template getLocalView<host_memory_space> ();
      auto X_lcl_1d = subview (X_lcl, pair_type (0, X.getLocalLength ()), 0);
      auto Y_lcl_1d = subview (Y_lcl, pair_type (0, Y.getLocalLength ()), 0);
      typedef typename decltype (X_lcl_1d)::const_type vec_view_type;
      typedef Impl::Idot<result_view_type, vec_view_type> impl_type;
      return impl_type::idot (result, X_lcl_1d, Y_lcl_1d, *comm);
    }
    else { // use device version
      auto X_lcl = X.template getLocalView<dev_memory_space> ();
      auto Y_lcl = Y.template getLocalView<dev_memory_space> ();
      auto X_lcl_1d = subview (X_lcl, pair_type (0, X.getLocalLength ()), 0);
      auto Y_lcl_1d = subview (Y_lcl, pair_type (0, Y.getLocalLength ()), 0);
      typedef typename decltype (X_lcl_1d)::const_type vec_view_type;
      typedef Impl::Idot<result_view_type, vec_view_type> impl_type;
      return impl_type::idot (result, X_lcl_1d, Y_lcl_1d, *comm);
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
  using ::Kokkos::subview;
  using ::Teuchos::Comm;
  using ::Teuchos::RCP;
  typedef ::Kokkos::pair<size_t, size_t> pair_type;
  typedef ::Tpetra::MultiVector<SC, LO, GO, NT> mv_type;
  typedef typename mv_type::device_type device_type;
  typedef typename device_type::memory_space dev_memory_space;
  typedef typename ::Tpetra::MultiVector<SC, LO, GO, NT>::dot_type dot_type;
  typedef ::Kokkos::View<dot_type*, device_type> result_view_type;
  // mfh 21 Dec 2016: It would be nicer to derive this from
  // Tpetra::MultiVector, but I'm not sure if Kokkos gives me a
  // convenient way to get the host memory space from a device_type,
  // rather than from a View.
  typedef typename result_view_type::host_mirror_space::memory_space
    host_memory_space;

  RCP<const typename mv_type::map_type> map = X.getMap ();
  RCP<const Comm<int> > comm = map.is_null () ? Teuchos::null : map->getComm ();
  if (! comm.is_null ()) {
    if (X.template need_sync<dev_memory_space> () &&
        ! X.template need_sync<host_memory_space> ()) { // use host version
      auto X_lcl = X.template getLocalView<host_memory_space> ();
      auto Y_lcl = Y.template getLocalView<host_memory_space> ();

      if (X.getNumVectors () == 1) {
        auto X_lcl_1d = subview (X_lcl, pair_type (0, X.getLocalLength ()), 0);
        auto Y_lcl_1d = subview (Y_lcl, pair_type (0, Y.getLocalLength ()), 0);
        typedef typename decltype (X_lcl_1d)::const_type vec_view_type;
        auto result_0d = subview (result, 0);
        typedef decltype (result_0d) result_0d_view_type;
        typedef Impl::Idot<result_0d_view_type, vec_view_type> impl_type;
        return impl_type::idot (result_0d, X_lcl_1d, Y_lcl_1d, *comm);
      }
      else {
        auto X_lcl_2d = subview (X_lcl, pair_type (0, X.getLocalLength ()),
                                 pair_type (0, X.getNumVectors ()));
        auto Y_lcl_2d = subview (Y_lcl, pair_type (0, Y.getLocalLength ()),
                                 pair_type (0, X.getNumVectors ()));
        typedef typename decltype (X_lcl_2d)::const_type vec_view_type;
        typedef Impl::Idot<result_view_type, vec_view_type> impl_type;
        return impl_type::idot (result, X_lcl_2d, Y_lcl_2d, *comm);
      }
    }
    else { // use device version
      auto X_lcl = X.template getLocalView<dev_memory_space> ();
      auto Y_lcl = Y.template getLocalView<dev_memory_space> ();

      if (X.getNumVectors () == 1) {
        auto X_lcl_1d = subview (X_lcl, pair_type (0, X.getLocalLength ()), 0);
        auto Y_lcl_1d = subview (Y_lcl, pair_type (0, Y.getLocalLength ()), 0);
        typedef typename decltype (X_lcl_1d)::const_type vec_view_type;
        auto result_0d = subview (result, 0);
        typedef decltype (result_0d) result_0d_view_type;
        typedef Impl::Idot<result_0d_view_type, vec_view_type> impl_type;
        return impl_type::idot (result_0d, X_lcl_1d, Y_lcl_1d, *comm);
      }
      else {
        auto X_lcl_2d = subview (X_lcl, pair_type (0, X.getLocalLength ()),
                                 pair_type (0, X.getNumVectors ()));
        auto Y_lcl_2d = subview (Y_lcl, pair_type (0, Y.getLocalLength ()),
                                 pair_type (0, X.getNumVectors ()));
        typedef typename decltype (X_lcl_2d)::const_type vec_view_type;
        typedef Impl::Idot<result_view_type, vec_view_type> impl_type;
        return impl_type::idot (result, X_lcl_2d, Y_lcl_2d, *comm);
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
