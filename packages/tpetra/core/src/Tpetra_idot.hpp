// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#include "Tpetra_iallreduce.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "KokkosBlas1_dot.hpp"
#include <stdexcept>
#include <sstream>

namespace Tpetra {
namespace Details {

// Temporary helper to get read-only view from multivector in requested space.
// TODO: when https://github.com/kokkos/kokkos/issues/3850 is resolved,
// remove this and just use templated getLocalView<Device>(ReadOnly).
template<typename MV>
struct GetReadOnly
{
  using DevView = typename MV::dual_view_type::t_dev::const_type;
  using HostView = typename MV::dual_view_type::t_host::const_type;

  template<typename exec_space>
  static DevView get(const MV& x, typename std::enable_if<std::is_same<exec_space, typename MV::execution_space>::value>::type* = nullptr)
  {
    return x.getLocalViewDevice(Tpetra::Access::ReadOnly);
  }

  template<typename exec_space>
  static HostView get(const MV& x, typename std::enable_if<!std::is_same<exec_space, typename MV::execution_space>::value>::type* = nullptr)
  {
    return x.getLocalViewHost(Tpetra::Access::ReadOnly);
  }
};

/// \brief Compute dot product locally. Where the kernel runs controlled by runOnDevice.

template<class MV, class ResultView, bool runOnDevice>
void idotLocal(const ResultView& localResult,
               const MV& X,
               const MV& Y)
{
  using pair_type = Kokkos::pair<size_t, size_t>;
  using exec_space = typename std::conditional<runOnDevice, typename MV::execution_space, Kokkos::DefaultHostExecutionSpace>::type;
  //if the execution space can access localResult, use it directly. Otherwise need to make a copy.
  static_assert(Kokkos::SpaceAccessibility<exec_space, typename ResultView::memory_space>::accessible,
      "idotLocal: Execution space must be able to access localResult");
  //If Dot executes on Serial, it requires the result to be HostSpace. If localResult is CudaUVMSpace,
  //we can just treat it like HostSpace.
  Kokkos::View<typename ResultView::data_type, typename exec_space::memory_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
    localResultUnmanaged(localResult.data(), localResult.extent(0));
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
  auto X_lcl = GetReadOnly<MV>::template get<exec_space>(X);
  auto Y_lcl = GetReadOnly<MV>::template get<exec_space>(Y);
  //Can the multivector dot kernel be used?
  bool useMVDot = X.isConstantStride() && Y.isConstantStride() && X_numVecs == Y_numVecs;
  if(useMVDot)
  {
    if (numVecs == size_t (1)) {
      auto X_lcl_1d = Kokkos::subview (X_lcl, rowRange, 0);
      auto Y_lcl_1d = Kokkos::subview (Y_lcl, rowRange, 0);
      auto result_0d = Kokkos::subview (localResultUnmanaged, 0);
      KokkosBlas::dot (result_0d, X_lcl_1d, Y_lcl_1d);
    }
    else {
      auto X_lcl_2d = Kokkos::subview (X_lcl, rowRange, pair_type (0, X_numVecs));
      auto Y_lcl_2d = Kokkos::subview (Y_lcl, rowRange, pair_type (0, Y_numVecs));
      KokkosBlas::dot (localResultUnmanaged, X_lcl_2d, Y_lcl_2d);
    }
  }
  else
  {
    auto XWhichVectors = Tpetra::getMultiVectorWhichVectors(X);
    auto YWhichVectors = Tpetra::getMultiVectorWhichVectors(Y);
    //Need to compute each dot individually, using 1D subviews of X_lcl and Y_lcl
    for(size_t vec = 0; vec < numVecs; vec++) {
      //Get the Vectors for each column of X and Y that will be dotted together.
      //Note: "which vectors" is not populated for constant stride MVs (but it's the identity mapping)
      size_t Xj = (numVecs == X_numVecs) ? vec : 0;
      Xj = X.isConstantStride() ? Xj : XWhichVectors[Xj];
      size_t Yj = (numVecs == Y_numVecs) ? vec : 0;
      Yj = Y.isConstantStride() ? Yj : YWhichVectors[Yj];
      auto Xcol = Kokkos::subview(X_lcl, rowRange, Xj);
      auto Ycol = Kokkos::subview(Y_lcl, rowRange, Yj);

      //Compute the rank-1 dot product, and place the result in an element of localResult
      KokkosBlas::dot(Kokkos::subview(localResultUnmanaged, vec), Xcol, Ycol);
    }
  }
}

//Helper to avoid extra instantiations of KokkosBlas::dot and iallreduce.
template<typename MV, typename ResultView>
struct IdotHelper
{
  using dot_type = typename MV::dot_type;

  //First version: use result directly
  template<typename exec_space>
  static std::shared_ptr< ::Tpetra::Details::CommRequest> run(
      const ResultView& globalResult, const MV& X, const MV& Y,
      typename std::enable_if<Kokkos::SpaceAccessibility<exec_space, typename ResultView::memory_space>::accessible>::type* = nullptr)

  {
    constexpr bool runOnDevice = std::is_same<exec_space, typename MV::execution_space>::value;
    idotLocal<MV, ResultView, runOnDevice>(globalResult, X, Y);
    //Fence because we're giving result of device kernel to MPI
    if(runOnDevice)
      exec_space().fence();
    auto comm = X.getMap()->getComm();
    return iallreduce(globalResult, globalResult, ::Teuchos::REDUCE_SUM, *comm);
  }

  //Second version: use a temporary local result, because exec_space can't write to globalResult
  template<typename exec_space>
  static std::shared_ptr< ::Tpetra::Details::CommRequest> run(
      const ResultView& globalResult, const MV& X, const MV& Y,
      typename std::enable_if<!Kokkos::SpaceAccessibility<exec_space, typename ResultView::memory_space>::accessible>::type* = nullptr)
  {
    constexpr bool runOnDevice = std::is_same<exec_space, typename MV::execution_space>::value;
    Kokkos::View<dot_type*, typename exec_space::memory_space> localResult(Kokkos::ViewAllocateWithoutInitializing("idot:localResult"), X.getNumVectors());
    idotLocal<MV, decltype(localResult), runOnDevice>(localResult, X, Y);
    //Fence because we're giving result of device kernel to MPI
    if(runOnDevice)
      exec_space().fence();
    auto comm = X.getMap()->getComm();
    return iallreduce(localResult, globalResult, ::Teuchos::REDUCE_SUM, *comm);
  }

};


/// \brief Internal (common) version of idot, a global dot product
/// that uses a non-blocking MPI reduction.
template<class MV, class ResultView>
std::shared_ptr< ::Tpetra::Details::CommRequest>
idotImpl(const ResultView& globalResult,
         const MV& X,
         const MV& Y)
{
  static_assert(std::is_same<typename ResultView::non_const_value_type, typename MV::dot_type>::value,
      "Tpetra::idot: result view's element type must match MV::dot_type");

  // Execution space to use for dot kernel(s) is whichever has access to up-to-date X.
  if(X.need_sync_device())
  {
    //run on host.
    return IdotHelper<MV, ResultView>::template run<Kokkos::DefaultHostExecutionSpace>(globalResult, X, Y);
  }
  else
  {
    //run on device.
    return IdotHelper<MV, ResultView>::template run<typename MV::execution_space>(globalResult, X, Y);
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
  Kokkos::View<dot_type*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
    resultView(resultRaw, numVecs);
  return Details::idotImpl(resultView, X, Y);
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
  return Details::idotImpl(result, X, Y);
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
  return Details::idotImpl(result1D, X, Y);
}

} // namespace Tpetra

#endif // TPETRA_DETAILS_IDOT_HPP
