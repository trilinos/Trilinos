// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_LOCALMAP_HPP
#define TPETRA_DETAILS_LOCALMAP_HPP

/// \file Tpetra_Details_LocalMap.hpp
/// \brief Declaration and definition of the Tpetra::Map class, an
///   implementation detail of Tpetra::Map.

#include "Tpetra_Details_FixedHashTable.hpp"
// #include "Tpetra_Details_OrdinalTraits.hpp" // comes in above
// #include "Kokkos_Core.hpp" // comes in above
#include "Tpetra_Details_LocalMap_fwd.hpp"

namespace Tpetra {
namespace Details {

/// \class LocalMap
/// \brief "Local" part of Map suitable for Kokkos kernels.
///
/// \warning This object's interface is not yet fixed.  We provide
///   this object currently only as a service to advanced users.
///
/// The "local" Map is suitable for use in Kokkos parallel operations
/// in the Map's native execution space, which is
/// <tt>Map::device_type::execution_space</tt>.
///
/// By "local," we mean that the object performs no MPI communication,
/// and can only access information that would never need MPI
/// communication, no matter what kind of Map this is.
template<class LocalOrdinal, class GlobalOrdinal, class DeviceType>
class LocalMap {
public:
  //! The type of local indices.
  using local_ordinal_type = LocalOrdinal;

  //! The type of global indices.
  using global_ordinal_type = GlobalOrdinal;

  //! The device type
  using device_type = DeviceType;

  //! The Kokkos execution space.
  using execution_space = typename device_type::execution_space;

  //! The Kokkos memory space.
  using memory_space = typename device_type::memory_space;

  //! Default constructor.
  KOKKOS_DEFAULTED_FUNCTION LocalMap() = default;

  LocalMap (const ::Tpetra::Details::FixedHashTable<GlobalOrdinal, LocalOrdinal, device_type>& glMap,
            const ::Kokkos::View<const GlobalOrdinal*, ::Kokkos::LayoutLeft, device_type>& lgMap,
            const GlobalOrdinal indexBase,
            const GlobalOrdinal myMinGid,
            const GlobalOrdinal myMaxGid,
            const GlobalOrdinal firstContiguousGid,
            const GlobalOrdinal lastContiguousGid,
            const LocalOrdinal numLocalElements,
            const bool contiguous) :
    glMap_ (glMap),
    lgMap_ (lgMap),
    indexBase_ (indexBase),
    myMinGid_ (myMinGid),
    myMaxGid_ (myMaxGid),
    firstContiguousGid_ (firstContiguousGid),
    lastContiguousGid_ (lastContiguousGid),
    numLocalElements_ (numLocalElements),
    contiguous_ (contiguous)
  {}

  //! The number of indices that live on the calling process.
  KOKKOS_INLINE_FUNCTION LocalOrdinal getLocalNumElements () const {
    return numLocalElements_;
  }

  //! The (global) index base.
  KOKKOS_INLINE_FUNCTION GlobalOrdinal getIndexBase () const {
    return indexBase_;
  }

  /// \brief Whether the Map is (locally) contiguous.
  ///
  /// This is conservative; a Map is "contiguous" if and only if
  /// it is stored that way.
  KOKKOS_INLINE_FUNCTION bool isContiguous () const {
    return contiguous_;
  }

  //! The minimum local index.
  KOKKOS_INLINE_FUNCTION LocalOrdinal getMinLocalIndex () const {
    return 0;
  }

  //! The maximum local index.
  KOKKOS_INLINE_FUNCTION LocalOrdinal
  getMaxLocalIndex () const
  {
    if (numLocalElements_ == 0) {
      return ::Tpetra::Details::OrdinalTraits<LocalOrdinal>::invalid ();
    } else { // Local indices are always zero-based.
      return static_cast<LocalOrdinal> (numLocalElements_ - 1);
    }
  }

  //! The minimum global index on the calling process.
  KOKKOS_INLINE_FUNCTION GlobalOrdinal getMinGlobalIndex () const {
    return myMinGid_;
  }

  //! The maximum global index on the calling process.
  KOKKOS_INLINE_FUNCTION GlobalOrdinal getMaxGlobalIndex () const {
    return myMaxGid_;
  }

  //! Get the local index corresponding to the given global index. (device only)
  KOKKOS_INLINE_FUNCTION LocalOrdinal
  getLocalElement (const GlobalOrdinal globalIndex) const
  {
    if (contiguous_) {
      if (globalIndex < myMinGid_ || globalIndex > myMaxGid_) {
        return ::Tpetra::Details::OrdinalTraits<LocalOrdinal>::invalid ();
      }
      return static_cast<LocalOrdinal> (globalIndex - myMinGid_);
    }
    else if (globalIndex >= firstContiguousGid_ &&
             globalIndex <= lastContiguousGid_) {
      return static_cast<LocalOrdinal> (globalIndex - firstContiguousGid_);
    }
    else {
      // If the given global index is not in the table, this returns
      // the same value as OrdinalTraits<LocalOrdinal>::invalid().
      return glMap_.get (globalIndex);
    }
  }

  //! Get the global index corresponding to the given local index. (device only)
  KOKKOS_INLINE_FUNCTION GlobalOrdinal
  getGlobalElement (const LocalOrdinal localIndex) const
  {
    if (localIndex < getMinLocalIndex () || localIndex > getMaxLocalIndex ()) {
      return ::Tpetra::Details::OrdinalTraits<GlobalOrdinal>::invalid ();
    }
    if (isContiguous ()) {
      return getMinGlobalIndex () + localIndex;
    }
    else {
      return lgMap_(localIndex);
    }
  }

private:
  //! Table that maps from global index to local index.
  ::Tpetra::Details::FixedHashTable<GlobalOrdinal, LocalOrdinal, device_type> glMap_;
  /// \brief Mapping from local indices to global indices.
  ///
  /// If this is empty, then it could be either that the Map is
  /// contiguous (meaning that we don't need to store all the
  /// global indices explicitly), or that the Map really does
  /// contain zero indices on the calling process.
  ///
  /// This has LayoutLeft so that we can call Kokkos::deep_copy to
  /// copy this between any two Kokkos Devices.  Otherwise, the
  /// Devices might have different default layouts, thus
  /// forbidding a deep_copy.  We use LayoutLeft instead of
  /// LayoutRight because LayoutRight is the default on non-CUDA
  /// Devices, and we want to make sure we catch assignment or
  /// copying from the default to the nondefault layout.
  ::Kokkos::View<const GlobalOrdinal*, ::Kokkos::LayoutLeft, device_type> lgMap_;

  GlobalOrdinal indexBase_          = 0;
  GlobalOrdinal myMinGid_           = Tpetra::Details::OrdinalTraits<GlobalOrdinal>::invalid();
  GlobalOrdinal myMaxGid_           = Tpetra::Details::OrdinalTraits<GlobalOrdinal>::invalid();
  GlobalOrdinal firstContiguousGid_ = Tpetra::Details::OrdinalTraits<GlobalOrdinal>::invalid();
  GlobalOrdinal lastContiguousGid_  = Tpetra::Details::OrdinalTraits<GlobalOrdinal>::invalid();
  LocalOrdinal numLocalElements_    = 0;
  bool contiguous_                  = false;
};

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_LOCALMAP_HPP

