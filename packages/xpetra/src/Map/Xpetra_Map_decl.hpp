// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_MAP_DECL_HPP
#define XPETRA_MAP_DECL_HPP

#include "Xpetra_ConfigDefs.hpp"

#include <Tpetra_KokkosCompat_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>

#ifdef HAVE_XPETRA_EPETRA
#include "Epetra_config.h"
#endif

#ifdef HAVE_XPETRA_TPETRA
#include <Tpetra_Map.hpp>
#endif

namespace Xpetra {

// TODO move this typedef to another place
// Node which is used for Epetra. This can be either the
// Serial node or OpenMP node (but not both)
#ifdef HAVE_XPETRA_EPETRA
#ifdef EPETRA_HAVE_OMP
typedef Tpetra::KokkosCompat::KokkosOpenMPWrapperNode EpetraNode;
#else
typedef Tpetra::KokkosCompat::KokkosSerialWrapperNode EpetraNode;
#endif
#endif

enum UnderlyingLib {
  UseEpetra,
  UseTpetra,
  NotSpecified
};

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class Map
  : public Teuchos::Describable {
 public:
  typedef LocalOrdinal local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef Node node_type;

  typedef Kokkos::View<const global_ordinal_type *, typename Node::device_type> global_indices_array_device_type;

  //! @name Constructor/Destructor Methods
  //@{

  //! Destructor.
  virtual ~Map();

  //@}

  //! @name Attributes
  //@{

  //! The number of elements in this Map.
  virtual global_size_t getGlobalNumElements() const = 0;

  //! The number of elements belonging to the calling process.
  virtual size_t getLocalNumElements() const = 0;

  //! The index base for this Map.
  virtual GlobalOrdinal getIndexBase() const = 0;

  //! The minimum local index.
  virtual LocalOrdinal getMinLocalIndex() const = 0;

  //! The maximum local index on the calling process.
  virtual LocalOrdinal getMaxLocalIndex() const = 0;

  //! The minimum global index owned by the calling process.
  virtual GlobalOrdinal getMinGlobalIndex() const = 0;

  //! The maximum global index owned by the calling process.
  virtual GlobalOrdinal getMaxGlobalIndex() const = 0;

  //! The minimum global index over all processes in the communicator.
  virtual GlobalOrdinal getMinAllGlobalIndex() const = 0;

  //! The maximum global index over all processes in the communicator.
  virtual GlobalOrdinal getMaxAllGlobalIndex() const = 0;

  //! The local index corresponding to the given global index.
  virtual LocalOrdinal getLocalElement(GlobalOrdinal globalIndex) const = 0;

  //! The global index corresponding to the given local index.
  virtual GlobalOrdinal getGlobalElement(LocalOrdinal localIndex) const = 0;

  //! Return the process ranks and corresponding local indices for the given global indices.
  virtual LookupStatus getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal> &GIDList,
                                          const Teuchos::ArrayView<int> &nodeIDList,
                                          const Teuchos::ArrayView<LocalOrdinal> &LIDList) const = 0;

  //! Return the process ranks for the given global indices.
  virtual LookupStatus getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal> &GIDList,
                                          const Teuchos::ArrayView<int> &nodeIDList) const = 0;

  //! Return a view of the global indices owned by this process.
  virtual Teuchos::ArrayView<const GlobalOrdinal> getLocalElementList() const = 0;

  //! Return a view of the global indices owned by this process on the Map's device.
  virtual global_indices_array_device_type getMyGlobalIndicesDevice() const = 0;

  //@}

  //! @name Boolean tests
  //@{

  //! Whether the given local index is valid for this Map on this process.
  virtual bool isNodeLocalElement(LocalOrdinal localIndex) const = 0;

  //! Whether the given global index is valid for this Map on this process.
  virtual bool isNodeGlobalElement(GlobalOrdinal globalIndex) const = 0;

  //! True if this Map is distributed contiguously, else false.
  virtual bool isContiguous() const = 0;

  //! Whether this Map is globally distributed or locally replicated.
  virtual bool isDistributed() const = 0;

  //! True if and only if map is compatible with this Map.
  virtual bool isCompatible(const Map<LocalOrdinal, GlobalOrdinal, Node> &map) const = 0;

  //! True if and only if map is identical to this Map.
  virtual bool isSameAs(const Map<LocalOrdinal, GlobalOrdinal, Node> &map) const = 0;

  //@}

  //! @name
  //@{

  //! Get this Map's Comm object.
  virtual Teuchos::RCP<const Teuchos::Comm<int> > getComm() const = 0;

  //@}

  //! @name
  //@{

  //! Return a simple one-line description of this object.
  virtual std::string description() const = 0;

  //! Print this object with the given verbosity level to the given Teuchos::FancyOStream.
  virtual void describe(Teuchos::FancyOStream &out,
                        const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const = 0;

  //@}

  //! @name
  //@{

  //! Return a new Map with processes with zero elements removed.
  virtual RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > removeEmptyProcesses() const = 0;

  //! Replace this Map's communicator with a subset communicator.
  virtual RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
  replaceCommWithSubset(const Teuchos::RCP<const Teuchos::Comm<int> > &newComm) const = 0;

  //@}

  //! @name Xpetra specific
  //@{

  //! Get the library used by this object (Tpetra or Epetra?)
  virtual UnderlyingLib lib() const = 0;

  // TODO: find a better solution for this hack
  // The problem is that EpetraMap, TpetraMap and StridedMap all inherit Map. To have proper toEpetra() we
  // need to understand the type of underlying matrix. But in src/Map we have no knowledge of StridedMaps, so
  // we cannot check for it by casting. This function allows us to avoid the restriction, as StridedMap redefines
  // it to return the base map.
  virtual RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getMap() const;

#ifdef HAVE_XPETRA_TPETRA
  typedef typename Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>::local_map_type local_map_type;

  /// \brief Get the local Map for Kokkos kernels.
  virtual local_map_type getLocalMap() const = 0;
#endif

  //@}

};  // Map class

}  // namespace Xpetra

#define XPETRA_MAP_SHORT
#endif  // XPETRA_MAP_DECL_HPP
