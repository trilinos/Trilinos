// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
/*
  Direct translation of parts of Galeri matrix generator.
*/
#ifndef GALERI_XPETRACARTESIAN_DEF_HPP
#define GALERI_XPETRACARTESIAN_DEF_HPP

#include "Galeri_XpetraCartesian_decl.hpp"

#include "Galeri_Exception.h"
#include "Galeri_MapTraits.hpp"
#include "Galeri_XpetraUtils.hpp"

#include "Tpetra_Details_initializeKokkos.hpp"

namespace Galeri::Xpetra::Maps {

using global_size_t = size_t;

template <class LocalOrdinal, class GlobalOrdinal, class Map>
Teuchos::RCP<Map> Cartesian1D(const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                              const GlobalOrdinal nx,
                              const GlobalOrdinal mx,
                              Teuchos::ParameterList& list) {
  if (nx <= 0 || mx <= 0 || (mx > nx))
    throw Exception(__FILE__, __LINE__,
                    "Incorrect input parameter to Maps::Cartesian1D()",
                    "nx = " + toString(nx) +
                        ", mx = " + toString(mx));

  using GO = GlobalOrdinal;
  using LO = LocalOrdinal;

  int myPID = comm->getRank();

  GO startx, endx;
  Utils::getSubdomainData<GO>(nx, mx, myPID, startx, endx);

  list.set("lnx", Teuchos::as<LO>(endx - startx));
  list.set("lny", -1);
  list.set("lnz", -1);

  size_t numMyElements = endx - startx;

  Tpetra::Details::initializeKokkos();

  using Node                             = typename Map::node_type;
  using exec_space                       = typename Node::execution_space;
  using range_type                       = Kokkos::RangePolicy<LocalOrdinal, exec_space>;
  using device_type                      = typename Node::device_type;
  using global_indices_array_device_type = Kokkos::View<const GlobalOrdinal*, device_type>;

  auto myGlobalElements = typename global_indices_array_device_type::non_const_type("global_indices", numMyElements);

  Kokkos::parallel_for(
      "Galeri::MapFill", range_type(0, endx - startx),
      KOKKOS_LAMBDA(const LocalOrdinal c0) {
        auto i               = startx + c0;
        myGlobalElements(c0) = i;
      });

  global_indices_array_device_type elementList = myGlobalElements;

  global_size_t numGlobalElements = nx;
  return MapTraits<GO, Map>::Build(numGlobalElements, elementList, 0 /*indexBase*/, comm);
}

template <class LocalOrdinal, class GlobalOrdinal, class Map>
Teuchos::RCP<Map> Cartesian2D(const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                              const GlobalOrdinal nx, const GlobalOrdinal ny,
                              const GlobalOrdinal mx, const GlobalOrdinal my,
                              Teuchos::ParameterList& list) {
  if (nx <= 0 || ny <= 0 || mx <= 0 || my <= 0 || (mx > nx) || (my > ny))
    throw(Exception(__FILE__, __LINE__,
                    "Incorrect input parameter to Maps::Cartesian2D()",
                    "nx = " + toString(nx) +
                        ", ny = " + toString(ny) +
                        ", mx = " + toString(mx) +
                        ", my = " + toString(my)));

  using GO = GlobalOrdinal;
  using LO = LocalOrdinal;

  int myPID = comm->getRank();

  GO startx, starty, endx, endy;
  Utils::getSubdomainData(nx, mx, myPID % mx, startx, endx);
  Utils::getSubdomainData(ny, my, myPID / mx, starty, endy);

  list.set("lnx", Teuchos::as<LO>(endx - startx));
  list.set("lny", Teuchos::as<LO>(endy - starty));
  list.set("lnz", -1);

  size_t numMyElements = (endx - startx) * (endy - starty);

  Tpetra::Details::initializeKokkos();

  using Node                             = typename Map::node_type;
  using exec_space                       = typename Node::execution_space;
  using range_type                       = Kokkos::MDRangePolicy<LocalOrdinal, exec_space, Kokkos::Rank<2>>;
  using device_type                      = typename Node::device_type;
  using global_indices_array_device_type = Kokkos::View<const GlobalOrdinal*, device_type>;

  auto myGlobalElements = typename global_indices_array_device_type::non_const_type("global_indices", numMyElements);

  Kokkos::parallel_for(
      "Galeri::MapFill", range_type({0, 0}, {endx - startx, endy - starty}),
      KOKKOS_LAMBDA(const LocalOrdinal c0, const LocalOrdinal c1) {
        auto i                  = startx + c0;
        auto j                  = starty + c1;
        auto count              = c1 * (endx - startx) + c0;
        myGlobalElements(count) = j * nx + i;
      });

  global_indices_array_device_type elementList = myGlobalElements;
  global_size_t numGlobalElements              = nx * ny;
  return MapTraits<GO, Map>::Build(numGlobalElements, elementList, 0 /*indexBase*/, comm);
}

template <class LocalOrdinal, class GlobalOrdinal, class Map>
Teuchos::RCP<Map> Cartesian3D(const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                              const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
                              const GlobalOrdinal mx, const GlobalOrdinal my, const GlobalOrdinal mz,
                              Teuchos::ParameterList& list) {
  if (nx <= 0 || ny <= 0 || nz <= 0 ||
      mx <= 0 || my <= 0 || mz <= 0 ||
      (mx > nx) || (my > ny) || (mz > nz))
    throw Exception(__FILE__, __LINE__,
                    "Incorrect input parameter to Maps::Cartesian3D()",
                    "nx = " + toString(nx) +
                        ", ny = " + toString(ny) +
                        ", nz = " + toString(nz) +
                        ", mx = " + toString(mx) +
                        ", my = " + toString(my) +
                        ", mz = " + toString(mz));

  using GO = GlobalOrdinal;
  using LO = LocalOrdinal;

  GO mxy = mx * my;

  int myPID = comm->getRank();

  GO startx, starty, startz, endx, endy, endz;
  Utils::getSubdomainData(nx, mx, (myPID % mxy) % mx, startx, endx);
  Utils::getSubdomainData(ny, my, (myPID % mxy) / mx, starty, endy);
  Utils::getSubdomainData(nz, mz, myPID / mxy, startz, endz);

  list.set("lnx", Teuchos::as<LO>(endx - startx));
  list.set("lny", Teuchos::as<LO>(endy - starty));
  list.set("lnz", Teuchos::as<LO>(endz - startz));

  size_t numMyElements = (endx - startx) * (endy - starty) * (endz - startz);

  Tpetra::Details::initializeKokkos();

  using Node                             = typename Map::node_type;
  using exec_space                       = typename Node::execution_space;
  using range_type                       = Kokkos::MDRangePolicy<LocalOrdinal, exec_space, Kokkos::Rank<3>>;
  using device_type                      = typename Node::device_type;
  using global_indices_array_device_type = Kokkos::View<const GlobalOrdinal*, device_type>;

  auto myGlobalElements = typename global_indices_array_device_type::non_const_type("global_indices", numMyElements);

  Kokkos::parallel_for(
      "Galeri::MapFill", range_type({0, 0, 0}, {endx - startx, endy - starty, endz - startz}),
      KOKKOS_LAMBDA(const LocalOrdinal c0, const LocalOrdinal c1, const LocalOrdinal c2) {
        auto i                  = startx + c0;
        auto j                  = starty + c1;
        auto k                  = startz + c2;
        auto count              = c2 * ((endx - startx) * (endy - starty)) + c1 * (endx - startx) + c0;
        myGlobalElements(count) = k * (nx * ny) + j * nx + i;
      });

  global_indices_array_device_type elementList = myGlobalElements;
  global_size_t numGlobalElements              = nx * ny * nz;
  return MapTraits<GO, Map>::Build(numGlobalElements, elementList, 0 /*indexBase*/, comm);
}

}  // namespace Galeri::Xpetra::Maps

#define GALERI_XPETRACARTESIAN_INSTANT_TPETRA(LO, GO, N)                                                                                                                                                                                \
  template Teuchos::RCP<Tpetra::Map<LO, GO, N>> Galeri::Xpetra::Maps::Cartesian1D<LO, GO, Tpetra::Map<LO, GO, N>>(const Teuchos::RCP<const Teuchos::Comm<int>>& comm, const GO, const GO, Teuchos::ParameterList&);                     \
  template Teuchos::RCP<Tpetra::Map<LO, GO, N>> Galeri::Xpetra::Maps::Cartesian2D<LO, GO, Tpetra::Map<LO, GO, N>>(const Teuchos::RCP<const Teuchos::Comm<int>>& comm, const GO, const GO, const GO, const GO, Teuchos::ParameterList&); \
  template Teuchos::RCP<Tpetra::Map<LO, GO, N>> Galeri::Xpetra::Maps::Cartesian3D<LO, GO, Tpetra::Map<LO, GO, N>>(const Teuchos::RCP<const Teuchos::Comm<int>>& comm, const GO, const GO, const GO, const GO, const GO, const GO, Teuchos::ParameterList&);

#define GALERI_XPETRACARTESIAN_INSTANT_XPETRA(LO, GO, N)                                                                                                                                                                                            \
  template Teuchos::RCP<Xpetra::TpetraMap<LO, GO, N>> Galeri::Xpetra::Maps::Cartesian1D<LO, GO, Xpetra::TpetraMap<LO, GO, N>>(const Teuchos::RCP<const Teuchos::Comm<int>>& comm, const GO, const GO, Teuchos::ParameterList&);                     \
  template Teuchos::RCP<Xpetra::TpetraMap<LO, GO, N>> Galeri::Xpetra::Maps::Cartesian2D<LO, GO, Xpetra::TpetraMap<LO, GO, N>>(const Teuchos::RCP<const Teuchos::Comm<int>>& comm, const GO, const GO, const GO, const GO, Teuchos::ParameterList&); \
  template Teuchos::RCP<Xpetra::TpetraMap<LO, GO, N>> Galeri::Xpetra::Maps::Cartesian3D<LO, GO, Xpetra::TpetraMap<LO, GO, N>>(const Teuchos::RCP<const Teuchos::Comm<int>>& comm, const GO, const GO, const GO, const GO, const GO, const GO, Teuchos::ParameterList&);

#ifdef HAVE_GALERI_XPETRA

#define GALERI_XPETRACARTESIAN_INSTANT(LO, GO, N)  \
  GALERI_XPETRACARTESIAN_INSTANT_TPETRA(LO, GO, N) \
  GALERI_XPETRACARTESIAN_INSTANT_XPETRA(LO, GO, N)

#else

#define GALERI_XPETRACARTESIAN_INSTANT(LO, GO, N) \
  GALERI_XPETRACARTESIAN_INSTANT_TPETRA(LO, GO, N)

#endif

#endif
