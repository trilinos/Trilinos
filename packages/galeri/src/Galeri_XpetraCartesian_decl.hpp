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
#ifndef GALERI_XPETRACARTESIAN2D_DECL_HPP
#define GALERI_XPETRACARTESIAN2D_DECL_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

namespace Galeri::Xpetra::Maps {

template <class LocalOrdinal, class GlobalOrdinal, class Map>
Teuchos::RCP<Map> Cartesian1D(const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                              const GlobalOrdinal nx,
                              const GlobalOrdinal mx,
                              Teuchos::ParameterList& list);

template <class LocalOrdinal, class GlobalOrdinal, class Map>
Teuchos::RCP<Map> Cartesian2D(const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                              const GlobalOrdinal nx, const GlobalOrdinal ny,
                              const GlobalOrdinal mx, const GlobalOrdinal my,
                              Teuchos::ParameterList& list);

template <class LocalOrdinal, class GlobalOrdinal, class Map>
Teuchos::RCP<Map> Cartesian3D(const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                              const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
                              const GlobalOrdinal mx, const GlobalOrdinal my, const GlobalOrdinal mz,
                              Teuchos::ParameterList& list);

}  // namespace Galeri::Xpetra::Maps

#endif
