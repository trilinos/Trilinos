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

  Differences with Galeri1:
   - This function only supports mapType=Cartesian2D and Cartesian3D
   - Parameters that are not set by user but computed inside of this function are saved on the parameter list. This allows users to retrieve these parameters after the creation of the map.
     (In Galeri1, such parameters was set to -1 instead)
*/

// TODO: Check is some GlobalOrdinal are avoidable.

#ifndef GALERI_XPETRAMAPS_DECL_HPP
#define GALERI_XPETRAMAPS_DECL_HPP

#include <string>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ParameterList.hpp>

#include "Galeri_ConfigDefs.h"
#include "Galeri_Exception.h"

#ifdef HAVE_GALERI_XPETRA
#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_Map.hpp>  // for enum UnderlyingLib
#endif  // HAVE_GALERI_XPETRA

namespace Galeri::Xpetra {

using Teuchos::RCP;

//! Map creation function (for Tpetra and Xpetra::TpetraMap)
template <class LocalOrdinal, class GlobalOrdinal, class Map>
RCP<Map> CreateMap(const std::string& mapType, const Teuchos::RCP<const Teuchos::Comm<int> >& comm, Teuchos::ParameterList& list);

#ifdef HAVE_GALERI_XPETRA
//! Map creation function (for Xpetra::Map with an UnderlyingLib parameter)
template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP< ::Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > CreateMap(::Xpetra::UnderlyingLib lib, const std::string& mapType, const Teuchos::RCP<const Teuchos::Comm<int> >& comm, Teuchos::ParameterList& list);
#endif

}  // namespace Galeri

#endif
