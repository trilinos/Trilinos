// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <string>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ParameterList.hpp>
#include <MueLu_ConfigDefs.hpp>

#include <Galeri_XpetraMaps.hpp>

#ifdef HAVE_GALERI_XPETRA

#include <TpetraCore_ETIHelperMacros.h>
TPETRA_ETI_MANGLING_TYPEDEFS()

namespace Galeri {
namespace Xpetra {

// Add other Galeri functions as needed
#define MUELU_ETI_GROUP(LO, GO, NO) \
  template Teuchos::RCP<::Xpetra::Map<LO, GO, NO>> CreateMap<LO, GO, NO>(::Xpetra::UnderlyingLib lib, const std::string& mapType, const Teuchos::RCP<const Teuchos::Comm<int>>& comm, Teuchos::ParameterList& list);

#include <MueLu_ETI_3arg.hpp>

}  // namespace Xpetra
}  // namespace Galeri

#endif  // ifdef HAVE_GALERI_XPETRA
//#endif //ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION
