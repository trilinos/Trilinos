// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Helper to get ride of template parameters

// This file can be use for two purpose:
// 1) As an header of a user program.
//    In this case, this file must be included *after* other headers
//    and the types Scalar, LocalOrdinal, GlobalOrdinal, Node must be defined.
//    Note also that there is no #ifndef/#endif to protect again the multiple inclusion of this file.
//    User should create is own header file including this one:
//
//    Example:
//     #ifndef MY_HEADER
//     #define MY_HEADER
//     #include <MueLu_UseDefaultTypes.hpp>
//     #include <MueLu_UseShortNames.hpp>
//     #endif
//
// 2) Inside of MueLu to enhance the readability.
//
// template <class Scalar = Xpetra::MultiVector<>::scalar_type,
//           class LocalOrdinal = typename Xpetra::MultiVector<Scalar>::local_ordinal_type,
//           class GlobalOrdinal = typename Xpetra::MultiVector<Scalar, LocalOrdinal>::global_ordinal_type,
//           class Node = typename Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
//  class TpetraMultiVector : public virtual Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
//
//  #include <MueLu_UseShortNames.hpp>
//
//  myMethod(RCP<const Map> & map) { [...] } // instead of myMethod(RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map)
//
//  [...]
//
// }
//

#include "MueLu_UseShortNamesOrdinal.hpp"
#include "MueLu_UseShortNamesScalar.hpp"

//! @file MueLu_UseShortNamesOrdinal.hpp

// TODO / NOTE: This file should not be included at the global scope (to avoid name collision)
