// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_VECTORFACTORY_DECL_HPP
#define XPETRA_VECTORFACTORY_DECL_HPP

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Vector.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraVector_decl.hpp"
#endif
#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraVector.hpp"
#include "Xpetra_EpetraIntVector.hpp"
#endif

#include "Xpetra_BlockedMap_decl.hpp"
#include "Xpetra_BlockedVector_decl.hpp"
#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

template <class Scalar /* = Vector<>::scalar_type*/,
          class LocalOrdinal /* = typename Vector<Scalar>::local_ordinal_type*/,
          class GlobalOrdinal /* = typename Vector<Scalar, LocalOrdinal>::local_ordinal_type*/,
          class Node /* = typename Vector<Scalar, LocalOrdinal, GlobalOrdinal>::node_type*/>
class VectorFactory {
#undef XPETRA_VECTORFACTORY_SHORT
#include "Xpetra_UseShortNames.hpp"

 private:
  //! Private constructor. This is a static class.
  VectorFactory() = default;

 public:
  //! Constructor specifying the number of non-zeros for all rows.
  static Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  Build(const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& map, bool zeroOut = true) {
    XPETRA_MONITOR("VectorFactory::Build");

    RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>
        bmap = Teuchos::rcp_dynamic_cast<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>(map);

    if (!bmap.is_null()) {
      return rcp(new Xpetra::BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(bmap, zeroOut));
    }

#ifdef HAVE_XPETRA_TPETRA
    if (map->lib() == UseTpetra) {
      return rcp(new TpetraVector(map, zeroOut));
    }
#endif

    XPETRA_FACTORY_ERROR_IF_EPETRA(map->lib());
    XPETRA_FACTORY_END;
    TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
  }

};  // class VectorFactory

#define XPETRA_VECTORFACTORY_SHORT

#if defined(HAVE_XPETRA_EPETRA)

// we need the Epetra specialization only if Epetra is enabled
#if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)

// Specialization for Scalar=double, LO=GO=int and EpetraNode node
// Used both for Epetra and Tpetra
// For any other node definition the general default implementation is used which allows Tpetra only
template <>
class VectorFactory<double, int, int, EpetraNode> {
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef EpetraNode Node;

#undef XPETRA_VECTORFACTORY_SHORT
#include "Xpetra_UseShortNames.hpp"

 private:
  //! Private constructor. This is a static class.
  VectorFactory() = default;

 public:
  static RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  Build(const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& map,
        bool zeroOut = true);
};
#endif  // #if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)

// Specialization for Scalar=double, LO=int, GO=long long and EpetraNode
// Used both for Epetra and Tpetra
// For any other node definition the general default implementation is used which allows Tpetra only
#if !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)

template <>
class VectorFactory<double, int, long long, EpetraNode> {
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef long long GlobalOrdinal;
  typedef EpetraNode Node;

#undef XPETRA_VECTORFACTORY_SHORT
#include "Xpetra_UseShortNames.hpp"

 private:
  //! Private constructor. This is a static class.
  VectorFactory() = default;

 public:
  static RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  Build(const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& map,
        bool zeroOut = true);
};
#endif  // #if !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)
#define XPETRA_VECTORFACTORY_SHORT

// we need the Epetra specialization only if Epetra is enabled
#if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)

// Specialization for Scalar=int, LO=GO=int and EpetraNode
// Used both for Epetra and Tpetra
// For any other node definition the general default implementation is used which allows Tpetra only
template <>
class VectorFactory<int, int, int, EpetraNode> {
  typedef int Scalar;
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef EpetraNode Node;

#undef XPETRA_VECTORFACTORY_SHORT
#include "Xpetra_UseShortNames.hpp"

 private:
  //! Private constructor. This is a static class.
  VectorFactory() = default;

 public:
  static RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  Build(const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& map,
        bool zeroOut = true);
};
#endif  // #if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)

// we need the Epetra specialization only if Epetra is enabled
#if !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)

// Specialization for Scalar=int, LO=int, GO=long long and Serial node
// Used both for Epetra and Tpetra
// For any other node definition the general default implementation is used which allows Tpetra only

template <>
class VectorFactory<int, int, long long, EpetraNode> {
  typedef int Scalar;
  typedef int LocalOrdinal;
  typedef long long GlobalOrdinal;
  typedef EpetraNode Node;

#undef XPETRA_VECTORFACTORY_SHORT
#include "Xpetra_UseShortNames.hpp"

 private:
  //! Private constructor. This is a static class.
  VectorFactory() = default;

 public:
  static RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  Build(const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& map,
        bool zeroOut = true);
};
#endif  // !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)

#endif  // #if defined(HAVE_XPETRA_EPETRA)

}  // namespace Xpetra

#define XPETRA_VECTORFACTORY_SHORT
#endif  // XPETRA_VECTORFACTORY_DECL_HPP
