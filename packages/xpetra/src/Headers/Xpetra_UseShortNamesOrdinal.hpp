// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Get rid of template parameters

// New definition of types using the types LocalOrdinal, GlobalOrdinal, Node of the current context.
#ifdef XPETRA_MAP_SHORT
using Map [[maybe_unused]] = Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_MAPUTILS_SHORT
using MapUtils [[maybe_unused]] = Xpetra::MapUtils<LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_MAPFACTORY_SHORT
using MapFactory [[maybe_unused]] = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_BLOCKEDMAP_SHORT
using BlockedMap [[maybe_unused]] = Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_CRSGRAPH_SHORT
using CrsGraph [[maybe_unused]] = Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_CRSGRAPHFACTORY_SHORT
using CrsGraphFactory [[maybe_unused]] = Xpetra::CrsGraphFactory<LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_VECTOR_SHORT
using LocalOrdinalVector [[maybe_unused]]  = Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;
using GlobalOrdinalVector [[maybe_unused]] = Xpetra::Vector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;
using LOVector [[maybe_unused]]            = LocalOrdinalVector;
using GOVector [[maybe_unused]]            = GlobalOrdinalVector;
#endif

#ifdef XPETRA_MULTIVECTOR_SHORT
using LocalOrdinalMultiVector [[maybe_unused]]  = Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;
using GlobalOrdinalMultiVector [[maybe_unused]] = Xpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;
using LOMultiVector [[maybe_unused]]            = LocalOrdinalMultiVector;
using GOMultiVector [[maybe_unused]]            = GlobalOrdinalMultiVector;
#endif

#ifdef XPETRA_VECTORFACTORY_SHORT
using LocalOrdinalVectorFactory [[maybe_unused]]  = Xpetra::VectorFactory<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;
using GlobalOrdinalVectorFactory [[maybe_unused]] = Xpetra::VectorFactory<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;
using LOVectorFactory [[maybe_unused]]            = LocalOrdinalVectorFactory;
using GOVectorFactory [[maybe_unused]]            = GlobalOrdinalVectorFactory;
#endif

#ifdef XPETRA_MULTIVECTORFACTORY_SHORT
using LocalOrdinalMultiVectorFactory [[maybe_unused]]  = Xpetra::MultiVectorFactory<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;
using GlobalOrdinalMultiVectorFactory [[maybe_unused]] = Xpetra::MultiVectorFactory<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;
using LOMultiVectorFactory [[maybe_unused]]            = LocalOrdinalMultiVectorFactory;
using GOMultiVectorFactory [[maybe_unused]]            = GlobalOrdinalMultiVectorFactory;
#endif

#ifdef XPETRA_IMPORT_SHORT
using Import [[maybe_unused]] = Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_EXPORT_SHORT
using Export [[maybe_unused]] = Xpetra::Export<LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_IMPORTFACTORY_SHORT
using ImportFactory [[maybe_unused]] = Xpetra::ImportFactory<LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_EXPORTFACTORY_SHORT
using ExportFactory [[maybe_unused]] = Xpetra::ExportFactory<LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_TPETRAMAP_SHORT
using TpetraMap [[maybe_unused]] = Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_TPETRACRSGRAPH_SHORT
using TpetraCrsGraph [[maybe_unused]] = Xpetra::TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_STRIDEDMAP_SHORT
using StridedMap [[maybe_unused]] = Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_STRIDEDMAPFACTORY_SHORT
using StridedMapFactory [[maybe_unused]] = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>;
#endif

// Note: There is no #ifndef/#define/#end in this header file because it can be included more than once (it can be included in methods templated by Scalar, LocalOrdinal, GlobalOrdinal, Node).

// TODO: add namespace {} for shortcut types

// Define convenient shortcut for data types
using LO [[maybe_unused]] = LocalOrdinal;
using GO [[maybe_unused]] = GlobalOrdinal;
using NO [[maybe_unused]] = Node;

// TODO: do the same for Epetra object (problem of namespace)
