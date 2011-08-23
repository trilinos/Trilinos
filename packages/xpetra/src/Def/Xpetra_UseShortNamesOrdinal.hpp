// Get rid of template parameters

// New definition of types using the types LocalOrdinal, GlobalOrdinal, Node, LocalMatOps of the current context.
#ifdef XPETRA_MAP_SHORT
typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
#endif

#ifdef XPETRA_MAPFACTORY_SHORT
typedef Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node> MapFactory;
#endif

#ifdef XPETRA_CRSGRAPH_SHORT
typedef Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsGraph;
#endif

#ifdef XPETRA_VECTOR_SHORT
typedef Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> LocalOrdinalVector;
typedef Xpetra::Vector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> GlobalOrdinalVector;
typedef LocalOrdinalVector  LOVector;
typedef GlobalOrdinalVector GOVector;
#endif

#ifdef XPETRA_MULTIVECTOR_SHORT
typedef Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> LocalOrdinalMultiVector;
typedef Xpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> GlobalOrdinalMultiVector;
typedef LocalOrdinalMultiVector  LOMultiVector;
typedef GlobalOrdinalMultiVector GOMultiVector;
#endif

#ifdef XPETRA_VECTORFACTORY_SHORT
typedef Xpetra::VectorFactory<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> LocalOrdinalVectorFactory;
typedef Xpetra::VectorFactory<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> GlobalOrdinalVectorFactory;
typedef LocalOrdinalVectorFactory  LOVectorFactory;
typedef GlobalOrdinalVectorFactory GOVectorFactory;
#endif

#ifdef XPETRA_MULTIVECTORFACTORY_SHORT
typedef Xpetra::MultiVectorFactory<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> LocalOrdinalMultiVectorFactory;
typedef Xpetra::MultiVectorFactory<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> GlobalOrdinalMultiVectorFactory;
typedef LocalOrdinalMultiVectorFactory  LOMultiVectorFactory;
typedef GlobalOrdinalMultiVectorFactory GOMultiVectorFactory;
#endif

#ifdef XPETRA_IMPORT_SHORT
typedef Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> Import;
#endif

#ifdef XPETRA_EXPORT_SHORT
typedef Xpetra::Export<LocalOrdinal, GlobalOrdinal, Node> Export;
#endif

#ifdef XPETRA_IMPORTFACTORY_SHORT
typedef Xpetra::ImportFactory<LocalOrdinal, GlobalOrdinal, Node> ImportFactory;
#endif

#ifdef XPETRA_EXPORTFACTORY_SHORT
typedef Xpetra::ExportFactory<LocalOrdinal, GlobalOrdinal, Node> ExportFactory;
#endif

#ifdef HAVE_XPETRA_TPETRA

#ifdef XPETRA_TPETRAMAP_SHORT
typedef Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMap;
#endif

#ifdef XPETRA_TPETRACRSGRAPH_SHORT
typedef Xpetra::TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> TpetraCrsGraph;
#endif

#endif

// Note: There is no #ifndef/#define/#end in this header file because it can be included more than once (it can be included in methods templated by Scalar, LocalOrdinal, GlobalOrdinal, Node).

// TODO: add namespace {} for shortcut types

// Define convenient shortcut for data types
typedef LocalOrdinal  LO;
typedef GlobalOrdinal GO;
typedef Node          NO;
typedef LocalMatOps   LMO;

// TODO: do the same for Epetra object (problem of namespace)
