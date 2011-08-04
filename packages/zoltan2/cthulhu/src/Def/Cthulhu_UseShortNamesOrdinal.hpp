#include "Cthulhu_ConfigDefs.hpp"

// Get rid of template parameters

// New definition of types using the types LocalOrdinal, GlobalOrdinal, Node, LocalMatOps of the current context.
#ifdef CTHULHU_MAP_SHORT
typedef Cthulhu::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
#endif

#ifdef CTHULHU_MAPFACTORY_SHORT
typedef Cthulhu::MapFactory<LocalOrdinal, GlobalOrdinal, Node> MapFactory;
#endif

#ifdef CTHULHU_CRSGRAPH_SHORT
typedef Cthulhu::CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsGraph;
#endif

#ifdef CTHULHU_VECTOR_SHORT
typedef Cthulhu::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> LocalOrdinalVector;
typedef Cthulhu::Vector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> GlobalOrdinalVector;
typedef LocalOrdinalVector  LOVector;
typedef GlobalOrdinalVector GOVector;
#endif

#ifdef CTHULHU_MULTIVECTOR_SHORT
typedef Cthulhu::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> LocalOrdinalMultiVector;
typedef Cthulhu::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> GlobalOrdinalMultiVector;
typedef LocalOrdinalMultiVector  LOMultiVector;
typedef GlobalOrdinalMultiVector GOMultiVector;
#endif

#ifdef CTHULHU_VECTORFACTORY_SHORT
typedef Cthulhu::VectorFactory<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> LocalOrdinalVectorFactory;
typedef Cthulhu::VectorFactory<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> GlobalOrdinalVectorFactory;
typedef LocalOrdinalVectorFactory  LOVectorFactory;
typedef GlobalOrdinalVectorFactory GOVectorFactory;
#endif

#ifdef CTHULHU_MULTIVECTORFACTORY_SHORT
typedef Cthulhu::MultiVectorFactory<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> LocalOrdinalMultiVectorFactory;
typedef Cthulhu::MultiVectorFactory<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> GlobalOrdinalMultiVectorFactory;
typedef LocalOrdinalMultiVectorFactory  LOMultiVectorFactory;
typedef GlobalOrdinalMultiVectorFactory GOMultiVectorFactory;
#endif

#ifdef CTHULHU_IMPORT_SHORT
typedef Cthulhu::Import<LocalOrdinal, GlobalOrdinal, Node> Import;
#endif

#ifdef CTHULHU_EXPORT_SHORT
typedef Cthulhu::Export<LocalOrdinal, GlobalOrdinal, Node> Export;
#endif

#ifdef CTHULHU_IMPORTFACTORY_SHORT
typedef Cthulhu::ImportFactory<LocalOrdinal, GlobalOrdinal, Node> ImportFactory;
#endif

#ifdef CTHULHU_EXPORTFACTORY_SHORT
typedef Cthulhu::ExportFactory<LocalOrdinal, GlobalOrdinal, Node> ExportFactory;
#endif

#ifdef HAVE_CTHULHU_TPETRA

#ifdef CTHULHU_TPETRAMAP_SHORT
typedef Cthulhu::TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMap;
#endif

#ifdef CTHULHU_TPETRACRSGRAPH_SHORT
typedef Cthulhu::TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> TpetraCrsGraph;
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
