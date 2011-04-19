#include "Cthulhu_ConfigDefs.hpp"

// Get rid of template parameters

// New definition of types using the types ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps of the current context.
#ifdef CTHULHU_MAP_SHORT
typedef Cthulhu::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
#endif

#ifdef CTHULHU_MAPFACTORY_SHORT
typedef Cthulhu::MapFactory<LocalOrdinal, GlobalOrdinal, Node> MapFactory;
#endif

#ifdef CTHULHU_CRSGRAPH_SHORT
typedef Cthulhu::CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsGraph;
#endif

#ifdef CTHULHU_CRSMATRIX_SHORT
typedef Cthulhu::CrsMatrix<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrix;
#endif

#ifdef CTHULHU_VECTOR_SHORT
typedef Cthulhu::Vector<ScalarType, LocalOrdinal, GlobalOrdinal, Node> Vector;
typedef Cthulhu::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> LocalOrdinalVector;
typedef Cthulhu::Vector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> GlobalOrdinalVector;
typedef LocalOrdinalVector  LOVector;
typedef GlobalOrdinalVector GOVector;
#endif

#ifdef CTHULHU_MULTIVECTOR_SHORT
typedef Cthulhu::MultiVector<ScalarType, LocalOrdinal, GlobalOrdinal, Node> MultiVector;
typedef Cthulhu::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> LocalOrdinalMultiVector;
typedef Cthulhu::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> GlobalOrdinalMultiVector;
typedef LocalOrdinalMultiVector  LOMultiVector;
typedef GlobalOrdinalMultiVector GOMultiVector;
#endif

#ifdef CTHULHU_OPERATOR_SHORT
typedef Cthulhu::Operator<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> Operator;
#endif

#ifdef CTHULHU_CRSOPERATOR_SHORT
typedef Cthulhu::CrsOperator<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsOperator;
#endif

#ifdef CTHULHU_VECTORFACTORY_SHORT
typedef Cthulhu::VectorFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node> VectorFactory;
typedef Cthulhu::VectorFactory<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> LocalOrdinalVectorFactory;
typedef Cthulhu::VectorFactory<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> GlobalOrdinalVectorFactory;
typedef LocalOrdinalVectorFactory  LOVectorFactory;
typedef GlobalOrdinalVectorFactory GOVectorFactory;
#endif

#ifdef CTHULHU_CRSMATRIXFACTORY_SHORT
typedef Cthulhu::CrsMatrixFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrixFactory;
#endif

#ifdef CTHULHU_MULTIVECTORFACTORY_SHORT
typedef Cthulhu::MultiVectorFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node> MultiVectorFactory;
typedef Cthulhu::MultiVectorFactory<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> LocalOrdinalMultiVectorFactory;
typedef Cthulhu::MultiVectorFactory<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> GlobalOrdinalMultiVectorFactory;
typedef LocalOrdinalMultiVectorFactory  LOMultiVectorFactory;
typedef GlobalOrdinalMultiVectorFactory GOMultiVectorFactory;
#endif

#ifdef CTHULHU_OPERATORFACTORY_SHORT
typedef Cthulhu::OperatorFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> OperatorFactory;
#endif

#ifdef CTHULHU_IMPORT_SHORT
typedef Cthulhu::Import<LocalOrdinal, GlobalOrdinal, Node> Import;
#endif

#ifdef CTHULHU_IMPORTFACTORY_SHORT
typedef Cthulhu::ImportFactory<LocalOrdinal, GlobalOrdinal, Node> ImportFactory;
#endif

#ifdef HAVE_CTHULHU_TPETRA

#ifdef CTHULHU_TPETRAMAP_SHORT
typedef Cthulhu::TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMap;
#endif

#ifdef CTHULHU_TPETRACRSGRAPH_SHORT
typedef Cthulhu::TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> TpetraCrsGraph;
#endif

#ifdef CTHULHU_TPETRACRSMATRIX_SHORT
typedef Cthulhu::TpetraCrsMatrix<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> TpetraCrsMatrix;
#endif

#ifdef CTHULHU_TPETRAMULTIVECTOR_SHORT
typedef Cthulhu::TpetraMultiVector<ScalarType, LocalOrdinal, GlobalOrdinal, Node> TpetraMultiVector;
#endif

#endif //

// Note: There is no #ifndef/#define/#end in this header file because it can be included more than once (it can be included in methods templated by ScalarType, LocalOrdinal, GlobalOrdinal, Node).

// TODO: add namespace {} for shortcut types

// Define convenient shortcut for data types
typedef ScalarType    SC;
typedef LocalOrdinal  LO;
typedef GlobalOrdinal GO;
typedef Node          NO;
typedef LocalMatOps   LMO;

// TODO: do the same for Epetra object (problem of namespace)
