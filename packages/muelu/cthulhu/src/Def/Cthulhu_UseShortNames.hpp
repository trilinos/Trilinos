#include "Cthulhu_ConfigDefs.hpp"

// Get ride of template parameters
// You have to include Cthulhu_Classes.hpp before this file.

// New definition of types using the types ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps of the current context.
#ifdef CTHULHU_MAP_HPP
typedef Cthulhu::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
#endif

#ifdef CTHULHU_CRSMATRIX_HPP
typedef Cthulhu::CrsMatrix<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrix;
#endif

#ifdef CTHULHU_VECTOR_HPP
typedef Cthulhu::Vector<ScalarType, LocalOrdinal, GlobalOrdinal, Node> Vector;
#endif

#ifdef CTHULHU_MULTIVECTOR_HPP
typedef Cthulhu::MultiVector<ScalarType, LocalOrdinal, GlobalOrdinal, Node> MultiVector;
#endif

#ifdef CTHULHU_OPERATOR_HPP
typedef Cthulhu::Operator<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> Operator;
#endif

#ifdef CTHULHU_CRSOPERATOR_HPP
typedef Cthulhu::CrsOperator<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsOperator;
#endif

#ifdef CTHULHU_VECTORFACTORY_HPP
typedef Cthulhu::VectorFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node> VectorFactory;
#endif

#ifdef CTHULHU_CRSMATRIXFACTORY_HPP
typedef Cthulhu::CrsMatrixFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrixFactory;
#endif
#ifdef CTHULHU_MULTIVECTORFACTORY_HPP
typedef Cthulhu::MultiVectorFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node> MultiVectorFactory;
#endif

#ifdef CTHULHU_OPERATORFACTORY_HPP
typedef Cthulhu::OperatorFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> OperatorFactory;
#endif

#ifdef HAVE_CTHULHU_TPETRA

#ifdef CTHULHU_TPETRAMAP_HPP
typedef Cthulhu::TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMap;
#endif

#ifdef CTHULHU_TPETRACRSMATRIX_HPP
typedef Cthulhu::TpetraCrsMatrix<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> TpetraCrsMatrix;
#endif
#ifdef CTHULHU_TPETRAMULTIVECTOR_HPP
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
