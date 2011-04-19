#include "Cthulhu_ConfigDefs.hpp"

// Get rid of template parameters

// New definition of types using the types ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps of the current context.

#ifdef CTHULHU_CRSMATRIX_SHORT
typedef Cthulhu::CrsMatrix<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrix;
#endif

#ifdef CTHULHU_VECTOR_SHORT
typedef Cthulhu::Vector<ScalarType, LocalOrdinal, GlobalOrdinal, Node> Vector;
#endif

#ifdef CTHULHU_MULTIVECTOR_SHORT
typedef Cthulhu::MultiVector<ScalarType, LocalOrdinal, GlobalOrdinal, Node> MultiVector;
#endif

#ifdef CTHULHU_OPERATOR_SHORT
typedef Cthulhu::Operator<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> Operator;
#endif

#ifdef CTHULHU_CRSOPERATOR_SHORT
typedef Cthulhu::CrsOperator<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsOperator;
#endif

#ifdef CTHULHU_VECTORFACTORY_SHORT
typedef Cthulhu::VectorFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node> VectorFactory;
#endif

#ifdef CTHULHU_CRSMATRIXFACTORY_SHORT
typedef Cthulhu::CrsMatrixFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrixFactory;
#endif

#ifdef CTHULHU_MULTIVECTORFACTORY_SHORT
typedef Cthulhu::MultiVectorFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node> MultiVectorFactory;
#endif

#ifdef CTHULHU_OPERATORFACTORY_SHORT
typedef Cthulhu::OperatorFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> OperatorFactory;
#endif

#ifdef HAVE_CTHULHU_TPETRA

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

// TODO: do the same for Epetra object (problem of namespace)
