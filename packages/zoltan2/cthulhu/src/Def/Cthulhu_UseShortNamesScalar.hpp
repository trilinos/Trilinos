#include "Cthulhu_ConfigDefs.hpp"

// Get rid of template parameters

// New definition of types using the types Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps of the current context.

// Note: There is no #ifndef/#define/#end in this header file because it can be included more than once (it can be included in methods templated by Scalar, LocalOrdinal, GlobalOrdinal, Node).

#ifdef CTHULHU_CRSMATRIX_SHORT
typedef Cthulhu::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrix;
#endif

#ifdef CTHULHU_VECTOR_SHORT
typedef Cthulhu::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> Vector;
#endif

#ifdef CTHULHU_MULTIVECTOR_SHORT
typedef Cthulhu::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVector;
#endif

#ifdef CTHULHU_OPERATOR_SHORT
typedef Cthulhu::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> Operator;
#endif

#ifdef CTHULHU_CRSOPERATOR_SHORT
typedef Cthulhu::CrsOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsOperator;
#endif

#ifdef CTHULHU_VECTORFACTORY_SHORT
typedef Cthulhu::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> VectorFactory;
#endif

#ifdef CTHULHU_CRSMATRIXFACTORY_SHORT
typedef Cthulhu::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrixFactory;
#endif

#ifdef CTHULHU_MULTIVECTORFACTORY_SHORT
typedef Cthulhu::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVectorFactory;
#endif

#ifdef CTHULHU_OPERATORFACTORY_SHORT
typedef Cthulhu::OperatorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> OperatorFactory;
#endif

#ifdef HAVE_CTHULHU_TPETRA

#ifdef CTHULHU_TPETRACRSMATRIX_SHORT
typedef Cthulhu::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> TpetraCrsMatrix;
#endif

#ifdef CTHULHU_TPETRAMULTIVECTOR_SHORT
typedef Cthulhu::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraMultiVector;
#endif

#endif //

// TODO: add namespace {} for shortcut types

// Define convenient shortcut for data types
typedef Scalar    SC;
typedef Teuchos::ScalarTraits<SC> ST;

// TODO: do the same for Epetra object (problem of namespace)
