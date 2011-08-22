#include "Xpetra_ConfigDefs.hpp"

// Get rid of template parameters

// New definition of types using the types Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps of the current context.

// Note: There is no #ifndef/#define/#end in this header file because it can be included more than once (it can be included in methods templated by Scalar, LocalOrdinal, GlobalOrdinal, Node).

#ifdef XPETRA_CRSMATRIX_SHORT
typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrix;
#endif

#ifdef XPETRA_VECTOR_SHORT
typedef Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> Vector;
#endif

#ifdef XPETRA_MULTIVECTOR_SHORT
typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVector;
#endif

#ifdef XPETRA_OPERATOR_SHORT
typedef Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> Operator;
#endif

#ifdef XPETRA_BLOCKEDCRSOPERATOR_SHORT
typedef Xpetra::BlockedCrsOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> BlockedCrsOperator;
#endif

#ifdef XPETRA_CRSOPERATOR_SHORT
typedef Xpetra::CrsOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsOperator;
#endif

#ifdef XPETRA_VECTORFACTORY_SHORT
typedef Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> VectorFactory;
#endif

#ifdef XPETRA_CRSMATRIXFACTORY_SHORT
typedef Xpetra::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrixFactory;
#endif

#ifdef XPETRA_MULTIVECTORFACTORY_SHORT
typedef Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVectorFactory;
#endif

#ifdef XPETRA_OPERATORFACTORY_SHORT
typedef Xpetra::OperatorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> OperatorFactory;
#endif

#ifdef HAVE_XPETRA_TPETRA

#ifdef XPETRA_TPETRACRSMATRIX_SHORT
typedef Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> TpetraCrsMatrix;
#endif

#ifdef XPETRA_TPETRAMULTIVECTOR_SHORT
typedef Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraMultiVector;
#endif

#endif //

// TODO: add namespace {} for shortcut types

// Define convenient shortcut for data types
typedef Scalar    SC;
typedef Teuchos::ScalarTraits<SC> ST;

// TODO: do the same for Epetra object (problem of namespace)
