// Get ride of template parameters
// You have to include Cthulhu_Classes.hpp before this file.

// New definition of types using the types ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps of the current context.
typedef Cthulhu::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
typedef Cthulhu::CrsMatrix<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrix;
typedef Cthulhu::MultiVector<ScalarType, LocalOrdinal, GlobalOrdinal, Node> MultiVector;

typedef Cthulhu::TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMap;
typedef Cthulhu::TpetraCrsMatrix<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> TpetraCrsMatrix;
typedef Cthulhu::TpetraMultiVector<ScalarType, LocalOrdinal, GlobalOrdinal, Node> TpetraMultiVector;

typedef Cthulhu::MultiVectorFactory<ScalarType, LocalOrdinal, GlobalOrdinal, Node> MultiVectorFactory;

// Note: There is no #ifndef/#define/#end in this header file because it can be included more than once (it can be included in methods templated by ScalarType, LocalOrdinal, GlobalOrdinal, Node).

// TODO: add namespace {} for shortcut types


// Define convenient shortcut for data types
typedef ScalarType    SC;
typedef LocalOrdinal  LO;
typedef GlobalOrdinal GO;
typedef Node          NO;
