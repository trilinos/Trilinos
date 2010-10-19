// Get ride of template parameters

// Declaration of Cthulhu interface classes.
namespace Cthulhu {
  template<class, class, class> class Map;
  template<class, class, class, class, class> class CrsMatrix;
  template<class, class, class, class> class MultiVector;
}

// New definition of types using the types ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps of the current context.
typedef Cthulhu::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
typedef Cthulhu::CrsMatrix<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrix;
typedef Cthulhu::MultiVector<ScalarType, LocalOrdinal, GlobalOrdinal, Node> MultiVector;

// Note: There is no #ifndef/#define/#end in this header file because it can be included more than once (it can be included in methods templated by ScalarType, LocalOrdinal, GlobalOrdinal, Node).

// TODO: add namespace {} for shortcut types
