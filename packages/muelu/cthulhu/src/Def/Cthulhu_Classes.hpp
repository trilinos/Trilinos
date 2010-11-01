#ifndef CTHULHU_CLASSES_DECL_HPP
#define CTHULHU_CLASSES_DECL_HPP

// Declaration of Cthulhu classes.
namespace Cthulhu {
  template<class, class, class> class Map;
  template<class, class, class, class, class> class CrsMatrix;
  template<class, class, class, class> class Vector;
  template<class, class, class, class> class MultiVector;

  template<class, class, class> class TpetraMap;
  template<class, class, class, class, class> class TpetraCrsMatrix;
  template<class, class, class, class> class MultiVector;
  template<class, class, class, class> class TpetraMultiVector;

  template<class, class, class, class, class> class CrsMatrixFactory;
  template<class, class, class, class> class VectorFactory;
  template<class, class, class, class> class MultiVectorFactory;

  class EpetraMap;
  //  class EpetraVector;
  class EpetraMultiVector;
  class EpetraCrsMatrix;
}

#endif
