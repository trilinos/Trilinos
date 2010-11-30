#ifndef CTHULHU_CLASSES_HPP
#define CTHULHU_CLASSES_HPP

#include "Cthulhu_ConfigDefs.hpp"

// Declaration of Cthulhu classes.
namespace Cthulhu {
  template<class, class, class> class Map;
  template<class, class, class, class, class> class CrsMatrix;
  template<class, class, class, class> class Vector;
  template<class, class, class, class> class MultiVector;

  template<class, class, class, class, class> class Operator;
  template<class, class, class, class, class> class CrsOperator;

#ifdef HAVE_CTHULHU_TPETRA
  template<class, class, class> class TpetraMap;
  template<class, class, class, class, class> class TpetraCrsMatrix;
  template<class, class, class, class> class MultiVector;
  template<class, class, class, class> class TpetraMultiVector;
#endif

#ifdef HAVE_CTHULHU_EPETRA
  class EpetraMap;
  class EpetraVector;
  class EpetraMultiVector;
  class EpetraCrsMatrix;
#endif

  template<class, class, class, class, class> class CrsMatrixFactory;
  template<class, class, class, class> class VectorFactory;
  template<class, class, class, class> class MultiVectorFactory;

  template<class, class, class, class, class> class OperatorFactory;
}

#endif
