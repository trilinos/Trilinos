
#ifndef NOX_EPETRA_TRANSFORM_H
#define NOX_EPETRA_TRANSFORM_H

#include <memory>

namespace Epetra_Transform {

template<class T,class U>
class Transform
{

 protected:

  typedef T originalType;
  typedef U newType;

 public:

  virtual ~Transform() {}

  virtual newType* operator()( const originalType & old ) = 0;

  virtual originalType* reverse( const newType & old ) { return (originalType*)(0); }

};

template<class T>
struct SameTypeTransform : public Transform<T,T> { virtual ~SameTypeTransform() {} };

} //namespace Epetra_Transform
  
#endif //NOX_EPETRA_TRANSFORM_H
