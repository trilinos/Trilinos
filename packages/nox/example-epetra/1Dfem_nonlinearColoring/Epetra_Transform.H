
#ifndef EPETRA_TRANSFORM_H
#define EPETRA_TRANSFORM_H

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

  virtual std::auto_ptr<newType> operator()( const originalType & old ) = 0;

  virtual std::auto_ptr<originalType> reverse( const newType & old ) { return std::auto_ptr<originalType>(0); }

};

template<class T>
struct SameTypeTransform : public Transform<T,T> { virtual ~SameTypeTransform() {} };

} //namespace Epetra_Transform
  
#endif //EPETRA_TRANSFORM_H
