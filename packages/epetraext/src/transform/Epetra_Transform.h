
#ifndef EPETRA_TRANSFORM_H
#define EPETRA_TRANSFORM_H

#ifdef HAVE_CONFIG_H
#include <EpetraExt_config.h>
#endif

namespace EpetraExt {

template<typename T, typename U>
struct Transform
{
  typedef T  OriginalType;
  typedef T* OriginalTypePtr;
  typedef T& OriginalTypeRef;

  typedef U  NewType;
  typedef U* NewTypePtr;
  typedef T& NewTypeRef;

  virtual ~Transform() {}

  virtual NewTypePtr operator()( OriginalTypeRef old ) = 0;

  virtual bool fwd() = 0;
  virtual bool rvs() = 0;
};

template<typename T, typename U>
struct StructuralTransform : public Transform<T,U>
{
  bool fwd() { return true; }
  bool rvs() { return true; }

  virtual ~StructuralTransform() {}
};

template<class T>
struct SameTypeTransform : public Transform<T,T>
{
  virtual ~SameTypeTransform() {}
};

template<typename T>
struct StructuralSameTypeTransform : public SameTypeTransform<T>
{
  bool fwd() { return true; }
  bool rvs() { return true; }

  virtual ~StructuralSameTypeTransform() {}
};

template<class T>
struct InPlaceTransform : public SameTypeTransform<T>
{
  NewTypePtr operator()( OriginalTypeRef old )
  { return NewTypePtr(0); }

  virtual ~InPlaceTransform() {}
};

template<class T>
struct ViewTransform : public SameTypeTransform<T>
{
  bool fwd() { return true; }
  bool rvs() { return true; }

  virtual ~ViewTransform() {}
};

} //namespace EpetraExt
  
#endif //EPETRA_TRANSFORM_H
