
#ifndef NOX_EPETRA_DISTTRANSFORM_H
#define NOX_EPETRA_DISTTRANSFORM_H

#include <nox_Epetra_Transform.H>

#include <Epetra_Import.h>
#include <Epetra_Export.h>

namespace Epetra_Transform {

template<class T,class U>
struct DistTransform : public Transform<T,U>
{
  virtual ~DistTransform() {}

  virtual Epetra_Import* importer() { return (Epetra_Import*)(0); }
  virtual Epetra_Export* exporter() { return (Epetra_Export*)(0); }

};

template<class T>
struct SameTypeDistTransform : public DistTransform<T,T> { virtual ~SameTypeDistTransform() {} };

} //namespace Epetra_Transform
  
#endif //NOX_EPETRA_DISTTRANSFORM_H
