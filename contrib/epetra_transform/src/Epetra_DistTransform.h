
#ifndef EPETRA_DISTTRANSFORM_H
#define EPETRA_DISTTRANSFORM_H

#include <Epetra_Transform.h>

#include <Epetra_Import.h>
#include <Epetra_Export.h>

namespace Epetra_Transform {

template<class T,class U>
struct DistTransform : public Transform<T,U>
{
  virtual ~DistTransform() {}

  virtual std::auto_ptr<Epetra_Import> importer() { return std::auto_ptr<Epetra_Import>(0); }
  virtual std::auto_ptr<Epetra_Export> exporter() { return std::auto_ptr<Epetra_Export>(0); }

};

template<class T>
struct SameTypeDistTransform : public DistTransform<T,T> { virtual ~SameTypeDistTransform() {} };

} //namespace Epetra_Transform
  
#endif //EPETRA_DISTTRANSFORM_H
