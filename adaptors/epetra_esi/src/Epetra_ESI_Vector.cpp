#ifndef _Epetra_ESI_Vector_cpp_
#define _Epetra_ESI_Vector_cpp_

#include "Epetra_ESI_Vector.h"

//------------------------------------------------------------------------------
template<class Scalar, class Ordinal>
epetra_esi::Vector<Scalar,Ordinal>::
Vector(epetra_esi::IndexSpace<Ordinal>& indexspace)
    : epetra_esi::Object(), Epetra_Vector(indexspace),
      ispace_(&indexspace),
      currentReads(0), currentReadWrites(0), coefPtr(NULL),
      whichConstructor(0)
{
   int err = addInterface("esi::Object", (void*)((esi::Object*)this) );
   if (err) msgAbort("epetra_esi::Vector ctor addInterface(esi::Object) error");

   err = ispace_->addReference();

   err += addInterface("esi::Vector", (void*)((esi::Vector<Scalar,Ordinal>*)this));
   if (err) msgAbort("epetra_esi::Vector ctor addInterface(esi::Vector) error");

   err = addInterface("epetra_esi::Vector",
                      (void*)((epetra_esi::Vector<Scalar,Ordinal>*)this));
   if (err) msgAbort("epetra_esi::Vector ctor addInterface(epetra_esi::Vector) error");

   err = addInterface("Epetra_Vector", (void*)((Epetra_Vector*)this) );
   if (err) msgAbort("epetra_esi::Vector ctor addInterface(Epetra_Vector) error");
}

//------------------------------------------------------------------------------
template<class Scalar, class Ordinal>
epetra_esi::Vector<Scalar,Ordinal>::
Vector(const Epetra_Vector& invec)
    : epetra_esi::Object(), Epetra_Vector(View, invec, 0),
      ispace_(NULL), currentReads(0), currentReadWrites(0), coefPtr(NULL),
      whichConstructor(1)
{
   int err = addInterface("esi::Object", (void*)((esi::Object*)this) );
   if (err) msgAbort("epetra_esi::Vector ctor addInterface(esi::Object) error");

   err = addInterface("esi::Vector", (void*)((esi::Vector<Scalar,Ordinal>*)this));
   if (err) msgAbort("epetra_esi::Vector ctor addInterface(esi::Vector) error");

   err = addInterface("epetra_esi::Vector",
                      (void*)((epetra_esi::Vector<Scalar,Ordinal>*)this));
   if (err) msgAbort("epetra_esi::Vector ctor addInterface(epetra_esi::Vector) error");

   err = addInterface("Epetra_Vector", (void*)((Epetra_Vector*)this) );
   if (err) msgAbort("epetra_esi::Vector ctor addInterface(Epetra_Vector) error");
}

//------------------------------------------------------------------------------
template<class Scalar, class Ordinal>
epetra_esi::Vector<Scalar,Ordinal>::~Vector()
{
  if (ispace_ != NULL) ispace_->deleteReference();
}

#endif

