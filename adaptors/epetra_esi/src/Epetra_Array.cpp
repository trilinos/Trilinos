#ifndef _Epetra_Array_cpp_
#define _Epetra_Array_cpp_

#include "Epetra_Array.h"

//------------------------------------------------------------------------------
template<class T>
Epetra_Array<T>::Epetra_Array(int length, int incr)
 : data_(NULL),
   length_(0),
   allocatedLength_(0),
   allocIncrement_(incr),
   userMemory_(false)
{
   int err = resize(length);
   if (err) {
      cerr << "ERROR in Epetra_Array construction: couldn't allocate data." << endl;
   }
}

//------------------------------------------------------------------------------
template<class T>
Epetra_Array<T>::Epetra_Array(int length, int allocLength, T* data)
  : data_(data),
    length_(length),
    allocatedLength_(allocLength),
    allocIncrement_(0),
    userMemory_(true)
{
}

//------------------------------------------------------------------------------
template<class T>
Epetra_Array<T>::Epetra_Array(const Epetra_Array<T>& src)
 : data_(NULL),
   length_(0),
   allocatedLength_(0),
   allocIncrement_(src.allocIncrement_),
   userMemory_(false)
{
  int err = resize(src.length_);
  if (err) cerr << "Epetra_Array::(copy constructor): ERROR, resize failed."<<endl;

  for(int i=0; i<src.length_; i++) data_[i] = src.data_[i];
}

#endif
