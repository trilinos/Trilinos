#ifndef EPETRA_VECTORHELPER_H
#define EPETRA_VECTORHELPER_H

#include "Python.h"

// Forward declarations
class Epetra_Vector  ;
class NumPyArrayBase ;
class Epetra_BlockMap;

// Class to hold static functions used by the swig Epetra interface
class Epetra_VectorHelper
{
public:
  static Epetra_Vector * new_Epetra_Vector(Epetra_BlockMap & map ,
					   PyObject        * p_pyObject);

  static void            loadViaCopy      (Epetra_Vector * p_epetraVector,
					   PyObject      * p_pyObject);
  static void            loadViaCopy      (Epetra_Vector & epetraVector,
					   const NumPyArrayBase & numPyArray);
  static void            unloadViaCopy    (const Epetra_Vector * p_epetraVector,
					   PyObject            * p_pyObject);
  static void            unloadViaCopy    (const Epetra_Vector & epetraVector,
					   NumPyArrayBase      & numPyArray);

protected:
  // Class should never be instantiated, it just holds
  // static functions
  Epetra_VectorHelper(PyObject * p_pyObject);
  ~Epetra_VectorHelper();

private:
  // Private and not implemented
  Epetra_VectorHelper();
  Epetra_VectorHelper(const Epetra_VectorHelper &);
  Epetra_VectorHelper & operator=(const Epetra_VectorHelper &);
};


#endif // EPETRA_VECTORHELPER_H
