#ifndef EPETRA_NUMPYVECTOR_H
#define EPETRA_NUMPYVECTOR_H

#include <Python.h>
#include <Numeric/arrayobject.h>

#include "Epetra_SerialComm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"

class Epetra_NumPyVector : public Epetra_Vector {

public:

  Epetra_NumPyVector(Epetra_BlockMap &);
  Epetra_NumPyVector(Epetra_BlockMap &,
		     PyObject        *);
  Epetra_NumPyVector(PyObject        *); 

  ~Epetra_NumPyVector();

  PyObject * getArray();

private:

  // Private methods thus not callable
  Epetra_NumPyVector();
  Epetra_NumPyVector(const Epetra_NumPyVector &);

  // Static helper functions
  static void         initialize();
  static Epetra_Map & getEpetraMap( PyObject *);
  static double     * getSourceData(PyObject *);

  // Static private data
  static       bool              initFlag;
  static const Epetra_SerialComm defaultComm;
  static       PyArrayObject *   tmp_array;
  static       Epetra_Map    *   tmp_map;

  // Private data
  Epetra_BlockMap * map;
  PyArrayObject   * array;

};

#endif
