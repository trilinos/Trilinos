#ifndef PYEPETRA_UTILS_H
#define PYEPETRA_UTILS_H

#include "Python.h"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include <vector>

PyObject * Epetra_RowMatrix_GetEntries(const Epetra_RowMatrix& Matrix,
                                       int GlobalRow) 
{
  if (!Matrix.Filled())
  {
    cout << "Matrix not FillComplete()'d" << endl;
    Py_INCREF(Py_None);
    return Py_None;
  }
  int ierr;
  PyObject* res;
  PyObject* PyIndices;
  PyObject* PyValues;
  int NumEntries, Length;
  Length = Matrix.MaxNumEntries();
  std::vector<int>    Indices(Length);
  std::vector<double> Values(Length);
  int MyRow = Matrix.RowMatrixRowMap().LID(GlobalRow);
  ierr = Matrix.ExtractMyRowCopy(MyRow, Length, NumEntries, &Values[0],
                                 &Indices[0]);
  if (ierr < 0)
  {
    Py_INCREF(Py_None);
    return Py_None;
  }

  PyIndices = PyList_New(NumEntries);
  PyValues  = PyList_New(NumEntries);

  // return global indices
  for (int i = 0 ; i < NumEntries ; ++i)
  {
    int GID = Matrix.RowMatrixColMap().GID(Indices[i]);
    PyList_SetItem(PyIndices, i, PyInt_FromLong(GID));
    PyList_SetItem(PyValues,  i, PyFloat_FromDouble(Values[i]));
  }
  res = PyTuple_New(2);
  PyTuple_SetItem(res, 0, PyIndices);
  PyTuple_SetItem(res, 1, PyValues);
  return(res);
}

PyObject* Epetra_RowMatrix_GetEntry(Epetra_RowMatrix& Matrix, 
                                    int GlobalRow, int GlobalCol)
{
  int NumEntries, Length, ierr;
  double val = 0.0;
  Length = Matrix.MaxNumEntries();
  std::vector<int>    Indices(Length);
  std::vector<double> Values(Length);
  int MyRow = Matrix.RowMatrixRowMap().LID(GlobalRow);
  int MyCol = Matrix.RowMatrixColMap().LID(GlobalCol);

  ierr = Matrix.ExtractMyRowCopy(MyRow, Length, NumEntries, &Values[0], 
                                 &Indices[0]);
  if (ierr) NumEntries = 0;
  for (int i = 0 ; i < Length ; ++i)
  {
    if (Indices[i] == MyCol)
    {
      val = Values[i];
      break;
    }
  }
  return(PyFloat_FromDouble(val));
}

#endif

