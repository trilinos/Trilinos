#include "Epetra_NumPyVector.h"

#define DEBUG 0
#if DEBUG
#include <iostream>
#include <string>
using namespace std;
#endif

// Static variables
bool                      Epetra_NumPyVector::initFlag    = false              ;
const Epetra_SerialComm   Epetra_NumPyVector::defaultComm = Epetra_SerialComm();
PyArrayObject           * Epetra_NumPyVector::tmp_array   = NULL               ;
Epetra_Map              * Epetra_NumPyVector::tmp_map     = NULL               ;

// Static helper functions
// =============================================================================
void Epetra_NumPyVector::initialize() {
  if ( !initFlag ) {
#if DEBUG
    cout << "Calling import_array() to initialize Numeric" << endl;
#endif
    import_array();
    initFlag = true;
  }
}

// =============================================================================
Epetra_Map & Epetra_NumPyVector::getEpetraMap(PyObject * pyObject)
{
#if DEBUG
  cout << "Inside Epetra_NumPyVector::getEpetraMap" << endl;
#endif
  getSourceData(pyObject);   // This creates the tmp_array
#if DEBUG
  cout << "  Getting PyArray size." << endl;
#endif
  const int totalLength = PyArray_Size((PyObject *)tmp_array);
  assert(NULL == tmp_map);
#if DEBUG
  cout << "  Constructing Epetra_Map" << endl;
#endif
  tmp_map = new Epetra_Map(totalLength,0,defaultComm);
  return *tmp_map;
}

// =============================================================================
double * Epetra_NumPyVector::getSourceData(PyObject * pyObject)
{
  initialize();
  assert(NULL == tmp_array);
  tmp_array = (PyArrayObject *) PyArray_ContiguousFromObject(pyObject, 'd', 0, 0);
  return (double *)tmp_array->data;
}
// =============================================================================

// Constructors
// =============================================================================
Epetra_NumPyVector::Epetra_NumPyVector(Epetra_BlockMap & blockMap):
  Epetra_Vector(blockMap, true)
{
#if DEBUG
  cout << "Inside Epetra_NumPyVector(Epetra_BlockMap &) constructor" << endl;
#endif
  // Create the array object
  initialize();
  int dims[ ] = { blockMap.NumMyElements() };
  double *v = NULL;
#if DEBUG
  cout << "  Calling ExtractView to get pointer to data" << endl;
#endif
  ExtractView(&v);
#if DEBUG
  cout << "  Creating PyArrayObject" << endl;
#endif
  array = (PyArrayObject *) PyArray_FromDimsAndData(1,dims,PyArray_DOUBLE,
						    (char *)v);

  // Copy the Epetra_BlockMap
  map = new Epetra_BlockMap(blockMap);
}

// =============================================================================
Epetra_NumPyVector::Epetra_NumPyVector(Epetra_BlockMap & blockMap,
                                       PyObject        * pyObject):
  Epetra_Vector(View, blockMap, getSourceData(pyObject))
{
  // Get the pointer to the array from static variable and clear
  assert(NULL != tmp_array);
  array     = tmp_array;
  tmp_array = NULL;

  // Copy the Epetra_BlockMap
  map = new Epetra_BlockMap(blockMap);
}

// =============================================================================
Epetra_NumPyVector::Epetra_NumPyVector(PyObject * pyObject):
  Epetra_Vector(View, getEpetraMap(pyObject), (double*)tmp_array->data) 
{
  // Store the pointer to the Epetra_Map
  assert(NULL != tmp_map);
  map     = tmp_map;
  tmp_map = NULL;

  // Store the pointer to the PyArrayObject
  assert(NULL != tmp_array);
  array     = tmp_array;
  tmp_array = NULL;
}

// =============================================================================
// Destrtuctor
Epetra_NumPyVector::~Epetra_NumPyVector()
{
  Py_XDECREF(array);
  delete map;
}

// =============================================================================
PyObject * Epetra_NumPyVector::getArray()
{
  return PyArray_Return(array);
}

