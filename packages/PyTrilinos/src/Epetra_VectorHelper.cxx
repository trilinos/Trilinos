#include "Epetra_VectorHelper.h"

#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#include "NumPyArray.h"
#include "NumPyWrapper.h"

Epetra_Vector * Epetra_VectorHelper::new_Epetra_Vector(Epetra_BlockMap & map,
							 PyObject * p_pyObject)
{
  assert (0 != p_pyObject && "PyObject pointer passed to function is NULL.");

  // Use NumPyArray for smart copying of non-contiguous vectors
  // NumPyArrayContiguous numPyArray(p_pyObject);
  NumPyArray numPyArray(p_pyObject);
  
  assert (map.NumMyElements() == numPyArray.getTotalLength() &&
          "BlockMap and NumPy Array objects passed to constructor "
          "are of different lengths"                               );
  
  Epetra_Vector * p_vector = 0;

  // Here we build and load the Epetra_Vector at the same time
  // if the array is contiguous.  Otherwise we let loadViaCopy()
  // take care of it.
  //   if (numPyArray.isContiguous()) {
  if (false) { // Just let loadViaCopy do it all
    p_vector = new Epetra_Vector(Copy                     ,
                                 map                      ,
                                 numPyArray.getDataArray() );
    assert (0 != p_vector && "New returned a NULL pointer");

  }
  else {
    p_vector = new Epetra_Vector(map);
    assert (0 != p_vector && "New returned a NULL pointer");
    
    loadViaCopy (*p_vector, numPyArray);
  }

  return p_vector;
}

void Epetra_VectorHelper::loadViaCopy (Epetra_Vector * p_epetraVector,
					 PyObject * p_pyObject)
{
  assert (0 != p_pyObject     && "PyObject pointer passed to function is NULL."     );
  assert (0 != p_epetraVector && "Epetra_Vector pointer passed to function is NULL.");

  const NumPyArrayContiguous numPyArray(p_pyObject);

  loadViaCopy (*p_epetraVector, numPyArray);
}


// Code common to constructor and other loadViaCopy
void Epetra_VectorHelper::loadViaCopy (Epetra_Vector & epetraVector, const NumPyArrayBase & numPyArray)
{
  assert (numPyArray.getTotalLength() == epetraVector.MyLength() 
          && "Source NumPy array and target Epetra Vector differ in length.");
  
  double * p_epetraView = NULL;
  
  epetraVector.ExtractView(&p_epetraView);
  
  assert (0 != p_epetraView && "p_epetraVector->ExtractView() failed to set pointer.");

  NumPyWrapper::CopyToExternalArray(NumPyWrapper::LOAD         ,   // Copy Action
                                    numPyArray.getArrayObject(),   // PyArrayObject
                                    numPyArray.getDataArray()  ,   // Source Array
                                    p_epetraView                ); // Target Array
}

void Epetra_VectorHelper::unloadViaCopy (const Epetra_Vector * p_epetraVector, PyObject * p_pyObject)
{
  assert (0 != p_pyObject     && "PyObject pointer passed to function is NULL."     );
  assert (0 != p_epetraVector && "Epetra_Vector pointer passed to function is NULL.");

  NumPyArray numPyArray(p_pyObject);

  unloadViaCopy (*p_epetraVector, numPyArray);
}


// Code common to constructor and other loadViaCopy
void Epetra_VectorHelper::unloadViaCopy (const Epetra_Vector & epetraVector, NumPyArrayBase & numPyArray)
{
  assert (numPyArray.getTotalLength() == epetraVector.MyLength() 
          && "Source NumPy array and target Epetra Vector differ in length.");
  
  double * p_epetraView = NULL;
  
  ((Epetra_Vector &)epetraVector).ExtractView(&p_epetraView);
  
  assert (0 != p_epetraView && "p_epetraVector->ExtractView() failed to set pointer.");

  NumPyWrapper::CopyToExternalArray(NumPyWrapper::UNLOAD       ,   // Copy Action
                                    numPyArray.getArrayObject(),   // PyArrayObject
                                    p_epetraView               ,   // Source Array
                                    numPyArray.getDataArray()   ); // Target Array
}
