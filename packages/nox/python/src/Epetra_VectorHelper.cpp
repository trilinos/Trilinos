// @HEADER
# ************************************************************************
# 
#            NOX: An Object-Oriented Nonlinear Solver Package
#                 Copyright (2002) Sandia Corporation
# 
#            LOCA: Library of Continuation Algorithms Package
#                 Copyright (2005) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# 
# Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
# Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
# ************************************************************************
#  CVS Information
#  $Source$
#  $Author$
#  $Date$
#  $Revision$
# ************************************************************************
// @HEADER

#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"

#include "Epetra_VectorHelper.h"

Epetra_Vector * Epetra_VectorHelper::new_Epetra_Vector(Epetra_BlockMap & map,
						       PyObject        * p_pyObject)
{
//   assert (0 != p_pyObject && "PyObject pointer passed to function is NULL.");

//   // Use NumPyArray for smart copying of non-contiguous vectors
//   // NumPyArrayContiguous numPyArray(p_pyObject);
//   NumPyArray numPyArray(p_pyObject);
  
//   assert (map.NumMyElements() == numPyArray.getTotalLength() &&
//           "BlockMap and NumPy Array objects passed to constructor "
//           "are of different lengths"                               );
  
//   Epetra_Vector * p_vector = 0;

//   // Here we build and load the Epetra_Vector at the same time
//   // if the array is contiguous.  Otherwise we let loadViaCopy()
//   // take care of it.
//   //   if (numPyArray.isContiguous()) {
//   if (false) { // Just let loadViaCopy do it all
//     p_vector = new Epetra_Vector(Copy                     ,
//                                  map                      ,
//                                  numPyArray.getDataArray() );
//     assert (0 != p_vector && "New returned a NULL pointer");

//   }
//   else {
//     p_vector = new Epetra_Vector(map);
//     assert (0 != p_vector && "New returned a NULL pointer");
    
//     loadViaCopy (*p_vector, numPyArray);
//   }

//   return p_vector;
  return NULL;
}

void Epetra_VectorHelper::loadViaCopy (Epetra_Vector * p_epetraVector,
				       PyObject      * p_pyObject     )
{
//   assert (0 != p_pyObject     && "PyObject pointer passed to function is NULL."     );
//   assert (0 != p_epetraVector && "Epetra_Vector pointer passed to function is NULL.");

//   const NumPyArrayContiguous numPyArray(p_pyObject);

//   loadViaCopy (*p_epetraVector, numPyArray);
}


// Code common to constructor and other loadViaCopy
// void Epetra_VectorHelper::loadViaCopy (      Epetra_Vector  & epetraVector,
// 				       const NumPyArrayBase & numPyArray   )
// {
//   assert (numPyArray.getTotalLength() == epetraVector.MyLength() 
//           && "Source NumPy array and target Epetra Vector differ in length.");
  
//   double * p_epetraView = NULL;
  
//   epetraVector.ExtractView(&p_epetraView);
  
//   assert (0 != p_epetraView && "p_epetraVector->ExtractView() failed to set pointer.");

//   NumPyWrapper::CopyToExternalArray(NumPyWrapper::LOAD         ,   // Copy Action
//                                     numPyArray.getArrayObject(),   // PyArrayObject
//                                     numPyArray.getDataArray()  ,   // Source Array
//                                     p_epetraView                ); // Target Array
// }

void Epetra_VectorHelper::unloadViaCopy (const Epetra_Vector * p_epetraVector, PyObject * p_pyObject)
{
  assert (0 != p_pyObject     && "PyObject pointer passed to function is NULL."     );
  assert (0 != p_epetraVector && "Epetra_Vector pointer passed to function is NULL.");

//   NumPyArray numPyArray(p_pyObject);

//   unloadViaCopy (*p_epetraVector, numPyArray);
}


// Code common to constructor and other loadViaCopy
// void Epetra_VectorHelper::unloadViaCopy (const Epetra_Vector  & epetraVector,
// 					       NumPyArrayBase & numPyArray   )
// {
//   assert (numPyArray.getTotalLength() == epetraVector.MyLength() 
//           && "Source NumPy array and target Epetra Vector differ in length.");
  
//   double * p_epetraView = NULL;
  
//   ((Epetra_Vector &)epetraVector).ExtractView(&p_epetraView);
  
//   assert (0 != p_epetraView && "p_epetraVector->ExtractView() failed to set pointer.");

//   NumPyWrapper::CopyToExternalArray(NumPyWrapper::UNLOAD       ,   // Copy Action
//                                     numPyArray.getArrayObject(),   // PyArrayObject
//                                     p_epetraView               ,   // Source Array
//                                     numPyArray.getDataArray()   ); // Target Array
// }
