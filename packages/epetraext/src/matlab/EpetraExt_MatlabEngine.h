//@HEADER
// ***********************************************************************
// 
//                       MATLAB Engine Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#ifndef MATLAB_ENGINE_H
#define MATLAB_ENGINE_H
#include <Epetra_ConfigDefs.h>

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Comm.h"

// the following deal with matlab provided headers:
#include "engine.h"
#include "mex.h"
#undef printf  // matlab has its own printf that we don't want to use

class Epetra_MultiVector;
class Epetra_RowMatrix;
class Epetra_CrsGraph;
class Epetra_SerialDenseMatrix;
class Epetra_BlockMap;
class Epetra_IntSerialDenseMatrix;
class Epetra_CrsMatrix;

namespace EpetraExt {

//! A class which provides data and command access to Matlab from Epetra.
/*! The EpetraExt_MatlabEngine class allows Epetra data objects to be
exported to Matlab and then operated on within Matlab using Matlab commands. 

When an EpetraExt_MatlabEngine object is constructed a new instance of the application Matlab
is started for EpetraExt_MatlabEngine to communicate with.  All communication between EpetraExt_MatlabEngine
and Matlab occurs on the root node (0) only.  For parallel environments all objects are
collected onto the root node before being put into Matlab.  Object data is put into a mxArray
which is then sent to the Matlab process.  All objects passed to Matlab
are copied into the Matlab memory space.  So at the point when Matlab receives its copy
of the data there is two copies of the mxArray in memory, one in the Epetra application
and one in Matlab.  Since Matlab has its own memory space the mxArray in the Epetra application
should be deleted in order to free up memory space.  All methods in EpetraExt_MatlabEngine that put Epetra
objects into Matlab delete the temporary mxArray as soon as it is put into Matlab.  If a user
desires to create his/her own mxArray's and then send them to Matlab the method PutIntoMatlab
can be used.  It is important to note that PutIntoMatlab does NOT delete the mxArray it is passed.
When the EpetraExt_MatlabEngine deconstructor is called the instance of Matlab that was started
during construction of the EpetraExt_MatlabEngine object exits.

<b>Error Codes</b>
<ul>
  <li> -1 engPutVariable (Matlab provided C function) returned a nonzero result
  <li> -2 internal gather of epetra object and copy to matlab object returned a nonzero result
  <li> -3 engEvalString (Matlab provided C function) returned a nonzero result
  <li> -4 engOutputBuffer (Matlab provided C function) returned a nonzero result
  <li> -5 GetmxArray returned a nonzero result, the variable requested probably does not exist in the Matlab workspace
  <li> -6 GetmxArrayDimensions determined that the mxArray from Matlab has greater than 2 dimensions
  <li> -7 a sparse Matlab mxArray was passed to a method that requires a dense mxArray
  <li> -8 a dense Matlab mxArray was passed to a method that requires a sparse mxArray
</ul>

<b>Build Instructions</b>
These instructions can be found in the file Trilinos/packages/epetraext/doc/matlab.README.

\verbinclude matlab.README

<b>Example Code</b>
The following example code generates simple Epetra objects and then
puts them into Matlab.

The point of this example is to illustrate the flow of calls when
using EpetraExt_MatlabEngine.  This example program can be found in the file
Trilinos/packages/epetraext/example/matlab/cxx_main.cpp.  

\verbinclude cxx_main.cpp

*/
class EpetraExt_MatlabEngine {

  public:

  //@{ \name Constructors/destructors.

  //! EpetraExt_MatlabEngine constructor which creates a MatlabEngine object with a connection to an instance of the application Matlab by starting a new Matlab process. 
  /*!
    \param Comm (In)
           An Epetra_Comm object.

    \return A MatlabEngine object
  */
  EpetraExt_MatlabEngine(const Epetra_Comm& Comm);

  //! EpetraExt_MatlabEngine destructor which closes the connection to Matlab which causes the Matlab process to also exit.
  ~EpetraExt_MatlabEngine();
  //@}
  
  //@{ \name General Epetra to Matlab Access Methods

  //! Sends a command to Matlab.
  /*! Any command that can normally be typed in at the Matlab command line can be passed in to EvalString(char* command).  
	    Commands such as 'help desk', 'edit', and 'plot(MATRIX)' will pop up an interactive window.

			\param command (In)
	       the matlab command to run

			\param outputBuffer (Out)
	       (Optional) a user preallocated buffer for Matlab text output

			\param outputBufferSize (In)
	       (Optional) the size of the outputBuffer
				 
			\return Error Codes, see Detailed Description for more information
   */
  int EvalString (char* command, char* outputBuffer = NULL, int outputBufferSize = -1);

  //@}

  //@{ \name Epetra to Matlab Data Transfer Methods

  //! Puts a copy of the serial or distributed MultiVector into the Matlab workspace.
  /*!
		\param A (In)
	  	the Epetra_MultiVector to put into Matlab
		  
		\param variableName (In)
			the variable name in the Matlab workspace of the Matlab double array (matrix)
		  that will contain the values of the MultiVector

		\return Error Codes, see Detailed Description for more information
  */
  int PutMultiVector(const Epetra_MultiVector& A, const char* variableName);

  //! Puts a copy of the serial or distributed RowMatrix into the Matlab workspace.
  /*!
		\param A (In)
	  	the Epetra_RowMatrix to put into Matlab
		  
		\param variableName (In)
	    the variable name in the Matlab workspace of the Matlab sparse double array (matrix)
		  that will contain the values of the RowMatrix

		\param transA (In)
	    if true then the transpose of A is put into Matlab
		  NOTE: It is faster to put the transpose of A into Matlab since Matlab stores matrices
		        in column-major form whereas Epetra stores them in row-major form.

		\return Error Codes, see Detailed Description for more information
  */
  int PutRowMatrix(const Epetra_RowMatrix& A, const char* variableName, bool transA);

  //! not implemented yet
  int PutCrsGraph(const Epetra_CrsGraph& A, const char* variableName, bool transA);

  //! Puts a copy of the SerialDenseMatrix into the Matlab workspace.
  /*!
		\param A (In)
	    the Epetra_SerialDenseMatrix to put into Matlab
		  
		\param variableName (In)
	    the variable name in the Matlab workspace of the Matlab double array (matrix)
		  that will contain the values of the SerialDenseMatrix

		\param proc (In)
	    for serial environment set to 0
		  for a parallel environment set to the process ID that owns the SerialDenseMatrix
		  \warning The same parameters must be passed to each process.

		\return Error Codes, see Detailed Description for more information
  */
  int PutSerialDenseMatrix(const Epetra_SerialDenseMatrix& A, const char* variableName, int proc=0);

  //! Puts a copy of the IntSerialDenseMatrix into the Matlab workspace.
  /*!
		\param A (In)
	    the Epetra_IntSerialDenseMatrix to put into Matlab
		  
		\param variableName (In)
	    the variable name in the Matlab workspace of the Matlab double array (matrix)
		  that will contain the values of the IntSerialDenseMatrix

		\param proc (In)
	    for serial environment set to 0
		  for a parallel environment set to the process ID that owns the IntSerialDenseMatrix
		  \warning The same parameters must be passed to each process.

		\return Error Codes, see Detailed Description for more information
  */
  int PutIntSerialDenseMatrix(const Epetra_IntSerialDenseMatrix& A, const char* variableName, int proc=0);

  //! Puts a copy of the BlockMap or Map into the Matlab workspace.
  /*!
	  \param blockMap (In)
	  	the Epetra_BlockMap to put into Matlab
		  
	  \param variableName (In)
	    the variable name in the Matlab workspace of the Matlab sparse double array (matrix)
		  that will contain the values of the BlockMap

	  \param transA (In)
	    if true then the transpose of blockMap is put into Matlab
		  NOTE: It is faster to put the transpose of blockMap into Matlab since Matlab stores matrices
		        in column-major form whereas Epetra stores them in row-major form.

	  \return Error Codes, see Detailed Description for more information
  */
  int PutBlockMap(const Epetra_BlockMap& blockMap, const char* variableName, bool transA);
  
  //! Puts a mxArray into Matlab.
  /*! The Matlab provided C library provides mxArray which is used to construct and fill a Matlab
      object before sending it to Matlab to be put into the Matlab workspace.  The mxArray is copied
      into the Matlab memory space and is not needed after it has been passed to Matlab.  mxArrays
      should be destroyed using mxDestoryArray(mxArray) when they are no longer needed by the C/C++
      program using them.  Objects in Matlab must be destroyed using EvalString(char* command) and
	    the appropriate Matlab command to destroy the object.  EpetraExt::MatlabEngine uses PutIntoMatlab
	    to pass all mxArrays it generates into Matlab.  However, a user can create, fill, and put his/her
	    own mxArrays into Matlab using this method.  To create a mxArray mex.h must be included.  For more
	    information on how to use mxArrays see Matlab's documentation (type helpdesk at the Matlab command prompt)
	    and see the External API Reference section.

	  \param variableName (In)
	    the name for the mxArray once it has been put into the Matlab workspace

	  \param matlabA (In)
	    the mxArray to put into the Matlab workspace

	  \return Matlab error code from engPutVariable for Matlab versions >= 6.5 or from engPutArray for all other versions
  */
  int PutIntoMatlab(const char* variableName, mxArray* matlabA);

  //@}

  //@{ \name Matlab to Epetra Data Transfer Methods

	//! Puts a Matlab variable into a MultiVector.
	/*! The values from the Matlab variable are exported to the MultiVector using an export object.  Therefore
	    the MultiVector must be prepared by the user just like any MultiVector would be before calling an export.
			
			\param variableName (In)
			   the name of the Matlab variable to be put into the given MultiVector
				
		  \param A (In)
			   the MultiVector to put the Matlab variable values into
							 
			\return Error Codes, see Detailed Description for more information
	*/
  int GetMultiVector(const char* variableName, Epetra_MultiVector& A);
	
	//! Puts a Matlab variable into a SerialDenseMatrix on the specified PE.
	/*! The SerialDenseMatrix must be constructed by the user and have the proper amount of space
	    to hold the values that will be copied from the given Matlab variable.
			
			\param variableName (In)
			   the name of the Matlab variable to be put into the given SerialDenseMatrix
				
		  \param A (In)
			   the SerialDenseMatrix to put the Matlab variable values into
	
	    \param proc (In)
			   the PE that will own the SerialDenseMatrix
							 
			\return Error Codes, see Detailed Description for more information
	*/
  int GetSerialDenseMatrix(const char* variableName, Epetra_SerialDenseMatrix& A, int proc=0);
	
  //! Puts a Matlab variable into a IntSerialDenseMatrix on the specified PE.
	/*! The IntSerialDenseMatrix must be constructed by the user and have the proper amount of space
	    to hold the values that will be copied from the given Matlab variable.
			
			\param variableName (In)
			   the name of the Matlab variable to be put into the given IntSerialDenseMatrix
				
		  \param A (In)
			   the IntSerialDenseMatrix to put the Matlab variable values into
	
	    \param proc (In)
			   the PE that will own the IntSerialDenseMatrix
							 
			\return Error Codes, see Detailed Description for more information
	*/
  int GetIntSerialDenseMatrix(const char* variableName, Epetra_IntSerialDenseMatrix& A, int proc=0);
	
	//! Puts a Matlab variable into a CrsMatrix.
	/*! The values from the Matlab variable are exported to the CrsMatrix using an export object.  Therefore
	    the CrsMatrix must be prepared by the user just like any CrsMatrix would be before calling an export.
			\warning Getting a CrsMatrix from Matlab will cause a temporary Matlab variable to be created and then
			         deleted named TRANS_variableName where variableName is the name of the Matlab variable.
			
			\param variableName (In)
			   the name of the Matlab variable to be put into the given CrsMatrix
				
		  \param A (In)
			   the CrsMatrix to put the Matlab variable values into
	
	    \param getTrans (In)
			   if false then a temporary Matlab variable is created TRANS_variableName
	       if true then the transpose of A copied from Matlab
		     NOTE: It is faster to copy the transpose of A from Matlab since Matlab stores matrices
		           in column-major form whereas Epetra stores them in row-major form.
							 
			\return Error Codes, see Detailed Description for more information
	*/
  int GetCrsMatrix(const char* variableName, Epetra_CrsMatrix& A, bool getTrans);
	
	//@}
		
  //@{ \name General Matlab to Epetra Access Methods
		
	//! Get general information about the mxArray.  For internal use but can be used by an advanced user.
	/*! Makes several Matlab function calls on the mxArray in order to determine the number of rows, columns,
	    nonzeros, and whether or not the mxArray is sparse.
			
			\param matlabA (In)
			  the mxArray to get information about
				
		  \param isSparse (Out)
			  true if the mxArray is sparse
				
			\param numRows (Out)
			  the number of rows in the mxArray
				
			\param numCols (Out)
			  the number of columns in the mxArray
				
			\param numNonZeros (Out)
			  the number of nonzeros in the mxArray
	
	    \return Error Codes, see Detailed Description for more information
	*/
  int GetmxArrayDimensions(mxArray* matlabA, bool& isSparse, int& numRows, int& numCols, int& numNonZeros);
  
  //! Get a mxArray from Matlab.  For internal use but can be used by an advanced user.
  /*! Calls the Matlab provided engGetVariable function which copies the specified variable from the Matlab workspace
	    into a mxArray.  Hence any changes to mxArray will not show up in Matlab, and any changes to the Matlab variable
			in Matlab will not show up in the mxArray.  When finished with the mxArray mxDestroyArray should be called on
			the mxArray.  Matlab appears to perform a copy to a new mxArray each time engGetVariable is called.
			
			\param variableName (In)
				the mxArray to get from the Matlab workspace
				
		  \param matlabA (Out)
			  a pointer to a mxArray* to put the mxArray that engGetVariable returns
			
      \return Error Codes, see Detailed Description for more information
  */
  int GetmxArray(const char* variableName, mxArray** matlabA);
  
	//@}

 private:
  // the Matlab Engine object provided by libeng and engine.h
  Engine* Engine_ ;
  // the Epetra_Comm object that will be used for all communication
  const Epetra_Comm& Comm_ ;

};
} // namespace EpetraExt

#endif /* MATLAB_ENGINE_H */
