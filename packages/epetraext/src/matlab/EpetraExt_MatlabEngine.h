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

//! MatlabEngine: 
/*! The EpetraExt_MatlabEngine class allows Epetra data objects to be
exported to Matlab and then operated on within Matlab using Matlab commands. 

When a MatlabEngine object is constructed a new instance of the application Matlab
is started for MatlabEngine to communicate with.  All communicatin between MatlabEngine
and Matlab occurs on the root node (0) only.  For parallel environments all objects are
collected onto the root node before being put into Matlab.  Object data is put into a mxArray
which is then sent to the Matlab process.  All objects passed to Matlab
are copied into the Matlab memory space.  So at the point when Matlab receives its copy
of the data there is two copies of the mxArray in memory, one in the Epetra application
and one in Matlab.  Since Matlab has its own memory space the mxArray in the Epetra application
should be deleted in order to free up memory space.  All methods in MatlabEngine that put Epetra
objects into Matlab delete the temporary mxArray as soon as it is put into Matlab.  If a user
desires to create his/her own mxArray's and then send them to Matlab the method PutIntoMatlab
can be used.  It is important to note that PutIntoMatlab does NOT delete the mxArray it is passed.
When the MatlabEngine deconstructor is called the instance of Matlab that was started
during construction of the MatlabEngine object exits.

<b>Error Codes</b>
<ul>
  <li> -1 engPutVariable returned a nonzero result
  <li> -2 internal gather of epetra object and copy to matlab object returned a nonzero result
  <li> -3 engEvalString returned a nonzero result
  <li> -4 engOutputBuffer returned a nonzero result
</ul>


*/

namespace EpetraExt {
//=========================================================================
class MatlabEngine {

  public:

  //@{ \name Constructors/destructors.
  /*! MatlabEngine constructor which creates a MatlabEngine object with a connection to an instance of the
 	  application Matlab by starting a new Matlab process. 

    \param Comm In
           An Epetra Communicator 

    \return  MatlabEngine object
  */
  MatlabEngine(const Epetra_Comm& Comm);
  //! MatlabEngine destructor which closes the connection to Matlab which causes the Matlab process to also exit.
  ~MatlabEngine();
  //@}
  
  //@{ \name General Matlab Access Methods

  //! EvalString method
  int EvalString (char* command);
  int EvalString (char* command, char* outputBuffer, int outputBufferSize);

  //@}

  //@{ \name Epetra to Matlab Data Transfer Methods
  //! Put a copy of the serial or distributed MultiVector into the Matlab workspace.
  /*!
	\param const Epetra_MultiVector& A
	       the Epetra_MultiVector to put into Matlab
		  
	\param const char* variableName
	       the variable name in the Matlab workspace of the Matlab double array (matrix)
		   that will contain the values of the MultiVector

	\return Error Codes, see Detailed Description for more information
  */
  int PutMultiVector(const Epetra_MultiVector& A, const char* variableName);

  //! Put a copy of the serial or distributed RowMatrix into the Matlab workspace.
  /*!
	\param const Epetra_RowMatrix& A
	       the 	\param const Epetra_RowMatrix& A to put into Matlab
		  
	\param const char* variableName
	       the variable name in the Matlab workspace of the Matlab sparse double array (matrix)
		   that will contain the values of the RowMatrix

	\param bool transA
	       if true then the transpose of A is put into Matlab
		   NOTE: It is faster to put the transpose of A into Matlab since Matlab stores matrices
		         in column-major form whereas Epetra stores them in row-major form.

	\return Error Codes, see Detailed Description for more information
  */
  int PutRowMatrix(const Epetra_RowMatrix& A, const char* variableName, bool transA);
  int PutCrsGraph(const Epetra_CrsGraph& A, const char* variableName, bool transA);

  //! Put a copy of the SerialDenseMatrix into the Matlab workspace.
  /*!
	\param const Epetra_SerialDenseMatrix& A
	       the Epetra_SerialDenseMatrix to put into Matlab
		  
	\param const char* variableName
	       the variable name in the Matlab workspace of the Matlab double array (matrix)
		   that will contain the values of the SerialDenseMatrix

	\param int proc
	       for serial environment set to 0
		   for a parallel environment set to the process ID that owns the SerialDenseMatrix
		   WARNING: the same parameters must be passed to each process

	\return Error Codes, see Detailed Description for more information
  */
  int PutSerialDenseMatrix(const Epetra_SerialDenseMatrix& A, const char* variableName, int proc);

  //! Put a copy of the IntSerialDenseMatrix into the Matlab workspace.
  /*!
	\param const Epetra_IntSerialDenseMatrix& A
	       the Epetra_IntSerialDenseMatrix to put into Matlab
		  
	\param const char* variableName
	       the variable name in the Matlab workspace of the Matlab double array (matrix)
		   that will contain the values of the IntSerialDenseMatrix

	\param int proc
	       for serial environment set to 0
		   for a parallel environment set to the process ID that owns the IntSerialDenseMatrix
		   WARNING: the same parameters must be passed to each process

	\return Error Codes, see Detailed Description for more information
  */
  int PutIntSerialDenseMatrix(const Epetra_IntSerialDenseMatrix& A, const char* variableName, int proc);

  int PutBlockMap(const Epetra_BlockMap& blockMap, const char* variableName);
  
  //*! Put a mxArray into Matlab.
  /*! The Matlab provided C library provides mxArray which is used to construct and fill a Matlab
      object before sending it to Matlab to be put into the Matlab workspace.  The mxArray is copied
      into the Matlab memory space and is not needed after it has been passed to Matlab.  mxArrays
      should be destoryed using mxDestoryArray(mxArray) when they are no longer needed by the C/C++
	  program using them.  Objects in Matlab must be destoryed using EvalString(char* command) and
	  the appropriate Matlab command to destory the object.  EpetraExt::MatlabEngine uses PutIntoMatlab
	  to pass all mxArrays it generates into Matlab.  However, a user can create, fill, and put his/her
	  own mxArrays into Matlab using this method.  To create a mxArray mex.h must be included.  For more
	  information on how to use mxArrays see Matlab's documentation (type helpdesk at the Matlab command prompt)
	  and see the External API Reference section.

	  \param const char* variableName
	      the name for the mxArray once it has been put into the Matlab workspace

	  \param mxArray* matlabA
	      the mxArray to put into the Matlab workspace

	  \return Matlab error code from engPutVariable for Matlab versions >= 6.5 or fromm engPutArray for all other versions
  */
  int PutIntoMatlab(const char* variableName, mxArray* matlabA);
  //@}

 private:
  // the Matlab Engine object provided by libeng and engine.h
  Engine* Engine_ ;
  // the Epetra_Comm object that will be used for all communication
  const Epetra_Comm& Comm_ ;

};
} // namespace EpetraExt

#endif /* MATLAB_ENGINE_H */
