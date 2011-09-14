/*
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef EPETRA_DATA_H
#define EPETRA_DATA_H

#include "Epetra_ConfigDefs.h"

//! Epetra_Data:  The Epetra Base Data Class.
/*! The Epetra_Data class is a base class for all Epetra Data Classes.
	  It provides a mechanism so that one data object can be shared by multiple
		class instances. However, it is meant only to be used internally by 
		another Epetra class. It does not provide smart pointer like capabilities.
		Incrementing and decrementing the reference count, and deleting the 
		data class instance (if necessary), are duties of the Epetra class 
		utilizing Epetra_Data.

		All of Epetra_Data's methods are protected. This is because Epetra_Data 
		should never be used directly. Rather, a class that derives from
		Epetra_Data should be used instead. For example, Epetra_MpiCommData or 
		Epetra_BlockMapData.

		DEVELOPER NOTES: 
		(1) Any class that inherits from Epetra_Data may need to define an 
		assignment operator, if it adds pointers. Epetra_Data doesn't have any, 
		and so the default (compiler-generated) assignment operator is good enough. 
		(2) The behavior of a derived class is left up to the 
		implementer(s) of that class. As such, it cannot be assumed that 
		just because a class inherits from Epetra_Data, that it supports copy 
		construction or assignment, or that it will perform as expected. 
*/

class EPETRA_LIB_DLL_EXPORT Epetra_Data {
 protected:
   //! @name Constructor/Destructor Methods
  //@{ 

  //! Epetra_Data Serial Constructor.
  Epetra_Data();

  //! Epetra_Data Copy Constructor.
  /*! Reference count will be set to 1 on new instance.*/
  Epetra_Data(const Epetra_Data & Data);

  //! Epetra_Data Destructor.
  virtual ~Epetra_Data();

  //@}

  //! @name Reference-Counting Methods
	//@{ 

	//! Increment reference count
	void IncrementReferenceCount();

	//! Decrement reference count
	void DecrementReferenceCount();

	//! Get reference count
	int ReferenceCount() const;

	//@}

	int ReferenceCount_;
  
};

#endif /* EPETRA_DATA_H */
