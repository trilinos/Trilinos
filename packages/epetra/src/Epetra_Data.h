
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef _EPETRA_DATA_H_
#define _EPETRA_DATA_H_

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

class Epetra_Data {
 protected:
  //@{ \name Constructor/Destructor Methods

  //! Epetra_Data Serial Constructor.
  Epetra_Data();

  //! Epetra_Data Copy Constructor.
  /*! Reference count will be set to 1 on new instance.*/
  Epetra_Data(const Epetra_Data & Data);

  //! Epetra_Data Destructor.
  virtual ~Epetra_Data();

  //@}

	//@{ \name Reference-Counting Methods

	//! Increment reference count
	void IncrementReferenceCount();

	//! Decrement reference count
	void DecrementReferenceCount();

	//! Get reference count
	int ReferenceCount();

	//@}

	int ReferenceCount_;
  
};

#endif /* _EPETRA_DATA_H_ */
