
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

#ifndef _EPETRA_SERIALCOMMDATA_H_
#define _EPETRA_SERIALCOMMDATA_H_

#include "Epetra_Data.h"

//! Epetra_SerialCommData:  The Epetra Serial Communication Data Class.
/*! The Epetra_SerialCommData class is an implementation detail of Epetra_SerialComm.
    It is reference-counted, and can be shared by multiple Epetra_SerialComm instances. 
		It derives from Epetra_Data, and inherits reference-counting from it.
*/

class Epetra_SerialCommData : public Epetra_Data {
	friend class Epetra_SerialComm;
 private:
  //@{ \name Constructor/Destructor Methods

  //! Epetra_SerialCommData Default Constructor.
  Epetra_SerialCommData();

  //! Epetra_SerialCommData Destructor.
  ~Epetra_SerialCommData();

  //@}

	int MyPID_;
	int NumProc_;

	// these are intentionally declared but not defined. See Epetra Developer's Guide for details.
  Epetra_SerialCommData(const Epetra_SerialCommData & CommData);
	Epetra_SerialCommData& operator=(const Epetra_SerialCommData & CommData);
  
};
#endif /* _EPETRA_SERIALCOMMDATA_H_ */
