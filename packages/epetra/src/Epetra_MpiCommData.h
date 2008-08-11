
//@HEADER
/*
************************************************************************

              Epetra: Linear Algebra Services Package 
                Copyright (2001) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov) 

************************************************************************
*/
//@HEADER

#ifndef EPETRA_MPICOMMDATA_H
#define EPETRA_MPICOMMDATA_H

#include "Epetra_Data.h"
#include <mpi.h>

//! Epetra_MpiCommData:  The Epetra Mpi Communication Data Class.
/*! The Epetra_MpiCommData class is an implementation detail of Epetra_MpiComm.
    It is reference-counted, and can be shared by multiple Epetra_MpiComm instances. 
		It derives from Epetra_Data, and inherits reference-counting from it.
*/

class Epetra_MpiCommData : public Epetra_Data {
	friend class Epetra_MpiComm;
 private:
  //! @name Constructor/Destructor Methods
  //@{ 

  //! Epetra_MpiCommData Default Constructor.
  Epetra_MpiCommData(MPI_Comm & Comm);

  //! Epetra_MpiCommData Destructor.
  ~Epetra_MpiCommData();

  //@}

	MPI_Comm Comm_; //!< \internal MPI_Comm variable.
  int rank_;
  int size_;
  enum {minTag_= 24050};
  enum {maxTag_= 24099};
  // Some day, when the Microsoft 6.0 C++ compiler disappears, we can use ANSI/ISO standard
  // declarations for minTag_ and maxTag_
  //static const int minTag_= 24050;
  //static const int maxTag_= 24099;

	mutable int curTag_;

	// these are intentionally declared but not defined. See Epetra Developer's Guide for details.
  Epetra_MpiCommData(const Epetra_MpiCommData & CommData);
	Epetra_MpiCommData& operator=(const Epetra_MpiCommData & CommData);
  
};
#endif /* EPETRA_MPICOMMDATA_H */
