
//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

#ifndef EPETRA_MPISMPCOMMDATA_H
#define EPETRA_MPISMPCOMMDATA_H

#include "Epetra_Data.h"
#include <mpi.h>

//! Epetra_MpiSmpCommData:  The Epetra Mpi Shared Memory ParallelCommunication Data Class.
/*! The Epetra_MpiSmpCommData class is an implementation detail of Epetra_MpiSmpComm.
    It is reference-counted, and can be shared by multiple Epetra_MpiSmpComm instances. 
		It derives from Epetra_Data, and inherits reference-counting from it.
*/

class Epetra_MpiSmpCommData : public Epetra_Data {
	friend class Epetra_MpiSmpComm;
 private:
  //@{ \name Constructor/Destructor Methods

  //! Epetra_MpiSmpCommData Default Constructor.
  Epetra_MpiSmpCommData(MPI_Comm & Comm);

  //! Epetra_MpiSmpCommData Destructor.
  ~Epetra_MpiSmpCommData();

  //@}

	MPI_Comm Comm_; //!< \internal MPI_Comm variable.
  int rank_;
  int size_;
  int minTag_;
  int maxTag_;
  mutable int curTag_;
  int ThreadID_;
  int NodeID_;
  int NumThreads_;

	// these are intentionally declared but not defined. See Epetra Developer's Guide for details.
  Epetra_MpiSmpCommData(const Epetra_MpiSmpCommData & CommData);
	Epetra_MpiSmpCommData& operator=(const Epetra_MpiSmpCommData & CommData);
  
};
#endif /* EPETRA_MPISMPCOMMDATA_H */
