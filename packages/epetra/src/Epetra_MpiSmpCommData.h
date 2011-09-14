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
  //! @name Constructor/Destructor Methods
  //@{ 

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
