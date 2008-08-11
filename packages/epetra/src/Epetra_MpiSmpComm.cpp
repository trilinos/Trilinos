
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
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
// ************************************************************************
//@HEADER

#include "Epetra_MpiSmpComm.h"
#include "Epetra_MpiComm.h"


//=============================================================================
Epetra_MpiSmpComm::Epetra_MpiSmpComm(MPI_Comm Comm) 
  : Epetra_Object("Epetra::MpiSmpComm"), 
    MpiSmpCommData_(Comm)
{}
//=============================================================================
Epetra_MpiSmpComm::Epetra_MpiSmpComm(const Epetra_MpiSmpComm& Comm) : 
  Epetra_Object(Comm.Label()),
  MpiSmpCommData_(Comm.MpiSmpCommData_)
{
	MpiSmpCommData_->IncrementReferenceCount();
}

//=============================================================================
void Epetra_MpiSmpComm::Barrier() const {
  MPI_Barrier(MpiSmpCommData->Comm_);
}
//=============================================================================
int Epetra_MpiSmpComm::Broadcast(double * Values, int Count, int Root) const {
  EPETRA_CHK_ERR(MPI_Bcast(Values, Count, MPI_DOUBLE, Root, MpiSmpCommData->Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiSmpComm::Broadcast(int * Values, int Count, int Root) const {
  EPETRA_CHK_ERR(MPI_Bcast(Values, Count, MPI_INT, Root, MpiSmpCommData->Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiSmpComm::Broadcast(long * Values, int Count, int Root) const {
  EPETRA_CHK_ERR(MPI_Bcast(Values, Count, MPI_LONG, Root, MpiSmpCommData->Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiSmpComm::Broadcast(char * Values, int Count, int Root) const {
  EPETRA_CHK_ERR(MPI_Bcast(Values, Count, MPI_CHAR, Root, MpiSmpCommData->Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiSmpComm::GatherAll(double * MyVals, double * AllVals, int Count) const {
  EPETRA_CHK_ERR(MPI_Allgather(MyVals, Count, MPI_DOUBLE, AllVals, Count, MPI_DOUBLE, MpiSmpCommData->Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiSmpComm::GatherAll(int * MyVals, int * AllVals, int Count) const {
  EPETRA_CHK_ERR(MPI_Allgather(MyVals, Count, MPI_INT, AllVals, Count, MPI_INT, MpiSmpCommData->Comm_)); 
  return(0);
}
//=============================================================================
int Epetra_MpiSmpComm::GatherAll(long * MyVals, long * AllVals, int Count) const {
  EPETRA_CHK_ERR(MPI_Allgather(MyVals, Count, MPI_LONG, AllVals, Count, MPI_LONG, MpiSmpCommData->Comm_)); 
  return(0);
}
//=============================================================================
int Epetra_MpiSmpComm::SumAll(double * PartialSums, double * GlobalSums, int Count) const {
  EPETRA_CHK_ERR(MPI_Allreduce(PartialSums, GlobalSums, Count, MPI_DOUBLE, MPI_SUM, MpiSmpCommData->Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiSmpComm::SumAll(int * PartialSums, int * GlobalSums, int Count) const {
  EPETRA_CHK_ERR(MPI_Allreduce(PartialSums, GlobalSums, Count, MPI_INT, MPI_SUM, MpiSmpCommData->Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiSmpComm::SumAll(long * PartialSums, long * GlobalSums, int Count) const {
  EPETRA_CHK_ERR(MPI_Allreduce(PartialSums, GlobalSums, Count, MPI_LONG, MPI_SUM, MpiSmpCommData->Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiSmpComm::MaxAll(double * PartialMaxs, double * GlobalMaxs, int Count) const {
  EPETRA_CHK_ERR(MPI_Allreduce(PartialMaxs, GlobalMaxs, Count, MPI_DOUBLE, MPI_MAX, MpiSmpCommData->Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiSmpComm::MaxAll(int * PartialMaxs, int * GlobalMaxs, int Count) const {
  EPETRA_CHK_ERR(MPI_Allreduce(PartialMaxs, GlobalMaxs, Count, MPI_INT, MPI_MAX, MpiSmpCommData->Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiSmpComm::MaxAll(long * PartialMaxs, long * GlobalMaxs, int Count) const {
  EPETRA_CHK_ERR(MPI_Allreduce(PartialMaxs, GlobalMaxs, Count, MPI_LONG, MPI_MAX, MpiSmpCommData->Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiSmpComm::MinAll(double * PartialMins, double * GlobalMins, int Count) const {
  EPETRA_CHK_ERR(MPI_Allreduce(PartialMins, GlobalMins, Count, MPI_DOUBLE, MPI_MIN, MpiSmpCommData->Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiSmpComm::MinAll(int * PartialMins, int * GlobalMins, int Count) const {
  EPETRA_CHK_ERR(MPI_Allreduce(PartialMins, GlobalMins, Count, MPI_INT, MPI_MIN, MpiSmpCommData->Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiSmpComm::MinAll(long * PartialMins, long * GlobalMins, int Count) const {
  EPETRA_CHK_ERR(MPI_Allreduce(PartialMins, GlobalMins, Count, MPI_LONG, MPI_MIN, MpiSmpCommData->Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiSmpComm::ScanSum(double * MyVals, double * ScanSums, int Count) const {
  EPETRA_CHK_ERR(MPI_Scan(MyVals, ScanSums, Count, MPI_DOUBLE, MPI_SUM, MpiSmpCommData->Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiSmpComm::ScanSum(int * MyVals, int * ScanSums, int Count) const {
  EPETRA_CHK_ERR(MPI_Scan(MyVals, ScanSums, Count, MPI_INT, MPI_SUM, MpiSmpCommData->Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiSmpComm::ScanSum(long * MyVals, long * ScanSums, int Count) const {
  EPETRA_CHK_ERR(MPI_Scan(MyVals, ScanSums, Count, MPI_LONG, MPI_SUM, MpiSmpCommData->Comm_));
  return(0);
}
//=============================================================================
Epetra_Distributor * Epetra_MpiSmpComm:: CreateDistributor() const {
  
  Epetra_MpiComm mpicomm(GetMpiComm());
  Epetra_Distributor * dist = dynamic_cast<Epetra_Distributor *>(new Epetra_MpiDistributor(mpicomm));
  return(dist);
}
//=============================================================================
Epetra_Directory * Epetra_SmpMpiComm:: CreateDirectory(const Epetra_BlockMap & map) const {

  Epetra_Directory * dir = dynamic_cast<Epetra_Directory *>(new Epetra_BasicDirectory(map));
  return(dir);
}
//=============================================================================
Epetra_MpiSmpComm::~Epetra_MpiSmpComm() {
	CleanupData();
}
//=============================================================================
void Epetra_MpiSmpComm::CleanupData() {
	if(MpiSmpCommData_ != 0) {
		MpiSmpCommData_->DecrementReferenceCount();
		if(MpiSmpCommData_->ReferenceCount() == 0) {
			delete MpiSmpCommData_;
			MpiSmpCommData_ = 0;
		}
	}
}
//=============================================================================
Epetra_MpiSmpComm& Epetra_MpiSmpComm::operator= (const Epetra_MpiSmpComm & Comm) {
	if((this != &Comm) && (MpiSmpCommData_ != Comm.MpiSmpCommData_)) {
		CleanupData();
		MpiSmpCommData_ = Comm.MpiSmpCommData_;
		MpiSmpCommData_->IncrementReferenceCount();
	}
	return(*this);
}

