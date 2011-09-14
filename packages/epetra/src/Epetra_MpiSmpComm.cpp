
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

