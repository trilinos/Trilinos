
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

