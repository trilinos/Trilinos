
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

#include "Epetra_MpiComm.h"


//=============================================================================
Epetra_MpiComm::Epetra_MpiComm(MPI_Comm Comm) 
  : Epetra_Object("Epetra::MpiComm"), 
    Comm_(Comm),
    curTag_(minTag_)
{
    MPI_Comm_size(Comm, &size_);
    MPI_Comm_rank(Comm, &rank_);
}
//=============================================================================
Epetra_MpiComm::Epetra_MpiComm(const Epetra_MpiComm& Comm) : 
  Epetra_Object(Comm.Label()),
  Comm_(Comm.Comm_),
  rank_(Comm.rank_),
  size_(Comm.size_),
  curTag_(Comm.curTag_)
{
}

//=============================================================================
void Epetra_MpiComm::Barrier() const {
  MPI_Barrier(Comm_);
}
//=============================================================================
int Epetra_MpiComm::Broadcast(double * Values, int Count, int Root) const {
  EPETRA_CHK_ERR(CheckInput(Values,Count));
  EPETRA_CHK_ERR(MPI_Bcast(Values, Count, MPI_DOUBLE, Root, Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiComm::Broadcast(int * Values, int Count, int Root) const {
  EPETRA_CHK_ERR(CheckInput(Values,Count));
  EPETRA_CHK_ERR(MPI_Bcast(Values, Count, MPI_INT, Root, Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiComm::GatherAll(double * MyVals, double * AllVals, int Count) const {
  EPETRA_CHK_ERR(CheckInput(MyVals,Count));
  EPETRA_CHK_ERR(CheckInput(AllVals,Count));
  EPETRA_CHK_ERR(MPI_Allgather(MyVals, Count, MPI_DOUBLE, AllVals, Count, MPI_DOUBLE, Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiComm::GatherAll(int * MyVals, int * AllVals, int Count) const {
  EPETRA_CHK_ERR(CheckInput(MyVals,Count));
  EPETRA_CHK_ERR(CheckInput(AllVals,Count));
  EPETRA_CHK_ERR(MPI_Allgather(MyVals, Count, MPI_INT, AllVals, Count, MPI_INT, Comm_)); 
  return(0);
}
//=============================================================================
int Epetra_MpiComm::SumAll(double * PartialSums, double * GlobalSums, int Count) const {
  EPETRA_CHK_ERR(CheckInput(PartialSums,Count));
  EPETRA_CHK_ERR(CheckInput(GlobalSums,Count));
  EPETRA_CHK_ERR(MPI_Allreduce(PartialSums, GlobalSums, Count, MPI_DOUBLE, MPI_SUM, Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiComm::SumAll(int * PartialSums, int * GlobalSums, int Count) const {
  EPETRA_CHK_ERR(CheckInput(PartialSums,Count));
  EPETRA_CHK_ERR(CheckInput(GlobalSums,Count));
  EPETRA_CHK_ERR(MPI_Allreduce(PartialSums, GlobalSums, Count, MPI_INT, MPI_SUM, Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiComm::MaxAll(double * PartialMaxs, double * GlobalMaxs, int Count) const {
  EPETRA_CHK_ERR(CheckInput(PartialMaxs,Count));
  EPETRA_CHK_ERR(CheckInput(GlobalMaxs,Count));
  EPETRA_CHK_ERR(MPI_Allreduce(PartialMaxs, GlobalMaxs, Count, MPI_DOUBLE, MPI_MAX, Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiComm::MaxAll(int * PartialMaxs, int * GlobalMaxs, int Count) const {
  EPETRA_CHK_ERR(CheckInput(PartialMaxs,Count));
  EPETRA_CHK_ERR(CheckInput(GlobalMaxs,Count));
  EPETRA_CHK_ERR(MPI_Allreduce(PartialMaxs, GlobalMaxs, Count, MPI_INT, MPI_MAX, Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiComm::MinAll(double * PartialMins, double * GlobalMins, int Count) const {
  EPETRA_CHK_ERR(CheckInput(PartialMins,Count));
  EPETRA_CHK_ERR(CheckInput(GlobalMins,Count));
  EPETRA_CHK_ERR(MPI_Allreduce(PartialMins, GlobalMins, Count, MPI_DOUBLE, MPI_MIN, Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiComm::MinAll(int * PartialMins, int * GlobalMins, int Count) const {
  EPETRA_CHK_ERR(CheckInput(PartialMins,Count));
  EPETRA_CHK_ERR(CheckInput(GlobalMins,Count));
  EPETRA_CHK_ERR(MPI_Allreduce(PartialMins, GlobalMins, Count, MPI_INT, MPI_MIN, Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiComm::ScanSum(double * MyVals, double * ScanSums, int Count) const {
  EPETRA_CHK_ERR(CheckInput(MyVals,Count));
  EPETRA_CHK_ERR(CheckInput(ScanSums,Count));
  EPETRA_CHK_ERR(MPI_Scan(MyVals, ScanSums, Count, MPI_DOUBLE, MPI_SUM, Comm_));
  return(0);
}
//=============================================================================
int Epetra_MpiComm::ScanSum(int * MyVals, int * ScanSums, int Count) const {
  EPETRA_CHK_ERR(CheckInput(MyVals,Count));
  EPETRA_CHK_ERR(CheckInput(ScanSums,Count));
  EPETRA_CHK_ERR(MPI_Scan(MyVals, ScanSums, Count, MPI_INT, MPI_SUM, Comm_));
  return(0);
}
//=============================================================================
Epetra_Distributor * Epetra_MpiComm:: CreateDistributor() const {

  Epetra_Distributor * dist = dynamic_cast<Epetra_Distributor *>(new Epetra_MpiDistributor(*this));
  return(dist);
}
//=============================================================================
Epetra_Directory * Epetra_MpiComm:: CreateDirectory(const Epetra_BlockMap & map) const {

  Epetra_Directory * dir = dynamic_cast<Epetra_Directory *>(new Epetra_BasicDirectory(map));
  return(dir);
}
//=============================================================================
Epetra_MpiComm::~Epetra_MpiComm()  {}
