
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

#include "Epetra_SerialComm.h"
#include "Epetra_BasicDirectory.h"

//=============================================================================
Epetra_SerialComm::Epetra_SerialComm()
  : Epetra_Object("Epetra::Comm"), 
		SerialCommData_(new Epetra_SerialCommData()) {}
 
//=============================================================================
Epetra_SerialComm::Epetra_SerialComm(const Epetra_SerialComm& Comm) 
	: Epetra_Object(Comm.Label()), 
		SerialCommData_(Comm.SerialCommData_) 
{
	SerialCommData_->IncrementReferenceCount();
}
//=============================================================================
void Epetra_SerialComm::Barrier() const {}
//=============================================================================
int Epetra_SerialComm::Broadcast(double * Values, int Count, int Root) const {
  (void)Values;
  (void)Count;
  (void)Root;
  return(0);
}
//=============================================================================
int Epetra_SerialComm::Broadcast(int * Values, int Count, int Root) const {
  (void)Values;
  (void)Count;
  (void)Root;
  return(0);
}
//=============================================================================
int Epetra_SerialComm::Broadcast(long * Values, int Count, int Root) const {
  (void)Values;
  (void)Count;
  (void)Root;
  return(0);
}
//=============================================================================
int Epetra_SerialComm::GatherAll(double * MyVals, double * AllVals, int Count) const {
  for (int i=0; i<Count; ++i) AllVals[i] = MyVals[i];
  return(0);
}
//=============================================================================
int Epetra_SerialComm::GatherAll(int * MyVals, int * AllVals, int Count) const {
  for (int i=0; i<Count; ++i) AllVals[i] = MyVals[i];
  return(0);
}
//=============================================================================
int Epetra_SerialComm::GatherAll(long * MyVals, long * AllVals, int Count) const {
  for (int i=0; i<Count; ++i) AllVals[i] = MyVals[i];
  return(0);
}
//=============================================================================
int Epetra_SerialComm::SumAll(double * PartialSums, double * GlobalSums, int Count) const {
  for (int i=0; i<Count; i++) GlobalSums[i] = PartialSums[i];
  return(0);
}
//=============================================================================
int Epetra_SerialComm::SumAll(int * PartialSums, int * GlobalSums, int Count) const {
  for (int i=0; i<Count; i++) GlobalSums[i] = PartialSums[i];
  return(0);
}
//=============================================================================
int Epetra_SerialComm::SumAll(long * PartialSums, long * GlobalSums, int Count) const {
  for (int i=0; i<Count; i++) GlobalSums[i] = PartialSums[i];
  return(0);
}
//=============================================================================
int Epetra_SerialComm::MaxAll(double * PartialMaxs, double * GlobalMaxs, int Count) const {
  for (int i=0; i<Count; i++) GlobalMaxs[i] = PartialMaxs[i];
  return(0);
}
//=============================================================================
int Epetra_SerialComm::MaxAll(int * PartialMaxs, int * GlobalMaxs, int Count) const {
  for (int i=0; i<Count; i++) GlobalMaxs[i] = PartialMaxs[i];
  return(0);
}
//=============================================================================
int Epetra_SerialComm::MaxAll(long * PartialMaxs, long * GlobalMaxs, int Count) const {
  for (int i=0; i<Count; i++) GlobalMaxs[i] = PartialMaxs[i];
  return(0);
}
//=============================================================================
int Epetra_SerialComm::MinAll(double * PartialMins, double * GlobalMins, int Count) const {
  for (int i=0; i<Count; i++) GlobalMins[i] = PartialMins[i];
  return(0);
}
//=============================================================================
int Epetra_SerialComm::MinAll(int * PartialMins, int * GlobalMins, int Count) const {
  for (int i=0; i<Count; i++) GlobalMins[i] = PartialMins[i];
  return(0);
}
//=============================================================================
int Epetra_SerialComm::MinAll(long * PartialMins, long * GlobalMins, int Count) const {
  for (int i=0; i<Count; i++) GlobalMins[i] = PartialMins[i];
  return(0);
}
//=============================================================================
int Epetra_SerialComm::ScanSum(double * MyVals, double * ScanSums, int Count) const {
  for (int i=0; i<Count; i++) ScanSums[i] = MyVals[i];
  return(0);
}
//=============================================================================
int Epetra_SerialComm::ScanSum(int * MyVals, int * ScanSums, int Count) const {
  for (int i=0; i<Count; i++) ScanSums[i] = MyVals[i];
  return(0);
}
//=============================================================================
int Epetra_SerialComm::ScanSum(long * MyVals, long * ScanSums, int Count) const {
  for (int i=0; i<Count; i++) ScanSums[i] = MyVals[i];
  return(0);
}
//=============================================================================
Epetra_Distributor * Epetra_SerialComm::CreateDistributor() const {

  Epetra_Distributor * dist = dynamic_cast<Epetra_Distributor *>(new Epetra_SerialDistributor(*this));
  return(dist);
}
//=============================================================================
Epetra_Directory * Epetra_SerialComm::CreateDirectory(const Epetra_BlockMap & map) const {

  Epetra_Directory * dir = dynamic_cast<Epetra_Directory *>(new Epetra_BasicDirectory(map));
  return(dir);
}
//=============================================================================
Epetra_SerialComm::~Epetra_SerialComm()  {
	CleanupData();
}
//=============================================================================
Epetra_SerialComm& Epetra_SerialComm::operator= (const Epetra_SerialComm & Comm) {
	if((this != &Comm) && (SerialCommData_ != Comm.SerialCommData_)) {
		CleanupData();
		SerialCommData_ = Comm.SerialCommData_;
		SerialCommData_->IncrementReferenceCount();
	}
	return(*this);
}
//=============================================================================
void Epetra_SerialComm::CleanupData() {
	if(SerialCommData_ != 0) {
		SerialCommData_->DecrementReferenceCount();
		if(SerialCommData_->ReferenceCount() == 0) {
			delete SerialCommData_;
			SerialCommData_ = 0;
		}
	}
}
//=============================================================================
int Epetra_SerialComm::ReferenceCount() const {
  return(SerialCommData_->ReferenceCount());
}
