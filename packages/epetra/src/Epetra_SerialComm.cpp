
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

#include "Epetra_SerialComm.h"
#include "Epetra_BasicDirectory.h"
#include "Epetra_SerialDistributor.h"

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
int Epetra_SerialComm::Broadcast(char * Values, int Count, int Root) const {
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
int Epetra_SerialComm::GatherAll(long long * MyVals, long long * AllVals, int Count) const {
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
int Epetra_SerialComm::SumAll(long long * PartialSums, long long * GlobalSums, int Count) const {
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
int Epetra_SerialComm::MaxAll(long long * PartialMaxs, long long * GlobalMaxs, int Count) const {
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
int Epetra_SerialComm::MinAll(long long * PartialMins, long long * GlobalMins, int Count) const {
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
int Epetra_SerialComm::ScanSum(long long * MyVals, long long * ScanSums, int Count) const {
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
