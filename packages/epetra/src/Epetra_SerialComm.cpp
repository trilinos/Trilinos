
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

#include "Epetra_SerialComm.h"
#include "Epetra_BasicDirectory.h"

//=============================================================================
Epetra_SerialComm::Epetra_SerialComm()
  : Epetra_Object("Epetra::Comm") {}
 
//=============================================================================
Epetra_SerialComm::Epetra_SerialComm(const Epetra_SerialComm& Comm) : 
  Epetra_Object(Comm.Label()){}

//=============================================================================
void Epetra_SerialComm::Barrier() const {}
//=============================================================================
int Epetra_SerialComm::Broadcast(double * Values, int Count, int Root) const {
  return(0);
}
//=============================================================================
int Epetra_SerialComm::Broadcast(int * Values, int Count, int Root) const {
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
Epetra_Distributor * Epetra_SerialComm::CreateDistributor() const {

  Epetra_Distributor * dist = dynamic_cast<Epetra_Distributor *>(new Epetra_SerialDistributor(*this));
  return(dist);
}
//=============================================================================
Epetra_Directory * Epetra_SerialComm:: CreateDirectory(const Epetra_BlockMap & map) const {

  Epetra_Directory * dir = dynamic_cast<Epetra_Directory *>(new Epetra_BasicDirectory(map));
  return(dir);
}
 //=============================================================================
Epetra_SerialComm::~Epetra_SerialComm()  {}
