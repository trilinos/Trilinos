//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// ***********************************************************************
//@HEADER

#ifndef EPETRAEXT_MULTISERIALCOMM_H
#define EPETRAEXT_MULTISERIALCOMM_H

#include "EpetraExt_ConfigDefs.h"
#include "EpetraExt_MultiComm.h" 
#include "Epetra_SerialComm.h" 
#include "Teuchos_RCP.hpp"

//! EpetraExt::MultiSerialComm is a class for keeping track of two levels of
//! parallelism. The class is an Epetra_SerialComm for the global problem 
//! and contains another Eptra_SerialComm of the split problem. Each processor
//! is part of the global problem, and part of a sub-domain.

/*! EpetraExt::MultSerialComm:   The class is an Epetra_SerialComm for the global problem
    and contains another Epetra_SerialComm of the split problem. Each processor
    is part of the global communicator, and a sub-domain communicator. 


<b>Constructing EpetraExt::MultSerialComm objects</b>

*/    

namespace EpetraExt {

class MultiSerialComm: public EpetraExt::MultiComm {
 public:

  //@{ \name Constructors/Destructor.
  //! MultiSerialComm constuctor
  /*! Creates a MultiSerialComm object and communicators for the global and sub- problems.
    
	\param In
	globalComm - MPI communciator (usually MPI_COMM_WORLD)
	\param In 
	subDomainProcss - number of processors in each subdomain. This must divide evenly into the total number of processors of the globalComm.
	\param In 
	numTimeSteps (Default=-1) - Piece of partitioning data needed specifically for parallel space-time project, corresponding to the total number of time steps.
  */
  MultiSerialComm(int numTimeSteps_=-1);
  
  //! Copy constructor.
  MultiSerialComm( const MultiSerialComm &MSC );

  //! Destructor
  virtual ~MultiSerialComm();
  //@}
  
  //! Get reference to split Communicator for sub-domain
  virtual Epetra_Comm& SubDomainComm() const {return *subComm;}

  //! Get reference to split Communicator for time domain
  virtual Epetra_Comm& TimeDomainComm() const { return *subComm; }

  //! Return number of sub-domains that the global problem is split into.
  virtual int NumSubDomains() const {return 1;}

  //! Return integer [0:numSubDomains-1} corresponding to this sub-domain's rank.
  virtual int SubDomainRank() const {return 0;}

  //! Return number of time domains that the global problem is split into.
  virtual int NumTimeDomains() const { return 1; }

  //! Return integer [0:numTimeDomains-1} corresponding to this time-domain's rank.
  virtual int TimeDomainRank() const { return 0; }

  //! Return number of time steps, first step number, on time domain.
  virtual int NumTimeStepsOnDomain() const {return numTimeSteps;}
  virtual int FirstTimeStepOnDomain() const {return 0;}

  //! Return total number of time steps.
  virtual int NumTimeSteps() const {return numTimeSteps;}

  //! Reset total number of time steps, allowing time steps per domain to
  //  be set later than the MultiLevel parallelism is set up.
  void ResetNumTimeSteps(int numTimeSteps);

  virtual Epetra_Comm * Clone() const  { return myComm->Clone(); };
  virtual void Barrier() const { myComm->Barrier(); };
  virtual int Broadcast(double * MyVals, int Count, int Root) const
          { return myComm->Broadcast( MyVals, Count, Root); };
  virtual int Broadcast(int * MyVals, int Count, int Root) const
          { return myComm->Broadcast( MyVals, Count, Root); };
  virtual int Broadcast(long * MyVals, int Count, int Root) const
          { return myComm->Broadcast( MyVals, Count, Root); };
  virtual int Broadcast(long long * MyVals, int Count, int Root) const
          { return myComm->Broadcast( MyVals, Count, Root); };
  virtual int Broadcast(char * MyVals, int Count, int Root) const
          { return myComm->Broadcast( MyVals, Count, Root); };
  virtual int GatherAll(double * MyVals, double * AllVals, int Count) const
          { return myComm->GatherAll( MyVals,  AllVals, Count); };
  virtual int GatherAll(int * MyVals, int * AllVals, int Count) const
          { return myComm->GatherAll( MyVals, AllVals, Count); };
  virtual int GatherAll(long * MyVals, long * AllVals, int Count) const
          { return myComm->GatherAll( MyVals,  AllVals, Count); };
  virtual int GatherAll(long long* MyVals, long long* AllVals, int Count) const
          { return myComm->GatherAll( MyVals,  AllVals, Count); };
  virtual int SumAll(double * PartialSums, double * GlobalSums, int Count) const
          { return myComm->SumAll( PartialSums,  GlobalSums, Count); };
  virtual int SumAll(int * PartialSums, int * GlobalSums, int Count) const
          { return myComm->SumAll( PartialSums,  GlobalSums, Count); };
  virtual int SumAll(long * PartialSums, long * GlobalSums, int Count) const
          { return myComm->SumAll( PartialSums,  GlobalSums, Count); };
  virtual int SumAll(long long* PartialSums, long long* GlobalSums, int Count) const
          { return myComm->SumAll( PartialSums,  GlobalSums, Count); };
  virtual int MaxAll(double * PartialMaxs, double * GlobalMaxs, int Count) const
          { return myComm->MaxAll( PartialMaxs,  GlobalMaxs, Count); };
  virtual int MaxAll(int * PartialMaxs, int * GlobalMaxs, int Count) const
          { return myComm->MaxAll( PartialMaxs,  GlobalMaxs, Count); };
  virtual int MaxAll(long * PartialMaxs, long * GlobalMaxs, int Count) const
          { return myComm->MaxAll( PartialMaxs, GlobalMaxs, Count); };
  virtual int MaxAll(long long* PartialMaxs, long long* GlobalMaxs, int Count) const
          { return myComm->MaxAll( PartialMaxs,  GlobalMaxs, Count); };
  virtual int MinAll(double * PartialMins, double * GlobalMins, int Count) const
          { return myComm->MinAll( PartialMins, GlobalMins, Count); };
  virtual int MinAll(int * PartialMins, int * GlobalMins, int Count) const
          { return myComm->MinAll( PartialMins, GlobalMins, Count); };
  virtual int MinAll(long * PartialMins, long * GlobalMins, int Count)const
          { return myComm->MinAll( PartialMins, GlobalMins, Count); };
  virtual int MinAll(long long* PartialMins, long long* GlobalMins, int Count) const
          { return myComm->MinAll( PartialMins, GlobalMins, Count); };
  virtual int ScanSum(double * MyVals, double * ScanSums, int Count)const
          { return myComm->ScanSum( MyVals,  ScanSums, Count); };
  virtual int ScanSum(int * MyVals, int * ScanSums, int Count) const
          { return myComm->ScanSum(MyVals, ScanSums, Count); };
  virtual int ScanSum(long * MyVals, long * ScanSums, int Count) const
          { return myComm->ScanSum(MyVals, ScanSums, Count); };
  virtual int ScanSum(long long* MyVals, long long* ScanSums, int Count) const
          { return myComm->ScanSum(MyVals, ScanSums, Count); };
  virtual int MyPID() const { return myComm->MyPID(); };
  virtual int NumProc() const { return myComm->NumProc(); };
  virtual Epetra_Distributor * CreateDistributor() const { return myComm->CreateDistributor(); };
  virtual Epetra_Directory * CreateDirectory(const Epetra_BlockMap & Map) const 
          { return myComm->CreateDirectory(Map); };
  virtual void PrintInfo(ostream & os) const { myComm->PrintInfo( os); };

 protected:

  Teuchos::RCP<Epetra_Comm> myComm;
  Epetra_SerialComm* subComm; 
  int numTimeSteps;
};

} //namespace EpetraExt

#endif /* EPETRAEXT_MULTISERIALCOMM_H */
