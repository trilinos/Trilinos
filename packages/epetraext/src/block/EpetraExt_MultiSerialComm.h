//@HEADER
/*
************************************************************************

              EpetraExt: Extended Linear Algebra Services Package 
                Copyright (2001) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov) 

************************************************************************
*/
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
  virtual int Broadcast(char * MyVals, int Count, int Root) const
          { return myComm->Broadcast( MyVals, Count, Root); };
  virtual int GatherAll(double * MyVals, double * AllVals, int Count) const
          { return myComm->GatherAll( MyVals,  AllVals, Count); };
  virtual int GatherAll(int * MyVals, int * AllVals, int Count) const
          { return myComm->GatherAll( MyVals, AllVals, Count); };
  virtual int GatherAll(long * MyVals, long * AllVals, int Count) const
          { return myComm->GatherAll( MyVals,  AllVals, Count); };
  virtual int SumAll(double * PartialSums, double * GlobalSums, int Count) const
          { return myComm->SumAll( PartialSums,  GlobalSums, Count); };
  virtual int SumAll(int * PartialSums, int * GlobalSums, int Count) const
          { return myComm->SumAll( PartialSums,  GlobalSums, Count); };
  virtual int SumAll(long * PartialSums, long * GlobalSums, int Count) const
          { return myComm->SumAll( PartialSums,  GlobalSums, Count); };
  virtual int MaxAll(double * PartialMaxs, double * GlobalMaxs, int Count) const
          { return myComm->MaxAll( PartialMaxs,  GlobalMaxs, Count); };
  virtual int MaxAll(int * PartialMaxs, int * GlobalMaxs, int Count) const
          { return myComm->MaxAll( PartialMaxs,  GlobalMaxs, Count); };
  virtual int MaxAll(long * PartialMaxs, long * GlobalMaxs, int Count) const
          { return myComm->MaxAll( PartialMaxs, GlobalMaxs, Count); };
  virtual int MinAll(double * PartialMins, double * GlobalMins, int Count) const
          { return myComm->MinAll( PartialMins, GlobalMins, Count); };
  virtual int MinAll(int * PartialMins, int * GlobalMins, int Count) const
          { return myComm->MinAll( PartialMins, GlobalMins, Count); };
  virtual int MinAll(long * PartialMins, long * GlobalMins, int Count)const
          { return myComm->MinAll( PartialMins, GlobalMins, Count); };
  virtual int ScanSum(double * MyVals, double * ScanSums, int Count)const
          { return myComm->ScanSum( MyVals,  ScanSums, Count); };
  virtual int ScanSum(int * MyVals, int * ScanSums, int Count) const
          { return myComm->ScanSum(MyVals, ScanSums, Count); };
  virtual int ScanSum(long * MyVals, long * ScanSums, int Count) const
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
