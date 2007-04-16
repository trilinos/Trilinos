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

#include "EpetraExt_MultiMpiComm.h" 

namespace EpetraExt {

MultiMpiComm::MultiMpiComm(MPI_Comm globalMpiComm, int subDomainProcs, int numTimeSteps_) :
        Epetra_MpiComm(globalMpiComm), subComm(0), numSubDomains(-1),
        subDomainRank(-1), numTimeSteps(numTimeSteps_),
	numTimeStepsOnDomain(-1), firstTimeStepOnDomain(-1)
{
  //Need to construct subComm for each sub domain, compute subDomainRank,
  //and check that all integer arithmatic works out correctly.
 
  int ierrmpi, size, rank;
  ierrmpi = MPI_Comm_size(globalMpiComm, &size);
  ierrmpi = MPI_Comm_rank(globalMpiComm, &rank);

  if (size % subDomainProcs != 0) {cout<<"ERROR: num subDomainProcs "<< subDomainProcs
     << " does not divide into num total procs " << size << endl; exit(-1);}

  numSubDomains = size / subDomainProcs;

  // Create split communicators, the size of subDomainProcs
  MPI_Comm split_MPI_Comm;
  subDomainRank = rank/subDomainProcs;
  ierrmpi =  MPI_Comm_split(globalMpiComm, subDomainRank, rank, &split_MPI_Comm);

  // Construct second epetra communicators
  subComm = new Epetra_MpiComm(split_MPI_Comm);

  // Compute number of time steps on this sub domain
  if (numTimeSteps > 0) {
    // Compute part for number of domains dividing evenly into number of steps
    numTimeStepsOnDomain = numTimeSteps / numSubDomains; 
    firstTimeStepOnDomain = numTimeStepsOnDomain * subDomainRank;

    // Dole out remainder
    int remainder = numTimeSteps % numSubDomains;
    if (subDomainRank < remainder) {
      numTimeStepsOnDomain++; 
      firstTimeStepOnDomain += subDomainRank; 
    }
    else firstTimeStepOnDomain += remainder; 
  }
  cout << "Processor " << rank << " is on subdomain " << subDomainRank 
       << " and owns " << numTimeStepsOnDomain << " time steps, starting with " 
       <<  firstTimeStepOnDomain << endl;
}

// This constructor is for just one subdomain, so only adds the info
// for multiple time steps on the domain. No two-level parallelism.
MultiMpiComm::MultiMpiComm(const Epetra_MpiComm& EpetraMpiComm_, int numTimeSteps_) :
        Epetra_MpiComm(EpetraMpiComm_), subComm(0), numSubDomains(1),
        subDomainRank(0), numTimeSteps(numTimeSteps_),
	numTimeStepsOnDomain(numTimeSteps_), firstTimeStepOnDomain(0)
{
   subComm = new Epetra_MpiComm(EpetraMpiComm_);
}
  
//Copy Constructor
MultiMpiComm::MultiMpiComm(const MultiMpiComm &MMC ) :
        Epetra_MpiComm(MMC), subComm(new Epetra_MpiComm(*(MMC.subComm))),
        numSubDomains(MMC.numSubDomains), subDomainRank(MMC.subDomainRank),
	numTimeSteps(MMC.numTimeSteps), 
	numTimeStepsOnDomain(MMC.numTimeStepsOnDomain), 
	firstTimeStepOnDomain(MMC.firstTimeStepOnDomain) 
{
}

MultiMpiComm::~MultiMpiComm()
{
  delete subComm;
}

} //namespace EpetraExt
