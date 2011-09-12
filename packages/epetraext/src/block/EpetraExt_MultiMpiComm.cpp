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

#include "EpetraExt_MultiMpiComm.h" 

namespace EpetraExt {

MultiMpiComm::MultiMpiComm(MPI_Comm globalMpiComm, int subDomainProcs, int numTimeSteps_) :
	Epetra_MpiComm(globalMpiComm),
	myComm(Teuchos::rcp(new Epetra_MpiComm(globalMpiComm))),
        subComm(0)
{
  //Need to construct subComm for each sub domain, compute subDomainRank,
  //and check that all integer arithmatic works out correctly.
 
  int ierrmpi, size, rank;
  ierrmpi = MPI_Comm_size(globalMpiComm, &size);
  ierrmpi = MPI_Comm_rank(globalMpiComm, &rank);

  if (size % subDomainProcs != 0) {
    cout << "ERROR: num subDomainProcs "<< subDomainProcs
	 << " does not divide into num total procs " << size << endl; 
    exit(-1);
  }

  numSubDomains = size / subDomainProcs;
  numTimeDomains = subDomainProcs;

  // Create split communicators
  MPI_Comm split_MPI_Comm;
  MPI_Comm time_split_MPI_Comm;
  subDomainRank = rank / subDomainProcs;
  timeDomainRank = rank % subDomainProcs;
  ierrmpi =  MPI_Comm_split(globalMpiComm, subDomainRank, rank, 
			    &split_MPI_Comm);
  ierrmpi =  MPI_Comm_split(globalMpiComm, timeDomainRank, rank, 
			    &time_split_MPI_Comm);

  // Construct second epetra communicators
  subComm = new Epetra_MpiComm(split_MPI_Comm);
  timeComm = new Epetra_MpiComm(time_split_MPI_Comm);

  // Compute number of time steps on this sub domain
  ResetNumTimeSteps(numTimeSteps_);

  if (numTimeSteps_ > 0)
    cout << "Processor " << rank << " is on subdomain " << subDomainRank 
         << " and owns " << numTimeStepsOnDomain << " time steps, starting with " 
         <<  firstTimeStepOnDomain << endl;
  else
    cout << "Processor " << rank << " is on subdomain " << subDomainRank << endl;
}

// This constructor is for just one subdomain, so only adds the info
// for multiple time steps on the domain. No two-level parallelism.
MultiMpiComm::MultiMpiComm(const Epetra_MpiComm& EpetraMpiComm_, int numTimeSteps_) :
	Epetra_MpiComm(EpetraMpiComm_),
	myComm(Teuchos::rcp(new Epetra_MpiComm(EpetraMpiComm_))),
        subComm(0)
{

  numSubDomains = 1;
  subDomainRank = 0;
  numTimeSteps = numTimeSteps_;
  numTimeStepsOnDomain = numTimeSteps_;
  firstTimeStepOnDomain = 0;
 
  subComm = new Epetra_MpiComm(EpetraMpiComm_);

  // Create split communicators for time domain
  MPI_Comm time_split_MPI_Comm;
  int rank = EpetraMpiComm_.MyPID();
  int ierrmpi =  MPI_Comm_split(EpetraMpiComm_.Comm(), rank, rank, 
				&time_split_MPI_Comm);
  timeComm = new Epetra_MpiComm(time_split_MPI_Comm);
  numTimeDomains = EpetraMpiComm_.NumProc();
  timeDomainRank = rank;
}
  
//Copy Constructor
MultiMpiComm::MultiMpiComm(const MultiMpiComm &MMC ) :
	Epetra_MpiComm(MMC),
	myComm(Teuchos::rcp(new Epetra_MpiComm(dynamic_cast<const Epetra_MpiComm&>(MMC)))),
        subComm(new Epetra_MpiComm(*MMC.subComm)),
	timeComm(new Epetra_MpiComm(*MMC.timeComm))
{
  numSubDomains = MMC.numSubDomains;
  numTimeDomains = MMC.numTimeDomains;
  subDomainRank = MMC.subDomainRank;
  timeDomainRank = MMC.timeDomainRank;

  numTimeSteps = MMC.numTimeSteps;
  numTimeStepsOnDomain = MMC.numTimeStepsOnDomain;
  firstTimeStepOnDomain = MMC.firstTimeStepOnDomain;
}

MultiMpiComm::~MultiMpiComm()
{
  delete subComm;
  delete timeComm;
}

void MultiMpiComm::ResetNumTimeSteps(int numTimeSteps_)
{
  numTimeSteps = numTimeSteps_;

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
  else {
    numTimeStepsOnDomain = -1;
    firstTimeStepOnDomain = -1;
  }
}

} //namespace EpetraExt
