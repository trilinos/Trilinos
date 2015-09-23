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

#ifndef EPETRAEXT_MULTICOMM_H
#define EPETRAEXT_MULTICOMM_H

#include "EpetraExt_ConfigDefs.h"
#include "Epetra_Comm.h" 

//! EpetraExt::MultiComm is an abstract class for keeping track of two levels of
//! parallelism. The class is an Epetra_Comm for the global problem 
//! and contains another Eptra_Comm of the split problem. Each processor
//! is part of the global problem, and part of a sub-domain.

/*! EpetraExt::MultComm:   The class is an Epetra_Comm for the global problem
    and contains another Epetra_Comm of the split problem. Each processor
    is part of the global communicator, and a sub-domain communicator. 


<b>Constructing EpetraExt::MultComm objects</b>

*/    

namespace EpetraExt {

class MultiComm: public virtual Epetra_Comm {
 public:

  //! Constructor
  MultiComm() {}

  //! Destructor
  virtual ~MultiComm() {};
  
  //! Get reference to split Communicator for sub-domain
  virtual Epetra_Comm& SubDomainComm() const = 0;

  //! Get reference to split Communicator for time domain
  virtual Epetra_Comm& TimeDomainComm() const = 0;

  //! Return number of sub-domains that the global problem is split into.
  virtual int NumSubDomains() const  = 0;

  //! Return integer [0:numSubDomains-1} corresponding to this sub-domain's rank.
  virtual int SubDomainRank() const  = 0;

  //! Return number of time domains that the global problem is split into.
  virtual int NumTimeDomains() const  = 0;

  //! Return integer [0:numTimeDomains-1} corresponding to this time-domain's rank.
  virtual int TimeDomainRank() const  = 0;

  //! Return number of time steps, first step number, on time domain.
  virtual int NumTimeStepsOnDomain() const = 0;
  virtual int FirstTimeStepOnDomain() const = 0;

  //! Return total number of time steps.
  virtual int NumTimeSteps() const = 0;

  //! Reset total number of time steps, allowing time steps per domain to
  //  be set later than the MultiLevel parallelism is set up.
  virtual void ResetNumTimeSteps(int numTimeSteps) = 0;

};

} //namespace EpetraExt

#endif /* EPETRAEXT_MULTICOMM_H */



