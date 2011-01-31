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



