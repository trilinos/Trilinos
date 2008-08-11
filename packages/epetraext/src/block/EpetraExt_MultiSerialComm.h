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
  Epetra_Comm& SubDomainComm() const {return *subComm;}

  //! Reset total number of time steps, allowing time steps per domain to
  //  be set later than the MultiLevel parallelism is set up.
  void ResetNumTimeSteps(int numTimeSteps);

 protected:

  Epetra_SerialComm* subComm; 
};

} //namespace EpetraExt

#endif /* EPETRAEXT_MULTISERIALCOMM_H */
