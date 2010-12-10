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

#include "EpetraExt_MultiSerialComm.h" 

namespace EpetraExt {

MultiSerialComm::MultiSerialComm(int numTimeSteps_) :
        myComm(Teuchos::rcp(new Epetra_SerialComm())),
        subComm(0)
{
  ResetNumTimeSteps(numTimeSteps_);

  // Construct second epetra communicators
  subComm = new Epetra_SerialComm();
}

//Copy Constructor
MultiSerialComm::MultiSerialComm(const MultiSerialComm &MSC ) :
        myComm(Teuchos::rcp(new Epetra_SerialComm())),
        subComm(new Epetra_SerialComm(*(MSC.subComm)))
{
  numTimeSteps = MSC.numTimeSteps;
}

MultiSerialComm::~MultiSerialComm()
{
  delete subComm;
}

void MultiSerialComm::ResetNumTimeSteps(int numTimeSteps_)
{
  numTimeSteps = numTimeSteps_;
}

} //namespace EpetraExt
