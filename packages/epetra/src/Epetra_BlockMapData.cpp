
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include "Epetra_BlockMapData.h"
#include "Epetra_HashTable.h"
#include "Epetra_Comm.h"
#include "Epetra_Directory.h"
//#include "Epetra_ConfigDefs.h" //DATA_DEBUG
// Use the new LID hash table approach by default
#define EPETRA_BLOCKMAP_NEW_LID

//=============================================================================
Epetra_BlockMapData::Epetra_BlockMapData(int NumGlobalElements, int ElementSize, int IndexBase, const Epetra_Comm & Comm) 
  : Comm_(Comm.Clone()),
    Directory_(0),
    NumGlobalElements_(NumGlobalElements),
    IndexBase_(IndexBase),
    ElementSize_(ElementSize),
    LastContiguousGID_(0),
    LastContiguousGIDLoc_(0),
    LIDHash_(0)
{
  //cout << "--BMD created, addr: " << this << endl; //DATA_DEBUG
}

//=============================================================================
Epetra_BlockMapData::~Epetra_BlockMapData()
{
  if(LIDHash_ != 0) {
    delete LIDHash_;
    LIDHash_ = 0;
  }

  if (Directory_ != 0) {
    delete Directory_;
    Directory_ = 0;
  }

  if(Comm_ != 0) {
    delete Comm_;
    Comm_ = 0;
  }
  //cout << "--BMD destroyed, addr: " << this << endl; //DATA_DEBUG
}
