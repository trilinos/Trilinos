
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

#include "Epetra_BlockMapData.h"
#include "Epetra_HashTable.h"
#include "Epetra_Directory.h"
#include "Epetra_Comm.h"
//#include "Epetra_ConfigDefs.h" //DATA_DEBUG

//=============================================================================
Epetra_BlockMapData::Epetra_BlockMapData(int NumGlobalElements, int ElementSize, int IndexBase, const Epetra_Comm & Comm) : 
  Comm_(Comm.Clone()),
  Directory_(0),
  NumGlobalElements_(NumGlobalElements),
  IndexBase_(IndexBase),
  ElementSize_(ElementSize),
  LastContiguousGIDLoc_(0),
  LIDHash_(0)
{
	//cout << "--BMD created, addr: " << (int) this << endl; //DATA_DEBUG
}

//=============================================================================
Epetra_BlockMapData::~Epetra_BlockMapData() {
  if(LIDHash_ != 0) {
    delete LIDHash_;
    LIDHash_ = 0;
  }
  if(Directory_ != 0) {
    delete Directory_;
    Directory_ = 0;
  }
	if(Comm_ != 0) {
		delete Comm_;
		Comm_ = 0;
	}
	//cout << "--BMD destroyed, addr: " << (int) this << endl; //DATA_DEBUG
}
