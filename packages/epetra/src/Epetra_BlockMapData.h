
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

#ifndef _EPETRA_BLOCKMAPDATA_H_
#define _EPETRA_BLOCKMAPDATA_H_

#include "Epetra_Data.h"
#include "Epetra_IntSerialDenseVector.h"

class Epetra_Comm;
class Epetra_Directory;
class Epetra_HashTable;

//! Epetra_BlockMapData:  The Epetra BlockMap Data Class.
/*! The Epetra_BlockMapData class is an implementation detail of Epetra_BlockMap.
    It is reference-counted, and can be shared by multiple Epetra_BlockMap instances. 
    It derives from Epetra_Data, and inherits reference-counting from it.
*/

class Epetra_BlockMapData : public Epetra_Data {
  friend class Epetra_BlockMap;

 private:

  //@{ \name Constructor/Destructor Methods

  //! Epetra_BlockMapData Default Constructor.
  Epetra_BlockMapData(int NumGlobalElements, int ElementSize, int IndexBase, const Epetra_Comm & Comm);

  //! Epetra_BlockMapData Destructor.
  ~Epetra_BlockMapData();

  //@}

  const Epetra_Comm * Comm_;
  Epetra_Directory * Directory_;
  
  Epetra_IntSerialDenseVector LID_;
  Epetra_IntSerialDenseVector MyGlobalElements_;
  Epetra_IntSerialDenseVector FirstPointInElementList_;
  Epetra_IntSerialDenseVector ElementSizeList_;
  Epetra_IntSerialDenseVector PointToElementList_;
  
  int NumGlobalElements_;
  int NumMyElements_;
  int IndexBase_;
  int ElementSize_;
  int MinMyElementSize_;
  int MaxMyElementSize_;
  int MinElementSize_;
  int MaxElementSize_;
  int MinAllGID_;
  int MaxAllGID_;
  int MinMyGID_;
  int MaxMyGID_;
  int MinLID_;
  int MaxLID_;
  int NumGlobalPoints_;
  int NumMyPoints_;
  
  bool ConstantElementSize_;
  bool LinearMap_;
  bool DistributedGlobal_;

  int LastContiguousGIDLoc_;
  Epetra_HashTable * LIDHash_;

	// these are intentionally declared but not defined. See Epetra Developer's Guide for details.
  Epetra_BlockMapData(const Epetra_BlockMapData & BlockMapData);
	Epetra_BlockMapData& operator=(const Epetra_BlockMapData & BlockMapData);

};
#endif /* _EPETRA_BLOCKMAPDATA_H_ */
