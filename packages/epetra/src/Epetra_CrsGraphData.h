
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

#ifndef EPETRA_CRSGRAPHDATA_H
#define EPETRA_CRSGRAPHDATA_H

#include "Epetra_Data.h"
#include "Epetra_DataAccess.h"
#include "Epetra_BlockMap.h"
#include "Epetra_IntSerialDenseVector.h"
class Epetra_Import;
class Epetra_Export;

//! Epetra_CrsGraphData:  The Epetra CrsGraph Data Class.
/*! The Epetra_CrsGraphData class is an implementation detail of Epetra_CrsGraph.
    It is reference-counted, and can be shared by multiple Epetra_CrsGraph instances. 
    It derives from Epetra_Data, and inherits reference-counting from it.
*/

class Epetra_CrsGraphData : public Epetra_Data {
  friend class Epetra_CrsGraph;

 private:

  //@{ \name Constructor/Destructor Methods

  //! Epetra_CrsGraphData Default Constructor.
  Epetra_CrsGraphData(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap);

  //! Epetra_CrsGraphData Constructor (user provided ColMap).
  Epetra_CrsGraphData(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, const Epetra_BlockMap& ColMap);

	//! Epetra_CrsGraphData copy constructor (not defined).
  Epetra_CrsGraphData(const Epetra_CrsGraphData& CrsGraphData);

  //! Epetra_CrsGraphData Destructor.
  ~Epetra_CrsGraphData();

  //@}

	//! Outputs state of almost all data members. (primarily used for testing purposes).
	/*! Output level: Uses same scheme as chmod. 4-bit = BlockMaps, 2-bit = Indices, 1-bit = Everything else.
		Default paramenter sets it to 3, which is everything but the BlockMaps. Commonly used options:
		1 = Everything except the BlockMaps & Indices_
		2 = Just Indices_
		3 = Everything except the BlockMaps
	*/
  void Print(ostream& os, int level = 3) const;

	//! Epetra_CrsGraphData assignment operator (not defined)
	Epetra_CrsGraphData& operator=(const Epetra_CrsGraphData& CrsGraphData);

	//@{ \name Helper methods called in CrsGraph. Mainly memory allocations and deallocations.

	//! called by FillComplete (and TransformToLocal)
	int MakeImportExport();

	//! called by PackAndPrepare
	int ReAllocateAndCast(char*& UserPtr, int& Length, const int IntPacketSizeTimesNumTrans);

	//! returns the IndicesSidekick array used by Indices() and the [] operators.
	int** Sidekick() const;
	void UpdateSidekick() const; // used internally by Graph to keep the Sidekick array updated when Indices_ changes.

	//@}

  // Defined by CrsGraph::FillComplete and related
	const Epetra_BlockMap RowMap_;
	Epetra_BlockMap ColMap_;
	Epetra_BlockMap DomainMap_;
	Epetra_BlockMap RangeMap_;

	const Epetra_Import* Importer_;
  const Epetra_Export* Exporter_;

	bool HaveColMap_;
  bool Filled_;
  bool Allocated_;
  bool Sorted_;
  bool StorageOptimized_;
  bool NoRedundancies_;
  bool IndicesAreGlobal_;
  bool IndicesAreLocal_;
  bool IndicesAreContiguous_;
  bool LowerTriangular_;
  bool UpperTriangular_;
  bool NoDiagonal_;
  bool GlobalConstantsComputed_;

  int IndexBase_;

  int NumGlobalEntries_;
  int NumGlobalBlockRows_;
  int NumGlobalBlockCols_;
  int NumGlobalBlockDiagonals_;
  int NumMyEntries_;
  int NumMyBlockRows_;
  int NumMyBlockCols_;
  int NumMyBlockDiagonals_;
  
  int MaxRowDim_;
  int MaxColDim_;
  int GlobalMaxRowDim_;
  int GlobalMaxColDim_;
  int MaxNumNonzeros_;
  int GlobalMaxNumNonzeros_;
  
  int NumGlobalNonzeros_;
  int NumGlobalRows_;
  int NumGlobalCols_;
  int NumGlobalDiagonals_;
  int NumMyNonzeros_;
  int NumMyRows_;
  int NumMyCols_;
  int NumMyDiagonals_;

	int MaxNumIndices_;
	int GlobalMaxNumIndices_;

	Epetra_IntSerialDenseVector* Indices_;
	mutable int** IndicesSidekick_;
	Epetra_IntSerialDenseVector NumAllocatedIndicesPerRow_;
	Epetra_IntSerialDenseVector NumIndicesPerRow_;
	Epetra_IntSerialDenseVector All_Indices_;
	Epetra_DataAccess CV_;

};

#endif /* EPETRA_CRSGRAPHDATA_H */
