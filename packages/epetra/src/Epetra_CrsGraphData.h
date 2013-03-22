/*
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
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
// ************************************************************************
//@HEADER
*/

#ifndef EPETRA_CRSGRAPHDATA_H
#define EPETRA_CRSGRAPHDATA_H

#include "Epetra_ConfigDefs.h"
#include "Epetra_Data.h"
#include "Epetra_DataAccess.h"
#include "Epetra_BlockMap.h"
#include "Epetra_IntSerialDenseVector.h"

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
#include "Epetra_LongLongSerialDenseVector.h"
#endif

// include STL vector
#include <vector>
class Epetra_Import;
class Epetra_Export;

//! Epetra_CrsGraphData:  The Epetra CrsGraph Data Class.
/*! The Epetra_CrsGraphData class is an implementation detail of Epetra_CrsGraph.
    It is reference-counted, and can be shared by multiple Epetra_CrsGraph instances. 
    It derives from Epetra_Data, and inherits reference-counting from it.
*/

class EPETRA_LIB_DLL_EXPORT Epetra_CrsGraphData : public Epetra_Data {
  friend class Epetra_CrsGraph;
  friend class Epetra_FECrsGraph;
  friend class Epetra_CrsMatrix;
 private:

  //! @name Constructor/Destructor Methods
  //@{ 

  //! Epetra_CrsGraphData Default Constructor.
  Epetra_CrsGraphData(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, bool StaticProfile);

  //! Epetra_CrsGraphData Constructor (user provided ColMap).
  Epetra_CrsGraphData(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, const Epetra_BlockMap& ColMap, bool StaticProfile);

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
  
  //! @name Helper methods called in CrsGraph. Mainly memory allocations and deallocations.
  //@{ 
                                      /**
                                       * Store some data for each row
                                       * describing which entries of
                                       * this row are nonzero. Data is
                                       * stored sorted in the @p
                                       * entries std::vector which is
                                       * kept sorted and without
                                       * duplicates.  The vector of
                                       * indices per row is dynamically
                                       * growing upon insertion.
                                       */
  template<typename int_type>
  struct EntriesInOneRow
  {
    public:
                /**
           * Storage for the column indices of
           * this row. This array is always
           * kept sorted.
           */
      std::vector<int_type> entries_;
       
                /**
           * Add the given column number to
           * this line.
           */
      void AddEntry (const int_type col_num);

              /**
         * Add many entries to one row.
         */
      void AddEntries (const int  n_cols,
          const int_type *col_nums);
  };
  
  //! called by FillComplete (and TransformToLocal)
  int MakeImportExport();
  
  //! called by PackAndPrepare
  int ReAllocateAndCast(char*& UserPtr, int& Length, const int IntPacketSizeTimesNumTrans);
  
  //@}
  
  // Defined by CrsGraph::FillComplete and related
  Epetra_BlockMap RowMap_;
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
  bool StaticProfile_;
  bool SortGhostsAssociatedWithEachProcessor_;

  long long IndexBase_;

  long long NumGlobalEntries_;
  long long NumGlobalBlockRows_;
  long long NumGlobalBlockCols_;
  long long NumGlobalBlockDiagonals_;
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
  
  long long NumGlobalNonzeros_;
  long long NumGlobalRows_;
  long long NumGlobalCols_;
  long long NumGlobalDiagonals_;
  int NumMyNonzeros_;
  int NumMyRows_;
  int NumMyCols_;
  int NumMyDiagonals_;

  int MaxNumIndices_;
  int GlobalMaxNumIndices_;

  int NumTempColIndices_;
  Epetra_IntSerialDenseVector NumAllocatedIndicesPerRow_;
  Epetra_IntSerialDenseVector NumIndicesPerRow_;
  Epetra_IntSerialDenseVector IndexOffset_;
  Epetra_DataAccess CV_;

  template<typename int_type>
  struct IndexData;

  IndexData<int>* data;
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  IndexData<long long>* LL_data;
#endif

  template<typename int_type>
  IndexData<int_type>& Data();
};

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES

template<>
struct Epetra_CrsGraphData::IndexData<long long>
{
  long long** Indices_;
  std::vector< EntriesInOneRow<long long> > SortedEntries_;
  long long* TempColIndices_;
  Epetra_LongLongSerialDenseVector All_Indices_;

  IndexData(int NumMyBlockRows, bool AllocSorted)
    :
    Indices_(NULL),
    SortedEntries_(),
    TempColIndices_(NULL),
    All_Indices_(0)
  {
    Allocate(NumMyBlockRows, AllocSorted);
  }

  void Allocate(int NumMyBlockRows, bool AllocSorted)
  {
    Deallocate();

    if(AllocSorted)
      SortedEntries_.resize(NumMyBlockRows);
    if(NumMyBlockRows > 0)
      Indices_ = new long long *[NumMyBlockRows];
  }

  void Deallocate()
  {
    delete[] Indices_;
    Indices_ = 0;

    std::vector< EntriesInOneRow<long long> > empty;
    SortedEntries_.swap(empty);

    delete [] TempColIndices_;
    TempColIndices_ = 0;

    All_Indices_.Resize(0);
  }
};

#endif

template<>
struct Epetra_CrsGraphData::IndexData<int>
{
  int** Indices_;
  std::vector< EntriesInOneRow<int> > SortedEntries_;
  int* TempColIndices_;
  Epetra_IntSerialDenseVector All_Indices_;

  IndexData(int NumMyBlockRows, bool AllocSorted)
    :
    Indices_(NULL),
    SortedEntries_(),
    TempColIndices_(NULL),
    All_Indices_(0)
  {
    Allocate(NumMyBlockRows, AllocSorted);
  }

  void Allocate(int NumMyBlockRows, bool AllocSorted)
  {
    Deallocate();

    if(AllocSorted)
      SortedEntries_.resize(NumMyBlockRows);

    if(NumMyBlockRows > 0)
      Indices_ = new int *[NumMyBlockRows];
  }

  void Deallocate()
  {
    delete[] Indices_;
    Indices_ = 0;

    std::vector< EntriesInOneRow<int> > empty;
    SortedEntries_.swap(empty);

    delete [] TempColIndices_;
    TempColIndices_ = 0;

    All_Indices_.Resize(0);
  }
};

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
template<>
inline Epetra_CrsGraphData::IndexData<long long>& Epetra_CrsGraphData::Data<long long>()
{
  if(RowMap_.GlobalIndicesLongLong() && !IndicesAreLocal_)
    return *LL_data;
  else
    throw "Epetra_CrsGraphData::Data<long long>: Map indices not long long or are local";
}
#endif

template<>
inline Epetra_CrsGraphData::IndexData<int>& Epetra_CrsGraphData::Data<int>()
{
  if(RowMap_.GlobalIndicesInt() || (RowMap_.GlobalIndicesLongLong() && !IndicesAreGlobal_))
    return *data;
  else
    throw "Epetra_CrsGraphData::Data<int>: Map indices not int or are global long long";
}


#endif /* EPETRA_CRSGRAPHDATA_H */
