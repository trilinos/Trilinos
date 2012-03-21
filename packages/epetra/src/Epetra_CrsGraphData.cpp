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

#include "Epetra_CrsGraphData.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
//#include "Epetra_ConfigDefs.h" //DATA_DEBUG

//=============================================================================
Epetra_CrsGraphData::Epetra_CrsGraphData(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, bool StaticProfile)
  // maps
  : RowMap_(RowMap),
    ColMap_(RowMap),
    DomainMap_(RowMap),
    RangeMap_(RowMap),
    // importer & exporter
    Importer_(0),
    Exporter_(0),
    // booleans
    HaveColMap_(false),
    Filled_(false),
    Allocated_(false),
    // for non-static profile, we insert always into sorted lists, so the
    // graph will always be sorted. The same holds for the redundancies.
    Sorted_(!StaticProfile),
    StorageOptimized_(false),
    NoRedundancies_(!StaticProfile),
    IndicesAreGlobal_(false),
    IndicesAreLocal_(false),
    IndicesAreContiguous_(false),
    LowerTriangular_(true),
    UpperTriangular_(true),
    NoDiagonal_(true),
    GlobalConstantsComputed_(false),
    StaticProfile_(StaticProfile),
    SortGhostsAssociatedWithEachProcessor_(false),

    // ints
    IndexBase_(RowMap.IndexBase()),
    NumGlobalEntries_(0),
    NumGlobalBlockRows_(RowMap.NumGlobalElements()),
    NumGlobalBlockCols_(NumGlobalBlockRows_),
    NumGlobalBlockDiagonals_(0),
    NumMyEntries_(0),
    NumMyBlockRows_(RowMap.NumMyElements()),
    NumMyBlockCols_(NumMyBlockRows_),
    NumMyBlockDiagonals_(0),
    MaxRowDim_(RowMap.MaxElementSize()),
    MaxColDim_(MaxRowDim_),
    GlobalMaxRowDim_(RowMap.MaxElementSize()),
    GlobalMaxColDim_(GlobalMaxRowDim_),
    MaxNumNonzeros_(0),
    GlobalMaxNumNonzeros_(0),
    NumGlobalNonzeros_(0),
    NumGlobalRows_(RowMap.NumGlobalPoints()),
    NumGlobalCols_(NumGlobalRows_),
    NumGlobalDiagonals_(0),
    NumMyNonzeros_(0),
    NumMyRows_(RowMap.NumMyPoints()),
    NumMyCols_(NumMyRows_),
    NumMyDiagonals_(0),	
    MaxNumIndices_(0),
    GlobalMaxNumIndices_(0),
    NumTempColIndices_(0),
    NumAllocatedIndicesPerRow_(0),
    NumIndicesPerRow_(0),
    IndexOffset_(0),
    CV_(CV),
    data(0),
	LL_data(0)
{
  if(RowMap.GlobalIndicesInt() == false && RowMap.GlobalIndicesLongLong() == false)
    throw "Epetra_CrsGraphData::Epetra_CrsGraphData: cannot be called without any index type for RowMap";

  data = new IndexData<int>(NumMyBlockRows_, ! StaticProfile);
  LL_data = new IndexData<long long>(RowMap.GlobalIndicesLongLong() ? NumMyBlockRows_ : 0, ! StaticProfile);

  //cout << "--CRSGD created(rowmap ctr), addr: " << this << endl; //DATA_DEBUG
}

//=============================================================================
Epetra_CrsGraphData::Epetra_CrsGraphData(Epetra_DataAccess CV, 
					 const Epetra_BlockMap& RowMap, 
					 const Epetra_BlockMap& ColMap, bool StaticProfile)
  // maps
  : RowMap_(RowMap),
    ColMap_(ColMap),
    DomainMap_(ColMap),
    RangeMap_(RowMap),
    // importer & exporter
    Importer_(0),
    Exporter_(0),
    // booleans
    HaveColMap_(true),
    Filled_(false),
    Allocated_(false),
    Sorted_(!StaticProfile),
    StorageOptimized_(false),
    NoRedundancies_(!StaticProfile),
    IndicesAreGlobal_(false),
    IndicesAreLocal_(false),
    IndicesAreContiguous_(false),
    LowerTriangular_(true),
    UpperTriangular_(true),
    NoDiagonal_(true),
    GlobalConstantsComputed_(false),
    StaticProfile_(StaticProfile),
    SortGhostsAssociatedWithEachProcessor_(false),
    // ints
    IndexBase_(RowMap.IndexBase()),
    NumGlobalEntries_(0),
    NumGlobalBlockRows_(RowMap.NumGlobalElements()),
    NumGlobalBlockCols_(ColMap.NumGlobalElements()),
    NumGlobalBlockDiagonals_(0),
    NumMyEntries_(0),
    NumMyBlockRows_(RowMap.NumMyElements()),
    NumMyBlockCols_(ColMap.NumMyElements()),
    NumMyBlockDiagonals_(0),
    MaxRowDim_(RowMap.MaxElementSize()),
    MaxColDim_(ColMap.MaxElementSize()),
    GlobalMaxRowDim_(RowMap.MaxElementSize()),
    GlobalMaxColDim_(ColMap.MaxElementSize()),
    MaxNumNonzeros_(0),
    GlobalMaxNumNonzeros_(0),
    NumGlobalNonzeros_(0),
    NumGlobalRows_(RowMap.NumGlobalPoints()),
    NumGlobalCols_(ColMap.NumGlobalPoints()),
    NumGlobalDiagonals_(0),
    NumMyNonzeros_(0),
    NumMyRows_(RowMap.NumMyPoints()),
    NumMyCols_(ColMap.NumMyPoints()),
    NumMyDiagonals_(0),	
    MaxNumIndices_(0),
    GlobalMaxNumIndices_(0),
    NumTempColIndices_(0),
    NumAllocatedIndicesPerRow_(0),
    NumIndicesPerRow_(0),
    IndexOffset_(0),
    CV_(CV),
    data(0),
	LL_data(0)
{
  if(RowMap.GlobalIndicesInt() == false && RowMap.GlobalIndicesLongLong() == false)
    throw "Epetra_CrsGraphData::Epetra_CrsGraphData: cannot be called without any index type for RowMap";

  if(!RowMap.GlobalIndicesTypeMatch(ColMap))
    throw "Epetra_CrsGraphData::Epetra_CrsGraphData: cannot be called with different indices types for RowMap and ColMap";

  data = new IndexData<int>(NumMyBlockRows_, ! StaticProfile);
  LL_data = new IndexData<long long>(RowMap.GlobalIndicesLongLong() ? NumMyBlockRows_ : 0, ! StaticProfile);

  //cout << "--CRSGD created(rowmap&colmap ctr), addr: " << this << endl; //DATA_DEBUG
}

//=============================================================================
Epetra_CrsGraphData::~Epetra_CrsGraphData() {

  if(data->Indices_ != 0 && !StorageOptimized_) {
    for (int i=0; i<NumMyBlockRows_; i++) {
      data->Indices_[i] = 0;
    } 
    delete[] data->Indices_;
    data->Indices_ = 0;
  }

  if(LL_data->Indices_ != 0 && !StorageOptimized_) {
    for (int i=0; i<NumMyBlockRows_; i++) {
      LL_data->Indices_[i] = 0;
    } 
    delete[] LL_data->Indices_;
    LL_data->Indices_ = 0;
  }

  if (data->TempColIndices_ != 0) {
    delete [] data->TempColIndices_;
    data->TempColIndices_ = 0;
  }

  if (LL_data->TempColIndices_ != 0) {
    delete [] LL_data->TempColIndices_;
    LL_data->TempColIndices_ = 0;
  }

  delete data;
  delete LL_data;

  if(Importer_ != 0) {
    delete Importer_;
    Importer_ = 0;
  }
  if(Exporter_ != 0) {
    delete Exporter_;
    Importer_ = 0;
  }

  NumMyBlockRows_ = 0;	// are these needed?
  Filled_ = false;      // they're about to go out of scope, after all
  Allocated_ = false;

  //cout << "--CRSGD destroyed, addr: " << this << endl; //DATA_DEBUG
}

//==========================================================================
int Epetra_CrsGraphData::MakeImportExport() {
  // Create Import object for use by matrix classes.    This is only needed if ColMap and DomainMap are different
  if (!ColMap_.SameAs(DomainMap_)) {
    if (Importer_ != 0) {
      delete Importer_;
      Importer_ = 0;
    }
    Importer_ = new Epetra_Import(ColMap_, DomainMap_);
  }
  
  // Now see if we need to define an export map.  This is only needed if RowMap and RangeMap are different
  if (!RowMap_.SameAs(RangeMap_)) {
    if (Exporter_ != 0) {
      delete Exporter_;
      Exporter_ = 0;
    }
    Exporter_ = new Epetra_Export(RowMap_, RangeMap_); // Create Export object. 
  }
   
  return(0);
}

//==========================================================================
int Epetra_CrsGraphData::ReAllocateAndCast(char*& UserPtr, int& Length, const int IntPacketSizeTimesNumTrans) {
  if(IntPacketSizeTimesNumTrans > Length) {
    if(Length > 0) 
      delete[] UserPtr;
    Length = IntPacketSizeTimesNumTrans;
    int* newPtr = new int[Length];
    UserPtr = reinterpret_cast<char*> (newPtr);
  }
  return(0);
}

//==========================================================================
void Epetra_CrsGraphData::Print(ostream& os, int level) const {
  bool four_bit = (level >= 4);      // 4-bit = BlockMaps
  bool two_bit = ((level % 4) >= 2); // 2-bit = Indices
  bool one_bit = ((level % 2) == 1); // 1-bit = Everything else

  os << "\n***** CrsGraphData (output level " << level << ") *****" << endl;

  if(four_bit) {
    os << "RowMap_:\n" << RowMap_ << endl;
    os << "ColMap_:\n" << ColMap_ << endl;
    os << "DomainMap_:\n" << DomainMap_ << endl;
    os << "RangeMap_:\n" << RangeMap_ << endl;
  }
	
  if(one_bit) {
    os.width(26); os << "HaveColMap_: "              << HaveColMap_;
    os.width(25); os << "Filled_: "                  << Filled_;
    os.width(25); os << "Allocated_: "               << Allocated_;
    os.width(25); os << "Sorted_: "                  << Sorted_ << endl;
    os.width(26); os << "StorageOptimized_: "        << StorageOptimized_;
    os.width(25); os << "SortGhostsAssociatedWithEachProcessor_: " << SortGhostsAssociatedWithEachProcessor_;
    os.width(25); os << "NoRedundancies_: "          << NoRedundancies_;
    os.width(25); os << "IndicesAreGlobal_: "        << IndicesAreGlobal_;
    os.width(25); os << "IndicesAreLocal_: "         << IndicesAreLocal_ << endl;
    os.width(26); os << "IndicesAreContiguous_: "    << IndicesAreContiguous_;
    os.width(25); os << "LowerTriangular_: "         << LowerTriangular_;
    os.width(25); os << "UpperTriangular_: "         << UpperTriangular_;
    os.width(25); os << "NoDiagonal_: "              << NoDiagonal_ << endl;
    os.width(25); os << "GlobalConstantsComputed_: " << GlobalConstantsComputed_ << endl;
    os.width(25); os << "StaticProfile_: " << StaticProfile_ << endl << endl;
		
    os.width(10); os << "NGBR_: " << NumGlobalBlockRows_;
    os.width(10); os << "NGBC_: " << NumGlobalBlockCols_;
    os.width(10); os << "NGBD_: " << NumGlobalBlockDiagonals_;
    os.width(10); os << "NGE_: "  << NumGlobalEntries_;
    os.width(10); os << "NGR_: "  << NumGlobalRows_;
    os.width(10); os << "NGC_: "  << NumGlobalCols_;
    os.width(10); os << "NGD_: "  << NumGlobalDiagonals_;
    os.width(10); os << "NGN_: "  << NumGlobalNonzeros_;
    os.width(10); os << "IB_: "   << IndexBase_ << endl;
    os.width(10); os << "GMRD_: " << GlobalMaxRowDim_;
    os.width(11); os << "GMCD_: " << GlobalMaxColDim_;
    os.width(11); os << "GMNI_: " << GlobalMaxNumIndices_;
    os.width(11); os << "NMBR_: " << NumMyBlockRows_;
    os.width(10); os << "NMBC_: " << NumMyBlockCols_;
    os.width(10); os << "NMBD_: " << NumMyBlockDiagonals_;
    os.width(10); os << "NME_: "  << NumMyEntries_;
    os.width(10); os << "NMR_: "  << NumMyRows_;
    os.width(10); os << "CV_: " << CV_ << endl;
    os.width(10); os << "NMC_: "  << NumMyCols_;
    os.width(10); os << "NMD_: "  << NumMyDiagonals_;
    os.width(10); os << "NMN_: "  << NumMyNonzeros_;
    os.width(10); os << "MRD_: "  << MaxRowDim_;
    os.width(11); os << "MCD_: "  << MaxColDim_;
    os.width(11); os << "MNI_: "  << MaxNumIndices_;
    os.width(11); os << "MNN_: "  << MaxNumNonzeros_;
    os.width(11); os << "GMNN_: " << GlobalMaxNumNonzeros_;
    os.width(11); os << "RC: " << ReferenceCount() << endl << endl;
		
    os << "NIPR_: " << NumIndicesPerRow_ << endl;
    os << "NAIPR_: " << NumAllocatedIndicesPerRow_ << endl;
    os << "IndexOffset_: " << IndexOffset_ << endl;
	if(RowMap_.GlobalIndicesInt() || (RowMap_.GlobalIndicesLongLong() && IndicesAreLocal_))
      os << "All_Indices_: " << data->All_Indices_ << endl;

	if(RowMap_.GlobalIndicesLongLong() && IndicesAreGlobal_)
      os << "All_Indices_: " << LL_data->All_Indices_ << endl;
  }
		
  if(two_bit) {
	if(RowMap_.GlobalIndicesInt() || (RowMap_.GlobalIndicesLongLong() && IndicesAreLocal_))
	{
      os << "Indices_: " << data->Indices_ << endl;
	  if(data->Indices_ != 0) {
        for(int i = 0; i < NumMyBlockRows_; i++) {
	  os << "Indices_[" << i << "]: (" << data->Indices_[i] << ") ";
	  if(data->Indices_[i] != 0) {
	    for(int j = 0; j < NumAllocatedIndicesPerRow_[i]; j++)
	      os << data->Indices_[i][j] << " ";
	  }
	  os << endl;
      }
	 }
	}

	if(RowMap_.GlobalIndicesLongLong() && IndicesAreGlobal_)
	{
      os << "Indices_: " << LL_data->Indices_ << endl;
	  if(LL_data->Indices_ != 0) {
        for(int i = 0; i < NumMyBlockRows_; i++) {
	  os << "Indices_[" << i << "]: (" << LL_data->Indices_[i] << ") ";
	  if(LL_data->Indices_[i] != 0) {
	    for(int j = 0; j < NumAllocatedIndicesPerRow_[i]; j++)
	      os << LL_data->Indices_[i][j] << " ";
	  }
	  os << endl;
      }
	}
   }
  }
	
  os << "***** End CrsGraphData *****" << endl;
}


//==========================================================================
template<typename int_type>
void
Epetra_CrsGraphData::EntriesInOneRow<int_type>::AddEntry (const int_type Col)
{
  // first check the last element (or if line is still empty)
  if ( (entries_.size()==0) || ( entries_.back() < Col) )
    {
      entries_.push_back(Col);
      return;
    }
 
  // do a binary search to find the place where to insert:
  typename std::vector<int_type>::iterator it = std::lower_bound(entries_.begin(),
                entries_.end(),
                Col);
 
  // If this entry is a duplicate, exit immediately
  if (*it == Col)
    return;
 
  // Insert at the right place in the vector. Vector grows automatically to
  // fit elements. Always doubles its size.
  entries_.insert(it, Col);
}


//==========================================================================
template<typename int_type>
void
Epetra_CrsGraphData::EntriesInOneRow<int_type>::AddEntries (const int  numCols,
             const int_type *Indices)
{
  if (numCols == 0)
    return;

  // Check whether the indices are sorted. Can do more efficient then.
  bool indicesAreSorted = true;
  for (int i=1; i<numCols; ++i)
    if (Indices[i] <= Indices[i-1]) {
      indicesAreSorted = false;
      break;
    }

  if (indicesAreSorted && numCols > 3) {
    const int_type * curInput = &Indices[0];
    int_type col = *curInput;
    const int_type * endInput = &Indices[numCols];

    // easy case: list of entries is empty or all entries are smaller than
    // the ones to be inserted
    if (entries_.size() == 0 || entries_.back() < col)
    {
      entries_.insert(entries_.end(), &Indices[0], &Indices[numCols]);
      return;
    }

    // find a possible insertion point for the first entry. check whether
    // the first entry is a duplicate before actually doing something.
    typename std::vector<int_type>::iterator it = 
      std::lower_bound(entries_.begin(), entries_.end(), col);
    while (*it == col) {
      ++curInput;
      if (curInput == endInput)
        break;
      col = *curInput;

      // check the very next entry in the current array
      ++it;
      if (it == entries_.end())
        break;
      if (*it > col)
        break;
      if (*it == col)
        continue;

      // ok, it wasn't the very next one, do a binary search to find the
      // insert point
      it = std::lower_bound(it, entries_.end(), col);
      if (it == entries_.end())
        break;
    }

    // all input entries were duplicates.
    if (curInput == endInput)
      return;

    // Resize vector by just inserting the list at the correct point. Note
    // that the list will not yet be sorted, but rather have the insert
    // indices in the middle and the old indices from the list on the
    // end. Next we will have to merge the two lists.
    const int pos1 = (int) (it - entries_.begin());
    entries_.insert (it, curInput, endInput);
    it = entries_.begin() + pos1;

    // Now merge the two lists...
    typename std::vector<int_type>::iterator it2 = it + (endInput - curInput);

    // As long as there are indices both in the end of the entries list and
    // in the input list, always continue with the smaller index.
    while (curInput != endInput && it2 != entries_.end())
    {
      if (*curInput < *it2)
        *it++ = *curInput++;
      else if (*curInput == *it2)
      {
        *it++ = *it2++;
        ++curInput;
      }
      else
        *it++ = *it2++;
    }
    // in case there are indices left in the input list or in the end of
    // entries, we have to use the entries that are left. Only one of the
    // two while loops will actually be entered.
    while (curInput != endInput)
      *it++ = *curInput++;

    while (it2 != entries_.end())
      *it++ = *it2++;

    // resize and return
    const int new_size = (int) (it - entries_.begin());
    entries_.resize (new_size);
    return;
  }

  // for unsorted or just a few indices, go to the one-by-one entry
  // function.
  for (int i=0; i<numCols; ++i)
    AddEntry(Indices[i]);
}

// explicit instantiation.
template void Epetra_CrsGraphData::EntriesInOneRow<int      >::AddEntry(const int Col);
template void Epetra_CrsGraphData::EntriesInOneRow<long long>::AddEntry(const long long Col);

template void Epetra_CrsGraphData::EntriesInOneRow<int      >::AddEntries(const int numCols, const int *Indices);
template void Epetra_CrsGraphData::EntriesInOneRow<long long>::AddEntries(const int numCols, const long long *Indices);
