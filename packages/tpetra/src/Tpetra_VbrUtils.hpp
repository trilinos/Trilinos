// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

#ifndef TPETRA_VBRUTILS_HPP
#define TPETRA_VBRUTILS_HPP

#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Tpetra_BlockMap.hpp"
#include "Tpetra_SrcDistObject.hpp"
#include "Tpetra_Packable.hpp"
#include <Tpetra_Distributor.hpp> // avoid error C2027: use of undefined type 'Tpetra::Distributor' at (void) distor below
#include <map>

/** \file Tpetra_VbrUtils.hpp

  Utilities for VbrMatrix INTERNAL USE ONLY.
*/
namespace Tpetra {

  // Forward declaration
  class Distributor;
  

/**
  Utilities for VbrMatrix INTERNAL USE ONLY.
*/
namespace VbrUtils {

template<typename LocalOrdinal, typename Scalar>
struct BlkInfo {
  LocalOrdinal numPtRows;
  LocalOrdinal numPtCols;
  Teuchos::ArrayRCP<Scalar> blkEntry;
};

template<typename LocalOrdinal, typename GlobalOrdinal, typename Scalar>
struct VbrData {

  typedef typename std::map<GlobalOrdinal,BlkInfo<LocalOrdinal,Scalar> > RowGlobalCols;

  VbrData() {}
  ~VbrData() {}

  void zeroEntries () {
    typedef typename Teuchos::Array<RowGlobalCols>::size_type size_type;
    typedef typename Teuchos::ScalarTraits<Scalar> STS;

    const size_type numEntries = this->data.size ();
    for (size_type i = 0; i < numEntries; ++i) {
      RowGlobalCols& rgc = (this->data)[i];
      typename RowGlobalCols::iterator rgc_it = rgc.begin(), rgc_end = rgc.end();
      for ( ; rgc_it != rgc_end; ++rgc_it) {
	BlkInfo<LocalOrdinal,Scalar>& blk = rgc_it->second;
	std::fill (blk.blkEntry.begin (), blk.blkEntry.end (), STS::zero ());
      }
    }
  }

  std::map<GlobalOrdinal,LocalOrdinal> row_map;
  Teuchos::Array<RowGlobalCols> data;
};

template<typename LocalOrdinal, typename GlobalOrdinal, typename Scalar>
void
getGlobalBlockEntryViewNonConst (VbrData<LocalOrdinal,GlobalOrdinal,Scalar>& vbrdata,
				 GlobalOrdinal globalBlockRow,
				 GlobalOrdinal globalBlockCol,
				 LocalOrdinal& numPtRows,
				 LocalOrdinal& numPtCols,
				 Teuchos::ArrayRCP<Scalar>& blockEntry)
{
  typename std::map<GlobalOrdinal,LocalOrdinal>::iterator miter =
      vbrdata.row_map.find(globalBlockRow);

  LocalOrdinal localBlockRow = 0;

  if (miter == vbrdata.row_map.end()) {
    localBlockRow = vbrdata.data.size ();
    vbrdata.row_map.insert(std::make_pair(globalBlockRow,localBlockRow));
    vbrdata.data.resize (localBlockRow+1);
  }
  else {
    localBlockRow = miter->second;
  }

  typedef typename VbrData<LocalOrdinal,GlobalOrdinal,Scalar>::RowGlobalCols RowGlobalCols;
  RowGlobalCols& blkrow = (vbrdata.data)[localBlockRow];
  typename RowGlobalCols::iterator iter = blkrow.find(globalBlockCol);
  if (iter == blkrow.end()) {
    BlkInfo<LocalOrdinal,Scalar> blk;
    blk.numPtRows = numPtRows;
    blk.numPtCols = numPtCols;
    size_t blockSize = numPtRows*numPtCols;
    blk.blkEntry = Teuchos::arcp(new Scalar[blockSize], 0, blockSize);
    std::fill(blk.blkEntry.begin(), blk.blkEntry.end(), (Scalar) 0);
    blkrow.insert(iter, std::make_pair(globalBlockCol, blk));
    blockEntry = blk.blkEntry;
  }
  else {
    BlkInfo<LocalOrdinal,Scalar>& blk = iter->second;
    numPtRows = blk.numPtRows;
    numPtCols = blk.numPtCols;
    blockEntry = blk.blkEntry;
  }
}

template<typename LocalOrdinal, typename GlobalOrdinal, typename Scalar, class Node>
Teuchos::RCP<const Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node> >
createOverlapMap(VbrData<LocalOrdinal,GlobalOrdinal,Scalar>& vbrdata,
                const Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node>& rowmap)
{
  Teuchos::Array<GlobalOrdinal> importBlocks;
  Teuchos::Array<int> src_procs;

  typename std::map<GlobalOrdinal,LocalOrdinal>::iterator
    iter = vbrdata.row_map.begin(),
    iter_end = vbrdata.row_map.end();
  for(; iter != iter_end; ++iter) {
    importBlocks.push_back(iter->first);
  }

  Teuchos::Array<GlobalOrdinal> importFirstPointInBlocks(importBlocks.size());
  Teuchos::Array<LocalOrdinal>  importBlockSizes(importBlocks.size());
  rowmap.getRemoteBlockInfo(importBlocks(), importFirstPointInBlocks(), importBlockSizes());

  return Teuchos::rcp(
     new Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node>(
           Teuchos::OrdinalTraits<global_size_t>::invalid(),
           importBlocks(),
           importFirstPointInBlocks(),
           importBlockSizes(),
           rowmap.getPointMap()->getIndexBase(),
           rowmap.getPointMap()->getComm(),
           rowmap.getPointMap()->getNode()
         )
  );
}

template<class LocalOrdinal, class GlobalOrdinal, class Scalar, class Node>
class VbrDataDist : public Tpetra::SrcDistObject, 
		    public Tpetra::Packable<char, LocalOrdinal> {
private:
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > pointMap_;
  VbrData<LocalOrdinal,GlobalOrdinal,Scalar>& vbrData_;

public:
  VbrDataDist (VbrData<LocalOrdinal,GlobalOrdinal,Scalar>& vbrData,
	       const Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node>& importMap)
    : pointMap_ (convertBlockMapToPointMap (importMap)),
      vbrData_ (vbrData)
  {}

  ~VbrDataDist() {}

  Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > 
  getMap () const { return pointMap_; }

  virtual void 
  pack (const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
	Teuchos::Array<char>& exports,
	const Teuchos::ArrayView<size_t>& numPacketsPerLID,
	size_t& constantNumPackets,
	Distributor& distor) const
  {
    (void) distor; // forestall "unused argument" compiler warning
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::av_reinterpret_cast;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef typename ArrayView<const LO>::size_type Tsize_t;
    typedef typename VbrData<LO,GO,Scalar>::RowGlobalCols RowGlobalCols;

    //We will pack each row's data into the exports buffer as follows:
    //[num-block-cols,numPtRows,{list-of-blk-cols},{list-of-ptColsPerBlkCol},{all vals}]
    //so the length of the char exports buffer for a row is:
    //sizeof(LO)*(2+2*(num-block-cols)) + sizeof(Scalar)*numPtRows*sum(numPtCols_i)

    //For each row corresponding to exportLIDs, accumulate the size that it will
    //occupy in the exports buffer:
    size_t total_exports_size = 0;
    for (Tsize_t i = 0; i < exportLIDs.size (); ++i) {
      const RowGlobalCols& rgc = (this->vbrData_.data)[exportLIDs[i]];
      typename RowGlobalCols::const_iterator rgc_it = rgc.begin();
      typename RowGlobalCols::const_iterator rgc_end = rgc.end();
      size_t numScalars = 0;
      for ( ; rgc_it != rgc_end; ++rgc_it) {
        const BlkInfo<LO,Scalar>& blk = rgc_it->second;
        numScalars += as<size_t> (blk.numPtRows) * as<size_t> (blk.numPtCols);
      }
      const size_t numBlkCols = as<size_t> (rgc.size ());
      const size_t size_for_this_row = 
	sizeof (GO) * (2 + 2*numBlkCols)
	+ sizeof (Scalar) * numScalars;
      numPacketsPerLID[i] = size_for_this_row;
      total_exports_size += size_for_this_row;
    }

    exports.resize (total_exports_size);

    ArrayView<char> avIndsC, avValsC;
    ArrayView<Scalar> avVals;

    Array<GO> blkCols;
    Array<LO> ptColsPerBlkCol;
    Array<Scalar> blkEntries;

    size_t offset = 0;
    for (Tsize_t i = 0; i < exportLIDs.size (); ++i) {
      blkCols.clear();
      ptColsPerBlkCol.clear();
      blkEntries.clear();

      const RowGlobalCols& rgc = (this->vbrData_.data)[exportLIDs[i]];
      typename RowGlobalCols::const_iterator rgc_it = rgc.begin ();
      typename RowGlobalCols::const_iterator rgc_end = rgc.end ();
      LO numPtRows = 0;
      for ( ; rgc_it != rgc_end; ++rgc_it) {
        blkCols.push_back (rgc_it->first);
        const BlkInfo<LO,Scalar>& blk = rgc_it->second;
        numPtRows = blk.numPtRows; // should be the same for all cols
        ptColsPerBlkCol.push_back (blk.numPtCols);
        for (LO j = 0; j < blk.numPtRows * blk.numPtCols; ++j) {
          blkEntries.push_back (blk.blkEntry[j]);
        }
      }

      LO numBlkCols = blkCols.size();
      LO numScalars = blkEntries.size();

      LO num_chars_for_ordinals = (2*numBlkCols + 2) * sizeof (GO);
      //get export views
      avIndsC = exports(offset, num_chars_for_ordinals);
      avValsC = exports(offset+ num_chars_for_ordinals, numScalars*sizeof(Scalar));
      ArrayView<GO> avInds = av_reinterpret_cast<GO>(avIndsC);
      typename ArrayView<GO>::iterator ind_it = avInds.begin ();

      *ind_it++ = numBlkCols;
      *ind_it++ = numPtRows;
      for (Tsize_t j = 0; j < blkCols.size (); ++j) {
        *ind_it++ = blkCols[j];
      }
      for (Tsize_t j = 0; j < ptColsPerBlkCol.size (); ++j) {
        *ind_it++ = ptColsPerBlkCol[j];
      }

      avVals = av_reinterpret_cast<Scalar> (avValsC);
      std::copy (blkEntries.begin (), blkEntries.end (), avVals.begin ());

      const size_t size_for_this_row = 
	sizeof (GO) * (2 + 2*numBlkCols) + 
	sizeof (Scalar) * numScalars;
      offset += size_for_this_row;
    }
    constantNumPackets = 0;
  }

};

}//namespace VbrUtils
}//namespace Tpetra

#endif

