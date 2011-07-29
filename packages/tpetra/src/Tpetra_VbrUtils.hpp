//@HEADER
// ************************************************************************
// 
//               Tpetra: Templated Linear Algebra Services Package 
//                 Copyright (2008) Sandia Corporation
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

#ifndef TPETRA_VBRUTILS_HPP
#define TPETRA_VBRUTILS_HPP

#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Tpetra_BlockMap.hpp"
#include <map>

/** \file Tpetra_VbrUtils.hpp

  Utilities for VbrMatrix INTERNAL USE ONLY.
*/
namespace Tpetra {
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

  VbrData()
   : row_map(),
     data(Teuchos::rcp(new Teuchos::Array<RowGlobalCols>))
  {}

  ~VbrData(){}

  std::map<GlobalOrdinal,LocalOrdinal> row_map;
  Teuchos::RCP<Teuchos::Array<RowGlobalCols> > data;
};

template<typename LocalOrdinal, typename GlobalOrdinal, typename Scalar>
inline
void zeroEntries(VbrData<LocalOrdinal,GlobalOrdinal,Scalar>& vbrdata)
{
  typedef typename VbrData<LocalOrdinal,GlobalOrdinal,Scalar>::RowGlobalCols RowGlobalCols;

  typedef typename Teuchos::Array<RowGlobalCols>::size_type Tsize_t;

  for(Tsize_t i=0; i<vbrdata.data->size(); ++i) {
    RowGlobalCols& rgc = (*vbrdata.data)[i];
    typename RowGlobalCols::iterator
      rgc_it = rgc.begin(), rgc_end = rgc.end();
    for(; rgc_it != rgc_end; ++rgc_it) {
      BlkInfo<LocalOrdinal,Scalar>& blk = rgc_it->second;
      std::fill(blk.blkEntry.begin(), blk.blkEntry.end(), 0.0);
    }    
  }
}

template<typename LocalOrdinal, typename GlobalOrdinal, typename Scalar>
inline
void getGlobalBlockEntryViewNonConst(
                   VbrData<LocalOrdinal,GlobalOrdinal,Scalar>& vbrdata,
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
    localBlockRow = vbrdata.data->size();
    vbrdata.row_map.insert(std::make_pair(globalBlockRow,localBlockRow));
    vbrdata.data->resize(localBlockRow+1);
  }
  else {
    localBlockRow = miter->second;
  }

  typedef typename VbrData<LocalOrdinal,GlobalOrdinal,Scalar>::RowGlobalCols RowGlobalCols;
  RowGlobalCols& blkrow = (*vbrdata.data)[localBlockRow];
  typename RowGlobalCols::iterator iter = blkrow.find(globalBlockCol);
  if (iter == blkrow.end()) {
    BlkInfo<LocalOrdinal,Scalar> blk;
    blk.numPtRows = numPtRows;
    blk.numPtCols = numPtCols;
    size_t blockSize = numPtRows*numPtCols;
    blk.blkEntry = Teuchos::arcp(new Scalar[blockSize], 0, blockSize);
    std::fill(blk.blkEntry.begin(), blk.blkEntry.end(), 0);
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

template<typename LocalOrdinal,typename GlobalOrdinal,typename Scalar,typename Node>
struct VbrDataDist : public Tpetra::DistObject<char,LocalOrdinal,GlobalOrdinal,Node> {
  VbrDataDist(VbrData<LocalOrdinal,GlobalOrdinal,Scalar>& vbr_data,
              const Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node>& importmap,
              const Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node>& rowmap)
   : Tpetra::DistObject<char,LocalOrdinal,GlobalOrdinal,Node>(
        convertBlockMapToPointMap(importmap)),
     vbrdata(vbr_data)
  {}

  ~VbrDataDist(){}

  VbrData<LocalOrdinal,GlobalOrdinal,Scalar>& vbrdata;
  Teuchos::RCP<Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node> > importMap;

  bool checkSizes(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node>& source)
  {
    bool ok = this->getMap()->getMinAllGlobalIndex() <= source.getMap()->getMinAllGlobalIndex();
    ok = ok && this->getMap()->getMaxAllGlobalIndex() >= source.getMap()->getMaxAllGlobalIndex();
    return ok;
  }

  void copyAndPermute(
   const DistObject<char, LocalOrdinal, GlobalOrdinal, Node>& source,
   size_t numSameIDs,
   const Teuchos::ArrayView<const LocalOrdinal>& permuteToLIDs,
   const Teuchos::ArrayView<const LocalOrdinal>& permuteFromLIDs)
  {
    throw std::runtime_error("VbrDataDist hasn't implemented copyAndPermute!!");
  }

  void packAndPrepare(
     const DistObject<char, LocalOrdinal, GlobalOrdinal, Node>& source,
     const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
     Teuchos::Array<char>& exports,
     const Teuchos::ArrayView<size_t>& numPacketsPerLID,
     size_t& constantNumPackets,
     Distributor& distor)
  {
    const VbrDataDist<LocalOrdinal,GlobalOrdinal,Scalar,Node>* vdd = dynamic_cast<const VbrDataDist<LocalOrdinal,GlobalOrdinal,Scalar,Node>*>(&source);
    if (vdd == NULL) {
      throw std::runtime_error("VbrDataDist::packAndPrepare ERROR, dynamic_cast failed.");
    }
    typedef typename Teuchos::ArrayView<const LocalOrdinal>::size_type Tsize_t;
    typedef typename VbrData<LocalOrdinal,GlobalOrdinal,Scalar>::RowGlobalCols RowGlobalCols;
  
    //We will pack each row's data into the exports buffer as follows:
    //[num-block-cols,numPtRows,{list-of-blk-cols},{list-of-ptColsPerBlkCol},{all vals}]
    //so the length of the char exports buffer for a row is:
    //sizeof(LocalOrdinal)*(2+2*(num-block-cols)) + sizeof(Scalar)*numPtRows*sum(numPtCols_i)
  
    //For each row corresponding to exportLIDs, accumulate the size that it will
    //occupy in the exports buffer:
    size_t total_exports_size = 0;
    for(Tsize_t i=0; i<exportLIDs.size(); ++i) {
      LocalOrdinal numBlkCols = 0;
      LocalOrdinal numScalars = 0;
  
      LocalOrdinal localIndex = exportLIDs[i];
      const RowGlobalCols& rgc = (*vdd->vbrdata.data)[localIndex];
      numBlkCols = rgc.size();
      typename RowGlobalCols::const_iterator rgc_it = rgc.begin(), rgc_end = rgc.end();
      for(; rgc_it!=rgc_end; ++rgc_it) {
        const BlkInfo<LocalOrdinal,Scalar>& blk = rgc_it->second;
        numScalars += blk.numPtRows*blk.numPtCols;
      }
  
      size_t size_for_this_row = sizeof(GlobalOrdinal)*(2+2*numBlkCols)
                   + sizeof(Scalar)*numScalars;
      numPacketsPerLID[i] = size_for_this_row;
      total_exports_size += size_for_this_row;
    }
  
    exports.resize(total_exports_size);
  
    ArrayView<char> avIndsC, avValsC;
    ArrayView<Scalar> avVals;
  
    Teuchos::Array<GlobalOrdinal> blkCols;
    Teuchos::Array<LocalOrdinal> ptColsPerBlkCol;
    Teuchos::Array<Scalar> blkEntries;
  
    size_t offset = 0;
    for(Tsize_t i=0; i<exportLIDs.size(); ++i) {
      blkCols.clear();
      ptColsPerBlkCol.clear();
      blkEntries.clear();
  
      LocalOrdinal numPtRows = 0;
      LocalOrdinal localIndex = exportLIDs[i];
      const RowGlobalCols& rgc = (*vdd->vbrdata.data)[localIndex];
      typename RowGlobalCols::const_iterator rgc_it = rgc.begin(), rgc_end = rgc.end();
      for(; rgc_it!=rgc_end; ++rgc_it) {
        blkCols.push_back(rgc_it->first);
        const BlkInfo<LocalOrdinal,Scalar>& blk = rgc_it->second;
        numPtRows = blk.numPtRows;//should be the same for all cols
        ptColsPerBlkCol.push_back(blk.numPtCols);
        for(LocalOrdinal j=0; j<blk.numPtRows*blk.numPtCols; ++j) {
          blkEntries.push_back(blk.blkEntry[j]);
        }
      }
  
      LocalOrdinal numBlkCols = blkCols.size();
      LocalOrdinal numScalars = blkEntries.size();
  
      LocalOrdinal num_chars_for_ordinals = (2*numBlkCols+2)*sizeof(GlobalOrdinal);
      //get export views
      avIndsC = exports(offset, num_chars_for_ordinals);
      avValsC = exports(offset+ num_chars_for_ordinals, numScalars*sizeof(Scalar));
      ArrayView<GlobalOrdinal> avInds = av_reinterpret_cast<GlobalOrdinal>(avIndsC);
      typename ArrayView<GlobalOrdinal>::iterator ind_it = avInds.begin();

      *ind_it++ = numBlkCols;
      *ind_it++ = numPtRows;
      for(Tsize_t j=0; j<blkCols.size(); ++j) {
        *ind_it++ = blkCols[j];
      }
      for(Tsize_t j=0; j<ptColsPerBlkCol.size(); ++j) {
        *ind_it++ = ptColsPerBlkCol[j];
      }

      avVals = av_reinterpret_cast<Scalar>(avValsC);
      std::copy(blkEntries.begin(), blkEntries.end(), avVals.begin());

      size_t size_for_this_row = sizeof(GlobalOrdinal)*(2+2*numBlkCols)
                     + sizeof(Scalar)*numScalars;
      offset += size_for_this_row;
    }
    constantNumPackets = 0;
  }

  void unpackAndCombine(
      const Teuchos::ArrayView<const LocalOrdinal>& importLIDs,
      const Teuchos::ArrayView<const char>& imports,
      const Teuchos::ArrayView<size_t>& numPacketsPerLID,
      size_t constantNumPackets,
      Distributor& distor,
      CombineMode CM)
  {
    throw std::runtime_error("VbrDataDist hasn't implemented unpackAndCombine!!");
  }

};

}//namespace VbrUtils
}//namespace Tpetra

#endif

