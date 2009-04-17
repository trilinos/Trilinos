//@HEADER
// ************************************************************************
// 
//               Tpetra: Templated Linear Algebra Services Package 
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

#ifndef TPETRA_CRSGRAPHDATA_HPP
#define TPETRA_CRSGRAPHDATA_HPP

#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"


namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration of CrsGraph, needed to prevent circular inclusions
  template<class LocalOrdinal, class GlobalOrdinal> class CrsGraph;
#endif

  template <class LocalOrdinal, class GlobalOrdinal>
  class CrsGraphData {
    public:
    friend class CrsGraph<LocalOrdinal,GlobalOrdinal>;

    ~CrsGraphData();

    private:

    // useful typedefs
    typedef Teuchos::OrdinalTraits<LocalOrdinal>    LOT;
    typedef Teuchos::OrdinalTraits<GlobalOrdinal>   GOT;

    CrsGraphData(const Map<LocalOrdinal,GlobalOrdinal> & RowMap);
    CrsGraphData(const Map<LocalOrdinal,GlobalOrdinal> & RowMap, const Map<LocalOrdinal,GlobalOrdinal> & ColMap);
    CrsGraphData(const CrsGraphData<LocalOrdinal,GlobalOrdinal> & CrsGraphData);  // disabled
    CrsGraphData& operator=(const CrsGraphData<LocalOrdinal,GlobalOrdinal> &rhs); // disabled

    Teuchos::RCP<const Teuchos::Comm<int> > comm_;
    Map<LocalOrdinal,GlobalOrdinal> rowMap_, colMap_;
    Map<LocalOrdinal,GlobalOrdinal> rangeMap_, domainMap_;
    GlobalOrdinal     numGlobalEntries_, numGlobalDiags_, globalMaxNumRowEntries_;
    Teuchos_Ordinal    numLocalEntries_,  numLocalDiags_,  localMaxNumEntries_;
    Teuchos::ArrayRCP<LocalOrdinal> rowNumToAlloc_; // exists only if indices are not allocated

    //
    // Unoptimized structure
    // these are allocated if storage is not optimized or allocation is not static
    //
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<GlobalOrdinal> > colGInds_;    // allocated only if indices are global
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<LocalOrdinal> >  colLInds_;    // allocated only if indices are local
    Teuchos::ArrayRCP<LocalOrdinal> rowNumEntries_;

    //
    // Optimized structure
    // structure used if allocation is static or after optimizeStorage()
    //
    Teuchos_Ordinal totalNumAllocated_;
    Teuchos::ArrayRCP<GlobalOrdinal> contigColGInds_;     // allocated only if indices are global
    Teuchos::ArrayRCP<LocalOrdinal>  contigColLInds_;     // allocated only if indices are local
    /* colIndsPtrs[j] is the begin() iterator from an ArrayView of 
       contigColInds_ corresponding to the proper row, of the appropriate length.
       In a debug build, it is an ArrayRCP, which does bounds checking. in an optimized
       build, it is a C pointer. colIndsPtrs is allocated to numLocalRows()+1; the span of the jth row begins with
       colIndsPtrs[j] and ends before colIndsPtrs[j+1] */
    Teuchos::ArrayRCP<typename Teuchos::ArrayRCP<GlobalOrdinal>::iterator> colGIndsPtrs_;  // allocated only if indicesAreGlobal()
    Teuchos::ArrayRCP<typename Teuchos::ArrayRCP<LocalOrdinal >::iterator> colLIndsPtrs_;  // allocated only if indicesAreLocal()

    // these are RCPs because they are optional
    // importer is needed if DomainMap is not sameas ColumnMap
    Teuchos::RCP<Import<LocalOrdinal,GlobalOrdinal> > importer_;
    // exporter is needed if RowMap is not sameas RangeMap
    Teuchos::RCP<Export<LocalOrdinal,GlobalOrdinal> > exporter_;

    bool fillComplete_;
    bool storageOptimized_;
    bool staticProfile_;
    bool hasColMap_;
    bool lowerTriangular_;
    bool upperTriangular_;
    bool indicesAreLocal_;
    bool indicesAreGlobal_;
    bool indicesAreSorted_;
    bool noRedundancies_;
    bool indicesAreAllocated_;
    bool haveGlobalConstants_;
  };

  template <class LocalOrdinal, class GlobalOrdinal> 
  CrsGraphData<LocalOrdinal,GlobalOrdinal>::CrsGraphData(const Map<LocalOrdinal,GlobalOrdinal> & rowMap)
  : comm_(rowMap.getComm())
  , rowMap_(rowMap)
  , colMap_(rowMap)     // 
  , rangeMap_(rowMap)   // these three must be set to something; we'll set them appropriately later
  , domainMap_(rowMap)  //
  , numGlobalEntries_(GOT::invalid())
  , numGlobalDiags_(GOT::zero())
  , globalMaxNumRowEntries_(GOT::zero())
  , numLocalEntries_(0)
  , numLocalDiags_(0)
  , localMaxNumEntries_(0)
  , colGInds_()
  , colLInds_()
  , totalNumAllocated_(0)
  , contigColGInds_(Teuchos::null)
  , contigColLInds_(Teuchos::null)
  , colGIndsPtrs_()
  , colLIndsPtrs_()
  , importer_(Teuchos::null)
  , exporter_(Teuchos::null)
  , fillComplete_(false)
  , storageOptimized_(false)
  , staticProfile_(false)
  , hasColMap_(false)
  , lowerTriangular_(true)
  , upperTriangular_(true)
  , indicesAreLocal_(false)
  , indicesAreGlobal_(false)
  , indicesAreSorted_(false)
  , noRedundancies_(false)
  , indicesAreAllocated_(false)
  , haveGlobalConstants_(false)
  {
    if (rowMap.getNumMyEntries() > 0) {
      rowNumToAlloc_ = Teuchos::arcp<LocalOrdinal>(rowMap.getNumMyEntries());
      rowNumEntries_ = Teuchos::arcp<LocalOrdinal>(rowMap.getNumMyEntries());
      std::fill(rowNumToAlloc_.begin(),rowNumToAlloc_.end(),0);
      std::fill(rowNumEntries_.begin(),rowNumEntries_.end(),0);
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal> 
  CrsGraphData<LocalOrdinal,GlobalOrdinal>::CrsGraphData(const Map<LocalOrdinal,GlobalOrdinal> & rowMap, const Map<LocalOrdinal,GlobalOrdinal> & colMap)
  : comm_(rowMap.getComm())
  , rowMap_(rowMap)
  , colMap_(colMap)     
  , rangeMap_(rowMap)   // these two must be set to something; we'll set them appropriately later
  , domainMap_(rowMap)  //
  , numGlobalEntries_(GOT::invalid())
  , numGlobalDiags_(GOT::zero())
  , globalMaxNumRowEntries_(GOT::zero())
  , numLocalEntries_(0)
  , numLocalDiags_(0)
  , localMaxNumEntries_(0)
  , colGInds_()
  , colLInds_()
  , totalNumAllocated_(0)
  , contigColGInds_(Teuchos::null)
  , contigColLInds_(Teuchos::null)
  , colGIndsPtrs_()
  , colLIndsPtrs_()
  , importer_(Teuchos::null)
  , exporter_(Teuchos::null)
  , fillComplete_(false)
  , storageOptimized_(false)
  , staticProfile_(false)
  , hasColMap_(true)
  , lowerTriangular_(true)
  , upperTriangular_(true)
  , indicesAreLocal_(false)
  , indicesAreGlobal_(false)
  , indicesAreSorted_(false)
  , noRedundancies_(false)
  , indicesAreAllocated_(false)
  , haveGlobalConstants_(false)
  {
    if (rowMap.getNumMyEntries() > 0) {
      rowNumToAlloc_ = Teuchos::arcp<LocalOrdinal>(rowMap.getNumMyEntries());
      rowNumEntries_ = Teuchos::arcp<LocalOrdinal>(rowMap.getNumMyEntries());
      std::fill(rowNumToAlloc_.begin(),rowNumToAlloc_.end(),0);
      std::fill(rowNumEntries_.begin(),rowNumEntries_.end(),0);
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal> 
  CrsGraphData<LocalOrdinal,GlobalOrdinal>::~CrsGraphData() 
  {}

} // namespace Tpetra

#endif /* TPETRA_CRSGRAPHDATA_HPP */

