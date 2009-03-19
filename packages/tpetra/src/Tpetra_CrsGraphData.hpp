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
    // structures used before optimizeStorage()
    Teuchos::Array<Teuchos::ArrayRCP<GlobalOrdinal> > colGInds_;                      // empty after fillComplete()
    Teuchos::Array<Teuchos::ArrayRCP<LocalOrdinal> >  colLInds_;                      // empty before fillComplete(), after optimizeStorage()
    Teuchos::Array<LocalOrdinal> rowNumEntries_, rowNumAlloc_;                        // empty after OptimizeStorage()
    // structure used after optimizeStorage()
    Teuchos::ArrayRCP<LocalOrdinal> contigColInds_;                                   // allocated during optimizeStorage()
    /* colIndsPtrs_[j] is always the begin() iterator from an ArrayView of the source or contigColInds_
       of the appropriate length. in a debug build, it is an ArrayRCP, which does bounds checking. in an optimized
       build, it is a C pointer. colIndsPtrs_ is allocated to numLocalRows()+1; the span of the jth row begins with
       colIndsPtrs_[j] and ends before colIndsPtrs_[j+1] */
    Teuchos::Array<typename Teuchos::ArrayRCP<LocalOrdinal>::iterator> colIndsPtrs_;
    // these are RCPs because they are optional
    // importer is needed if DomainMap is not sameas ColumnMap
    Teuchos::RCP<Import<LocalOrdinal,GlobalOrdinal> > importer_;
    // exporter is needed if RowMap is not sameas RangeMap
    Teuchos::RCP<Export<LocalOrdinal,GlobalOrdinal> > exporter_;

    bool fillComplete_;
    bool hasColMap_;
    bool lowerTriangular_;
    bool upperTriangular_;
    bool noDiagonal_;
    bool indicesAreLocal_;
    bool indicesAreGlobal_;
    bool indicesAreSorted_;
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
  , colGInds_(0)
  , colLInds_(0)
  , rowNumEntries_(rowMap.getNumMyEntries(),0)
  , rowNumAlloc_(rowMap.getNumMyEntries(),0)
  , contigColInds_(Teuchos::null)
  , colIndsPtrs_(0)
  , importer_(Teuchos::null)
  , exporter_(Teuchos::null)
  , fillComplete_(false)
  , hasColMap_(false)
  , lowerTriangular_(true)
  , upperTriangular_(true)
  , noDiagonal_(true)
  , indicesAreLocal_(false)
  , indicesAreGlobal_(false)
  , indicesAreSorted_(false)
  , indicesAreAllocated_(false)
  , haveGlobalConstants_(false)
  {}

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
  , colGInds_(0)
  , colLInds_(0)
  , rowNumEntries_(rowMap.getNumMyEntries(),0)
  , rowNumAlloc_(rowMap.getNumMyEntries(),0)
  , contigColInds_(Teuchos::null)
  , colIndsPtrs_(0)
  , importer_(Teuchos::null)
  , exporter_(Teuchos::null)
  , fillComplete_(false)
  , hasColMap_(true)
  , lowerTriangular_(true)
  , upperTriangular_(true)
  , noDiagonal_(true)
  , indicesAreLocal_(false)
  , indicesAreGlobal_(false)
  , indicesAreSorted_(false)
  , indicesAreAllocated_(false)
  , haveGlobalConstants_(false)
  {}

  template <class LocalOrdinal, class GlobalOrdinal> 
  CrsGraphData<LocalOrdinal,GlobalOrdinal>::~CrsGraphData() 
  {}

} // namespace Tpetra

#endif /* TPETRA_CRSGRAPHDATA_HPP */

