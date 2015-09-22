// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Panzer_OrientationContainer_impl_hpp__
#define __Panzer_OrientationContainer_impl_hpp__

#include "Panzer_ConfigDefs.hpp"

#include "Panzer_UniqueGlobalIndexer.hpp"

namespace panzer {

template <typename Scalar,typename Array,typename LocalOrdinal,typename GlobalOrdinal>
OrientationContainer<Scalar,Array,LocalOrdinal,GlobalOrdinal>::
OrientationContainer(const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinal,GlobalOrdinal> > & globalIndexer,
                     const std::string & fieldName)
  : globalIndexer_(globalIndexer)
  , fieldName_(fieldName)
{
}

template <typename Scalar,typename Array,typename LocalOrdinal,typename GlobalOrdinal>
void
OrientationContainer<Scalar,Array,LocalOrdinal,GlobalOrdinal>::
getOrientations(const std::string & blockId,
                const std::vector<std::size_t> & cell_local_ids,
                Array & orientationsArray) const
{
  int fieldNum = globalIndexer_->getFieldNum(fieldName_);
  const std::vector<int> & elmtOffset = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);
 
  // gather operation for each cell in workset
  for(std::size_t cellIndex=0;cellIndex<cell_local_ids.size();++cellIndex) {
    std::vector<double> orientation;
    std::size_t cellLocalId = cell_local_ids[cellIndex];
 
    globalIndexer_->getElementOrientation(cellLocalId,orientation); 

    // loop over basis functions and fill the fields
    for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
      int offset = elmtOffset[basis];
      orientationsArray(cellIndex,basis) = orientation[offset];
    }
  }
}

template <typename Scalar,typename Array>
Teuchos::RCP<const OrientationContainerBase<Scalar,Array> > 
buildOrientationContainer(const Teuchos::RCP<const UniqueGlobalIndexerBase> & globalIndexer,
                          const std::string & fieldName)
{
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::RCP;
  using Teuchos::rcp;
  
  // int, int 
  {
    typedef int LO;
    typedef int GO;

    RCP<const UniqueGlobalIndexer<LO,GO> > ugi 
        = rcp_dynamic_cast<const UniqueGlobalIndexer<LO,GO> >(globalIndexer);
    if(ugi!=Teuchos::null)
      return rcp(new OrientationContainer<Scalar,Array,LO,GO>(ugi,fieldName));
  }

  // int, Ordinal64
  {
    typedef int LO;
    typedef Ordinal64 GO;

    RCP<const UniqueGlobalIndexer<LO,GO> > ugi 
        = rcp_dynamic_cast<const UniqueGlobalIndexer<LO,GO> >(globalIndexer);
    if(ugi!=Teuchos::null)
      return rcp(new OrientationContainer<Scalar,Array,LO,GO>(ugi,fieldName));
  }

  // int, pair<int,int>
  {
    typedef int LO;
    typedef std::pair<int,int> GO;

    RCP<const UniqueGlobalIndexer<LO,GO> > ugi 
        = rcp_dynamic_cast<const UniqueGlobalIndexer<LO,GO> >(globalIndexer);
    if(ugi!=Teuchos::null)
      return rcp(new OrientationContainer<Scalar,Array,LO,GO>(ugi,fieldName));
  }

  // int, pair<int,Ordinal64>
  {
    typedef int LO;
    typedef std::pair<int,int> GO;

    RCP<const UniqueGlobalIndexer<LO,GO> > ugi 
        = rcp_dynamic_cast<const UniqueGlobalIndexer<LO,GO> >(globalIndexer);
    if(ugi!=Teuchos::null)
      return rcp(new OrientationContainer<Scalar,Array,LO,GO>(ugi,fieldName));
  }

  TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                             "panzer::buildOrientationContainer: Could not cast UniqueGlobalIndexerBase");
}

} // end namespace panzer

#endif
