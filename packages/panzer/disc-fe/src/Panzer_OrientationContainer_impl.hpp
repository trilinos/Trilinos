// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_OrientationContainer_impl_hpp__
#define __Panzer_OrientationContainer_impl_hpp__

#include "PanzerDiscFE_config.hpp"

#include "Panzer_GlobalIndexer.hpp"

namespace panzer {

template <typename Scalar,typename Array,typename LocalOrdinal,typename GlobalOrdinal>
OrientationContainer<Scalar,Array,LocalOrdinal,GlobalOrdinal>::
OrientationContainer(const Teuchos::RCP<const GlobalIndexer> & globalIndexer,
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
buildOrientationContainer(const Teuchos::RCP<const GlobalIndexer> & globalIndexer,
                          const std::string & fieldName)
{
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::RCP;
  using Teuchos::rcp;
  
  // int, int 
  {
    typedef int LO;
    typedef int GO;

    RCP<const GlobalIndexer> ugi 
        = rcp_dynamic_cast<const GlobalIndexer>(globalIndexer);
    if(ugi!=Teuchos::null)
      return rcp(new OrientationContainer<Scalar,Array,LO,GO>(ugi,fieldName));
  }

  // int, panzer::GlobalOrdinal
  {
    typedef int LO;
    typedef panzer::GlobalOrdinal GO;

    RCP<const GlobalIndexer> ugi 
        = rcp_dynamic_cast<const GlobalIndexer>(globalIndexer);
    if(ugi!=Teuchos::null)
      return rcp(new OrientationContainer<Scalar,Array,LO,GO>(ugi,fieldName));
  }

  // int, pair<int,int>
  {
    typedef int LO;
    typedef std::pair<int,int> GO;

    RCP<const GlobalIndexer> ugi 
        = rcp_dynamic_cast<const GlobalIndexer>(globalIndexer);
    if(ugi!=Teuchos::null)
      return rcp(new OrientationContainer<Scalar,Array,LO,GO>(ugi,fieldName));
  }

  // int, pair<int,panzer::GlobalOrdinal>
  {
    typedef int LO;
    typedef std::pair<int,panzer::GlobalOrdinal> GO;

    RCP<const GlobalIndexer> ugi 
        = rcp_dynamic_cast<const GlobalIndexer>(globalIndexer);
    if(ugi!=Teuchos::null)
      return rcp(new OrientationContainer<Scalar,Array,LO,GO>(ugi,fieldName));
  }

  TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                             "panzer::buildOrientationContainer: Could not cast GlobalIndexer");
}

} // end namespace panzer

#endif
