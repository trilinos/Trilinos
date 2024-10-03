// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Panzer_SetupPartitionedWorksetUtilities.hpp>
#include "Panzer_STK_WorksetFactory.hpp"

#include "Panzer_LocalMeshInfo.hpp"

#include "Panzer_WorksetFactoryBase.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_LocalMeshUtilities.hpp"



namespace panzer_stk {

/** Set mesh
  */
void WorksetFactory::setMesh(const Teuchos::RCP<const panzer_stk::STK_Interface> & mesh)
{
   mesh_ = mesh;
}

Teuchos::RCP<std::map<unsigned,panzer::Workset> > WorksetFactory::
getSideWorksets(const panzer::WorksetDescriptor & desc,
                const panzer::WorksetNeeds & needs) const
{
  TEUCHOS_ASSERT(desc.useSideset());

  return panzer_stk::buildBCWorksets(*mesh_,needs,desc.getElementBlock(0),desc.getSideset());
}

Teuchos::RCP<std::map<unsigned,panzer::Workset> > WorksetFactory::
getSideWorksets(const panzer::WorksetDescriptor & desc,
                const panzer::WorksetNeeds & needs_a,
                const panzer::WorksetNeeds & needs_b) const
{
  // ensure that this is a interface descriptor
  TEUCHOS_ASSERT(desc.connectsElementBlocks());
  TEUCHOS_ASSERT(desc.getSideset(0)==desc.getSideset(1));
  return panzer_stk::buildBCWorksets(*mesh_, needs_a, desc.getElementBlock(0),
                                             needs_b, desc.getElementBlock(1),
                                             desc.getSideset(0));
}

Teuchos::RCP<std::vector<panzer::Workset> > WorksetFactory::
getWorksets(const panzer::WorksetDescriptor & worksetDesc,
            const panzer::WorksetNeeds & needs) const
{
  PANZER_FUNC_TIME_MONITOR("panzer_stk::WorksetFactory::getWorksets");
  if(worksetDesc.requiresPartitioning()){

    // Generate the local mesh info if it doesn't already exist
    if(mesh_info_ == Teuchos::null){
      TEUCHOS_ASSERT(mesh_ != Teuchos::null);
      mesh_info_ = panzer_stk::generateLocalMeshInfo(*mesh_);
    }

    auto worksets = panzer::buildPartitionedWorksets(*mesh_info_, worksetDesc, this->getOrientationsInterface());

    // Fill in whatever is in the needs object
    // FIXME: This will just get optimized out... Adding volatile to the calls makes the worksets pretty ugly
    for(auto & workset : *worksets){

      // Initialize IntegrationValues from integration descriptors
      for(const auto & id : needs.getIntegrators())
        workset.getIntegrationValues(id);

      // Initialize PointValues from point descriptors
      for(const auto & pd : needs.getPoints())
        workset.getPointValues(pd);

      // Initialize BasisValues
      for(const auto & bd : needs.getBases()){

        // Initialize BasisValues from integrators
        for(const auto & id : needs.getIntegrators())
          workset.getBasisValues(bd,id);

        // Initialize BasisValues from points
        for(const auto & pd : needs.getPoints())
          workset.getBasisValues(bd,pd);
      }
    }

    return worksets;

  } else if(!worksetDesc.useSideset()) {
    // The non-partitioned case always creates worksets with just the
    // owned elements.  CLASSIC_MODE gets the workset size directly
    // from needs that is populated externally. As we transition away
    // from classic mode, we need to create a copy of needs and
    // override the workset size with values from WorksetDescription.
    if (worksetDesc.getWorksetSize() == panzer::WorksetSizeType::CLASSIC_MODE)
      return panzer_stk::buildWorksets(*mesh_,worksetDesc.getElementBlock(), needs);
    else {
      int worksetSize = worksetDesc.getWorksetSize();
      if (worksetSize == panzer::WorksetSizeType::ALL_ELEMENTS) {
        std::vector<stk::mesh::Entity> elements;
        mesh_->getMyElements(worksetDesc.getElementBlock(),elements);
        worksetSize = elements.size();
      }
      panzer::WorksetNeeds tmpNeeds(needs);
      tmpNeeds.cellData = panzer::CellData(worksetSize,needs.cellData.getCellTopology());
      return panzer_stk::buildWorksets(*mesh_,worksetDesc.getElementBlock(), tmpNeeds);
    }
  }
  else if(worksetDesc.useSideset() && worksetDesc.sideAssembly()) {
    // uses cascade by default, each subcell has its own workset
    return panzer_stk::buildWorksets(*mesh_,needs,worksetDesc.getSideset(),worksetDesc.getElementBlock(),true);
  }
  else {
    TEUCHOS_ASSERT(false);
  }
}

}
