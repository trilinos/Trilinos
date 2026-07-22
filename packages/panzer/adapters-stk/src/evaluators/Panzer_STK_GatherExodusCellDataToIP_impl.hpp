// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_STK_GATHER_EXODUS_CELL_DATA_TO_IP_IMPL_HPP
#define PANZER_STK_GATHER_EXODUS_CELL_DATA_TO_IP_IMPL_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Teuchos_FancyOStream.hpp"
#include <sstream>

// **********************************************************************
// Specialization: Residual
// **********************************************************************

template<typename EvalT, typename Traits>
panzer_stk::GatherExodusCellDataToIP<EvalT, Traits>::
GatherExodusCellDataToIP(const Teuchos::RCP<const STK_Interface> & mesh,
                         const std::vector<std::string>& fieldNames,
                         const std::vector<std::string>& exodusNames,
                         const Teuchos::RCP<panzer::IntegrationRule>& integrationRule)
  : mesh_(mesh),
    exodusNames_(exodusNames)
{
  using panzer::Cell;
  using panzer::IP;

  TEUCHOS_ASSERT(fieldNames.size() == exodusNames.size());
  TEUCHOS_ASSERT(nonnull(mesh));
  TEUCHOS_ASSERT(nonnull(integrationRule));
  
  gatherFields_.resize(fieldNames.size());
  stkFields_.resize(fieldNames.size());
  for (std::size_t fd = 0; fd < fieldNames.size(); ++fd) {
    gatherFields_[fd] = PHX::MDField<ScalarT,Cell,IP>(fieldNames[fd],integrationRule->dl_scalar);
    this->addEvaluatedField(gatherFields_[fd]);
  }

  std::ostringstream os;
  for (size_t i=0; i < fieldNames.size(); ++i) {
    os << fieldNames[i];
    if (i < fieldNames.size()-1)
      os << ",";
  }

  this->setName("GatherExodusCellDataToIP: "+os.str());
}

// **********************************************************************
template<typename EvalT, typename Traits>
void panzer_stk::GatherExodusCellDataToIP<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData /* d */,
		      PHX::FieldManager<Traits>& /* fm */)
{
  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    std::string fieldName = gatherFields_[fd].fieldTag().name();

    stkFields_[fd] = mesh_->getMetaData()->template get_field<double>(stk::topology::ELEMENT_RANK, exodusNames_[fd]);

    if(stkFields_[fd]==0) {
      std::stringstream ss;
      mesh_->printMetaData(ss);
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                                 "panzer_stk::GatherExodusCellDataToIP: STK field " << "\"" << fieldName << "\" "
                                 "not found.\n STK meta data follows: \n\n" << ss.str());
    }
  }
}

// **********************************************************************
template<typename EvalT, typename Traits>
void panzer_stk::GatherExodusCellDataToIP<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  const std::vector<stk::mesh::Entity>& localElementEntities = *mesh_->getElementsOrderedByLID();
  const std::vector<std::size_t>& localWorksetCellIds = this->wda(workset).cell_local_ids;

  auto host_mirror = Kokkos::create_mirror_view(Kokkos::WithoutInitializing,gatherFields_[0].get_static_view());

  for (std::size_t fieldIndex=0; fieldIndex < gatherFields_.size(); ++fieldIndex) {
    VariableField* field = stkFields_[fieldIndex];
    const std::size_t numQuadPoints = gatherFields_[fieldIndex].extent(1);

    for(std::size_t worksetCellIndex=0;worksetCellIndex<localWorksetCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localWorksetCellIds[worksetCellIndex];
      for(std::size_t qp=0; qp < numQuadPoints; ++qp) {
        double value = *stk::mesh::field_data(*field, localElementEntities[cellLocalId]);
        host_mirror(worksetCellIndex,qp) = value;
      }
    }

    Kokkos::deep_copy(gatherFields_[fieldIndex].get_static_view(),host_mirror);
  }
}
#endif
