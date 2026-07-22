// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __Panzer_GatherTangent_BlockedEpetra_impl_hpp__
#define   __Panzer_GatherTangent_BlockedEpetra_impl_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Epetra
#include "Epetra_Map.h"

// Panzer
#include "Panzer_BlockedVector_ReadOnly_GlobalEvaluationData.hpp"
#include "Panzer_GlobalEvaluationData.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_GlobalIndexer_Utilities.hpp"

// Phalanx
#include "Phalanx_DataLayout.hpp"

// Teuchos
#include "Teuchos_Assert.hpp"
#include "Teuchos_FancyOStream.hpp"

// Thyra
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_SpmdVectorBase.hpp"

///////////////////////////////////////////////////////////////////////////////
//
//  Initializing Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename TRAITS, typename LO, typename GO>
panzer::GatherTangent_BlockedEpetra<EvalT, TRAITS, LO, GO>::
GatherTangent_BlockedEpetra(
  const std::vector<Teuchos::RCP<const GlobalIndexer>>&
    indexers,
  const Teuchos::ParameterList& p)
  :
  indexers_(indexers),
  useTimeDerivativeSolutionVector_(false),
  globalDataKey_("Tangent Gather Container")
{
  using panzer::PureBasis;
  using PHX::MDField;
  using PHX::print;
  using std::size_t;
  using std::string;
  using std::vector;
  using Teuchos::RCP;

  // Get the necessary information from the ParameterList.
  const vector<string>& names = *(p.get<RCP<vector<string>>>("DOF Names"));
  indexerNames_ = p.get<RCP<vector<string>>>("Indexer Names");
  RCP<PureBasis> basis = p.get<RCP<PureBasis>>("Basis");
  if (p.isType<bool>("Use Time Derivative Solution Vector"))
    useTimeDerivativeSolutionVector_ =
      p.get<bool>("Use Time Derivative Solution Vector");
  if (p.isType<string>("Global Data Key"))
    globalDataKey_ = p.get<string>("Global Data Key");

  // Allocate the fields.
  int numFields(names.size());
  gatherFields_.resize(numFields);
  for (int fd(0); fd < numFields; ++fd)
  {
    gatherFields_[fd] =
      MDField<ScalarT, Cell, NODE>(names[fd], basis->functional);
    this->addEvaluatedField(gatherFields_[fd]);
  } // end loop over names

  // Figure out what the first active name is.
  string firstName("<none>");
  if (numFields > 0)
    firstName = names[0];
  string n("GatherTangent (Blocked Epetra):  " + firstName + " (" +
    print<EvalT>() + ")");
  this->setName(n);
} // end of Initializing Constructor

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename TRAITS, typename LO, typename GO>
void
panzer::GatherTangent_BlockedEpetra<EvalT, TRAITS, LO, GO>::
postRegistrationSetup(
  typename TRAITS::SetupData /* d  */,
  PHX::FieldManager<TRAITS>& /* fm */)
{
  using std::size_t;
  using std::string;
  using Teuchos::null;
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_->size());
  int numFields(gatherFields_.size());
  indexerIds_.resize(numFields);
  subFieldIds_.resize(numFields);
  for (int fd(0); fd < numFields; ++fd)
  {
    // Get the field ID from the DOF manager.
    const string& fieldName((*indexerNames_)[fd]);
    indexerIds_[fd]  = getFieldBlock(fieldName, indexers_);
    subFieldIds_[fd] = indexers_[indexerIds_[fd]]->getFieldNum(fieldName);
    TEUCHOS_ASSERT(indexerIds_[fd] >= 0);
  } // end loop over gatherFields_
  indexerNames_ = null;
} // end of postRegistrationSetup()

///////////////////////////////////////////////////////////////////////////////
//
//  preEvaluate()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename TRAITS, typename LO, typename GO>
void
panzer::GatherTangent_BlockedEpetra<EvalT, TRAITS, LO, GO>::
preEvaluate(
  typename TRAITS::PreEvalData d)
{
  using std::logic_error;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::typeName;
  using Thyra::ProductVectorBase;
  using BVROGED = BlockedVector_ReadOnly_GlobalEvaluationData;
  using GED     = GlobalEvaluationData;
  if (d.gedc->containsDataObject(globalDataKey_))
  {
    RCP<GED> ged = d.gedc->getDataObject(globalDataKey_);
    xBvRoGed_    = rcp_dynamic_cast<BVROGED>(ged, true);
  } // end if (d.gedc.containsDataObject(globalDataKey_))
} // end of preEvaluate()

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename TRAITS, typename LO, typename GO>
void
panzer::GatherTangent_BlockedEpetra<EvalT, TRAITS, LO, GO>::
evaluateFields(
  typename TRAITS::EvalData workset)
{
  using PHX::MDField;
  using std::size_t;
  using std::string;
  using std::vector;
  using Teuchos::ArrayRCP;
  using Teuchos::ptrFromRef;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::VectorBase;
  using Thyra::SpmdVectorBase;

  // If no global evaluation data container was set, then this evaluator
  // becomes a no-op.
  if (xBvRoGed_.is_null())
    return;

  // For convenience, pull out some objects from the workset.
  string blockId(this->wda(workset).block_id);
  const vector<size_t>& localCellIds = this->wda(workset).cell_local_ids;
  int numFields(gatherFields_.size()), numCells(localCellIds.size());

  // Loop over the fields to be gathered.
  for (int fieldIndex(0); fieldIndex < numFields; ++fieldIndex)
  {
    MDField<ScalarT, Cell, NODE>& field = gatherFields_[fieldIndex];
    int indexerId(indexerIds_[fieldIndex]),
      subFieldNum(subFieldIds_[fieldIndex]);

    // Grab the local data for inputing.
    auto xEvRoGed = xBvRoGed_->getGEDBlock(indexerId);
    auto subRowIndexer = indexers_[indexerId];
    const vector<int>& elmtOffset =
      subRowIndexer->getGIDFieldOffsets(blockId, subFieldNum);
    int numBases(elmtOffset.size());

    // Gather operation for each cell in the workset.
    for (int cell(0); cell < numCells; ++cell)
    {
      LO cellLocalId = localCellIds[cell];
      auto LIDs = subRowIndexer->getElementLIDs(cellLocalId);

      // Loop over the basis functions and fill the fields.
      for (int basis(0); basis < numBases; ++basis)
      {
        int offset(elmtOffset[basis]), lid(LIDs[offset]);
        field(cell, basis) = (*xEvRoGed)[lid];
      } // end loop over the basis functions
    } // end loop over localCellIds
  } // end loop over the fields to be gathered
} // end of evaluateFields()

#endif // __Panzer_GatherTangent_BlockedEpetra_impl_hpp__
