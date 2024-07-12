// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __Panzer_GatherTangent_Epetra_impl_hpp__
#define   __Panzer_GatherTangent_Epetra_impl_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Epetra
#include "Epetra_Vector.h"

// Panzer
#include "Panzer_EpetraVector_ReadOnly_GlobalEvaluationData.hpp"
#include "Panzer_GlobalEvaluationData.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"

// Teuchos
#include "Teuchos_Assert.hpp"

// Thyra
#include "Thyra_SpmdVectorBase.hpp"

///////////////////////////////////////////////////////////////////////////////
//
//  Initializing Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename TRAITS, typename LO, typename GO>
panzer::GatherTangent_Epetra<EvalT, TRAITS, LO, GO>::
GatherTangent_Epetra(
  const Teuchos::RCP<const panzer::GlobalIndexer>& indexer,
  const Teuchos::ParameterList& p)
  :
  globalIndexer_(indexer),
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
  RCP<const PureBasis> basis;
  if (p.isType<RCP<PureBasis>>("Basis"))
    basis = p.get<RCP<PureBasis>>("Basis");
  else // if (not p.isType<RCP<PureBasis>>("Basis"))
    basis = p.get<RCP<const PureBasis>>("Basis");
  if (p.isType<bool>("Use Time Derivative Solution Vector"))
    useTimeDerivativeSolutionVector_ =
      p.get<bool>("Use Time Derivative Solution Vector");
  if (p.isType<string>("Global Data Key"))
    globalDataKey_ = p.get<string>("Global Data Key");

  // Allocate fields.
  int numFields(names.size());
  gatherFields_.resize(numFields);
  for (int fd(0); fd < numFields; ++fd)
  {
    gatherFields_[fd] =
      MDField<ScalarT, Cell, NODE>(names[fd], basis->functional);
    this->addEvaluatedField(gatherFields_[fd]);

    // This fixes the case of dxdpEvRoGed_ being null and no
    // operations performed during evalaute. Keeps the field with
    // initial zero state.
    this->addUnsharedField(gatherFields_[fd].fieldTag().clone());
  } // end loop over names

  // Figure out what the first active name is.
  string firstName("<none>");
  if (numFields > 0)
    firstName = names[0];
  string n("GatherTangent (Epetra):  " + firstName + " (" +
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
panzer::GatherTangent_Epetra<EvalT, TRAITS, LO, GO>::
postRegistrationSetup(
  typename TRAITS::SetupData /* d  */,
  PHX::FieldManager<TRAITS>& /* fm */)
{
  using std::logic_error;
  using std::size_t;
  using std::string;
  using Teuchos::null;
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_->size());
  int numFields(gatherFields_.size());
  fieldIds_.resize(numFields);
  for (int fd(0); fd < numFields; ++fd)
  {
    // Get the field ID from the DOF manager.
    const string& fieldName((*indexerNames_)[fd]);
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    // This is the error return code; raise the alarm.
    TEUCHOS_TEST_FOR_EXCEPTION(fieldIds_[fd] == -1, logic_error,
      "GatherTangent_Epetra<Residual>: Could not find field \"" + fieldName +
      "\" in the global indexer. ");
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
panzer::GatherTangent_Epetra<EvalT, TRAITS, LO, GO>::
preEvaluate(
  typename TRAITS::PreEvalData d)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using EVROGED = EpetraVector_ReadOnly_GlobalEvaluationData;
  using GED     = GlobalEvaluationData;
  if (d.gedc->containsDataObject(globalDataKey_))
  {
    RCP<GED> ged = d.gedc->getDataObject(globalDataKey_);
    dxdpEvRoGed_ = rcp_dynamic_cast<EVROGED>(ged, true);
  }
} // end of preEvaluate()

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename TRAITS, typename LO, typename GO>
void
panzer::GatherTangent_Epetra<EvalT, TRAITS, LO, GO>::
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
  using Thyra::SpmdVectorBase;

  // If no global evaluation data container was set, then this evaluator
  // becomes a no-op.
  if (dxdpEvRoGed_.is_null())
    return;

  // For convenience, pull out some objects from the workset.
  string blockId(this->wda(workset).block_id);
  const vector<size_t>& localCellIds = this->wda(workset).cell_local_ids;
  int numCells(localCellIds.size()), numFields(gatherFields_.size());

  // NOTE:  A reordering of these loops will likely improve performance.  The
  //        "getGIDFieldOffsets may be expensive.  However the "getElementGIDs"
  //        can be cheaper.  However the lookup for LIDs may be more expensive!

  // Gather operation for each cell in the workset.
  auto LIDs = globalIndexer_->getLIDs();
  auto LIDs_h = Kokkos::create_mirror_view(LIDs);
  Kokkos::deep_copy(LIDs_h, LIDs);
  // Loop over the fields to be gathered.
  for (int fieldIndex(0); fieldIndex < numFields; ++fieldIndex)
  {
    MDField<ScalarT, Cell, NODE>& field = gatherFields_[fieldIndex];
    auto field_h = Kokkos::create_mirror_view(field.get_static_view());
    for (int cell(0); cell < numCells; ++cell)
    {
      size_t cellLocalId(localCellIds[cell]);
      int fieldNum(fieldIds_[fieldIndex]);
      const vector<int>& elmtOffset =
        globalIndexer_->getGIDFieldOffsets(blockId, fieldNum);
      int numBases(elmtOffset.size());

      // Loop over the basis functions and fill the fields.
      for (int basis(0); basis < numBases; ++basis)
      {
        int offset(elmtOffset[basis]), lid(LIDs_h(cellLocalId, offset));
        field_h(cell, basis) = (*dxdpEvRoGed_)[lid];
      } // end loop over the basis functions
    } // end loop over the cells in the workset
    Kokkos::deep_copy(field.get_static_view(), field_h);
  } // end loop over the fields to be gathered
} // end of evaluateFields()

#endif // __Panzer_GatherTangent_Epetra_impl_hpp__
