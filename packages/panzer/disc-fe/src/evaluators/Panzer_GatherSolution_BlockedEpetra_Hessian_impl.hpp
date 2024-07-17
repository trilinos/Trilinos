// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __Panzer_GatherSolution_BlockedEpetra_Hessian_impl_hpp__
#define   __Panzer_GatherSolution_BlockedEpetra_Hessian_impl_hpp__

// Only do this if required by the user.
#ifdef    Panzer_BUILD_HESSIAN_SUPPORT

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Panzer
#include "Panzer_BlockedEpetraLinearObjContainer.hpp"
#include "Panzer_BlockedVector_ReadOnly_GlobalEvaluationData.hpp"
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_ParameterList_GlobalEvaluationData.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_GlobalIndexer.hpp"

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
//  Initializing Constructor (Hessian Specialization)
//
///////////////////////////////////////////////////////////////////////////////
template<typename TRAITS, typename LO, typename GO>
panzer::GatherSolution_BlockedEpetra<panzer::Traits::Hessian, TRAITS, LO, GO>::
GatherSolution_BlockedEpetra(
  const std::vector<Teuchos::RCP<const GlobalIndexer<LO, int>>>&
    indexers,
  const Teuchos::ParameterList& p)
  :
  indexers_(indexers)
{
  using panzer::PureBasis;
  using PHX::MDField;
  using PHX::print;
  using std::size_t;
  using std::string;
  using std::vector;
  using Teuchos::RCP;
  GatherSolution_Input input;
  input.setParameterList(p);
  const vector<string>& names = input.getDofNames();
  RCP<const PureBasis>  basis = input.getBasis();
  indexerNames_                    = input.getIndexerNames();
  useTimeDerivativeSolutionVector_ = input.useTimeDerivativeSolutionVector();
  globalDataKey_                   = input.getGlobalDataKey();
  gatherSeedIndex_                 = input.getGatherSeedIndex();
  sensitivitiesName_               = input.getSensitivitiesName();
  firstSensitivitiesAvailable_     = input.firstSensitivitiesAvailable();
  secondSensitivitiesAvailable_    = input.secondSensitivitiesAvailable();
  sensitivities2ndPrefix_          = input.getSecondSensitivityDataKeyPrefix();

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
  string firstName("<none>"), n("GatherSolution (BlockedEpetra");
  if (numFields > 0)
    firstName = names[0];
  if (not firstSensitivitiesAvailable_)
    n += ", No First Sensitivities";
  n += "):  " + firstName + " (" + print<EvalT>() + ")";
  this->setName(n);
} // end of Initializing Constructor (Hessian Specialization)

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup() (Hessian Specialization)
//
///////////////////////////////////////////////////////////////////////////////
template<typename TRAITS, typename LO, typename GO>
void
panzer::GatherSolution_BlockedEpetra<panzer::Traits::Hessian, TRAITS, LO, GO>::
postRegistrationSetup(
  typename TRAITS::SetupData /* d  */,
  PHX::FieldManager<TRAITS>& /* fm */)
{
  using std::size_t;
  using std::string;
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_.size());
  int numFields(gatherFields_.size());
  indexerIds_.resize(numFields);
  subFieldIds_.resize(numFields);
  for (int fd(0); fd < numFields; ++fd)
  {
    // Get the field ID from the DOF manager.
    const string& fieldName(indexerNames_[fd]);
    indexerIds_[fd]  = getFieldBlock(fieldName, indexers_);
    subFieldIds_[fd] = indexers_[indexerIds_[fd]]->getFieldNum(fieldName);
    TEUCHOS_ASSERT(indexerIds_[fd] >= 0);
  } // end loop over gatherFields_
  indexerNames_.clear();
} // end of postRegistrationSetup() (Hessian Specialization)

///////////////////////////////////////////////////////////////////////////////
//
//  preEvaluate() (Hessian Specialization)
//
///////////////////////////////////////////////////////////////////////////////
template<typename TRAITS, typename LO, typename GO>
void
panzer::GatherSolution_BlockedEpetra<panzer::Traits::Hessian, TRAITS, LO, GO>::
preEvaluate(
  typename TRAITS::PreEvalData d)
{
  using std::endl;
  using std::logic_error;
  using std::string;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::typeName;
  using Thyra::ProductVectorBase;
  using BELOC   = BlockedEpetraLinearObjContainer;
  using BVROGED = BlockedVector_ReadOnly_GlobalEvaluationData;

  // Manage sensitivities.
  firstApplySensitivities_ = false;
  if ((firstSensitivitiesAvailable_                    )  and
      (d.first_sensitivities_name == sensitivitiesName_))
    firstApplySensitivities_ = true;
  secondApplySensitivities_ = false;
  if ((secondSensitivitiesAvailable_                    )  and
      (d.second_sensitivities_name == sensitivitiesName_))
    secondApplySensitivities_ = true;

  // First try the refactored ReadOnly container.
  RCP<GlobalEvaluationData> ged;
  string post(useTimeDerivativeSolutionVector_ ? " - Xdot" : " - X");
  if (d.gedc->containsDataObject(globalDataKey_ + post))
  {
    ged       = d.gedc->getDataObject(globalDataKey_ + post);
    xBvRoGed_ = rcp_dynamic_cast<BVROGED>(ged, true);
  }

  // Otherwise, try the old path.
  else if (d.gedc->containsDataObject(globalDataKey_))
  {
    ged = d.gedc->getDataObject(globalDataKey_);

    // Try to extract the linear object container.
    auto roGed            = rcp_dynamic_cast<const BVROGED>(ged);
    auto blockedContainer = rcp_dynamic_cast<const BELOC>(ged);
    if (not roGed.is_null())
      xBvRoGed_ = rcp_dynamic_cast<BVROGED>(ged, true);
    else if (not blockedContainer.is_null())
    {
      if (useTimeDerivativeSolutionVector_)
        x_ = rcp_dynamic_cast<ProductVectorBase<double>>
          (blockedContainer->get_dxdt());
      else // if (not useTimeDerivativeSolutionVector_)
        x_ = rcp_dynamic_cast<ProductVectorBase<double>>
          (blockedContainer->get_x());
    } // end if roGed or blockedContainer is non-null
  } // end if we're doing things the new or old way

  // Ensure that we actually have something.
  TEUCHOS_ASSERT((not x_.is_null()) or (not xBvRoGed_.is_null()))

  // Don't try to extract dx if it's not required.
  if (not secondApplySensitivities_)
    return;

  // Now parse the second derivative direction.
  if (d.gedc->containsDataObject(sensitivities2ndPrefix_ + globalDataKey_))
  {
    ged = d.gedc->getDataObject(sensitivities2ndPrefix_ + globalDataKey_);
    dxBvRoGed_ = rcp_dynamic_cast<BVROGED>(ged, true);
  } // end if (d.gedc->containsDataObject(...))

  // Ensure that we actually have something.
  TEUCHOS_TEST_FOR_EXCEPTION(dxBvRoGed_.is_null(), logic_error,
    "Cannot find sensitivity vector associated with \"" +
    sensitivities2ndPrefix_ + globalDataKey_ + "\" and \"" + post + "\".");
} // end of preEvaluate() (Hessian Specialization)

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields() (Hessian Specialization)
//
///////////////////////////////////////////////////////////////////////////////
template<typename TRAITS, typename LO, typename GO>
void
panzer::GatherSolution_BlockedEpetra<panzer::Traits::Hessian, TRAITS, LO, GO>::
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
  using Thyra::ProductVectorBase;
  using Thyra::SpmdVectorBase;
  using Thyra::VectorBase;

  // For convenience, pull out some objects from the workset.
  string blockId(this->wda(workset).block_id);
  const vector<size_t>& localCellIds = this->wda(workset).cell_local_ids;
  int numFields(gatherFields_.size()), numCells(localCellIds.size());
  double seedValue(0);
  if (firstApplySensitivities_)
  {
    if ((useTimeDerivativeSolutionVector_) and (gatherSeedIndex_ < 0))
      seedValue = workset.alpha;
    else if (gatherSeedIndex_ < 0)
      seedValue = workset.beta;
    else if (not useTimeDerivativeSolutionVector_)
      seedValue = workset.gather_seeds[gatherSeedIndex_];
    else
      TEUCHOS_ASSERT(false);
  } // end if (firstApplySensitivities_)

  // Turn off sensitivies.  This may be faster if we don't expand the term, but
  // I suspect not, because anywhere it is used the full complement of
  // sensitivies will be needed anyway.
  if (not firstApplySensitivities_)
     seedValue = 0.0;

  vector<int> blockOffsets;
  computeBlockOffsets(blockId, indexers_, blockOffsets);
  if (x_.is_null())
  {
    // Loop over the fields to be gathered.
    for (int fieldInd(0); fieldInd < numFields; ++fieldInd)
    {
      MDField<ScalarT, Cell, NODE>& field = gatherFields_[fieldInd];
      int indexerId(indexerIds_[fieldInd]),
        subFieldNum(subFieldIds_[fieldInd]);

      // Grab the local data for inputing.
      auto xEvRoGed = xBvRoGed_->getGEDBlock(indexerId);
      auto subRowIndexer = indexers_[indexerId];
      const vector<int>& elmtOffset =
        subRowIndexer->getGIDFieldOffsets(blockId, subFieldNum);
      int numBases(elmtOffset.size());

      // Gather operation for each cell in the workset.
      for (int cell(0); cell < numCells; ++cell)
      {
        size_t cellLocalId(localCellIds[cell]);
        auto LIDs = subRowIndexer->getElementLIDs(cellLocalId);

        // Loop over the basis functions and fill the fields.
        for (int basis(0); basis < numBases; ++basis)
        {
          int offset(elmtOffset[basis]), lid(LIDs[offset]);
          field(cell, basis) = (*xEvRoGed)[lid];
        } // end loop over the basis functions
      } // end loop over localCellIds
    } // end loop over the fields to be gathered
  }
  else // if (not x_.is_null())
  {
    // Loop over the fields to be gathered.
    for (int fieldInd(0); fieldInd < numFields; ++fieldInd)
    {
      MDField<ScalarT, Cell, NODE>& field = gatherFields_[fieldInd];
      int indexerId(indexerIds_[fieldInd]),
        subFieldNum(subFieldIds_[fieldInd]);

      // Grab the local data for inputing.
      ArrayRCP<const double> x;
      rcp_dynamic_cast<SpmdVectorBase<double>>(x_->
        getNonconstVectorBlock(indexerId))->getLocalData(ptrFromRef(x));
      auto subRowIndexer = indexers_[indexerId];
      const vector<int>& elmtOffset =
        subRowIndexer->getGIDFieldOffsets(blockId, subFieldNum);
      int numBases(elmtOffset.size());

      // Gather operation for each cell in the workset.
      for (int cell(0); cell < numCells; ++cell)
      {
        size_t cellLocalId(localCellIds[cell]);
        auto LIDs = subRowIndexer->getElementLIDs(cellLocalId);

        // Loop over the basis functions and fill the fields.
        for (int basis(0); basis < numBases; ++basis)
        {
          int offset(elmtOffset[basis]), lid(LIDs[offset]);
          field(cell, basis) = x[lid];
        } // end loop over the basis functions
      } // end loop over localCellIds
    } // end loop over the fields to be gathered
  } // end if (x_.is_null()) or not

  // Deal with the first sensitivities.
  if (firstApplySensitivities_)
  {
    // Loop over the fields to be gathered.
    for (int fieldInd(0); fieldInd < numFields; ++fieldInd)
    {
      MDField<ScalarT, Cell, NODE>& field = gatherFields_[fieldInd];
      int indexerId(indexerIds_[fieldInd]),
        subFieldNum(subFieldIds_[fieldInd]);
      auto subRowIndexer = indexers_[indexerId];
      const vector<int>& elmtOffset =
        subRowIndexer->getGIDFieldOffsets(blockId, subFieldNum);
      int startBlkOffset(blockOffsets[indexerId]),
        numBases(elmtOffset.size());

      // Gather operation for each cell in the workset.
      for (int cell(0); cell < numCells; ++cell)
      {
        // Loop over the basis functions and fill the fields.
        for (int basis(0); basis < numBases; ++basis)
        {
          int offset(elmtOffset[basis]);
          field(cell, basis).fastAccessDx(startBlkOffset + offset) = seedValue;
        } // end loop over the basis functions
      } // end loop over localCellIds
    } // end loop over the fields to be gathered
  } // end if (firstApplySensitivities_)

  // Deal with the second sensitivities.
  if (secondApplySensitivities_)
  {
    // Loop over the fields to be gathered.
    for (int fieldInd(0); fieldInd < numFields; ++fieldInd)
    {
      MDField<ScalarT, Cell, NODE>& field = gatherFields_[fieldInd];
      int indexerId(indexerIds_[fieldInd]),
        subFieldNum(subFieldIds_[fieldInd]);

      // Grab the local data for inputing.
      auto dxEvRoGed = dxBvRoGed_->getGEDBlock(indexerId);
      auto subRowIndexer = indexers_[indexerId];
      const vector<int>& elmtOffset =
        subRowIndexer->getGIDFieldOffsets(blockId, subFieldNum);
      int numBases(elmtOffset.size());

      // Gather operation for each cell in the workset.
      for (int cell(0); cell < numCells; ++cell)
      {
        size_t cellLocalId(localCellIds[cell]);
        auto LIDs = subRowIndexer->getElementLIDs(cellLocalId);

        // Loop over the basis functions and fill the fields.
        for (int basis(0); basis < numBases; ++basis)
        {
          int offset(elmtOffset[basis]), lid(LIDs[offset]);
          field(cell, basis).val().fastAccessDx(0) = (*dxEvRoGed)[lid];
        } // end loop over the basis functions
      } // end loop over localCellIds
    } // end loop over the fields to be gathered
  } // end if (secondApplySensitivities_)
} // end of evaluateFields() (Hessian Specialization)

#endif // Panzer_BUILD_HESSIAN_SUPPORT

#endif // __Panzer_GatherSolution_BlockedEpetra_Hessian_impl_hpp__
