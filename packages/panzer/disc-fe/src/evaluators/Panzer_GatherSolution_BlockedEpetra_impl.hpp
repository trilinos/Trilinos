// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __Panzer_GatherSolution_BlockedEpetra_impl_hpp__
#define   __Panzer_GatherSolution_BlockedEpetra_impl_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Panzer
#include "Panzer_BlockedEpetraLinearObjContainer.hpp"
#include "Panzer_BlockedVector_ReadOnly_GlobalEvaluationData.hpp"
#include "Panzer_GatherSolution_Input.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_GlobalIndexer_Utilities.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"

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
//  Initializing Constructor (Residual Specialization)
//
///////////////////////////////////////////////////////////////////////////////
template<typename TRAITS, typename LO, typename GO>
panzer::
GatherSolution_BlockedEpetra<panzer::Traits::Residual, TRAITS, LO, GO>::
GatherSolution_BlockedEpetra(
  const std::vector<Teuchos::RCP<const GlobalIndexer>>&
    indexers,
  const Teuchos::ParameterList& p)
  :
  indexers_(indexers),
  hasTangentFields_(false)
{
  using panzer::PureBasis;
  using PHX::MDField;
  using PHX::print;
  using std::size_t;
  using std::string;
  using std::vector;
  using Teuchos::RCP;
  using vvstring = std::vector<std::vector<std::string>>;
  GatherSolution_Input input;
  input.setParameterList(p);
  const vector<string>& names             = input.getDofNames();
  RCP<const PureBasis>  basis             = input.getBasis();
  const vvstring&       tangentFieldNames = input.getTangentNames();
  indexerNames_                    = input.getIndexerNames();
  useTimeDerivativeSolutionVector_ = input.useTimeDerivativeSolutionVector();
  globalDataKey_                   = input.getGlobalDataKey();

  // Allocate the fields.
  int numFields(names.size());
  gatherFields_.resize(numFields);
  for (int fd(0); fd < numFields; ++fd)
  {
    gatherFields_[fd] =
      MDField<ScalarT, Cell, NODE>(names[fd], basis->functional);
    this->addEvaluatedField(gatherFields_[fd]);
  } // end loop over names

  // Setup the dependent tangent fields, if requested.
  if (tangentFieldNames.size() > 0)
  {
    TEUCHOS_ASSERT(gatherFields_.size() == tangentFieldNames.size());
    hasTangentFields_ = true;
    tangentFields_.resize(numFields);
    for (int fd(0); fd < numFields; ++fd)
    {
      int numTangentFields(tangentFieldNames[fd].size());
      tangentFields_[fd].resize(numTangentFields);
      for (int i(0); i < numTangentFields; ++i)
      {
        tangentFields_[fd][i] =
          MDField<const ScalarT, Cell, NODE>(tangentFieldNames[fd][i],
          basis->functional);
        this->addDependentField(tangentFields_[fd][i]);
      } // end loop over tangentFieldNames[fd]
    } // end loop over gatherFields_
  } // end if we have tangent fields

  // Figure out what the first active name is.
  string firstName("<none>");
  if (numFields > 0)
    firstName = names[0];
  string n("GatherSolution (BlockedEpetra): " + firstName + " (" +
    print<EvalT>() + ")");
  this->setName(n);
} // end of Initializing Constructor (Residual Specialization)

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup() (Residual Specialization)
//
///////////////////////////////////////////////////////////////////////////////
template<typename TRAITS, typename LO, typename GO>
void
panzer::
GatherSolution_BlockedEpetra<panzer::Traits::Residual, TRAITS, LO, GO>::
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
} // end of postRegistrationSetup() (Residual Specialization)

///////////////////////////////////////////////////////////////////////////////
//
//  preEvaluate() (Residual Specialization)
//
///////////////////////////////////////////////////////////////////////////////
template<typename TRAITS, typename LO, typename GO>
void
panzer::
GatherSolution_BlockedEpetra<panzer::Traits::Residual, TRAITS, LO, GO>::
preEvaluate(
  typename TRAITS::PreEvalData d)
{
  using std::logic_error;
  using std::string;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::typeName;
  using Thyra::ProductVectorBase;
  using BELOC   = panzer::BlockedEpetraLinearObjContainer;
  using BVROGED = panzer::BlockedVector_ReadOnly_GlobalEvaluationData;
  using GED     = panzer::GlobalEvaluationData;

  // First try the refactored ReadOnly container.
  RCP<GED> ged;
  string post(useTimeDerivativeSolutionVector_ ? " - Xdot" : " - X");
  if (d.gedc->containsDataObject(globalDataKey_ + post))
  {
    ged       = d.gedc->getDataObject(globalDataKey_ + post);
    xBvRoGed_ = rcp_dynamic_cast<BVROGED>(ged, true);
    return;
  } // end of the refactored ReadOnly way

  // Now try the old path.
  {
    ged = d.gedc->getDataObject(globalDataKey_);

    // Extract the linear object container.
    auto roGed = rcp_dynamic_cast<const BVROGED>(ged);
    auto beLoc = rcp_dynamic_cast<const BELOC>(ged);
    if (not roGed.is_null())
      xBvRoGed_ = rcp_dynamic_cast<BVROGED>(ged, true);
    else if (not beLoc.is_null())
    {
      if (useTimeDerivativeSolutionVector_)
        x_ = rcp_dynamic_cast<ProductVectorBase<double>>(beLoc->get_dxdt());
      else // if (not useTimeDerivativeSolutionVector_)
        x_ = rcp_dynamic_cast<ProductVectorBase<double>>(beLoc->get_x());
      TEUCHOS_TEST_FOR_EXCEPTION(x_.is_null(), logic_error, "Gather "         \
        "Residual:  Can't find the x_ vector in GEDC \"" << globalDataKey_ <<
        "\" (" << post << ").  A cast failed for " << ged << ".  Type is " <<
        typeName(*ged));
    } // end if we have a roGed or beLoc
  } // end of the old path
} // end of preEvaluate() (Residual Specialization)

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields() (Residual Specialization)
//
///////////////////////////////////////////////////////////////////////////////
template<typename TRAITS, typename LO, typename GO>
void
panzer::
GatherSolution_BlockedEpetra<panzer::Traits::Residual, TRAITS, LO, GO>::
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

  // Loop over the fields to be gathered.
  for (int fieldInd(0); fieldInd < numFields; ++fieldInd)
  {
    MDField<ScalarT, Cell, NODE>& field = gatherFields_[fieldInd];
    auto field_h = Kokkos::create_mirror_view(field.get_static_view());

    int indexerId(indexerIds_[fieldInd]),
      subFieldNum(subFieldIds_[fieldInd]);

    // Grab the local data for inputing.
    ArrayRCP<const double> x;
    Teuchos::RCP<const ReadOnlyVector_GlobalEvaluationData> xEvRoGed;

    if(not x_.is_null()) {
      rcp_dynamic_cast<SpmdVectorBase<double>>(x_->
        getNonconstVectorBlock(indexerId))->getLocalData(ptrFromRef(x));
    }
    else {
      xEvRoGed = xBvRoGed_->getGEDBlock(indexerId);
    }

    auto subRowIndexer = indexers_[indexerId];
    const vector<int>& elmtOffset =
      subRowIndexer->getGIDFieldOffsets(blockId, subFieldNum);
    int numBases(elmtOffset.size());

    auto LIDs = subRowIndexer->getLIDs();
    auto LIDs_h = Kokkos::create_mirror_view(LIDs);
    Kokkos::deep_copy(LIDs_h, LIDs);

    // Gather operation for each cell in the workset.
    for (int cell(0); cell < numCells; ++cell)
    {
      LO cellLocalId = localCellIds[cell];

      // Loop over the basis functions and fill the fields.
      for (int basis(0); basis < numBases; ++basis)
      {
        int offset(elmtOffset[basis]), lid(LIDs_h(cellLocalId, offset));
        if(x_.is_null())
          field_h(cell, basis) = (*xEvRoGed)[lid];
        else
          field_h(cell, basis) = x[lid];
      } // end loop over the basis functions
    } // end loop over localCellIds
    Kokkos::deep_copy(field.get_static_view(), field_h);
  } // end loop over the fields to be gathered
} // end of evaluateFields() (Residual Specialization)

///////////////////////////////////////////////////////////////////////////////
//
//  Initializing Constructor (Tangent Specialization)
//
///////////////////////////////////////////////////////////////////////////////
template<typename TRAITS, typename LO, typename GO>
panzer::GatherSolution_BlockedEpetra<panzer::Traits::Tangent, TRAITS, LO, GO>::
GatherSolution_BlockedEpetra(
  const std::vector<Teuchos::RCP<const GlobalIndexer>>&
    indexers,
  const Teuchos::ParameterList& p)
  :
  indexers_(indexers),
  hasTangentFields_(false)
{
  using panzer::PureBasis;
  using PHX::MDField;
  using PHX::print;
  using std::size_t;
  using std::string;
  using std::vector;
  using Teuchos::RCP;
  using vvstring = std::vector<std::vector<std::string>>;
  GatherSolution_Input input;
  input.setParameterList(p);
  const vector<string>& names             = input.getDofNames();
  RCP<const PureBasis>  basis             = input.getBasis();
  const vvstring&       tangentFieldNames = input.getTangentNames();
  indexerNames_                    = input.getIndexerNames();
  useTimeDerivativeSolutionVector_ = input.useTimeDerivativeSolutionVector();
  globalDataKey_                   = input.getGlobalDataKey();

  // Allocate the fields.
  int numFields(names.size());
  gatherFields_.resize(numFields);
  for (int fd(0); fd < numFields; ++fd)
  {
    gatherFields_[fd] =
      MDField<ScalarT, Cell, NODE>(names[fd], basis->functional);
    this->addEvaluatedField(gatherFields_[fd]);
  } // end loop over names

  // Set up the dependent tangent fields, if requested.
  if (tangentFieldNames.size() > 0)
  {
    TEUCHOS_ASSERT(gatherFields_.size() == tangentFieldNames.size());
    hasTangentFields_ = true;
    tangentFields_.resize(numFields);
    for (int fd(0); fd < numFields; ++fd)
    {
      int numTangentFields(tangentFieldNames[fd].size());
      tangentFields_[fd].resize(numTangentFields);
      for (int i(0); i < numTangentFields; ++i)
      {
        tangentFields_[fd][i] =
          MDField<const ScalarT, Cell, NODE>(tangentFieldNames[fd][i],
          basis->functional);
        this->addDependentField(tangentFields_[fd][i]);
      } // end loop over tangentFieldNames
    } // end loop over gatherFields_
  } // end if we have tangent fields

  // Figure out what the first active name is.
  string firstName("<none>");
  if (numFields > 0)
    firstName = names[0];
  string n("GatherSolution Tangent (BlockedEpetra): " + firstName + " (" +
    print<EvalT>() + ")");
  this->setName(n);
} // end of Initializing Constructor (Tangent Specialization)

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup() (Tangent Specialization)
//
///////////////////////////////////////////////////////////////////////////////
template<typename TRAITS, typename LO, typename GO>
void
panzer::GatherSolution_BlockedEpetra<panzer::Traits::Tangent, TRAITS, LO, GO>::
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
} // end of postRegistrationSetup() (Tangent Specialization)

///////////////////////////////////////////////////////////////////////////////
//
//  preEvaluate() (Tangent Specialization)
//
///////////////////////////////////////////////////////////////////////////////
template<typename TRAITS, typename LO, typename GO>
void
panzer::GatherSolution_BlockedEpetra<panzer::Traits::Tangent, TRAITS, LO, GO>::
preEvaluate(
  typename TRAITS::PreEvalData d)
{
  using std::logic_error;
  using std::string;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::typeName;
  using Thyra::ProductVectorBase;
  using BELOC   = panzer::BlockedEpetraLinearObjContainer;
  using BVROGED = panzer::BlockedVector_ReadOnly_GlobalEvaluationData;
  using GED     = panzer::GlobalEvaluationData;

  // First try the refactored ReadOnly container.
  RCP<GED> ged;
  string post(useTimeDerivativeSolutionVector_ ? " - Xdot" : " - X");
  if (d.gedc->containsDataObject(globalDataKey_ + post))
  {
    ged       = d.gedc->getDataObject(globalDataKey_ + post);
    xBvRoGed_ = rcp_dynamic_cast<BVROGED>(ged, true);
    return;
  } // end of the refactored ReadOnly way

  // Now try the old path.
  {
    ged = d.gedc->getDataObject(globalDataKey_);

    // Extract the linear object container.
    auto roGed = rcp_dynamic_cast<const BVROGED>(ged);
    auto beLoc = rcp_dynamic_cast<const BELOC>(ged);
    if (not roGed.is_null())
      xBvRoGed_ = rcp_dynamic_cast<BVROGED>(ged, true);
    else if (not beLoc.is_null())
    {
      if (useTimeDerivativeSolutionVector_)
        x_ = rcp_dynamic_cast<ProductVectorBase<double>>(beLoc->get_dxdt());
      else // if (not useTimeDerivativeSolutionVector_)
        x_ = rcp_dynamic_cast<ProductVectorBase<double>>(beLoc->get_x());
      TEUCHOS_TEST_FOR_EXCEPTION(x_.is_null(), logic_error, "Gather "         \
        "Tangent:  Can't find the x_ vector in GEDC \"" << globalDataKey_ <<
        "\" (" << post << ").  A cast failed for " << ged << ".  Type is " <<
        typeName(*ged));
    } // end if we have a roGed or beLoc
  } // end of the old path
} // end of preEvaluate() (Tangent Specialization)

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields() (Tangent Specialization)
//
///////////////////////////////////////////////////////////////////////////////
template<typename TRAITS, typename LO, typename GO>
void
panzer::GatherSolution_BlockedEpetra<panzer::Traits::Tangent, TRAITS, LO, GO>::
evaluateFields(
  typename TRAITS::EvalData workset)
{
  using PHX::MDField;
  using std::pair;
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
  vector<pair<int, GO>> GIDs;
  string blockId(this->wda(workset).block_id);
  const vector<size_t>& localCellIds = this->wda(workset).cell_local_ids;
  int numFields(gatherFields_.size()), numCells(localCellIds.size());

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
        LO cellLocalId = localCellIds[cell];
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

  // Deal with the tangent fields.
  if (hasTangentFields_)
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
      int numBases(elmtOffset.size());

      // Gather operation for each cell in the workset.
      for (int cell(0); cell < numCells; ++cell)
      {
        // Loop over the basis functions and fill the fields.
        for (int basis(0); basis < numBases; ++basis)
        {
          int numTangentFields(tangentFields_[fieldInd].size());
          for (int i(0); i < numTangentFields; ++i)
            field(cell, basis).fastAccessDx(i) =
              tangentFields_[fieldInd][i](cell, basis).val();
        } // end loop over the basis functions
      } // end loop over localCellIds
    } // end loop over the fields to be gathered
  } // end if (hasTangentFields_)
} // end of evaluateFields() (Tangent Specialization)

///////////////////////////////////////////////////////////////////////////////
//
//  Initializing Constructor (Jacobian Specialization)
//
///////////////////////////////////////////////////////////////////////////////
template<typename TRAITS, typename LO, typename GO>
panzer::
GatherSolution_BlockedEpetra<panzer::Traits::Jacobian, TRAITS, LO, GO>::
GatherSolution_BlockedEpetra(
  const std::vector<Teuchos::RCP<const GlobalIndexer>>&
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
  disableSensitivities_            = not input.firstSensitivitiesAvailable();

  // Allocate the fields.
  int numFields(names.size());
  gatherFields_.resize(numFields);
  for (int fd(0); fd < numFields; ++fd)
  {
    MDField<ScalarT, Cell, NODE> f(names[fd], basis->functional);
    gatherFields_[fd] = f;
    this->addEvaluatedField(gatherFields_[fd]);
  } // end loop over names

  // Figure out what the first active name is.
  string firstName("<none>"), n("GatherSolution (BlockedEpetra");
  if (numFields > 0)
    firstName = names[0];
  if (disableSensitivities_)
    n += ", No Sensitivities";
  n += "):  " + firstName + " (" + print<EvalT>() + ")";
  this->setName(n);
} // end of Initializing Constructor (Jacobian Specialization)

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup() (Jacobian Specialization)
//
///////////////////////////////////////////////////////////////////////////////
template<typename TRAITS, typename LO, typename GO>
void
panzer::
GatherSolution_BlockedEpetra<panzer::Traits::Jacobian, TRAITS, LO, GO>::
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
} // end of postRegistrationSetup() (Jacobian Specialization)

///////////////////////////////////////////////////////////////////////////////
//
//  preEvaluate() (Jacobian Specialization)
//
///////////////////////////////////////////////////////////////////////////////
template<typename TRAITS, typename LO, typename GO>
void
panzer::
GatherSolution_BlockedEpetra<panzer::Traits::Jacobian, TRAITS, LO, GO>::
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
  using BELOC   = panzer::BlockedEpetraLinearObjContainer;
  using BVROGED = panzer::BlockedVector_ReadOnly_GlobalEvaluationData;
  using GED     = panzer::GlobalEvaluationData;

  // Manage sensitivities.
  applySensitivities_ = false;
  if ((not disableSensitivities_                       )  and
      (d.first_sensitivities_name == sensitivitiesName_))
    applySensitivities_ = true;

  // First try the refactored ReadOnly container.
  RCP<GED> ged;
  string post(useTimeDerivativeSolutionVector_ ? " - Xdot" : " - X");
  if (d.gedc->containsDataObject(globalDataKey_ + post))
  {
    ged       = d.gedc->getDataObject(globalDataKey_ + post);
    xBvRoGed_ = rcp_dynamic_cast<BVROGED>(ged, true);
    return;
  } // end of the refactored ReadOnly way

  // Now try the old path.
  {
    ged = d.gedc->getDataObject(globalDataKey_);

    // Extract the linear object container.
    auto roGed = rcp_dynamic_cast<const BVROGED>(ged);
    auto beLoc = rcp_dynamic_cast<const BELOC>(ged);
    if (not roGed.is_null())
      xBvRoGed_ = rcp_dynamic_cast<BVROGED>(ged, true);
    else if (not beLoc.is_null())
    {
      if (useTimeDerivativeSolutionVector_)
        x_ = rcp_dynamic_cast<ProductVectorBase<double>>(beLoc->get_dxdt());
      else // if (not useTimeDerivativeSolutionVector_)
        x_ = rcp_dynamic_cast<ProductVectorBase<double>>(beLoc->get_x());
      TEUCHOS_TEST_FOR_EXCEPTION(x_.is_null(), logic_error, "Gather "         \
        "Jacobian:  Can't find x vector in GEDC \"" << globalDataKey_ <<
        "\" (" << post << ").  A cast failed for " << ged << ".  Type is " <<
        typeName(*ged));
    } // end if we have a roGed or beLoc
  } // end of the old path
} // end of preEvaluate() (Jacobian Specialization)

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields() (Jacobian Specialization)
//
///////////////////////////////////////////////////////////////////////////////
template<typename TRAITS, typename LO, typename GO>
void
panzer::
GatherSolution_BlockedEpetra<panzer::Traits::Jacobian, TRAITS, LO, GO>::
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
  if (applySensitivities_)
  {
    if ((useTimeDerivativeSolutionVector_) and (gatherSeedIndex_ < 0))
      seedValue = workset.alpha;
    else if (gatherSeedIndex_ < 0)
      seedValue = workset.beta;
    else if (not useTimeDerivativeSolutionVector_)
      seedValue = workset.gather_seeds[gatherSeedIndex_];
    else
      TEUCHOS_ASSERT(false);
  } // end if (applySensitivities_)

  // Turn off sensitivies:  This may be faster if we don't expand the term, but
  // I suspect not, because anywhere it is used the full complement of
  // sensitivies will be needed anyway.
  if (not applySensitivities_)
    seedValue = 0.0;

  vector<int> blockOffsets;
  computeBlockOffsets(blockId, indexers_, blockOffsets);
  int numDerivs(blockOffsets[blockOffsets.size() - 1]);

  // NOTE:  A reordering of these loops will likely improve performance.  The
  //        "getGIDFieldOffsets may be expensive.  However the "getElementGIDs"
  //        can be cheaper.  However the lookup for LIDs may be more expensive!

   // Loop over the fields to be gathered.
  for (int fieldInd(0); fieldInd < numFields; ++fieldInd)
  {
    MDField<ScalarT, Cell, NODE>& field = gatherFields_[fieldInd];
    auto field_h = Kokkos::create_mirror_view(field.get_view());

    int indexerId(indexerIds_[fieldInd]), subFieldNum(subFieldIds_[fieldInd]);

    // Grab the local data for inputing.
    ArrayRCP<const double> x;
    Teuchos::RCP<const ReadOnlyVector_GlobalEvaluationData> xEvRoGed;
    if(not x_.is_null()) {
      rcp_dynamic_cast<SpmdVectorBase<double>>(x_->getNonconstVectorBlock(indexerId))->getLocalData(ptrFromRef(x));
    }else {
      xEvRoGed = xBvRoGed_->getGEDBlock(indexerId);
    }

    auto subRowIndexer = indexers_[indexerId];
    const vector<int>& elmtOffset =
      subRowIndexer->getGIDFieldOffsets(blockId, subFieldNum);
    int startBlkOffset(blockOffsets[indexerId]), numBases(elmtOffset.size());

    auto LIDs = subRowIndexer->getLIDs();
    auto LIDs_h = Kokkos::create_mirror_view(LIDs);
    Kokkos::deep_copy(LIDs_h, LIDs);

    // Gather operation for each cell in the workset.
    for (int cell(0); cell < numCells; ++cell)
    {
      LO cellLocalId = localCellIds[cell];

      // Loop over the basis functions and fill the fields.
      for (int basis(0); basis < numBases; ++basis)
      {
        // Set the value and seed the FAD object.
        int offset(elmtOffset[basis]), lid(LIDs_h(cellLocalId, offset));
        if(x_.is_null())
          field_h(cell, basis) = ScalarT(numDerivs, (*xEvRoGed)[lid]);
        else
          field_h(cell, basis) = ScalarT(numDerivs, x[lid]);

        field_h(cell, basis).fastAccessDx(startBlkOffset + offset) = seedValue;
      } // end loop over the basis functions
    } // end loop over localCellIds
    Kokkos::deep_copy(field.get_static_view(), field_h);
  } // end loop over the fields to be gathered
} // end of evaluateFields() (Jacobian Specialization)

#endif // __Panzer_GatherSolution_BlockedEpetra_impl_hpp__
