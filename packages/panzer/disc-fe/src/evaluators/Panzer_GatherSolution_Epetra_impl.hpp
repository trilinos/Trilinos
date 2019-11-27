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

#ifndef   __Panzer_GatherSolution_Epetra_impl_hpp__
#define   __Panzer_GatherSolution_Epetra_impl_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Epetra
#include "Epetra_Vector.h"

// Panzer
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Panzer_EpetraVector_ReadOnly_GlobalEvaluationData.hpp"
#include "Panzer_GatherSolution_Input.hpp"
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"

// Teuchos
#include "Teuchos_Assert.hpp"

// Thyra
#include "Thyra_SpmdVectorBase.hpp"

///////////////////////////////////////////////////////////////////////////////
//
//  Initializing Constructor (Residual Specialization)
//
///////////////////////////////////////////////////////////////////////////////
template<typename TRAITS, typename LO, typename GO>
panzer::GatherSolution_Epetra<panzer::Traits::Residual, TRAITS, LO, GO>::
GatherSolution_Epetra(
  const Teuchos::RCP<const panzer::GlobalIndexer>& indexer,
  const Teuchos::ParameterList& p)
  :
  globalIndexer_(indexer),
  hasTangentFields_(false)
{
  using panzer::PureBasis;
  using PHX::MDField;
  using PHX::print;
  using std::size_t;
  using std::vector;
  using std::string;
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
    TEUCHOS_ASSERT(gatherFields_.size() == tangentFieldNames.size())
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
    } // end loop over gatherFields
  } // end if (tangentFieldNames.size() > 0)

  // Figure out what the first active name is.
  string firstName("<none>");
  if (numFields > 0)
    firstName = names[0];
  string n("GatherSolution (Epetra): " + firstName + " (" +
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
panzer::GatherSolution_Epetra<panzer::Traits::Residual, TRAITS, LO, GO>::
postRegistrationSetup(
  typename TRAITS::SetupData /* d  */,
  PHX::FieldManager<TRAITS>& /* fm */)
{
  using std::logic_error;
  using std::size_t;
  using std::string;
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_.size())
  int numFields(gatherFields_.size());
  fieldIds_.resize(numFields);
  for (int fd(0); fd < numFields; ++fd)
  {
    // Get the field ID from the DOF manager.
    const string& fieldName(indexerNames_[fd]);
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    // This is the error return code; raise the alarm.
    TEUCHOS_TEST_FOR_EXCEPTION(fieldIds_[fd] == -1, logic_error,
      "GatherSolution_Epetra<Residual>: Could not find field \"" +
      fieldName + "\" in the global indexer. ")
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
panzer::GatherSolution_Epetra<panzer::Traits::Residual, TRAITS, LO, GO>::
preEvaluate(
  typename TRAITS::PreEvalData d)
{
  using std::string;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using ELOC    = panzer::EpetraLinearObjContainer;
  using EVROGED = panzer::EpetraVector_ReadOnly_GlobalEvaluationData;
  using GED     = panzer::GlobalEvaluationData;
  using LOC     = panzer::LinearObjContainer;
  using LPGED   = panzer::LOCPair_GlobalEvaluationData;

  // First try the refactored ReadOnly container.
  RCP<GED> ged;
  string post(useTimeDerivativeSolutionVector_ ? " - Xdot" : " - X");
  if (d.gedc->containsDataObject(globalDataKey_ + post))
  {
    ged       = d.gedc->getDataObject(globalDataKey_ + post);
    xEvRoGed_ = rcp_dynamic_cast<EVROGED>(ged, true);
    return;
  } // end of the refactored ReadOnly way

  // Now try the old path.
  ged = d.gedc->getDataObject(globalDataKey_);
  {
    // Try to extract the linear object container.
    auto epetraContainer = rcp_dynamic_cast<ELOC>(ged);
    auto locPair         = rcp_dynamic_cast<LPGED>(ged);
    if (not locPair.is_null())
    {
      RCP<LOC> loc = locPair->getGhostedLOC();
      epetraContainer = rcp_dynamic_cast<ELOC>(loc);
    } // end if (not locPair.is_null())
    if (not epetraContainer.is_null())
    {
      if (useTimeDerivativeSolutionVector_)
        x_ = epetraContainer->get_dxdt();
      else // if (not useTimeDerivativeSolutionVector_)
        x_ = epetraContainer->get_x();
      return;
    } // end if (not epetraContainer.is_null())
  } // end of the old path

  // As a last resort, try to extract an EpetraVector_ReadOnly object.  This
  // will throw an exception if it doesn't work.
  xEvRoGed_ = rcp_dynamic_cast<EVROGED>(ged, true);
} // end of preEvaluate() (Residual Specialization)

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields() (Residual Specialization)
//
///////////////////////////////////////////////////////////////////////////////
template<typename TRAITS, typename LO, typename GO>
void
panzer::GatherSolution_Epetra<panzer::Traits::Residual, TRAITS, LO, GO>::
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

  // For convenience, pull out some objects from the workset.
  string blockId(this->wda(workset).block_id);
  const vector<size_t>& localCellIds = this->wda(workset).cell_local_ids;
  int numCells(localCellIds.size()), numFields(gatherFields_.size());

  // NOTE:  A reordering of these loops will likely improve performance.  The
  //        "getGIDFieldOffsets may be expensive.  However the "getElementGIDs"
  //        can be cheaper.  However the lookup for LIDs may be more expensive!
  if (x_.is_null())
  {
    // Gather operation for each cell in the workset.
    for (int cell(0); cell < numCells; ++cell)
    {
      size_t cellLocalId(localCellIds[cell]);
      auto LIDs = globalIndexer_->getElementLIDs(cellLocalId);

      // Loop over the fields to be gathered.
      for (int fieldInd(0); fieldInd < numFields; ++fieldInd)
      {
        MDField<ScalarT, Cell, NODE>& field = gatherFields_[fieldInd];
        int fieldNum(fieldIds_[fieldInd]);
        const vector<int>& elmtOffset =
          globalIndexer_->getGIDFieldOffsets(blockId, fieldNum);
        int numBases(elmtOffset.size());

        // Loop over the basis functions and fill the fields.
        for (int basis(0); basis < numBases; ++basis)
        {
          int offset(elmtOffset[basis]), lid(LIDs[offset]);
          field(cell, basis) = (*xEvRoGed_)[lid];
        } // end loop over the basis functions
      } // end loop over the fields to be gathered
    } // end loop over localCellIds
  }
  else // if (not x_.is_null())
  {
    // Gather operation for each cell in the workset.
    for (int cell(0); cell < numCells; ++cell)
    {
      size_t cellLocalId(localCellIds[cell]);
      auto LIDs = globalIndexer_->getElementLIDs(cellLocalId);

      // Loop over the fields to be gathered.
      for (int fieldInd(0); fieldInd < numFields; ++fieldInd)
      {
        MDField<ScalarT, Cell, NODE>& field = gatherFields_[fieldInd];
        int fieldNum(fieldIds_[fieldInd]);
        const vector<int>& elmtOffset =
          globalIndexer_->getGIDFieldOffsets(blockId, fieldNum);
        int numBases(elmtOffset.size());

        // Loop over the basis functions and fill the fields.
        for (int basis(0); basis < numBases; ++basis)
        {
          int offset(elmtOffset[basis]), lid(LIDs[offset]);
          field(cell, basis) = (*x_)[lid];
        } // end loop over the basis functions
      } // end loop over the fields to be gathered
    } // end loop over localCellIds
  } // end if (x_.is_null()) or not
} // end of evaluateFields() (Residual Specialization)

///////////////////////////////////////////////////////////////////////////////
//
//  Initializing Constructor (Tangent Specialization)
//
///////////////////////////////////////////////////////////////////////////////
template<typename TRAITS, typename LO, typename GO>
panzer::GatherSolution_Epetra<panzer::Traits::Tangent, TRAITS, LO, GO>::
GatherSolution_Epetra(
  const Teuchos::RCP<const panzer::GlobalIndexer>& indexer,
  const Teuchos::ParameterList& p)
  :
  globalIndexer_(indexer),
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
    TEUCHOS_ASSERT(gatherFields_.size() == tangentFieldNames.size())
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
      } // end loop over tangeng_field_names[fd]
    } // end loop over gatherFields_
  } // end if (tangentFieldNames.size() > 0)

  // Figure out what the first active name is.
  string firstName("<none>");
  if (numFields > 0)
    firstName = names[0];
  string n("GatherSolution (Epetra): " + firstName + " (" +
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
panzer::GatherSolution_Epetra<panzer::Traits::Tangent, TRAITS, LO, GO>::
postRegistrationSetup(
  typename TRAITS::SetupData /* d  */,
  PHX::FieldManager<TRAITS>& /* fm */)
{
  using std::logic_error;
  using std::size_t;
  using std::string;
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_.size())
  int numFields(gatherFields_.size());
  fieldIds_.resize(numFields);
  for (int fd(0); fd < numFields; ++fd)
  {
    // Get the field ID from the DOF manager.
    const string& fieldName(indexerNames_[fd]);
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    // This is the error return code; raise the alarm.
    TEUCHOS_TEST_FOR_EXCEPTION(fieldIds_[fd] == -1, logic_error,
      "GatherSolution_Epetra<Tangent>:  Could not find field \"" + fieldName
      + "\" in the global indexer. ")
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
panzer::GatherSolution_Epetra<panzer::Traits::Tangent, TRAITS, LO, GO>::
preEvaluate(
  typename TRAITS::PreEvalData d)
{
  using std::string;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using ELOC    = panzer::EpetraLinearObjContainer;
  using EVROGED = panzer::EpetraVector_ReadOnly_GlobalEvaluationData;
  using GED     = panzer::GlobalEvaluationData;
  using LOC     = panzer::LinearObjContainer;
  using LPGED   = panzer::LOCPair_GlobalEvaluationData;

  // First try the refactored ReadOnly container.
  RCP<GED> ged;
  string post(useTimeDerivativeSolutionVector_ ? " - Xdot" : " - X");
  if (d.gedc->containsDataObject(globalDataKey_ + post))
  {
    ged       = d.gedc->getDataObject(globalDataKey_ + post);
    xEvRoGed_ = rcp_dynamic_cast<EVROGED>(ged, true);
    return;
  } // end of the refactored ReadOnly way

  // Now try the old path.
  ged = d.gedc->getDataObject(globalDataKey_);
  {
    // Try to extract the linear object container.
    auto epetraContainer = rcp_dynamic_cast<ELOC>(ged);
    auto locPair         = rcp_dynamic_cast<LPGED>(ged);
    if (not locPair.is_null())
    {
      RCP<LOC> loc = locPair->getGhostedLOC();
      epetraContainer = rcp_dynamic_cast<ELOC>(loc);
    } // end if (not locPair.is_null())
    if (not epetraContainer.is_null())
    {
      if (useTimeDerivativeSolutionVector_)
        x_ = epetraContainer->get_dxdt();
      else // if (not useTimeDerivativeSolutionVector_)
        x_ = epetraContainer->get_x();
      return;
    } // end if (not epetraContainer.is_null())
  } // end of the old path

  // As a last resort, try to extract an EpetraVector_ReadOnly object.  This
  // will throw an exception if it doesn't work.
  xEvRoGed_ = rcp_dynamic_cast<EVROGED>(ged, true);
} // end of preEvaluate() (Tangent Specialization)

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields() (Tangent Specialization)
//
///////////////////////////////////////////////////////////////////////////////
template<typename TRAITS, typename LO, typename GO>
void
panzer::GatherSolution_Epetra<panzer::Traits::Tangent, TRAITS, LO, GO>::
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

  // For convenience, pull out some objects from the workset.
  string blockId(this->wda(workset).block_id);
  const vector<size_t>& localCellIds = this->wda(workset).cell_local_ids;
  int numCells(localCellIds.size()), numFields(gatherFields_.size());

  // NOTE:  A reordering of these loops will likely improve performance.  The
  //        "getGIDFieldOffsets may be expensive.  However the "getElementGIDs"
  //        can be cheaper.  However the lookup for LIDs may be more expensive!
  if (x_.is_null())
  {
    // Gather operation for each cell in the workset.
    for (int cell(0); cell < numCells; ++cell)
    {
      size_t cellLocalId(localCellIds[cell]);
      auto LIDs = globalIndexer_->getElementLIDs(cellLocalId);

      // Loop over the fields to be gathered.
      for (int fieldInd(0); fieldInd < numFields; ++fieldInd)
      {
        MDField<ScalarT, Cell, NODE>& field = gatherFields_[fieldInd];
        int fieldNum(fieldIds_[fieldInd]);
        const vector<int>& elmtOffset =
          globalIndexer_->getGIDFieldOffsets(blockId, fieldNum);
        int numBases(elmtOffset.size());

        // Loop over the basis functions and fill the fields.
        for (int basis(0); basis < numBases; ++basis)
        {
          int offset(elmtOffset[basis]), lid(LIDs[offset]);
          field(cell, basis) = (*xEvRoGed_)[lid];
        } // end loop over the basis functions
      } // end loop over the fields to be gathered
    } // end loop over localCellIds
  }
  else // if (not x_.is_null())
  {
    // Gather operation for each cell in the workset.
    for (int cell(0); cell < numCells; ++cell)
    {
      size_t cellLocalId(localCellIds[cell]);
      auto LIDs = globalIndexer_->getElementLIDs(cellLocalId);

      // Loop over the fields to be gathered.
      for (int fieldInd(0); fieldInd < numFields; ++fieldInd)
      {
        MDField<ScalarT, Cell, NODE>& field = gatherFields_[fieldInd];
        int fieldNum(fieldIds_[fieldInd]);
        const vector<int>& elmtOffset =
          globalIndexer_->getGIDFieldOffsets(blockId, fieldNum);
        int numBases(elmtOffset.size());

        // Loop over the basis functions and fill the fields.
        for (int basis(0); basis < numBases; ++basis)
        {
          int offset(elmtOffset[basis]), lid(LIDs[offset]);
          field(cell, basis) = (*x_)[lid];
        } // end loop over the basis functions
      } // end loop over the fields to be gathered
    } // end loop over localCellIds
  } // end if (x_.is_null()) or not

  // Deal with the tangent fields.
  if (hasTangentFields_)
  {
    // Loop over the cells in the workset.
    for (int cell(0); cell < numCells; ++cell)
    {
      // Loop over the fields to be gathered.
      for (int fieldInd(0); fieldInd < numFields; ++fieldInd)
      {
        MDField<ScalarT, Cell, NODE>& field = gatherFields_[fieldInd];
        int fieldNum(fieldIds_[fieldInd]);
        const vector<int>& elmtOffset =
          globalIndexer_->getGIDFieldOffsets(blockId, fieldNum);
        int numBases(elmtOffset.size());

        // Loop over the basis functions.
        for (int basis(0); basis < numBases; ++basis)
        {
          // Loop over the tangent fields and fill them in.
          int numTangentFields(tangentFields_[fieldInd].size());
          for (int i(0); i < numTangentFields; ++i)
            field(cell, basis).fastAccessDx(i) =
              tangentFields_[fieldInd][i](cell, basis).val();
        } // end loop over the basis functions
      } // end loop over the fields to be gathered
    } // end loop over the cells in the workset
  } // end if (hasTangentFields_)
} // end of evaluateFields() (Tangent Specialization)

///////////////////////////////////////////////////////////////////////////////
//
//  Initializing Constructor (Jacobian Specialization)
//
///////////////////////////////////////////////////////////////////////////////
template<typename TRAITS, typename LO, typename GO>
panzer::GatherSolution_Epetra<panzer::Traits::Jacobian, TRAITS, LO, GO>::
GatherSolution_Epetra(
  const Teuchos::RCP<const panzer::GlobalIndexer>& indexer,
  const Teuchos::ParameterList& p)
  :
  globalIndexer_(indexer)
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
  string firstName("<none>"), n("GatherSolution (Epetra");
  if (numFields > 0)
    firstName = names[0];
  if (disableSensitivities_)
    n += ", No Sensitivities";
  n += "): " + firstName + " (" + print<EvalT>() + ")";
  this->setName(n);
} // end of Initializing Constructor (Jacobian Specialization)

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup() (Jacobian Specialization)
//
///////////////////////////////////////////////////////////////////////////////
template<typename TRAITS, typename LO, typename GO>
void
panzer::GatherSolution_Epetra<panzer::Traits::Jacobian, TRAITS, LO, GO>::
postRegistrationSetup(
  typename TRAITS::SetupData /* d  */,
  PHX::FieldManager<TRAITS>& /* fm */)
{
  using std::logic_error;
  using std::size_t;
  using std::string;
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_.size())
  int numFields(gatherFields_.size());
  fieldIds_.resize(numFields);
  for (int fd(0); fd < numFields; ++fd)
  {
    // Get the field ID from the DOF manager.
    const string& fieldName(indexerNames_[fd]);
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    // This is the error return code; raise the alarm.
    TEUCHOS_TEST_FOR_EXCEPTION(fieldIds_[fd] == -1, logic_error,
      "GatherSolution_Epetra<Jacobian>: Could not find field \"" +
      fieldName + "\" in the global indexer. ")
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
panzer::GatherSolution_Epetra<panzer::Traits::Jacobian, TRAITS, LO, GO>::
preEvaluate(
  typename TRAITS::PreEvalData d)
{
  using std::string;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using ELOC    = panzer::EpetraLinearObjContainer;
  using EVROGED = panzer::EpetraVector_ReadOnly_GlobalEvaluationData;
  using GED     = panzer::GlobalEvaluationData;
  using LOC     = panzer::LinearObjContainer;
  using LPGED   = panzer::LOCPair_GlobalEvaluationData;

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
    xEvRoGed_ = rcp_dynamic_cast<EVROGED>(ged, true);
    return;
  } // end of the refactored ReadOnly way

  // Now try the old path.
  ged = d.gedc->getDataObject(globalDataKey_);
  {
    // Try to extract the linear object container.
    auto epetraContainer = rcp_dynamic_cast<ELOC>(ged);
    auto locPair         = rcp_dynamic_cast<LPGED>(ged);
    if (not locPair.is_null())
    {
      RCP<LOC> loc = locPair->getGhostedLOC();
      epetraContainer = rcp_dynamic_cast<ELOC>(loc);
    } // end if (not locPair.is_null())
    if (not epetraContainer.is_null())
    {
      if (useTimeDerivativeSolutionVector_)
        x_ = epetraContainer->get_dxdt();
      else // if (not useTimeDerivativeSolutionVector_)
        x_ = epetraContainer->get_x();
      return;
    } // end if (not epetraContainer.is_null())
  } // end of the old path

  // As a last resort, try to extract an EpetraVector_ReadOnly object.  This
  // will throw an exception if it doesn't work.
  xEvRoGed_ = rcp_dynamic_cast<EVROGED>(ged, true);
} // end of preEvaluate() (Jacobian Specialization)

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields() (Jacobian Specialization)
//
///////////////////////////////////////////////////////////////////////////////
template<typename TRAITS, typename LO, typename GO>
void
panzer::GatherSolution_Epetra<panzer::Traits::Jacobian, TRAITS, LO, GO>::
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

  // For convenience, pull out some objects from the workset.
  string blockId(this->wda(workset).block_id);
  const vector<size_t>& localCellIds = this->wda(workset).cell_local_ids;
  int numFields(gatherFields_.size()), numCells(localCellIds.size());

  // Set a sensitivity seed value.
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
      TEUCHOS_ASSERT(false)
  } // end if (applySensitivities_)

  // Turn off sensitivies.  This may be faster if we don't expand the term, but
  // I suspect not, because anywhere it is used the full complement of
  // sensitivies will be needed anyway.
  if (not applySensitivities_)
    seedValue = 0.0;

  // Interface worksets handle DOFs from two element blocks.  The derivative
  // offset for the other element block must be shifted by the derivative side
  // of my element block.
  int dos(0);
  if (this->wda.getDetailsIndex() == 1)
  {
    // Get the DOF count for my element block.
    dos = globalIndexer_->getElementBlockGIDCount(workset.details(0).block_id);
  } // end if (this->wda.getDetailsIndex() == 1)

  // NOTE:  A reordering of these loops will likely improve performance.  The
  //        "getGIDFieldOffsets may be expensive.  However the "getElementGIDs"
  //        can be cheaper.  However the lookup for LIDs may be more expensive!
  if (x_.is_null())
  {
    // Loop over the fields to be gathered.
    for (int fieldInd(0); fieldInd < numFields; ++fieldInd)
    {
      MDField<ScalarT, Cell, NODE>& field = gatherFields_[fieldInd];
      int fieldNum(fieldIds_[fieldInd]);
      const vector<int>& elmtOffset =
        globalIndexer_->getGIDFieldOffsets(blockId, fieldNum);
      int numBases(elmtOffset.size());

      // Gather operation for each cell in the workset.
      for (int cell(0); cell < numCells; ++cell)
      {
        size_t cellLocalId(localCellIds[cell]);
	auto LIDs = globalIndexer_->getElementLIDs(cellLocalId);

        // Loop over the basis functions and fill the fields.
        for (int basis(0); basis < numBases; ++basis)
        {
          int offset(elmtOffset[basis]), lid(LIDs[offset]);
          field(cell, basis) = (*xEvRoGed_)[lid];
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
      int fieldNum(fieldIds_[fieldInd]);
      const vector<int>& elmtOffset =
        globalIndexer_->getGIDFieldOffsets(blockId, fieldNum);
      int numBases(elmtOffset.size());

      // Gather operation for each cell in the workset.
      for (int cell(0); cell < numCells; ++cell)
      {
        size_t cellLocalId(localCellIds[cell]);
	auto LIDs = globalIndexer_->getElementLIDs(cellLocalId);

        // Loop over the basis functions and fill the fields.
        for (int basis(0); basis < numBases; ++basis)
        {
          int offset(elmtOffset[basis]), lid(LIDs[offset]);
          field(cell, basis) = (*x_)[lid];
        } // end loop over the basis functions
      } // end loop over localCellIds
    } // end loop over the fields to be gathered
  } // end if (x_.is_null()) or not

  // Deal with the sensitivities.
  if (applySensitivities_)
  {
    // Loop over the fields to be gathered.
    for (int fieldInd(0); fieldInd < numFields; ++fieldInd)
    {
      MDField<ScalarT, Cell, NODE>& field = gatherFields_[fieldInd];
      int fieldNum(fieldIds_[fieldInd]);
      const vector<int>& elmtOffset =
        globalIndexer_->getGIDFieldOffsets(blockId, fieldNum);
      int numBases(elmtOffset.size());

      // Gather operation for each cell in the workset.
      for (int cell(0); cell < numCells; ++cell)
      {
        // Loop over the basis functions.
        for (int basis(0); basis < numBases; ++basis)
        {
          // Seed the FAD object.
          int offset(elmtOffset[basis]);
          field(cell, basis).fastAccessDx(dos + offset) = seedValue;
        } // end loop over the basis functions
      } // end loop over localCellIds
    } // end loop over the fields to be gathered
  } // end if (applySensitivities_)
} // end of evaluateFields() (Jacobian Specialization)

#endif // __Panzer_GatherSolution_Epetra_impl_hpp__
