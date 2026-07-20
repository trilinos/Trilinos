/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_transfer_util/EntityCentroidRecoverField.hpp"
#include <stk_transfer_util/RecoverField.hpp>     // for RecoverField, RecoverFiel...
#include <ostream>                          // for operator<<, basic_ostream...
#include <string>                           // for char_traits, operator<<
#include <vector>                           // for vector, vector<>::const_i...

#include "stk_search_util/ObjectCoordinates.hpp"               // for compute_entity_centroid
#include "stk_mesh/base/Entity.hpp"         // for Entity
#include "stk_mesh/base/MetaData.hpp"       // for MetaData
#include "stk_util/util/ReportHandler.hpp"  // for STK_ThrowRequireMsg
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace transfer {

EntityCentroidLinearRecoverField::EntityCentroidLinearRecoverField(RecoverField::RecoveryType recType,
                                                                   const std::vector<std::shared_ptr<stk::search::CachedFieldDataBase>>& recVars,
                                                                   const std::shared_ptr<stk::search::CachedFieldDataBase>& recNodeVar, int nSampEntity,
                                                                   stk::mesh::Entity entity, FieldTransform transform)
  : RecoverField(recType, nSampEntity)
  , m_recoverVars(recVars)
  , m_nodeVar(recNodeVar)
  , m_transform(transform)
{
  m_totalVarComponents = 0;

  for(auto const& var : recVars) {
    m_totalVarComponents += stk::mesh::field_scalars_per_entity(*var->m_field, entity);
  }
}

EntityCentroidLinearRecoverField::~EntityCentroidLinearRecoverField() {}

void EntityCentroidLinearRecoverField::sample_patch(const std::vector<stk::mesh::Entity>& patch, int nSampPatch,
                                                    double* fieldSample, double* basisSample) const
{
  // Function meant for piecewise linear variables with sampling at centroid of mesh entity
  STK_ThrowRequireMsg(m_nSampleElements == 1, "The EntityCentroidLinear method only works with piecewise"
                                      << " linear interpolation with sampling at the centroid of a"
                                      << " mesh object.  The number of samples per entity should"
                                      << " be 1.  Value found:" << m_nSampleElements);

  std::vector<stk::mesh::Entity>::const_iterator entityEnd = patch.end(), entityBeg = patch.begin(), entityIter;

  stk::search::CachedEntityFieldData entityFieldData;

  for(auto& var : m_recoverVars) {
    auto field = var->m_field;
    int nCompField = 0;
    for(entityIter = entityBeg; entityIter != entityEnd; ++entityIter) {
      var->populate_entity_data(*entityIter, entityFieldData);

      int scalarsPerEntity = entityFieldData.numComponents * entityFieldData.numCopies;
      if (nCompField == 0) {
        nCompField = scalarsPerEntity;
      }
      else {
        STK_ThrowRequireMsg(scalarsPerEntity == nCompField, "Error in EntityCentroidLinear, field "<<field->name()
                                   <<" has two different sizes ("<<scalarsPerEntity<<" and "<<nCompField
                                   <<") within the patch.");
      }
    }

    STK_ThrowRequireMsg(nCompField != 0, "Variable to interpolate has zero intrinsic_length."
                                          << " It is like a vector of length zero."
                                          << " Field name:" << field->name());

    STK_ThrowRequireMsg(nCompField == m_totalVarComponents,
                    "Number of field values found: " << nCompField
                                                     << " is incompatible with expected value: " << m_totalVarComponents);


    int entityCount = 0;
    for(entityIter = entityBeg; entityIter != entityEnd; ++entityIter, ++entityCount) {
      stk::mesh::Entity entity = *entityIter;
      var->populate_entity_data(entity, entityFieldData);
      const double* entityData = entityFieldData.constPointer;

      double data = 0.0;
      int entityFieldCopyComponent = 0;

      if(nullptr == entityData) {
        for(int copy = 0; copy < entityFieldData.numCopies; copy++) {
          for(int component = 0; component < entityFieldData.numComponents; component++) {
            fieldSample[entityFieldCopyComponent * nSampPatch + entityCount] = data;
            entityFieldCopyComponent++;
          }
        }
      } else {
        for(int copy = 0; copy < entityFieldData.numCopies; copy++) {
          for(int component = 0; component < entityFieldData.numComponents; component++) {
            data = m_transform(entityData[component*entityFieldData.componentStride + copy*entityFieldData.copyStride]);
            fieldSample[entityFieldCopyComponent * nSampPatch + entityCount] = data;
            entityFieldCopyComponent++;
          }
        }
      }
    }
  }

  STK_ThrowRequireMsg(m_recoveryType == RecoverField::TRILINEAR,
                  "The only recover type supported by the EntityCentroidLinear"
                      << " algorithm is stk::transfer::RecoverField::TRILINEAR. \n"
                      << "stk::transfer::RecoverField::TRILINEAR is enumeration number " << RecoverField::TRILINEAR << "\n"
                      << "The recovery type found was: " << m_recoveryType);

  STK_ThrowRequireMsg((int)m_recoveryType == 8, "The only recover type supported by the Element_Centroid_Linear"
                                              << " algorithm is stk::transfer::RecoverField::TRILINEAR.\n"
                                              << "stk::transfer::RecoverField::TRILINEAR should be enumeration number 8.\n"
                                              << "stk::transfer::RecoverField::TRILINEAR now appears to be enumeration number "
                                              << RecoverField::TRILINEAR
                                              << " The recovery type found was: " << m_recoveryType);

  const unsigned ndim = m_nodeVar->m_field->mesh_meta_data().spatial_dimension();

  int sampCount = 0;

  for(entityIter = entityBeg; entityIter != entityEnd; ++entityIter, ++sampCount) {
    stk::mesh::Entity entity = *entityIter;

    m_scratchSpace.clear();
    stk::search::determine_centroid(ndim, entity, m_nodeVar, m_scratchSpace);

    double x = m_scratchSpace[0];
    double y = 1 < ndim ? m_scratchSpace[1] : 0;
    double z = 2 < ndim ? m_scratchSpace[2] : 0;

    evaluate_trilinear_basis(x, y, z, basisSample + sampCount, nSampPatch);
  }

  STK_ThrowRequireMsg(sampCount == nSampPatch, "Internal programming error in EntityCentroidLinear algorithm."
                                                 << " Sample count does not match expected count.\n"
                                                 << "Expected count was: " << sampCount << ".\n"
                                                 << "Actual count found was: " << nSampPatch);
}

EntityCentroidQuadraticRecoverField::EntityCentroidQuadraticRecoverField(stk::transfer::RecoverField::RecoveryType recType,
                                                                     const std::vector<std::shared_ptr<stk::search::CachedFieldDataBase>>& recVars,
                                                                     const std::shared_ptr<stk::search::CachedFieldDataBase>& recNodeVar, int nsampElem,
                                                                     stk::mesh::Entity entity, FieldTransform transform)
  : RecoverField(recType, nsampElem)
  , m_recoverVars(recVars)
  , m_nodeVar(recNodeVar)
  , m_transform(transform)
{
  m_totalVarComponents = 0;

  for(auto const& v_i : recVars) {
    m_totalVarComponents += stk::mesh::field_scalars_per_entity(*v_i->m_field, entity);
  }
}

EntityCentroidQuadraticRecoverField::~EntityCentroidQuadraticRecoverField() {}

void EntityCentroidQuadraticRecoverField::sample_patch(const std::vector<stk::mesh::Entity>& patch, int nSampPatch,
                                                       double* fieldSample, double* basisSample) const
{
  // Function meant for piecewise linear variables with sampling at centroid of mesh entity
  STK_ThrowRequireMsg(m_nSampleElements == 1, "The EntityCentroidQuadratic method only works with piecewise"
                                      << " quadratic interpolation with sampling at the centroid of a"
                                      << " mesh object.  The number of samples per entity should"
                                      << " be 1.  Value found:" << m_nSampleElements);

  std::vector<stk::mesh::Entity>::const_iterator entityEnd = patch.end(), entityBeg = patch.begin(), entityIter;

  stk::search::CachedEntityFieldData entityFieldData;

  for(auto& var : m_recoverVars) {
    auto field = var->m_field;
    int nCompField = 0;
    for(entityIter = entityBeg; entityIter != entityEnd; ++entityIter) {
      var->populate_entity_data(*entityIter, entityFieldData);

      int scalarsPerEntity = entityFieldData.numComponents * entityFieldData.numCopies;
      if(nCompField == 0) {
        nCompField = scalarsPerEntity;
      }
      else {
        STK_ThrowRequireMsg(scalarsPerEntity == nCompField, "Error in EntityCentroidQuadratic, field "
                                                             << field->name() << " has two different sizes ("
                                                             << scalarsPerEntity << " and " << nCompField
                                                             << ") within the patch.");
      }
    }

    STK_ThrowRequireMsg(nCompField != 0, "Variable to interpolate has zero intrinsic_length."
                                          << " It is like a vector of length zero."
                                          << " Field name:" << field->name());

    STK_ThrowRequireMsg(nCompField == m_totalVarComponents,
                    "Number of field values found: " << nCompField
                                                     << " is incompatible with expected value: " << m_totalVarComponents);

    int entityCount = 0;
    for(entityIter = entityBeg; entityIter != entityEnd; ++entityIter, ++entityCount) {
      stk::mesh::Entity entity = *entityIter;
      var->populate_entity_data(entity, entityFieldData);
      const double* entityData = entityFieldData.constPointer;

      double data = 0.0;
      int entityFieldCopyComponent = 0;

      if(nullptr == entityData) {
        for(int copy = 0; copy < entityFieldData.numCopies; copy++) {
          for(int component = 0; component < entityFieldData.numComponents; component++) {
            fieldSample[entityFieldCopyComponent * nSampPatch + entityCount] = data;
            entityFieldCopyComponent++;
          }
        }
      } else {
        for(int copy = 0; copy < entityFieldData.numCopies; copy++) {
          for(int component = 0; component < entityFieldData.numComponents; component++) {
            data = m_transform(entityData[component*entityFieldData.componentStride + copy*entityFieldData.copyStride]);
            fieldSample[entityFieldCopyComponent * nSampPatch + entityCount] = data;
            entityFieldCopyComponent++;
          }
        }
      }
    }
  }

  STK_ThrowRequireMsg(m_recoveryType == RecoverField::TRIQUADRATIC,
                  "The only recover type supported by the EntityCentroidQuadratic"
                      << " algorithm is stk::transfer::RecoverField::TRIQUADRATIC. \n"
                      << "stk::transfer::RecoverField::TRIQUADRATIC is enumeration number " << RecoverField::TRIQUADRATIC
                      << "\n"
                      << "The recovery type found was: " << m_recoveryType);

  STK_ThrowRequireMsg((int)m_recoveryType == 27,
                  "The only recover type supported by the EntityCentroidQuadratic"
                      << " algorithm is stk::transfer::RecoverField::TRIQUADRATIC.\n"
                      << "stk::transfer::RecoverField::TRIQUADRATIC should be enumeration number 27.\n"
                      << "stk::transfer::RecoverField::TRIQUADRATIC now appears to be enumeration number "
                      << RecoverField::TRIQUADRATIC << " The recovery type found was: " << m_recoveryType);

  const int ndim = m_nodeVar->m_field->mesh_meta_data().spatial_dimension();

  int sampCount = 0;

  for(entityIter = entityBeg; entityIter != entityEnd; ++entityIter, ++sampCount) {
    stk::mesh::Entity entity = *entityIter;

    m_scratchSpace.clear();
    stk::search::determine_centroid(ndim, entity, m_nodeVar, m_scratchSpace);

    double x = m_scratchSpace[0];
    double y = 1 < ndim ? m_scratchSpace[1] : 0;
    double z = 2 < ndim ? m_scratchSpace[2] : 0;

    evaluate_triquadratic_basis(x, y, z, basisSample + sampCount, nSampPatch);
  }

  STK_ThrowRequireMsg(sampCount == nSampPatch, "Internal programming error in EntityCentroidQuadratic algorithm."
                                                 << " Sample count does not match expected count.\n"
                                                 << "Expected count was: " << sampCount << ".\n"
                                                 << "Actual count found was: " << nSampPatch);
}

EntityCentroidCubicRecoverField::EntityCentroidCubicRecoverField(stk::transfer::RecoverField::RecoveryType recType,
                                                                 const std::vector<std::shared_ptr<stk::search::CachedFieldDataBase>>& recVars,
                                                                 const std::shared_ptr<stk::search::CachedFieldDataBase>& recNodeVar, int nsampElem,
                                                                 stk::mesh::Entity entity, FieldTransform transform)
  : RecoverField(recType, nsampElem)
  , m_recoverVars(recVars)
  , m_nodeVar(recNodeVar)
  , m_transform(transform)
{
  m_totalVarComponents = 0;

  for(auto const& v_i : recVars) {
    m_totalVarComponents += stk::mesh::field_scalars_per_entity(*v_i->m_field, entity);
  }
}

EntityCentroidCubicRecoverField::~EntityCentroidCubicRecoverField() {}

void EntityCentroidCubicRecoverField::sample_patch(const std::vector<stk::mesh::Entity>& patch, int nSampPatch,
                                                   double* fieldSample, double* basisSample) const
{
  // Function meant for piecewise linear variables with sampling at centroid of mesh entity
  STK_ThrowRequireMsg(m_nSampleElements == 1, "The EntityCentroidCubic method only works with piecewise"
                                      << " cubic interpolation with sampling at the centroid of a"
                                      << " mesh object.  The number of samples per entity should"
                                      << " be 1.  Value found:" << m_nSampleElements);

  std::vector<stk::mesh::Entity>::const_iterator entityEnd = patch.end(), entityBeg = patch.begin(), entityIter;

  stk::search::CachedEntityFieldData entityFieldData;

  for(auto& var : m_recoverVars) {
    auto field = var->m_field;
    int nCompField = 0;
    for(entityIter = entityBeg; entityIter != entityEnd; ++entityIter) {
      var->populate_entity_data(*entityIter, entityFieldData);

      int scalarsPerEntity = entityFieldData.numComponents * entityFieldData.numCopies;
      if(nCompField == 0) {
        nCompField = scalarsPerEntity;
      }
      else {
        STK_ThrowRequireMsg(scalarsPerEntity == nCompField, "Error in EntityCentroidCubic, field "
                                                             << field->name() << " has two different sizes ("
                                                             << scalarsPerEntity << " and " << nCompField
                                                             << ") within the patch.");
      }
    }

    STK_ThrowRequireMsg(nCompField != 0, "Variable to interpolate has zero intrinsic_length."
                                          << " It is like a vector of length zero."
                                          << " Field name:" << var->name());

    STK_ThrowRequireMsg(nCompField == m_totalVarComponents,
                    "Number of field values found: " << nCompField
                                                     << " is incompatible with expected value: " << m_totalVarComponents);

    int entityCount = 0;
    for(entityIter = entityBeg; entityIter != entityEnd; ++entityIter, ++entityCount) {
      stk::mesh::Entity entity = *entityIter;
      var->populate_entity_data(entity, entityFieldData);
      const double* entityData = entityFieldData.constPointer;

      double data = 0.0;
      int entityFieldCopyComponent = 0;

      if(nullptr == entityData) {
        for(int copy = 0; copy < entityFieldData.numCopies; copy++) {
          for(int component = 0; component < entityFieldData.numComponents; component++) {
            fieldSample[entityFieldCopyComponent * nSampPatch + entityCount] = data;
            entityFieldCopyComponent++;
          }
        }
      } else {
        for(int copy = 0; copy < entityFieldData.numCopies; copy++) {
          for(int component = 0; component < entityFieldData.numComponents; component++) {
            data = m_transform(entityData[component*entityFieldData.componentStride + copy*entityFieldData.copyStride]);
            fieldSample[entityFieldCopyComponent * nSampPatch + entityCount] = data;
            entityFieldCopyComponent++;
          }
        }
      }
    }
  }

  STK_ThrowRequireMsg(m_recoveryType == RecoverField::TRICUBIC,
                  "The only recover type supported by the EntityCentroidCubic"
                      << " algorithm is stk::transfer::RecoverField::TRICUBIC. \n"
                      << "stk::transfer::RecoverField::TRICUBIC is enumeration number " << RecoverField::TRICUBIC
                      << "\n"
                      << "The recovery type found was: " << m_recoveryType);

  STK_ThrowRequireMsg((int)m_recoveryType == 64,
                  "The only recover type supported by the EntityCentroidCubic"
                      << " algorithm is stk::transfer::RecoverField::TRICUBIC.\n"
                      << "stk::transfer::RecoverField::TRICUBIC should be enumeration number 64.\n"
                      << "stk::transfer::RecoverField::TRICUBIC now appears to be enumeration number "
                      << RecoverField::TRICUBIC << " The recovery type found was: " << m_recoveryType);

  const int ndim = m_nodeVar->m_field->mesh_meta_data().spatial_dimension();

  int sampCount = 0;

  for(entityIter = entityBeg; entityIter != entityEnd; ++entityIter, ++sampCount) {
    stk::mesh::Entity entity = *entityIter;

    m_scratchSpace.clear();
    stk::search::determine_centroid(ndim, entity, m_nodeVar, m_scratchSpace);

    double x = m_scratchSpace[0];
    double y = 1 < ndim ? m_scratchSpace[1] : 0;
    double z = 2 < ndim ? m_scratchSpace[2] : 0;

    evaluate_tricubic_basis(x, y, z, basisSample + sampCount, nSampPatch);
  }

  STK_ThrowRequireMsg(sampCount == nSampPatch, "Internal programming error in EntityCentroidCubic algorithm."
                                                 << " Sample count does not match expected count.\n"
                                                 << "Expected count was: " << sampCount << ".\n"
                                                 << "Actual count found was: " << nSampPatch);
}

} // namespace transfer
} // namespace stk
