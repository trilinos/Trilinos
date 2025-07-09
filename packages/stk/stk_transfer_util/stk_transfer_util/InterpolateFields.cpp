// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_mesh/base/Bucket.hpp"                         // for Bucket
#include <stk_mesh/base/BulkData.hpp>                       // for BulkData
#include <stk_mesh/base/Entity.hpp>                         // for Entity
#include "stk_mesh/base/FieldBase.hpp"                      // for FieldBase
#include "stk_transfer/TransferInterface.hpp"                   // for ProvideMa...
#include "stk_transfer_util/InterpolateFields.hpp"
#include "stk_transfer_util/LeastSquaresInterpolation.hpp"  // for EntityInt...
#include "stk_transfer_util/PointInterpolation.hpp"
#include "stk_topology/topology.hpp"                        // for operator<<
#include "stk_util/util/ReportHandler.hpp"

#include <algorithm>                                        // for copy, max
#include <memory>                                           // for allocator...
#include <ostream>                                          // for operator<<
#include <string>                                           // for operator<<
namespace stk::mesh { class Part; }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace transfer {

InterpolateFieldsInterface::InterpolateFieldsInterface(stk::mesh::BulkData& bulk)
  : m_bulk(bulk)
  , m_meta(bulk.mesh_meta_data())
  , m_fieldsAreSet(false)
{
}

void InterpolateFieldsInterface::set_fields(const FieldSpecVector& fieldSpecs)
{
  m_fieldVec = stk::transfer::get_fields(m_meta, fieldSpecs);
  m_lowerBound = stk::transfer::get_lower_bounds(fieldSpecs);
  m_upperBound = stk::transfer::get_upper_bounds(fieldSpecs);
  m_preTransform = stk::transfer::get_pre_transforms(fieldSpecs);
  m_postTransform = stk::transfer::get_post_transforms(fieldSpecs);
  m_defaultFieldValue = stk::transfer::get_default_field_values(fieldSpecs);

  m_fieldsAreSet = true;
}

void InterpolateFieldsInterface::apply_bounds(const unsigned index, const unsigned length, double* fieldData) const
{
  stk::transfer::apply_bounds(length, fieldData, m_lowerBound[index], m_upperBound[index]);
}

MasterElementFieldInterpolator::MasterElementFieldInterpolator(
    stk::mesh::BulkData& bulk, std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider)
  : InterpolateFieldsInterface(bulk)
  , m_masterElemProvider(masterElemProvider)
{
}

void MasterElementFieldInterpolator::interpolate_fields(const stk::search::spmd::EntityKeyPair& key,
                                                        const std::vector<double>& /*evalPoint*/,
                                                        const std::vector<double>& parametricCoords,
                                                        InterpolationData& data) const
{
  const stk::mesh::Bucket& bucket = m_bulk.bucket(key);
  const stk::topology& topo = bucket.topology();

  const auto num_nodes = m_bulk.num_nodes(key);

  unsigned numFieldComponents;
  unsigned numNodes;

  for(unsigned k = 0; k < data.nFields; ++k) {
    unsigned n = data.fieldKey[k];

    FieldTransform preTransform = m_preTransform[n];
    FieldTransform postTransform = m_postTransform[n];

    double* toField = data.fieldPtr[n];
    unsigned sizeOfRecvField = data.fieldSize[n];
    unsigned sizeOfSendField = stk::transfer::get_size_of_field(m_fieldVec[n].field);

    if(nullptr == toField) {
      continue;
    }

    double defaultValue = m_defaultFieldValue[n];

    unsigned sendIndex = m_fieldVec[n].index;
    unsigned recvIndex = data.fieldDataIndex[n];

    m_fieldDataBuffer.resize(num_nodes * sizeOfSendField);
    m_resultBuffer.resize(sizeOfSendField);

    auto meTopo = stk::search::SearchTopology(topo, key, &bucket);
    auto meField = stk::search::SearchField(m_fieldVec[n].field, preTransform, defaultValue);

    m_masterElemProvider->nodal_field_data(key, meField, numFieldComponents, numNodes, m_fieldDataBuffer);
    m_masterElemProvider->evaluate_field(meTopo, parametricCoords, sizeOfSendField, m_fieldDataBuffer, m_resultBuffer);

    const unsigned r = recvIndex;
    const unsigned s = sendIndex;
    if(r || s) {
      STK_ThrowAssertMsg(r <= sizeOfRecvField, "Error: Index too large for receiving field: r = "
                                                << r << " sizeOfRecvField = " << sizeOfRecvField
                                                << " field name = " << m_fieldVec[n].field->name());
      STK_ThrowAssertMsg(s <= sizeOfSendField, "Error: Index too large for sending field: s = "
                                                << s << " sizeOfRecvField = " << sizeOfSendField
                                                << " field name = " << m_fieldVec[n].field->name());
      std::fill(toField, toField + sizeOfRecvField, defaultValue);
      const unsigned i = r ? r - 1 : 0;
      const unsigned j = s ? s - 1 : 0;
      toField[i] = postTransform(m_resultBuffer[j]);
      apply_bounds(n, 1, &toField[i]);
    }
    else {
      STK_ThrowAssertMsg(sizeOfRecvField <= sizeOfSendField, "Error: recv field too long on entity "
                                                              << key << ". Length(send,recv) = (" << sizeOfSendField
                                                              << "," << sizeOfRecvField
                                                              << ") send field = " << m_fieldVec[n].field->name());
      for(unsigned i(0); i < sizeOfRecvField; ++i) {
        toField[i] = postTransform(m_resultBuffer[i]);
      }
      apply_bounds(n, sizeOfRecvField, toField);
    }
  }
}

PatchRecoveryFieldInterpolator::PatchRecoveryFieldInterpolator(stk::mesh::BulkData& bulk,
                                                               const stk::mesh::FieldBase* coordinates,
                                                               const stk::mesh::Part* meshPart,
                                                               const stk::mesh::Selector* meshSelector,
                                                               const stk::transfer::PatchRecoveryEvaluationType type)
  : InterpolateFieldsInterface(bulk)
  , m_coordinates(coordinates)
  , m_meshPart(meshPart)
  , m_meshSelector(nullptr != meshSelector ? new stk::mesh::Selector(*meshSelector) : nullptr)
  , m_activeSelector(stk::mesh::Selector().complement())
  , m_patchFilter(m_meshPart, m_meshSelector)
{
  set_interpolator(type);
}

PatchRecoveryFieldInterpolator::PatchRecoveryFieldInterpolator(stk::mesh::BulkData& bulk,
                                                               const stk::mesh::FieldBase* coordinates,
                                                               const stk::mesh::Selector& activeSelector,
                                                               const stk::mesh::Part* meshPart,
                                                               const stk::mesh::Selector* meshSelector,
                                                               const stk::transfer::PatchRecoveryEvaluationType type)
  : InterpolateFieldsInterface(bulk)
  , m_coordinates(coordinates)
  , m_meshPart(meshPart)
  , m_meshSelector(nullptr != meshSelector ? new stk::mesh::Selector(*meshSelector) : nullptr)
  , m_activeSelector(activeSelector)
  , m_patchFilter(m_meshPart, m_meshSelector)
{
  set_interpolator(type);
}

void PatchRecoveryFieldInterpolator::set_interpolator(const stk::transfer::PatchRecoveryEvaluationType type)
{
  if(type == stk::transfer::LINEAR_LEAST_SQUARES) {
    m_leastSquaresInterpolator = [this] (const EntityInterpolationData& interpData,
                                         const std::vector<double>& evalPoint, unsigned fieldSize, double* fieldPtr) {
                                           m_patchFilter.initialize(interpData.sendPart, interpData.sendSelector);
                                           return element_centroid_linear<EntityPatchFilter>(
                                             m_patchFilter, m_activeSelector, interpData, evalPoint, fieldSize, fieldPtr);
                                        };
  } else if(type == stk::transfer::LINEAR_MOVING_LEAST_SQUARES) {
    m_leastSquaresInterpolator = [this] (const EntityInterpolationData& interpData,
                                         const std::vector<double>& evalPoint, unsigned fieldSize, double* fieldPtr) {
                                           m_patchFilter.initialize(interpData.sendPart, interpData.sendSelector);
                                           return element_centroid_linear_mls<EntityPatchFilter>(
                                             m_patchFilter, m_activeSelector, interpData, evalPoint, fieldSize, fieldPtr);
                                        };
  } else if(type == stk::transfer::QUADRATIC_LEAST_SQUARES) {
    m_leastSquaresInterpolator = [this] (const EntityInterpolationData& interpData,
                                         const std::vector<double>& evalPoint, unsigned fieldSize, double* fieldPtr) {
                                           m_patchFilter.initialize(interpData.sendPart, interpData.sendSelector);
                                           return element_centroid_quadratic<stk::transfer::EntityPatchFilter>(
                                             m_patchFilter, m_activeSelector, interpData, evalPoint, fieldSize, fieldPtr);
                                        };
  } else if(type == stk::transfer::QUADRATIC_MOVING_LEAST_SQUARES) {
    m_leastSquaresInterpolator = [this] (const EntityInterpolationData& interpData,
                                         const std::vector<double>& evalPoint, unsigned fieldSize, double* fieldPtr) {
                                           m_patchFilter.initialize(interpData.sendPart, interpData.sendSelector);
                                           return element_centroid_quadratic_mls<stk::transfer::EntityPatchFilter>(
                                             m_patchFilter, m_activeSelector, interpData, evalPoint, fieldSize, fieldPtr);
                                        };
  } else if(type == stk::transfer::CUBIC_LEAST_SQUARES) {
    m_leastSquaresInterpolator = [this] (const EntityInterpolationData& interpData,
                                         const std::vector<double>& evalPoint, unsigned fieldSize, double* fieldPtr) {
                                           m_patchFilter.initialize(interpData.sendPart, interpData.sendSelector);
                                           return element_centroid_cubic<stk::transfer::EntityPatchFilter>(
                                             m_patchFilter, m_activeSelector, interpData, evalPoint, fieldSize, fieldPtr);
                                        };
  } else if(type == stk::transfer::CUBIC_MOVING_LEAST_SQUARES) {
    m_leastSquaresInterpolator = [this] (const EntityInterpolationData& interpData,
                                         const std::vector<double>& evalPoint, unsigned fieldSize, double* fieldPtr) {
                                           m_patchFilter.initialize(interpData.sendPart, interpData.sendSelector);
                                           return element_centroid_cubic_mls<stk::transfer::EntityPatchFilter>(
                                             m_patchFilter, m_activeSelector, interpData, evalPoint, fieldSize, fieldPtr);
                                        };
  } else {
    STK_ThrowRequireMsg(false, "Invalid enum: " << type << " for PatchRecoveryFieldInterpolator interpolation type");
  }
}

void PatchRecoveryFieldInterpolator::interpolate_fields(const stk::search::spmd::EntityKeyPair& key,
                                                        const std::vector<double>& evalPoint,
                                                        const std::vector<double>& /*parametricCoords*/,
                                                        InterpolationData& data) const
{
  for(unsigned i = 0; i < data.nFields; ++i) {
    unsigned fieldIndex = data.fieldKey[i];

    FieldTransform preTransform = m_preTransform[fieldIndex];
    FieldTransform postTransform = m_postTransform[fieldIndex];

    const unsigned sendDataIndex = m_fieldVec[fieldIndex].index;
    const unsigned recvDataIndex = data.fieldDataIndex[fieldIndex];

    stk::transfer::EntityInterpolationData interpData(m_bulk, m_fieldVec[fieldIndex].field,
                                                      sendDataIndex, recvDataIndex,
                                                      m_coordinates, key, m_meshPart, m_meshSelector,
                                                      preTransform, postTransform, m_defaultFieldValue[fieldIndex]);

    bool solvable = m_leastSquaresInterpolator(interpData, evalPoint, data.fieldSize[fieldIndex], data.fieldPtr[fieldIndex]);

    if(false == solvable) {
      std::stringstream oss;

      oss << "Could not construct least squares patch for field " << m_fieldVec[i].field->name()
          << " at coordinate: " << evalPoint[0] << "," << evalPoint[1]
          << ((evalPoint.size() > 2) ? "," : "");

      if(evalPoint.size() > 2) oss << evalPoint[2];

      STK_ThrowRequireMsg(solvable, oss.str());
    }
    else {
      const int r = recvDataIndex;
      if(r) {
        STK_ThrowAssertMsg(r <= data.fieldSize[fieldIndex], "Error: Index too large for receiving field.");
        const int ri = r - 1;
        apply_bounds(fieldIndex, 1, &data.fieldPtr[fieldIndex][ri]);
      }
      else {
        apply_bounds(fieldIndex, data.fieldSize[fieldIndex], data.fieldPtr[fieldIndex]);
      }
    }
  }
}

CopyFieldInterpolator::CopyFieldInterpolator(stk::mesh::BulkData& bulk)
  : InterpolateFieldsInterface(bulk)
{
}

void CopyFieldInterpolator::interpolate_fields(const stk::search::spmd::EntityKeyPair& key,
                                               const std::vector<double>& /*evalPoint*/,
                                               const std::vector<double>& /*parametricCoords*/,
                                               InterpolationData& data) const
{
  for(unsigned i = 0; i < data.nFields; ++i) {

    unsigned fieldIndex = data.fieldKey[i];

    FieldTransform preTransform = m_preTransform[fieldIndex];
    FieldTransform postTransform = m_postTransform[fieldIndex];

    unsigned sendIndex = m_fieldVec[fieldIndex].index;
    unsigned recvIndex = data.fieldDataIndex[fieldIndex];

    const unsigned sendLength = stk::mesh::field_scalars_per_entity(*m_fieldVec[fieldIndex].field, key);
    [[maybe_unused]] const unsigned recvLength = data.fieldSize[fieldIndex];

    const double* sendField = static_cast<const double *>(stk::mesh::field_data(*m_fieldVec[fieldIndex].field, key));
    double* recvField = data.fieldPtr[fieldIndex];

    if (nullptr == recvField){
        continue;
    }

    double defaultValue = m_defaultFieldValue[fieldIndex];

    if(sendLength) {
      STK_ThrowAssertMsg(sendField, "The field does not exist, nullptr pointer found.\n"
                                      << "Field name: " << m_fieldVec[fieldIndex].field->name() << "\n"
                                      << "Entity: " << key);

      STK_ThrowAssertMsg(recvField, "Recv field does not exist, nullptr pointer found.\n"
                                      << "Send field Name: " << m_fieldVec[fieldIndex].field->name() << "\n"
                                      << "Send entity: " << key);

      STK_ThrowAssertMsg(m_bulk.bucket(key).field_data_is_allocated(*m_fieldVec[fieldIndex].field),
                          "The field does not exist on send entity. Field name: "
                          << m_fieldVec[fieldIndex].field->name() << "\n"
                          << "Entity: " << key);

      if(stk::topology::NODE_RANK == m_fieldVec[fieldIndex].field->entity_rank() &&
         stk::topology::NODE_RANK != m_bulk.entity_rank(key)) {
        std::ostringstream oss;
        oss << "The INTERPOLATION FUNCTION COPY can not be used for NODE transfers.\n"
            << "Attempt to use the COPY interpolation to interpolate field name: "
            << m_fieldVec[fieldIndex].field->name() << "\n"
            << "(which is of type node) "
            << "\n"
            << "Remove INTERPOLATION FUNCTION COPY and use standard interpolation instead.\n\n";

        throw std::runtime_error(oss.str());
      }

      const unsigned r = recvIndex;
      const unsigned s = sendIndex;
      if(r || s) {
        STK_ThrowAssertMsg(r <= recvLength, "Error: Index too large for receiving field.");
        STK_ThrowAssertMsg(s <= sendLength, "Error: Index too large for sending field.");
        std::fill(recvField, recvField + recvLength, defaultValue);
        const unsigned ri = r ? r - 1 : 0;
        const unsigned si = s ? s - 1 : 0;
        recvField[ri] = postTransform( preTransform(sendField[si]) );
        apply_bounds(fieldIndex, 1, &recvField[ri]);
      }
      else {
        STK_ThrowAssertMsg(sendLength <= recvLength, "Error: Send field too long on entity "
                                                      << key << ". Length(send,recv) = (" << sendLength << ","
                                                      << recvLength
                                                      << ") send field = " << m_fieldVec[fieldIndex].field->name());
        std::copy(sendField, sendField + sendLength, recvField);

        for(unsigned j(0); j < sendLength; ++j) {
          recvField[j] = postTransform( preTransform(recvField[j]) );
        }

        apply_bounds(fieldIndex, sendLength, recvField);
      }
    }
  }
}

SumFieldInterpolator::SumFieldInterpolator(stk::mesh::BulkData& bulk)
  : InterpolateFieldsInterface(bulk)
{
}

void SumFieldInterpolator::interpolate_fields(const stk::search::spmd::EntityKeyPair& key,
                                              const std::vector<double>& /*evalPoint*/,
                                              const std::vector<double>& /*parametricCoords*/,
                                              InterpolationData& data) const
{
  for(unsigned i = 0; i < data.nFields; ++i) {
    unsigned fieldIndex = data.fieldKey[i];

    FieldTransform preTransform = m_preTransform[fieldIndex];
    FieldTransform postTransform = m_postTransform[fieldIndex];

    unsigned sendIndex = m_fieldVec[fieldIndex].index;
    unsigned recvIndex = data.fieldDataIndex[fieldIndex];

    const unsigned sendLength = stk::mesh::field_scalars_per_entity(*m_fieldVec[fieldIndex].field, key);
    const unsigned recvLength = data.fieldSize[fieldIndex];

    const double* sendField = static_cast<const double *>(stk::mesh::field_data(*m_fieldVec[fieldIndex].field, key));
    double* recvField = data.fieldPtr[fieldIndex];

    double defaultValue = m_defaultFieldValue[fieldIndex];

    if(sendLength) {
      STK_ThrowAssertMsg(sendField, "The field does not exist, nullptr pointer found.\n"
                                      << "Field name: " << m_fieldVec[fieldIndex].field->name() << "\n"
                                      << "Entity: " << key);

      STK_ThrowAssertMsg(recvField, "Recv field does not exist, nullptr pointer found.\n"
                                      << "Send field name: " << m_fieldVec[fieldIndex].field->name() << "\n"
                                      << "Send entity: " << key);

      STK_ThrowAssertMsg(m_bulk.bucket(key).field_data_is_allocated(*m_fieldVec[fieldIndex].field),
                          "The field does not exist on send entity. Field name: "
                          << m_fieldVec[fieldIndex].field->name() << "\n"
                          << "Entity: " << key);

      if(stk::topology::NODE_RANK == m_fieldVec[fieldIndex].field->entity_rank() &&
         stk::topology::NODE_RANK != m_bulk.entity_rank(key)) {
        std::ostringstream oss;
        oss << "The INTERPOLATION FUNCTION COPY can not be used for NODE transfers.\n"
            << "Attempt to use the COPY interpolation to interpolate field Name:"
            << m_fieldVec[fieldIndex].field->name() << "\n"
            << "(which is of type node) "
            << "\n"
            << "Remove INTERPOLATION FUNCTION COPY and use standard interpolation instead.\n\n";
        throw std::runtime_error(oss.str());
      }

      const unsigned r = recvIndex;
      const unsigned s = sendIndex;
      if(r || s) {
        STK_ThrowAssertMsg(r <= recvLength, "Error: Index too large for receiving field.");
        STK_ThrowAssertMsg(s <= sendLength, "Error: Index too large for sending field.");
        std::fill(recvField, recvField + recvLength, defaultValue);
        const unsigned ri = r ? r - 1 : 0;
        const unsigned si = s ? s - 1 : 0;
        recvField[ri] += postTransform( preTransform(sendField[si]) );
        apply_bounds(fieldIndex, 1, &recvField[ri]);
      }
      else {
        STK_ThrowAssertMsg(recvLength <= sendLength, "Error: Recv field too long on entity "
                                                      << key << ". Length(send,recv) = (" << sendLength << ","
                                                      << recvLength
                                                      << ") send field = " << m_fieldVec[fieldIndex].field->name());
        for(unsigned idx(0); idx < recvLength; ++idx) {
          recvField[idx] += postTransform( preTransform(sendField[idx]) );
        }
        apply_bounds(fieldIndex, recvLength, recvField);
      }
    }
  }
}

} // namespace transfer
} // namespace stk



