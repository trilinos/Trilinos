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


#ifndef STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_LEASTSQUARESINTERPOLATION_HPP_
#define STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_LEASTSQUARESINTERPOLATION_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <limits>                        // for numeric_limits
#include <string>                        // for string, operator==, basic_st...
#include <utility>                       // for pair
#include <vector>                        // for vector
#include <functional>

#include "stk_mesh/base/Types.hpp"       // for EntityRank
#include <stk_transfer_util/LeastSquares.hpp>
#include "stk_mesh/base/BulkData.hpp"             // for BulkData
#include "stk_mesh/base/Entity.hpp"               // for Entity
#include "stk_mesh/base/FieldBase.hpp"            // for FieldBase, field_data
#include "stk_mesh/base/FieldRestriction.hpp"     // for FieldRestriction
#include "stk_mesh/base/MetaData.hpp"             // for get_field_by_name
#include "stk_topology/topology.hpp"              // for topology
#include "stk_transfer_util/RecoverField.hpp"           // for RecoverField, Recov...
#include "stk_transfer_util/EntityCentroidRecoverField.hpp"
#include "stk_transfer_util/Patch.hpp"
#include "stk_util/util/ReportHandler.hpp"        // for ThrowRequireMsg

namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class Selector; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace transfer {

using FieldTransform = std::function<double(double)>;

struct EntityInterpolationData {
  const stk::mesh::BulkData& sendBulk;
  const stk::mesh::FieldBase* sendField;
  const unsigned sendDataIndex;
  const unsigned recvDataIndex;
  const stk::mesh::FieldBase* sendNodalCoordField;
  stk::mesh::Entity sendEntity;
  const stk::mesh::Part* sendPart;
  const stk::mesh::Selector* sendSelector;
  FieldTransform  preTransform;
  FieldTransform postTransform;

  EntityInterpolationData(const stk::mesh::BulkData& sendBulk_, const stk::mesh::FieldBase* sendField_,
                          const unsigned sendDataIndex_, const unsigned recvDataIndex_,
                          const stk::mesh::FieldBase* sendNodalCoordField_, stk::mesh::Entity sendEntity_,
                          const stk::mesh::Part* sendPart_, const stk::mesh::Selector* sendSelector_,
                          FieldTransform preTransform_ = [](double value) {return value;},
                          FieldTransform postTransform_ = [](double value) {return value;})
  : sendBulk(sendBulk_),
    sendField(sendField_),
    sendDataIndex(sendDataIndex_),
    recvDataIndex(recvDataIndex_),
    sendNodalCoordField(sendNodalCoordField_),
    sendEntity(sendEntity_),
    sendPart(sendPart_),
    sendSelector(sendSelector_),
    preTransform(preTransform_),
    postTransform(postTransform_)
  { }

private:
  EntityInterpolationData();
  EntityInterpolationData(const EntityInterpolationData&);
};

template <typename PATCHTYPE>
bool least_squares_linear_interpolation(stk::transfer::Patch<PATCHTYPE>& patch,
                                        stk::transfer::LeastSquares& leastSquaresCalculator,
                                        const std::vector<double>& evalPoint,
                                        const EntityInterpolationData& interpData,
                                        unsigned toFieldSize, double* toFieldPtr)
{
  const stk::mesh::FieldBase* sendField = interpData.sendField;
  const unsigned sendDataIndex = interpData.sendDataIndex;
  const unsigned recvDataIndex = interpData.recvDataIndex;
  const stk::mesh::FieldBase* sendCoordField = interpData.sendNodalCoordField;

  bool solvable = false;

  stk::mesh::Entity sendEntity = patch.get_patch_seed();
  const stk::mesh::BulkData& bulk = patch.get_bulk_data();
  const std::vector<stk::mesh::Entity>& patchEntities = patch.get_patch_entities();

  const std::vector<const stk::mesh::FieldBase*> sendFieldList(1, sendField);
  stk::transfer::EntityCentroidLinearRecoverField fieldInfo(stk::transfer::RecoverField::TRILINEAR, sendFieldList,
                                                            *sendCoordField, 1, sendEntity, interpData.preTransform);

  const auto basisSize = (unsigned)stk::transfer::RecoverField::TRILINEAR;
  const unsigned totalLength = stk::mesh::field_scalars_per_entity(*sendField, sendEntity);
  const int numCoeff = basisSize * totalLength;
  STK_ThrowRequireMsg(totalLength == (unsigned)fieldInfo.total_recov_var_comp(), " Can only handle one variable at a time.");

  std::vector<double> recoveredCoeff(2 * numCoeff, 0.0); // Include scratch space

  // Filter objects so that only those participating in the
  // transfer (if supplied) are included in the patch
  const double diagonal_tol = 1.0e-10;

  solvable = (bool)fieldInfo.recover_coeff(bulk, sendEntity, patchEntities, leastSquaresCalculator,
                                           recoveredCoeff.data(), diagonal_tol);

  // Interpolate new value for receive object

  double* result = toFieldPtr;
  if(!result) {
    STK_ThrowErrorMsg( "Data not found on object.  nullptr pointer found instead"
                        << " during interpolation phase of transfer.\n"
                        << " Data field to be accessed from'" << sendField->name() << "'"
                        << "\n");
  }
  else {
    int ndim = evalPoint.size();

    STK_ThrowRequireMsg(ndim == 2 || ndim == 3, " Transfer only handles two or three dimensions.");

    double x = evalPoint[0];
    double y = evalPoint[1];
    double z = (ndim == 2) ? 0.0 : evalPoint[2];
    const unsigned s = sendDataIndex;
    const unsigned r = recvDataIndex;
    if(r || s) {
      STK_ThrowRequireMsg(r <= toFieldSize, "Error: Index too large for receiving field.");
      STK_ThrowRequireMsg(s <= totalLength, "Error: Index too large for sending field.");
      const unsigned i = r ? r - 1 : 0;
      const unsigned j = s ? s - 1 : 0;
      if(solvable) {
        double* basisEval = &recoveredCoeff[numCoeff];
        double* coeff_comp = &recoveredCoeff[basisSize * j];
        stk::transfer::evaluate_trilinear_basis(x, y, z, basisEval, 1);
        result[i] = 0.0;
        for(unsigned k = 0; k < basisSize; k++) {
          result[i] += coeff_comp[k] * basisEval[k];
        }
        interpData.postTransform(result[i]);
      }
      else {
        result[i] = 0;
      }
    }
    else {
      if(solvable) {
        double* basisEval = &recoveredCoeff[numCoeff];
        stk::transfer::evaluate_trilinear_basis(x, y, z, basisEval, 1);
        for(unsigned i_comp = 0; i_comp < totalLength; ++i_comp) {
          double* coeff_comp = &recoveredCoeff[basisSize * i_comp];
          result[i_comp] = 0.0;
          for(unsigned k = 0; k < basisSize; k++) {
            result[i_comp] += coeff_comp[k] * basisEval[k];
          }
          interpData.postTransform(result[i_comp]);
        }
      }
      else {
        for(unsigned i_comp = 0; i_comp < totalLength; ++i_comp) {
          result[i_comp] = 0;
        }
      }
    }
  }

  return solvable;
}



template <typename PATCHTYPE>
bool least_squares_quadratic_interpolation(stk::transfer::Patch<PATCHTYPE>& patch,
                                           stk::transfer::LeastSquares& leastSquaresCalculator,
                                           const std::vector<double>& evalPoint,
                                           const EntityInterpolationData& interpData,
                                           unsigned toFieldSize, double* toFieldPtr)
{
  const stk::mesh::FieldBase* sendField = interpData.sendField;
  const unsigned sendDataIndex = interpData.sendDataIndex;
  const unsigned recvDataIndex = interpData.recvDataIndex;
  const stk::mesh::FieldBase* sendCoordField = interpData.sendNodalCoordField;

  bool solvable = false;

  stk::mesh::Entity sendEntity = patch.get_patch_seed();
  const stk::mesh::BulkData& bulk = patch.get_bulk_data();
  const std::vector<stk::mesh::Entity>& patchEntities = patch.get_patch_entities();

  const std::vector<const stk::mesh::FieldBase*> sendFieldList(1, sendField);
  stk::transfer::EntityCentroidQuadraticRecoverField fieldInfo(stk::transfer::RecoverField::TRIQUADRATIC, sendFieldList,
                                                               *sendCoordField, 1, sendEntity, interpData.preTransform);

  const auto basisSize = (unsigned)stk::transfer::RecoverField::TRIQUADRATIC;
  const unsigned totalLength = stk::mesh::field_scalars_per_entity(*sendField, sendEntity);
  const int numCoeff = basisSize * totalLength;
  STK_ThrowRequireMsg(totalLength == (unsigned)fieldInfo.total_recov_var_comp(), " Can only handle one variable at a time.");

  std::vector<double> recoveredCoeff(2 * numCoeff, 0.0); // Include scratch space

  // Filter objects so that only those participating in the
  // transfer (if supplied) are included in the patch
  const double diagonal_tol = 1.0e-10;

  solvable = (bool)fieldInfo.recover_coeff(bulk, sendEntity, patchEntities, leastSquaresCalculator,
                                           recoveredCoeff.data(), diagonal_tol);

  // Interpolate new value for receive object

  double* result = toFieldPtr;
  if(!result) {
    STK_ThrowErrorMsg( "Data not found on object.  nullptr pointer found instead"
                        << " during interpolation phase of transfer.\n"
                        << " Data field to be accessed from'" << sendField->name() << "'"
                        << "\n");
  }
  else {
    int ndim = evalPoint.size();

    STK_ThrowRequireMsg(ndim == 2 || ndim == 3, " Transfer only handles two or three dimensions.");

    double x = evalPoint[0];
    double y = evalPoint[1];
    double z = (ndim == 2) ? 0.0 : evalPoint[2];
    const unsigned s = sendDataIndex;
    const unsigned r = recvDataIndex;
    if(r || s) {
      STK_ThrowRequireMsg(r <= toFieldSize, "Error: Index too large for receiving field.");
      STK_ThrowRequireMsg(s <= totalLength, "Error: Index too large for sending field.");
      const unsigned i = r ? r - 1 : 0;
      const unsigned j = s ? s - 1 : 0;
      if(solvable) {
        double* coeff_comp = &recoveredCoeff[basisSize * j];
        double* basisEval = &recoveredCoeff[numCoeff];
        stk::transfer::evaluate_triquadratic_basis(x, y, z, basisEval, 1);
        result[i] = 0.0;
        for(unsigned k = 0; k < basisSize; k++) {
          result[i] += coeff_comp[k] * basisEval[k];
        }
        interpData.postTransform(result[i]);
      }
      else {
        result[i] = 0;
      }
    }
    else {
      if(solvable) {
        double* basisEval = &recoveredCoeff[numCoeff];
        stk::transfer::evaluate_triquadratic_basis(x, y, z, basisEval, 1);
        for(unsigned i_comp = 0; i_comp < totalLength; ++i_comp) {
          double* coeff_comp = &recoveredCoeff[basisSize * i_comp];
          result[i_comp] = 0.0;
          for(unsigned k = 0; k < basisSize; k++) {
            result[i_comp] += coeff_comp[k] * basisEval[k];
          }
          interpData.postTransform(result[i_comp]);
        }
      }
      else {
        for(unsigned i_comp = 0; i_comp < totalLength; ++i_comp) {
          result[i_comp] = 0;
        }
      }
    }
  }

  return solvable;
}

template <typename PATCHTYPE>
bool least_squares_cubic_interpolation(stk::transfer::Patch<PATCHTYPE>& patch,
                                       stk::transfer::LeastSquares& leastSquaresCalculator,
                                       const std::vector<double>& evalPoint,
                                       const EntityInterpolationData& interpData,
                                       unsigned toFieldSize, double* toFieldPtr)
{
  const stk::mesh::FieldBase* sendField = interpData.sendField;
  const unsigned sendDataIndex = interpData.sendDataIndex;
  const unsigned recvDataIndex = interpData.recvDataIndex;
  const stk::mesh::FieldBase* sendCoordField = interpData.sendNodalCoordField;

  bool solvable = false;

  stk::mesh::Entity sendEntity = patch.get_patch_seed();
  const stk::mesh::BulkData& bulk = patch.get_bulk_data();
  const std::vector<stk::mesh::Entity>& patchEntities = patch.get_patch_entities();

  const std::vector<const stk::mesh::FieldBase*> sendFieldList(1, sendField);
  stk::transfer::EntityCentroidCubicRecoverField fieldInfo(stk::transfer::RecoverField::TRICUBIC, sendFieldList,
                                                           *sendCoordField, 1, sendEntity, interpData.preTransform);

  const auto basisSize = (unsigned)stk::transfer::RecoverField::TRICUBIC;
  const unsigned totalLength = stk::mesh::field_scalars_per_entity(*sendField, sendEntity);
  const int numCoeff = basisSize * totalLength;
  STK_ThrowRequireMsg(totalLength == (unsigned)fieldInfo.total_recov_var_comp(), " Can only handle one variable at a time.");

  std::vector<double> recoveredCoeff(2 * numCoeff, 0.0); // Include scratch space

  // Filter objects so that only those participating in the
  // transfer (if supplied) are included in the patch
  const double diagonal_tol = 1.0e-10;

  solvable = (bool)fieldInfo.recover_coeff(bulk, sendEntity, patchEntities, leastSquaresCalculator,
                                           recoveredCoeff.data(), diagonal_tol);

  // Interpolate new value for receive object

  double* result = toFieldPtr;
  if(!result) {
    STK_ThrowErrorMsg( "Data not found on object.  nullptr pointer found instead"
                        << " during interpolation phase of transfer.\n"
                        << " Data field to be accessed from'" << sendField->name() << "'"
                        << "\n");
  }
  else {
    int ndim = evalPoint.size();

    STK_ThrowRequireMsg(ndim == 2 || ndim == 3, " Transfer only handles two or three dimensions.");

    double x = evalPoint[0];
    double y = evalPoint[1];
    double z = (ndim == 2) ? 0.0 : evalPoint[2];
    const unsigned s = sendDataIndex;
    const unsigned r = recvDataIndex;
    if(r || s) {
      STK_ThrowRequireMsg(r <= toFieldSize, "Error: Index too large for receiving field.");
      STK_ThrowRequireMsg(s <= totalLength, "Error: Index too large for sending field.");
      const unsigned i = r ? r - 1 : 0;
      const unsigned j = s ? s - 1 : 0;
      if(solvable) {
        double* coeff_comp = &recoveredCoeff[basisSize * j];
        double* basisEval = &recoveredCoeff[numCoeff];
        stk::transfer::evaluate_tricubic_basis(x, y, z, basisEval, 1);
        result[i] = 0.0;
        for(unsigned k = 0; k < basisSize; k++) {
          result[i] += coeff_comp[k] * basisEval[k];
        }
        interpData.postTransform(result[i]);
      }
      else {
        result[i] = 0;
      }
    }
    else {
      if(solvable) {
        double* basisEval = &recoveredCoeff[numCoeff];
        stk::transfer::evaluate_tricubic_basis(x, y, z, basisEval, 1);
        for(unsigned i_comp = 0; i_comp < totalLength; ++i_comp) {
          double* coeff_comp = &recoveredCoeff[basisSize * i_comp];
          result[i_comp] = 0.0;
          for(unsigned k = 0; k < basisSize; k++) {
            result[i_comp] += coeff_comp[k] * basisEval[k];
          }
          interpData.postTransform(result[i_comp]);
        }
      }
      else {
        for(unsigned i_comp = 0; i_comp < totalLength; ++i_comp) {
          result[i_comp] = 0;
        }
      }
    }
  }

  return solvable;
}

} // namespace transfer
} // namespace stk


#endif /* STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_LEASTSQUARESINTERPOLATION_HPP_ */
