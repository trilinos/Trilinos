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

#ifndef STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_INTERPOLATEFIELDS_HPP_
#define STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_INTERPOLATEFIELDS_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_mesh/base/EntityKey.hpp"                // for EntityKey
#include "stk_mesh/base/Selector.hpp"                 // for Selector
#include "stk_search_util/CachedEntity.hpp"
#include "stk_search/SearchInterface.hpp"             // for ExternalPointHa...
#include "stk_search_util/MasterElementProvider.hpp"  // for ProvideMaste...
#include "stk_search_util/spmd/EntityKeyPair.hpp"
#include "stk_transfer/TransferInterface.hpp"
#include "stk_transfer/TransferTypes.hpp"
#include "stk_transfer_util/LeastSquares.hpp"         // for LINEAR_LEAST_SQ...
#include "stk_transfer_util/PointInterpolation.hpp"

#include <memory>                                     // for shared_ptr
#include <vector>                                     // for vector
#include <functional>                                 // for function
namespace stk::mesh { class BulkData; }
namespace stk::mesh { class MetaData; }
namespace stk::mesh { class FieldBase; }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace transfer {

using InterpolateFieldsInterfaceBase = stk::transfer::FieldInterpolatorInterface<stk::search::spmd::EntityKeyPair>;

class InterpolateFieldsInterface : public InterpolateFieldsInterfaceBase {
 public:
  InterpolateFieldsInterface(stk::mesh::BulkData& bulk);

  void set_fields(const FieldSpecVector& fieldSpecs);

  const std::vector<IndexedField>& get_fields() const { return m_fieldVec; }

  const std::vector<double>& get_lower_bounds() const { return m_lowerBound; }
  const std::vector<double>& get_upper_bounds() const { return m_upperBound; }

  const std::vector<FieldTransform>& get_pre_transforms() const { return m_preTransform; }
  const std::vector<FieldTransform>& get_post_transforms() const { return m_postTransform; }

  const std::vector<double>& get_default_field_values() const { return m_defaultFieldValue; }

  bool fields_are_set() const { return m_fieldsAreSet; }

  virtual void apply_bounds(const unsigned index, const unsigned length, double* fieldData) const;

  virtual ~InterpolateFieldsInterface() {}

 protected:
  stk::mesh::BulkData& m_bulk;
  const stk::mesh::MetaData& m_meta;
  std::vector<IndexedField> m_fieldVec;
  std::vector<double> m_lowerBound;
  std::vector<double> m_upperBound;
  std::vector<FieldTransform> m_preTransform;
  std::vector<FieldTransform> m_postTransform;
  std::vector<double> m_defaultFieldValue;

  bool m_fieldsAreSet{false};

  InterpolateFieldsInterface(const InterpolateFieldsInterface&) = delete;
  const InterpolateFieldsInterface& operator()(const InterpolateFieldsInterface&) = delete;
};

class MasterElementFieldInterpolator : public InterpolateFieldsInterface {
 public:
  MasterElementFieldInterpolator(stk::mesh::BulkData& bulk,
                                 std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider);

  void interpolate_fields(const stk::search::spmd::EntityKeyPair& k, const std::vector<double>& evalPoint,
                          const std::vector<double>& parametricCoords, InterpolationData& data) const override;

  ~MasterElementFieldInterpolator() = default;

 private:
  mutable std::vector<double> m_resultBuffer;
  mutable std::vector<double> m_fieldDataBuffer;
  std::shared_ptr<stk::search::MasterElementProviderInterface> m_masterElemProvider;

  MasterElementFieldInterpolator(const MasterElementFieldInterpolator&) = delete;
  const MasterElementFieldInterpolator& operator()(const MasterElementFieldInterpolator&) = delete;
};

class PatchRecoveryFieldInterpolator : public InterpolateFieldsInterface {
 public:

  PatchRecoveryFieldInterpolator(stk::mesh::BulkData& bulk, const stk::mesh::FieldBase* coordinates,
                                 const stk::mesh::Part* meshPart = nullptr,
                                 const stk::mesh::Selector* meshSelector = nullptr,
                                 const stk::transfer::PatchRecoveryEvaluationType type = stk::transfer::LINEAR_LEAST_SQUARES);

  PatchRecoveryFieldInterpolator(stk::mesh::BulkData& bulk, const stk::mesh::FieldBase* coordinates,
                                 const stk::mesh::Selector& activeSelector,
                                 const stk::mesh::Part* meshPart = nullptr,
                                 const stk::mesh::Selector* meshSelector = nullptr,
                                 const stk::transfer::PatchRecoveryEvaluationType type = stk::transfer::LINEAR_LEAST_SQUARES);

  void interpolate_fields(const stk::search::spmd::EntityKeyPair& k, const std::vector<double>& evalPoint,
                          const std::vector<double>& parametricCoords, InterpolationData& data) const override;

  virtual ~PatchRecoveryFieldInterpolator() { delete m_meshSelector; }

 protected:
  const stk::mesh::FieldBase* m_coordinates{nullptr};
  const stk::mesh::Part* m_meshPart{nullptr};
  stk::mesh::Selector* m_meshSelector{nullptr};
  const stk::mesh::Selector m_activeSelector;
  std::function<bool(const stk::transfer::EntityInterpolationData&,
                     const std::vector<double>&, unsigned, double*)> m_leastSquaresInterpolator;

  virtual void set_interpolator(const stk::transfer::PatchRecoveryEvaluationType type);

 private:
  EntityPatchFilter m_patchFilter;

  PatchRecoveryFieldInterpolator(const PatchRecoveryFieldInterpolator&) = delete;
  const PatchRecoveryFieldInterpolator& operator()(const PatchRecoveryFieldInterpolator&) = delete;
};

class CopyFieldInterpolator : public InterpolateFieldsInterface {
 public:
  CopyFieldInterpolator(stk::mesh::BulkData& bulk);

  void interpolate_fields(const stk::search::spmd::EntityKeyPair& k, const std::vector<double>& evalPoint,
                          const std::vector<double>& parametricCoords, InterpolationData& data) const override;

  ~CopyFieldInterpolator() {}

 private:
  CopyFieldInterpolator(const CopyFieldInterpolator&) = delete;
  const CopyFieldInterpolator& operator()(const CopyFieldInterpolator&) = delete;
};

class SumFieldInterpolator : public InterpolateFieldsInterface {
 public:
  SumFieldInterpolator(stk::mesh::BulkData& bulk);

  void interpolate_fields(const stk::search::spmd::EntityKeyPair& k, const std::vector<double>& evalPoint,
                          const std::vector<double>& parametricCoords, InterpolationData& data) const override;

  ~SumFieldInterpolator() {}

 private:
  SumFieldInterpolator(const SumFieldInterpolator&) = delete;
  const SumFieldInterpolator& operator()(const SumFieldInterpolator&) = delete;
};

}
}

#endif /* STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_INTERPOLATEFIELDS_HPP_ */
