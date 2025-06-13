/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_TRANSFER_UTIL_ENTITY_CENTROID_RECOVER_FIELD_HPP
#define STK_TRANSFER_UTIL_ENTITY_CENTROID_RECOVER_FIELD_HPP

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <stk_transfer_util/RecoverField.hpp>  // for RecoverField, RecoverField::...
#include "stk_transfer/TransferTypes.hpp"
#include <vector>                        // for vector
#include <functional>
#include "stk_mesh/base/Bucket.hpp"      // for Bucket
#include "stk_mesh/base/BulkData.hpp"    // for BulkData
#include "stk_mesh/base/Entity.hpp"      // for Entity
#include "stk_mesh/base/FieldBase.hpp"   // for FieldBase, field_data
#include "stk_mesh/base/Selector.hpp"    // for Selector
namespace stk { namespace mesh { class Part; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

// namespace stk

namespace stk {
namespace transfer {

class EntityCentroidLinearRecoverField : public RecoverField {
 public:
  EntityCentroidLinearRecoverField(RecoverField::RecoveryType recType,
                                   const std::vector<const stk::mesh::FieldBase*>& recVars,
                                   const stk::mesh::FieldBase& recNodeVar, int nSampElem,
                                   stk::mesh::Entity entity,
                                   FieldTransform  transform = [](double value) {return value;});

  ~EntityCentroidLinearRecoverField() override;

  void sample_patch(const std::vector<stk::mesh::Entity>& patch, int nSampPatch, double* fieldSample,
                    double* basisSample) const override;

 protected:
  std::vector<const stk::mesh::FieldBase*> m_recoverVars;
  const stk::mesh::FieldBase& m_nodeVar;
  FieldTransform  m_transform;
  mutable std::vector<double> m_scratchSpace;
};

class EntityCentroidQuadraticRecoverField : public RecoverField {
 public:
  EntityCentroidQuadraticRecoverField(RecoverField::RecoveryType recType,
                                      const std::vector<const stk::mesh::FieldBase*>& recVars,
                                      const stk::mesh::FieldBase& recNodeVar, int nSampElem,
                                      stk::mesh::Entity entity,
                                      FieldTransform  transform = [](double value) {return value;});

  ~EntityCentroidQuadraticRecoverField() override;

  void sample_patch(const std::vector<stk::mesh::Entity>& patch, int nSampPatch, double* fieldSample,
                    double* basisSample) const override;

 protected:
  std::vector<const stk::mesh::FieldBase*> m_recoverVars;
  const stk::mesh::FieldBase& m_nodeVar;
  FieldTransform  m_transform;
  mutable std::vector<double> m_scratchSpace;
};

class EntityCentroidCubicRecoverField : public RecoverField {
 public:
  EntityCentroidCubicRecoverField(RecoverField::RecoveryType recType,
                                  const std::vector<const stk::mesh::FieldBase*>& recVars,
                                  const stk::mesh::FieldBase& recNodeVar, int nSampElem,
                                  stk::mesh::Entity entity,
                                  FieldTransform  transform = [](double value) {return value;});

  ~EntityCentroidCubicRecoverField() override;

  void sample_patch(const std::vector<stk::mesh::Entity>& patch, int nSampPatch, double* fieldSample,
                    double* basisSample) const override;

 protected:
  std::vector<const stk::mesh::FieldBase*> m_recoverVars;
  const stk::mesh::FieldBase& m_nodeVar;
  FieldTransform  m_transform;
  mutable std::vector<double> m_scratchSpace;
};

struct EntityPatchFilter {
  EntityPatchFilter(const stk::mesh::Part* part,
                    const stk::mesh::FieldBase* field = nullptr,
                    const stk::mesh::Selector* selector = nullptr)
    : m_meshPart(part)
    , m_field(field)
    , m_selector(selector)
  {
  }

  EntityPatchFilter(const stk::mesh::Part* part,
                    const stk::mesh::Selector* selector = nullptr)
    : m_meshPart(part)
    , m_field(nullptr)
    , m_selector(selector)
  {
  }

  virtual ~EntityPatchFilter() {}

  /// Return true if entity participates or if none are supplied. Return false otherwise.
  virtual bool pass(stk::mesh::Entity entity, const stk::mesh::BulkData& bulkData) const
  {
    bool isContainedInMeshPart = true;
    if(nullptr != m_meshPart) {
      isContainedInMeshPart = bulkData.bucket(entity).member(*m_meshPart);
    }

    if(isContainedInMeshPart && (nullptr != m_field)) {
      double* result = static_cast<double *>(stk::mesh::field_data(*m_field, entity));
      isContainedInMeshPart = (nullptr != result) ? true : false;
    }

    if(isContainedInMeshPart && (nullptr != m_selector)) {
      isContainedInMeshPart = (*m_selector)(bulkData.bucket(entity));
    }

    return isContainedInMeshPart;
  }

  void initialize(const stk::mesh::Part* part,
                  const stk::mesh::FieldBase* field = nullptr,
                  const stk::mesh::Selector* selector = nullptr)
  {
    m_meshPart = part;
    m_field = field;
    m_selector = selector;
  }

  void initialize(const stk::mesh::Part* part,
                  const stk::mesh::Selector* selector = nullptr)
  {
    m_meshPart = part;
    m_selector = selector;
  }

 private:
  const stk::mesh::Part* m_meshPart{nullptr};
  const stk::mesh::FieldBase* m_field{nullptr};
  const stk::mesh::Selector* m_selector{nullptr};
};

} // namespace transfer
} // namespace stk

#endif // STK_TRANSFER_UTIL_ENTITY_CENTROID_RECOVER_FIELD_HPP
