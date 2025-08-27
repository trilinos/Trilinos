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


#ifndef STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_FIELDEVALUATOR_HPP_
#define STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_FIELDEVALUATOR_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <limits>                            // for numeric_limits
#include <string>                            // for string
#include <vector>                            // for vector
#include "stk_mesh/base/Entity.hpp"          // for Entity
#include "stk_mesh/base/Selector.hpp"        // for Selector
#include "stk_mesh/base/Types.hpp"           // for EntityRank, ConstPartVector
#include <stk_mesh/base/BulkData.hpp>        // for BulkData
#include <stk_mesh/base/MetaData.hpp>        // for MetaData
#include "stk_mesh/base/FieldBase.hpp"       // for field_data, Fie...
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_topology/topology.hpp"         // for topology, topology::NODE...
#include "stk_util/util/ReportHandler.hpp"   // for ThrowRequireMsg

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk
{
namespace unit_test_util
{
inline double linear_function(double coeff_c,
                              double coeff_x, double coeff_y, double coeff_z,
                              double       x, double       y, double       z)
{
  return coeff_c + coeff_x * x + coeff_y * y + coeff_z * z;
}

struct FieldEvaluator {
  static constexpr unsigned InvalidIndex = std::numeric_limits<unsigned>::max();
  virtual double operator()(stk::mesh::Entity entity, const double x, const double y, const double z,
                            const unsigned index = InvalidIndex) const = 0;
  virtual ~FieldEvaluator(){};
};

struct ConstantFieldEvaluator : public FieldEvaluator {
  double operator()(stk::mesh::Entity /*entity*/, const double /*x*/, const double /*y*/, const double /*z*/,
                    [[maybe_unused]] const unsigned index = InvalidIndex) const
  {
    return value;
  }
  ConstantFieldEvaluator(double constantValue)
    : value(constantValue)
  {
  }
  ~ConstantFieldEvaluator() {}

 private:
  ConstantFieldEvaluator();
  ConstantFieldEvaluator(const ConstantFieldEvaluator&);
  const ConstantFieldEvaluator& operator()(const ConstantFieldEvaluator&);
  double value;
};

struct ConstantFieldIndexEvaluator : public FieldEvaluator {
  double operator()(stk::mesh::Entity /*entity*/, const double /*x*/, const double /*y*/, const double /*z*/,
                    const unsigned index = InvalidIndex) const
  {
    if(fieldIndex == 0 || index == fieldIndex) {
      return value;
    }

    return 0.0;
  }
  ConstantFieldIndexEvaluator(double constantValue, unsigned index = 0)
    : value(constantValue)
    , fieldIndex(index)
  {
  }
  ~ConstantFieldIndexEvaluator() {}

 private:
  ConstantFieldIndexEvaluator();
  ConstantFieldIndexEvaluator(const ConstantFieldIndexEvaluator&);
  const ConstantFieldIndexEvaluator& operator()(const ConstantFieldIndexEvaluator&);
  double value;
  unsigned fieldIndex;
};

struct PlanarLinearFieldEvaluator : public FieldEvaluator {
  virtual double operator()(stk::mesh::Entity /*entity*/, const double /*x*/, const double y, const double z,
                            [[maybe_unused]] const unsigned index = InvalidIndex) const
  {
    return m_spatialDimension == 2 ? linear_function(m_coeffC, m_coeffX, m_coeffY, m_coeffZ, 0, y, 0)
                                   : linear_function(m_coeffC, m_coeffX, m_coeffY, m_coeffZ, 0, y, z);
  }
  PlanarLinearFieldEvaluator(const unsigned spatialDimension,
                             double coeffC = 1.0, double coeffX = 1.0, double coeffY = 1.0, double coeffZ = 1.0)
    : m_spatialDimension(spatialDimension)
    , m_coeffC(coeffC)
    , m_coeffX(coeffX)
    , m_coeffY(coeffY)
    , m_coeffZ(coeffZ)
  {
    STK_ThrowRequireMsg(3 == m_spatialDimension || 2 == m_spatialDimension,
                    "Invalid spatial dimension" << m_spatialDimension);
  }

  ~PlanarLinearFieldEvaluator() {}

 protected:
  PlanarLinearFieldEvaluator(const PlanarLinearFieldEvaluator&);
  const PlanarLinearFieldEvaluator& operator()(const PlanarLinearFieldEvaluator&);
  PlanarLinearFieldEvaluator() {}

  unsigned m_spatialDimension = 3;
  double m_coeffC = 1.0;
  double m_coeffX = 1.0;
  double m_coeffY = 1.0;
  double m_coeffZ = 1.0;
};

struct LinearFieldEvaluator : public FieldEvaluator {
  virtual double operator()(stk::mesh::Entity /*entity*/, const double /*x*/, const double y, const double z,
                            [[maybe_unused]] const unsigned index = InvalidIndex) const
  {
    return 2 == m_spatialDimension ? linear_function(m_coeffC, m_coeffX, m_coeffY, m_coeffZ, 0, y, 0)
                                   : linear_function(m_coeffC, m_coeffX, m_coeffY, m_coeffZ, 0, y, z);
  }
  LinearFieldEvaluator(const unsigned spatialDimension,
                       double coeffC = 1.0, double coeffX = 1.0, double coeffY = 1.0, double coeffZ = 1.0)
    : m_spatialDimension(spatialDimension)
    , m_coeffC(coeffC)
    , m_coeffX(coeffX)
    , m_coeffY(coeffY)
    , m_coeffZ(coeffZ)
  {
    STK_ThrowRequireMsg((2 == m_spatialDimension) || (3 == m_spatialDimension),
                    "Invalid spatial dimension" << m_spatialDimension);
  }

  ~LinearFieldEvaluator() {}

 protected:
  LinearFieldEvaluator(const LinearFieldEvaluator&);
  const LinearFieldEvaluator& operator()(const LinearFieldEvaluator&);
  LinearFieldEvaluator() {}

  unsigned m_spatialDimension = 3;
  double m_coeffC = 1.0;
  double m_coeffX = 1.0;
  double m_coeffY = 1.0;
  double m_coeffZ = 1.0;
};

struct BoundedLinearFieldEvaluator : public LinearFieldEvaluator {
  void apply_bounds(const unsigned length, double* fieldData, const double lowerBound, const double upperBound) const
  {
    if(std::numeric_limits<double>::lowest() != lowerBound) {
      for(unsigned i(0); i < length; ++i) {
        if(fieldData[i] < lowerBound) {
          fieldData[i] = lowerBound;
        }
      }
    }
    if(std::numeric_limits<double>::max() != upperBound) {
      for(unsigned i(0); i < length; ++i) {
        if(upperBound < fieldData[i]) {
          fieldData[i] = upperBound;
        }
      }
    }
  }

  double operator()(stk::mesh::Entity entity, const double x, const double y, const double z,
                    [[maybe_unused]] const unsigned index = InvalidIndex) const
  {
    double value = LinearFieldEvaluator::operator()(entity, x, y, z);
    apply_bounds(1u, &value, m_lowerBound, m_upperBound);
    return value;
  }
  BoundedLinearFieldEvaluator(const unsigned spatialDimension, const double lowerBound, const double upperBound)
    : LinearFieldEvaluator(spatialDimension)
    , m_lowerBound(lowerBound)
    , m_upperBound(upperBound)
  {
    STK_ThrowRequireMsg(lowerBound <= upperBound, "Invalid bounds");
  }

  ~BoundedLinearFieldEvaluator() {}

 private:
  BoundedLinearFieldEvaluator(const BoundedLinearFieldEvaluator&);
  const BoundedLinearFieldEvaluator& operator()(const BoundedLinearFieldEvaluator&);
  BoundedLinearFieldEvaluator() {}

  double m_lowerBound = std::numeric_limits<double>::lowest();
  double m_upperBound = std::numeric_limits<double>::max();
};

struct EntityIdFieldEvaluator : public FieldEvaluator {
  double operator()(stk::mesh::Entity entity, const double /*x*/, const double /*y*/, const double /*z*/,
                    [[maybe_unused]] const unsigned index = InvalidIndex) const override
  {
    return static_cast<double>(m_bulk.identifier(entity));
  }

  EntityIdFieldEvaluator(const stk::mesh::BulkData& bulk)
    : m_bulk(bulk)
  {
  }
  ~EntityIdFieldEvaluator() {}

 private:
  EntityIdFieldEvaluator(const EntityIdFieldEvaluator&);
  const EntityIdFieldEvaluator& operator()(const EntityIdFieldEvaluator&);
  const stk::mesh::BulkData& m_bulk;
};

struct TriPolynomialFieldEvaluator : public FieldEvaluator {
  double evaluate_polynomial(double x, const std::vector<double>& coeffs) const
  {
    double value = 0.0;

    for(unsigned i = 0; i < m_dim; ++i) {
      value += coeffs[i] * std::pow(x, i);
    }

    return value;
  }

  double operator()(stk::mesh::Entity /*entity*/, const double x, const double y, const double z,
                    [[maybe_unused]] const unsigned index = InvalidIndex) const override
  {
    double xValue = evaluate_polynomial(x, m_xcoeffs);
    double yValue = evaluate_polynomial(y, m_ycoeffs);
    double zValue = evaluate_polynomial(z, m_zcoeffs);

    return xValue * yValue * zValue;
  }

  TriPolynomialFieldEvaluator(unsigned dim, const std::vector<double>& xcoeffs, const std::vector<double>& ycoeffs,
                              const std::vector<double>& zcoeffs)
    : m_dim(dim)
    , m_xcoeffs(xcoeffs)
    , m_ycoeffs(ycoeffs)
    , m_zcoeffs(zcoeffs)
  {
    if(m_xcoeffs.size() < dim) {
      m_xcoeffs.resize(dim, 0.0);
    }

    if(m_ycoeffs.size() < dim) {
      m_ycoeffs.resize(dim, 0.0);
    }

    if(m_zcoeffs.size() < dim) {
      m_zcoeffs.resize(dim, 0.0);
    }
  }

  ~TriPolynomialFieldEvaluator() {}

 protected:
  TriPolynomialFieldEvaluator(const TriPolynomialFieldEvaluator&);
  const TriPolynomialFieldEvaluator& operator()(const TriPolynomialFieldEvaluator&);
  TriPolynomialFieldEvaluator() {}

  unsigned m_dim = 0;

  std::vector<double> m_xcoeffs;
  std::vector<double> m_ycoeffs;
  std::vector<double> m_zcoeffs;
};

struct LinearFieldEvaluatorWithCoefficients : public TriPolynomialFieldEvaluator {
  LinearFieldEvaluatorWithCoefficients(const std::vector<double>& xcoeffs, const std::vector<double>& ycoeffs,
                                       const std::vector<double>& zcoeffs)
    : TriPolynomialFieldEvaluator(2, xcoeffs, ycoeffs, zcoeffs)
  {
  }

  ~LinearFieldEvaluatorWithCoefficients() {}

  using TriPolynomialFieldEvaluator::operator();

 protected:
  LinearFieldEvaluatorWithCoefficients(const LinearFieldEvaluatorWithCoefficients&);
  const LinearFieldEvaluatorWithCoefficients& operator()(const LinearFieldEvaluatorWithCoefficients&);
  LinearFieldEvaluatorWithCoefficients() {}
};

struct QuadraticFieldEvaluatorWithCoefficients : public TriPolynomialFieldEvaluator {
  QuadraticFieldEvaluatorWithCoefficients(const std::vector<double>& xcoeffs, const std::vector<double>& ycoeffs,
                                          const std::vector<double>& zcoeffs)
    : TriPolynomialFieldEvaluator(3, xcoeffs, ycoeffs, zcoeffs)
  {
  }

  ~QuadraticFieldEvaluatorWithCoefficients() {}

  using TriPolynomialFieldEvaluator::operator();

 protected:
  QuadraticFieldEvaluatorWithCoefficients(const QuadraticFieldEvaluatorWithCoefficients&);
  const QuadraticFieldEvaluatorWithCoefficients& operator()(const QuadraticFieldEvaluatorWithCoefficients&);
  QuadraticFieldEvaluatorWithCoefficients() {}
};

struct CubicFieldEvaluatorWithCoefficients : public TriPolynomialFieldEvaluator {
  CubicFieldEvaluatorWithCoefficients(const std::vector<double>& xcoeffs, const std::vector<double>& ycoeffs,
                                      const std::vector<double>& zcoeffs)
    : TriPolynomialFieldEvaluator(4, xcoeffs, ycoeffs, zcoeffs)
  {
  }

  ~CubicFieldEvaluatorWithCoefficients() {}

  using TriPolynomialFieldEvaluator::operator();

 protected:
  CubicFieldEvaluatorWithCoefficients(const CubicFieldEvaluatorWithCoefficients&);
  const CubicFieldEvaluatorWithCoefficients& operator()(const CubicFieldEvaluatorWithCoefficients&);
  CubicFieldEvaluatorWithCoefficients() {}
};

inline void determine_centroid(const int spatialDimension, stk::mesh::Entity entity,
                               const stk::mesh::FieldBase& nodalCoordField, double* centroid)
{
  const stk::mesh::BulkData& bulkData = nodalCoordField.get_mesh();
  const stk::mesh::Entity* const nodes = bulkData.begin_nodes(entity);
  const unsigned numNodes = bulkData.num_nodes(entity);

  stk::mesh::field_data_execute<double, stk::mesh::ReadOnly>(nodalCoordField,
    [&](auto& coordData) {
      for(stk::mesh::ComponentIdx i = 0_comp; i < spatialDimension; ++i) {
        for(unsigned iNode = 0; iNode < numNodes; ++iNode) {
          auto coords = coordData.entity_values(nodes[iNode]);
          centroid[i] += coords(i);
        }
      }
      for(int i = 0; i < spatialDimension; ++i) {
        centroid[i] /= numNodes;
      }
    }
  );
}

inline void determine_centroid(const unsigned spatialDimension, stk::mesh::Entity entity,
                               const stk::mesh::FieldBase& nodalCoordField, std::vector<double>& centroid)
{
  // Collect up the coordinate values.
  centroid.clear();
  centroid.resize(spatialDimension);

  determine_centroid(spatialDimension, entity, nodalCoordField, centroid.data());
}

inline void set_entity_field(const stk::mesh::BulkData& bulk,
                             const stk::mesh::FieldBase& stkField,
                             const FieldEvaluator& eval)
{
  stk::mesh::EntityRank rank = stkField.entity_rank();
  STK_ThrowRequireMsg(stk::topology::NODE_RANK != rank, "Input entity rank cannot be NODE_RANK");

  stk::mesh::EntityVector entities;
  stk::mesh::Selector selector(stkField);
  stk::mesh::get_selected_entities(selector, bulk.buckets(rank), entities);

  const unsigned spatialDimension = bulk.mesh_meta_data().spatial_dimension();
  stk::mesh::FieldBase const* coord = bulk.mesh_meta_data().coordinate_field();

  stk::mesh::field_data_execute<double, stk::mesh::ReadWrite>(stkField,
    [&](auto& stkFieldData) {
      for(stk::mesh::Entity entity : entities) {
        std::vector<double> centroid;
        determine_centroid(spatialDimension, entity, *coord, centroid);

        double x = centroid[0];
        double y = centroid[1];
        double z = (spatialDimension == 2 ? 0.0 : centroid[2]);

        auto data = stkFieldData.entity_values(entity);
        for(stk::mesh::ComponentIdx i : data.components()) {
          data(i) = eval(entity, x, y, z, i + 1);
        }
      }
    }
  );
}

inline void set_node_field(const stk::mesh::BulkData& bulk,
                           const stk::mesh::FieldBase& stkField,
                           const FieldEvaluator& eval)
{
  stk::mesh::EntityRank rank = stkField.entity_rank();
  STK_ThrowRequireMsg(stk::topology::NODE_RANK == rank, "Input entity rank must be NODE_RANK");

  stk::mesh::EntityVector entities;
  stk::mesh::Selector selector(stkField);
  stk::mesh::get_selected_entities(selector, bulk.buckets(rank), entities);

  const unsigned spatialDimension = bulk.mesh_meta_data().spatial_dimension();
  const auto& coordField = *bulk.mesh_meta_data().coordinate_field();

  stk::mesh::field_data_execute<double, double, stk::mesh::ReadOnly, stk::mesh::ReadWrite>(coordField, stkField,
    [&](auto& coordFieldData, auto& stkFieldData) {
      for(stk::mesh::Entity entity : entities) {
        auto coordData = coordFieldData.entity_values(entity);
        auto fieldData = stkFieldData.entity_values(entity);
        double x = coordData(0_comp);
        double y = coordData(1_comp);
        double z = (spatialDimension == 2 ? 0.0 : coordData(2_comp));

        for(stk::mesh::ComponentIdx i : fieldData.components()) {
          fieldData(i) = eval(entity, x, y, z, i + 1);
        }
      }
    }
  );
}


}
}



#endif /* STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_FIELDEVALUATOR_HPP_ */
