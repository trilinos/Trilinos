// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef FieldTypes_hpp
#define FieldTypes_hpp

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldTraits.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stdint.h>

namespace percept {

  using GenericFieldType = stk::mesh::Field<double>;

  using ScalarFieldType         = stk::mesh::Field<double>;
  using ScalarIntFieldType      = stk::mesh::Field<int>;
  using ScalarEntityIdFieldType = stk::mesh::Field<int64_t>;
  using ScalarUnsignedFieldType = stk::mesh::Field<unsigned>;
  using VectorFieldType         = stk::mesh::Field<double>;
  using VectorIntFieldType      = stk::mesh::Field<int>;
  using ArrayDimType            = stk::mesh::Field<double>;
  using ArrayEntityIdDimType    = stk::mesh::Field<int64_t>;

  using RefineFieldType                 = ScalarIntFieldType;
  using RefineLevelType                 = ScalarIntFieldType;
  using TransitionElementType           = ScalarIntFieldType;
  using ParentElementType               = ScalarFieldType;
  using NewNodesType                    = ScalarIntFieldType;
  using ErrorFieldType                  = ScalarFieldType;
  using WeightsFieldType                = ScalarFieldType;
  using GregoryControlPointsType        = ScalarFieldType;
  using CoordinatesFieldType            = VectorFieldType;
  using UnprojectedCoordinatesFieldType = ArrayDimType;
  using NormalsFieldType                = CoordinatesFieldType;
  using WallDistanceFieldType           = ScalarIntFieldType;
  using NodeRegistryFieldType           = GenericFieldType;
  using SGCellNodeIdsType               = ArrayEntityIdDimType;

  using RefineFieldType_type                 = stk::mesh::FieldTraits<RefineFieldType>::data_type;
  using RefineLevelType_type                 = stk::mesh::FieldTraits<RefineLevelType>::data_type;
  using TransitionElementType_type           = stk::mesh::FieldTraits<TransitionElementType>::data_type;
  using ParentElementType_type               = stk::mesh::FieldTraits<ParentElementType>::data_type;
  using NewNodesType_type                    = stk::mesh::FieldTraits<NewNodesType>::data_type;
  using ErrorFieldType_type                  = stk::mesh::FieldTraits<ErrorFieldType>::data_type;
  using WallDistanceFieldType_type           = stk::mesh::FieldTraits<WallDistanceFieldType>::data_type;
  using CoordinatesFieldType_type            = stk::mesh::FieldTraits<CoordinatesFieldType>::data_type;
  using NodeRegistryFieldType_type           = stk::mesh::FieldTraits<NodeRegistryFieldType>::data_type;
  using SGCellNodeIdsType_type               = stk::mesh::FieldTraits<SGCellNodeIdsType>::data_type;
  using WeightsFieldType_type                = stk::mesh::FieldTraits<WeightsFieldType>::data_type;
  using UnprojectedCoordinatesFieldType_type = stk::mesh::FieldTraits<UnprojectedCoordinatesFieldType>::data_type;
  using GregoryControlPointsType_type        = stk::mesh::FieldTraits<GregoryControlPointsType>::data_type;
  using NormalsFieldType_type                = stk::mesh::FieldTraits<NormalsFieldType>::data_type;
  using ScalarFieldType_type                 = stk::mesh::FieldTraits<ScalarFieldType>::data_type;
}

#endif
