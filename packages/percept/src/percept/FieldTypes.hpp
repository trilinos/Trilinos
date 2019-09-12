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

  typedef stk::mesh::Field<double, stk::mesh::SimpleArrayTag>  GenericFieldType;

  typedef stk::mesh::Field<double>                          ScalarFieldType ;
  typedef stk::mesh::Field<int>                             ScalarIntFieldType ;
  typedef stk::mesh::Field<int64_t>                         ScalarEntityIdFieldType ;
  typedef stk::mesh::Field<unsigned>                        ScalarUnsignedFieldType ;
  typedef stk::mesh::Field<double, stk::mesh::Cartesian>    VectorFieldType ;
  typedef stk::mesh::Field<int, stk::mesh::Cartesian>       VectorIntFieldType ;
  typedef stk::mesh::Field<double, shards::ArrayDimension>  ArrayDimType;
  typedef stk::mesh::Field<int64_t, shards::ArrayDimension>  ArrayEntityIdDimType;

  typedef ScalarIntFieldType      RefineFieldType;
  typedef ScalarIntFieldType      RefineLevelType;
  typedef ScalarIntFieldType      TransitionElementType;
  typedef ScalarFieldType         ParentElementType;
  typedef ScalarIntFieldType      NewNodesType;
  typedef GenericFieldType        ErrorFieldType;
  typedef ScalarFieldType         WeightsFieldType;
  typedef ArrayDimType            GregoryControlPointsType;
  typedef VectorFieldType         CoordinatesFieldType;
  typedef ArrayDimType            UnprojectedCoordinatesFieldType;
  typedef CoordinatesFieldType    NormalsFieldType;
  typedef ScalarIntFieldType      WallDistanceFieldType;
  typedef GenericFieldType        NodeRegistryFieldType;
  typedef ArrayEntityIdDimType    SGCellNodeIdsType;

  typedef stk::mesh::FieldTraits<RefineFieldType>::data_type       RefineFieldType_type;
  typedef stk::mesh::FieldTraits<RefineLevelType>::data_type       RefineLevelType_type;
  typedef stk::mesh::FieldTraits<TransitionElementType>::data_type TransitionElementType_type;
  typedef stk::mesh::FieldTraits<ParentElementType>::data_type     ParentElementType_type;
  typedef stk::mesh::FieldTraits<NewNodesType>::data_type          NewNodesType_type;
  typedef stk::mesh::FieldTraits<ErrorFieldType>::data_type        ErrorFieldType_type;
  typedef stk::mesh::FieldTraits<WallDistanceFieldType>::data_type WallDistanceFieldType_type;
  typedef stk::mesh::FieldTraits<CoordinatesFieldType>::data_type  CoordinatesFieldType_type;
  typedef stk::mesh::FieldTraits<NodeRegistryFieldType>::data_type NodeRegistryFieldType_type;
  typedef stk::mesh::FieldTraits<SGCellNodeIdsType>::data_type     SGCellNodeIdsType_type;

}

#endif
