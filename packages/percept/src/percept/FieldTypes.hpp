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

}

#endif
