// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef TransferHelper_h
#define TransferHelper_h

#include <stk_mesh/base/Field.hpp>

namespace shards {
  class CellTopology;
}

namespace stk {
  namespace mesh {
    struct Entity;
  }  
}

namespace percept {

enum TransferType {THREED_TO_THREED=0, 
		   TWOD_TO_TWOD,
		   TWOD_AXI_TO_THREED};

  enum SrcFieldType {SRC_FIELD=0, // generic field type
                     SRC_RZN_FIELD, // normal (out of plane) component for 2D axisymmetric -> 3D
                     SRC_RZP_FIELD}; // in plane components (2) for 2D axisymmetric -> 3D

double 
parametricDistanceToEntity(const double*                 point,
			   const shards::CellTopology &  cellTopo);
  
void compute_element_centroid(const stk::mesh::FieldBase & coordinates,
                              const stk::mesh::Entity & entity, double * coords);
void compute_nodal_coords(    const stk::mesh::FieldBase & coordinates,
                              const stk::mesh::Entity & entity, double * coords);
}

#endif
