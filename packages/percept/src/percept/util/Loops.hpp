// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_Loops_hpp
#define percept_Loops_hpp

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

#include <percept/function/ElementOp.hpp>
#include <percept/function/BucketOp.hpp>
#include <percept/function/internal/GenericFunction.hpp>

namespace percept {

  /** \brief Loop over all buckets and apply \param bucketOp
      passing in the argument \param field to \param bucketOp
  */
  void bucketOpLoop(stk::mesh::BulkData& bulkData,
                    BucketOp& bucketOp,
                    stk::mesh::FieldBase *field = 0,
                    stk::mesh::Part *part = 0);

  void bucketOpLoop(stk::mesh::BulkData& bulkData,
                    BucketOp& bucketOp,
                    stk::mesh::FieldBase *field,
                    stk::mesh::Selector *selector,
                    bool is_surface_norm = false);

  void elementOpLoop(stk::mesh::BulkData& bulkData,
                     ElementOp& elementOp,
                     stk::mesh::FieldBase *field=0,
                     stk::mesh::Part *part = 0);

  void elementOpLoop(stk::mesh::BulkData& bulkData,
                     ElementOp& elementOp,
                     stk::mesh::FieldBase *field,
                     stk::mesh::Selector* selector,
                     bool is_surface_norm=false );

  void nodalOpLoop(stk::mesh::BulkData& bulkData,
                   GenericFunction& nodalOp,
                   stk::mesh::FieldBase *field=0,
                   stk::mesh::Selector* selector=0);

} // namespace percept

#endif
