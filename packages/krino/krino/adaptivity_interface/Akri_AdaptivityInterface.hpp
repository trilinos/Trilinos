// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_ADAPTIVITYINTERFACE_H_
#define KRINO_INCLUDE_AKRI_ADAPTIVITYINTERFACE_H_

#include <string>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace diag { class Timer; } }

namespace krino {

namespace HAdapt
{
  void setup(stk::mesh::MetaData & meta, stk::mesh::Part & active_part, stk::diag::Timer & root_timer);
  void do_adaptive_refinement(stk::mesh::MetaData & meta, const std::string & marker_field_name);
  void do_uniform_refinement(stk::mesh::MetaData & meta, const int num_levels);
  void mark_based_on_indicator_field(const stk::mesh::BulkData & bulk,
    const std::string & marker_field_name,
    const std::string & indicator_field_name,
    const int max_refinement_level,
    const int current_refinement_level,
    const uint64_t target_elem_count);
}

} // namespace krino

#endif /* KRINO_INCLUDE_AKRI_ADAPTIVITYINTERFACE_H_ */
