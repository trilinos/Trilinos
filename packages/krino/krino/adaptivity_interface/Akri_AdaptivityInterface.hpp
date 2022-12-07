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

class PerceptRefinement;
class RefinementInterface;

PerceptRefinement & create_percept_refinement(stk::mesh::MetaData & meta, stk::diag::Timer & parentTimer);
RefinementInterface & create_refinement(stk::mesh::MetaData & meta, const bool usePercept, stk::diag::Timer & parentTimer);

namespace HAdapt
{
  void setup(stk::mesh::MetaData & meta, stk::mesh::Part & active_part, stk::diag::Timer & root_timer);
  void do_adaptive_refinement(stk::mesh::MetaData & meta, const std::string & marker_field_name);
  void do_uniform_refinement(stk::mesh::MetaData & meta, const int num_levels);
}

} // namespace krino

#endif /* KRINO_INCLUDE_AKRI_ADAPTIVITYINTERFACE_H_ */
