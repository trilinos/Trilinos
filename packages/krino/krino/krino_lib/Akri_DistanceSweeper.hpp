// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_DistanceSweeper_h
#define Akri_DistanceSweeper_h

namespace stk { namespace mesh { class BulkData; } }
namespace krino { class FieldRef; }

namespace krino {
namespace DistanceSweeper {

  void fix_sign_by_sweeping(const stk::mesh::BulkData& mesh, const FieldRef distance_field, const double signed_narrow_band);

}}

#endif // Akri_DistanceSweeper_h
