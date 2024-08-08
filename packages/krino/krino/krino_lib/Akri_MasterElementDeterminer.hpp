// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_MasterElementDeterminer_h
#define Akri_MasterElementDeterminer_h

#include <Akri_MasterElement.hpp>

namespace krino { class FieldRef; }
namespace stk { namespace mesh { class Bucket; } }
namespace stk { struct topology; }

namespace krino {

class MasterElementDeterminer
{
public:
  static const MasterElement& getMasterElement(stk::mesh::Bucket & bucket, FieldRef field);
  static const MasterElement& getMasterElement(stk::topology topology);
  static stk::topology get_field_topology(const stk::mesh::Bucket & b, const FieldRef field);
  static void clear_master_elements();
private:
  static std::vector<std::unique_ptr<MasterElement>> theMasterElements;
};

}  // end namespace krino

#endif // Akri_MasterElementDeterminer_h
