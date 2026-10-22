// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef AKRI_ORDEREDIDPAIR_H_
#define AKRI_ORDEREDIDPAIR_H_

#include <utility>
#include <stk_mesh/base/Types.hpp>

namespace krino {

class OrderedIdPair {
public:
  OrderedIdPair(stk::mesh::EntityId id0, stk::mesh::EntityId id1)
  : my_ids(id0 < id1 ? id0 : id1, id0 < id1 ? id1 : id0) {}

  stk::mesh::EntityId first() const { return my_ids.first; }
  stk::mesh::EntityId second() const { return my_ids.second; }

private:
  std::pair<stk::mesh::EntityId, stk::mesh::EntityId> my_ids;
};

inline bool operator<(const OrderedIdPair lhs, const OrderedIdPair rhs)
{
  if(lhs.first() < rhs.first()) return true;
  if(rhs.first() < lhs.first()) return false;
  if(lhs.second() < rhs.second()) return true;
  return false;
}

inline bool operator==(const OrderedIdPair lhs, const OrderedIdPair rhs)
{
  return (lhs.first() == rhs.first()) && (lhs.second() == rhs.second());
}

inline bool operator!=(const OrderedIdPair lhs, const OrderedIdPair rhs)
{
  return !(lhs == rhs);
}

}


#endif /* AKRI_ORDEREDIDPAIR_H_ */
