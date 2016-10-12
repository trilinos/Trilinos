
#ifndef SIDESETENTRY_HPP_
#define SIDESETENTRY_HPP_

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace stk
{
namespace mesh
{

struct SideSetEntry
{
  SideSetEntry() : element(stk::mesh::Entity()), side(stk::mesh::INVALID_CONNECTIVITY_ORDINAL){};
  SideSetEntry(stk::mesh::Entity in_element, stk::mesh::ConnectivityOrdinal in_side)
    : element(in_element),
      side(in_side)
  {  }
  SideSetEntry(stk::mesh::Entity in_element, int in_side)
    : SideSetEntry(in_element, static_cast<stk::mesh::ConnectivityOrdinal>(in_side))
  {  }

  stk::mesh::Entity element;
  stk::mesh::ConnectivityOrdinal side;
};

typedef std::vector<SideSetEntry> StkSideSet;

}
}

#endif /* SIDESETENTRY_HPP_ */
