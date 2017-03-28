
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

  bool operator==(const SideSetEntry &rhs) const
  {
      return ((element == rhs.element) && (side == rhs.side));
  }

  bool operator<(const SideSetEntry &rhs) const
  {
      if(element < rhs.element)
          return true;
      else if (element == rhs.element && side < rhs.side)
          return true;
      else return false;
  }

  stk::mesh::Entity element;
  stk::mesh::ConnectivityOrdinal side;
};

typedef std::vector<SideSetEntry> SideSet;

}
}

#endif /* SIDESETENTRY_HPP_ */
