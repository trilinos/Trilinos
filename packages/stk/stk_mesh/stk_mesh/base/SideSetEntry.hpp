
#ifndef SIDESETENTRY_HPP_
#define SIDESETENTRY_HPP_

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>
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

  stk::mesh::Entity element;
  stk::mesh::ConnectivityOrdinal side;
};


class SideSetEntryLess
{
public:
    SideSetEntryLess(const BulkData& mesh);
    bool operator()(const SideSetEntry& lhs, const SideSetEntry& rhs) const;
private:
  const BulkData& m_mesh;
};

class SideSetEntryEquals
{
public:
    SideSetEntryEquals(const BulkData& mesh);
    bool operator()(const SideSetEntry& lhs, const SideSetEntry& rhs) const;
private:
  const BulkData& m_mesh;
};

//////////////

inline
SideSetEntryLess::SideSetEntryLess(const BulkData& mesh) : m_mesh(mesh){}

inline
bool SideSetEntryLess::operator()(const SideSetEntry& lhs, const SideSetEntry& rhs) const
{
    if(m_mesh.identifier(lhs.element) < m_mesh.identifier(rhs.element))
        return true;
    else if(m_mesh.identifier(lhs.element) > m_mesh.identifier(rhs.element))
        return false;
    else
    {
        if(lhs.side<rhs.side)
            return true;
        else
            return false;
    }
    return false;
}

//////////////
inline
SideSetEntryEquals::SideSetEntryEquals(const BulkData& mesh) : m_mesh(mesh){}

inline
bool SideSetEntryEquals::operator()(const SideSetEntry& lhs, const SideSetEntry& rhs) const
{
    if(m_mesh.identifier(lhs.element) == m_mesh.identifier(rhs.element) &&
            lhs.side == rhs.side)
        return true;
    return false;
}

}
}

#endif /* SIDESETENTRY_HPP_ */
