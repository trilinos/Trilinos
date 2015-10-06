#ifndef DeletedElementInfo_hpp
#define DeletedElementInfo_hpp

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace stk
{
namespace mesh
{
namespace impl
{

struct DeletedElementInfo
{
    stk::mesh::Entity entity;
    stk::mesh::EntityId identifier;
    bool isShell;
};

typedef std::vector<DeletedElementInfo> DeletedElementInfoVector;

}
}
}

#endif
