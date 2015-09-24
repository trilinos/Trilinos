#ifndef ModificationObserver_hpp
#define ModificationObserver_hpp

#include <stk_mesh/base/Types.hpp>

namespace stk { namespace mesh { class Entity; } }

namespace stk
{
namespace mesh
{

class ModificationObserver
{
public:
    virtual ~ModificationObserver()
    {
    }

    virtual void entity_added(stk::mesh::Entity entity)
    {
    }

    virtual void entity_deleted(stk::mesh::Entity entity)
    {
    }

    virtual void started_modification_end_notification()
    {
    }

    virtual void finished_modification_end_notification()
    {
    }

    virtual void elements_about_to_move_procs_notification(const stk::mesh::EntityProcVec &elemProcPairsToMove)
    {
    }

    virtual void elements_moved_procs_notification(const stk::mesh::EntityProcVec &elemProcPairsToMove)
    {
    }
};

}
}

#endif
