#ifndef ModificationNotifier_hpp
#define ModificationNotifier_hpp

#include <stk_mesh/base/ModificationObserver.hpp>
#include <stk_mesh/base/Types.hpp>

namespace stk { namespace mesh { class Entity; } }

namespace stk
{
namespace mesh
{

class ModificationNotifier
{
public:
    ModificationNotifier()
    {
    }

    void register_observer(ModificationObserver *observer)
    {
        observers.push_back(observer);
    }

    void entity_added(stk::mesh::Entity entity)
    {
        for(ModificationObserver *observer : observers)
        {
            observer->entity_added(entity);
        }
    }

    void entity_deleted(stk::mesh::Entity entity)
    {
        for(ModificationObserver *observer : observers)
        {
            observer->entity_deleted(entity);
        }
    }

    void notify_started_modification_end()
    {
        for(ModificationObserver *observer : observers)
        {
            observer->started_modification_end_notification();
        }
    }

    void notify_finished_modification_end()
    {
        for(ModificationObserver *observer : observers)
        {
            observer->finished_modification_end_notification();
        }
    }

    void notify_elements_about_to_move_procs(const stk::mesh::EntityProcVec &elemProcPairsToMove)
    {
        for(ModificationObserver *observer : observers)
        {
            observer->elements_about_to_move_procs_notification(elemProcPairsToMove);
        }
    }

    void notify_elements_moved_procs(const stk::mesh::EntityProcVec &elemProcPairsToMove)
    {
        for(ModificationObserver *observer : observers)
        {
            observer->elements_moved_procs_notification(elemProcPairsToMove);
        }
    }
private:
    std::vector<ModificationObserver *> observers;
};

}
}

#endif
