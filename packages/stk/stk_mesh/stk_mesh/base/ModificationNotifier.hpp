#ifndef ModificationNotifier_hpp
#define ModificationNotifier_hpp

#include <stk_mesh/base/ModificationObserver.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>

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

    void notify_entity_added(stk::mesh::Entity entity)
    {
        for(ModificationObserver *observer : observers)
        {
            observer->entity_added(entity);
        }
    }

    void notify_entity_deleted(stk::mesh::Entity entity)
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

    void notify_finished_modification_end(MPI_Comm communicator)
    {
        reduce_values_for_observers(communicator);
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

    void notify_local_entities_created_or_deleted(stk::mesh::EntityRank rank)
    {
        for(ModificationObserver *observer : observers)
        {
            observer->local_entities_created_or_deleted_notification(rank);
        }
    }

    void notify_local_entity_comm_info_changed(stk::mesh::EntityRank rank)
    {
        for(ModificationObserver *observer : observers)
        {
            observer->local_entity_comm_info_changed_notification(rank);
        }
    }

    void notify_local_buckets_changed(stk::mesh::EntityRank rank)
    {
        for(ModificationObserver *observer : observers)
        {
            observer->local_buckets_changed_notification(rank);
        }
    }

private:
    std::vector<ModificationObserver *> observers;

    void reduce_values_for_observers(MPI_Comm communicator);
    void get_values_to_reduce_from_observers(std::vector<size_t> &allObserverValues, std::vector<std::vector<size_t> >& observerValues);
    void set_max_values_on_observers(const std::vector<size_t> &maxValues, std::vector<std::vector<size_t> >& observerValues);
};

}
}

#endif
