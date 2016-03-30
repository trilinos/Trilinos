#ifndef ModificationNotifier_hpp
#define ModificationNotifier_hpp

#include <stk_mesh/base/ModificationObserver.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

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
        std::vector<size_t> allObserverValues;
        std::vector<std::vector<size_t>> observerValues(observers.size());
        for(size_t i=0; i<observers.size(); ++i)
        {
            observers[i]->fill_values_to_reduce(observerValues[i]);
            allObserverValues.insert(allObserverValues.end(), observerValues[i].begin(), observerValues[i].end());
        }
        std::vector<size_t> maxValues(allObserverValues.size());
        stk::all_reduce_max(communicator, allObserverValues.data(), maxValues.data(), maxValues.size());
        size_t offset = 0;
        for(size_t i=0; i<observers.size(); ++i)
        {
            observerValues[i].assign(maxValues.begin()+offset, maxValues.begin()+offset+observerValues[i].size());
            offset += observerValues[i].size();
            observers[i]->set_reduced_values(observerValues[i]);
        }

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
};

}
}

#endif
