#ifndef ModificationObserver_hpp
#define ModificationObserver_hpp

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>

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

    virtual void local_entities_created_or_deleted_notification(stk::mesh::EntityRank rank)
    {
    }

    virtual void local_entity_comm_info_changed_notification(stk::mesh::EntityRank rank)
    {
    }

    virtual void local_buckets_changed_notification(stk::mesh::EntityRank rank)
    {
    }

    virtual void fill_values_to_reduce(std::vector<size_t> &valuesToReduce)
    {
        valuesToReduce.clear();
    }

    virtual void set_reduced_values(const std::vector<size_t> &reducedValues)
    {
    }
};

}
}

#endif
