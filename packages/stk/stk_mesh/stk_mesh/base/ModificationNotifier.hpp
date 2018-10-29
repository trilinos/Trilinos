// Copyright (c) 2013, Sandia Corporation.
 // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 // the U.S. Government retains certain rights in this software.
 // 
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 // 
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 // 
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 // 
 //     * Neither the name of Sandia Corporation nor the names of its
 //       contributors may be used to endorse or promote products derived
 //       from this software without specific prior written permission.
 // 
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef ModificationNotifier_hpp
#define ModificationNotifier_hpp

#include <stk_mesh/base/ModificationObserver.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <memory>

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

    ~ModificationNotifier()
    {

    }

    void register_observer(std::shared_ptr<ModificationObserver> observer)
    {
        observers.push_back(observer);
    }

    void unregister_observer(std::shared_ptr<ModificationObserver> observer)
    {
        auto iter = std::find(observers.begin(), observers.end(), observer);
        if(iter != observers.end())
            observers.erase(iter);
    }

    void notify_entity_added(stk::mesh::Entity entity)
    {
        for(std::shared_ptr<ModificationObserver>& observer : observers)
        {
            observer->entity_added(entity);
        }
    }

    void notify_entity_deleted(stk::mesh::Entity entity)
    {
        for(std::shared_ptr<ModificationObserver>& observer : observers)
        {
            observer->entity_deleted(entity);
        }
    }

    void notify_entity_parts_added(stk::mesh::Entity entity, const stk::mesh::OrdinalVector& parts)
    {
        for(std::shared_ptr<ModificationObserver>& observer : observers)
        {
            observer->entity_parts_added(entity, parts);
        }
    }

    void notify_entity_parts_removed(stk::mesh::Entity entity, const stk::mesh::OrdinalVector& parts)
    {
        for(std::shared_ptr<ModificationObserver>& observer : observers)
        {
            observer->entity_parts_removed(entity, parts);
        }
    }

    void notify_started_modification_end()
    {
        for(std::shared_ptr<ModificationObserver>& observer : observers)
        {
            observer->started_modification_end_notification();
        }
    }

    void notify_finished_modification_end(MPI_Comm communicator)
    {
        reduce_values_for_observers(communicator);
        for(std::shared_ptr<ModificationObserver>& observer : observers)
        {
            observer->finished_modification_end_notification();
        }
    }

    void notify_elements_about_to_move_procs(const stk::mesh::EntityProcVec &elemProcPairsToMove)
    {
        for(std::shared_ptr<ModificationObserver>& observer : observers)
        {
            observer->elements_about_to_move_procs_notification(elemProcPairsToMove);
        }
    }

    void notify_elements_moved_procs(const stk::mesh::EntityProcVec &elemProcPairsToMove)
    {
        for(std::shared_ptr<ModificationObserver>& observer : observers)
        {
            observer->elements_moved_procs_notification(elemProcPairsToMove);
        }
    }

    void notify_local_entities_created_or_deleted(stk::mesh::EntityRank rank)
    {
        for(std::shared_ptr<ModificationObserver>& observer : observers)
        {
            observer->local_entities_created_or_deleted_notification(rank);
        }
    }

    void notify_local_entity_comm_info_changed(stk::mesh::EntityRank rank)
    {
        for(std::shared_ptr<ModificationObserver>& observer : observers)
        {
            observer->local_entity_comm_info_changed_notification(rank);
        }
    }

    void notify_local_buckets_changed(stk::mesh::EntityRank rank)
    {
        for(std::shared_ptr<ModificationObserver>& observer : observers)
        {
            observer->local_buckets_changed_notification(rank);
        }
    }

    template<typename ObserverType>
    bool has_observer_type() const
    {
        bool found = false;

        for(const std::shared_ptr<ModificationObserver>& observer : observers)
        {
            if (dynamic_cast<const ObserverType*>(observer.get()) != nullptr)
            {
                found = true;
                break;
            }
        }
        return found;
    }

    template<typename ObserverType>
    std::vector<ObserverType*> get_observer_type() const
    {
        std::vector<ObserverType*> typed_observers;

        for(const std::shared_ptr<ModificationObserver> &observer : observers)
        {
            ObserverType* typed_observer = dynamic_cast<ObserverType*>(observer.get());
            if (typed_observer != nullptr)
            {
                typed_observers.push_back(typed_observer);
            }
        }

        return typed_observers;
    }

private:
    std::vector< std::shared_ptr<ModificationObserver> > observers;

    void reduce_values_for_observers(MPI_Comm communicator);
    void get_values_to_reduce_from_observers(std::vector<size_t> &allObserverValues, std::vector<std::vector<size_t> >& observerValues);
    void set_max_values_on_observers(const std::vector<size_t> &maxValues, std::vector<std::vector<size_t> >& observerValues);
};

}
}

#endif
