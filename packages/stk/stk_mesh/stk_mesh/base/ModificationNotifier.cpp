#include <stk_mesh/base/ModificationNotifier.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

namespace stk
{
namespace mesh
{

std::vector<size_t> get_global_max(const std::vector<size_t> &allObserverValues, MPI_Comm communicator)
{
    std::vector<size_t> maxValues(allObserverValues.size());
    stk::all_reduce_max(communicator, allObserverValues.data(), maxValues.data(), maxValues.size());
    return maxValues;
}

void ModificationNotifier::reduce_values_for_observers(MPI_Comm communicator)
{
    std::vector<size_t> allObserverValues;
    std::vector<std::vector<size_t> > observerValues(observers.size());
    get_values_to_reduce_from_observers(allObserverValues, observerValues);
    std::vector<size_t> maxValues = get_global_max(allObserverValues, communicator);
    set_max_values_on_observers(maxValues, observerValues);
}

void ModificationNotifier::get_values_to_reduce_from_observers(std::vector<size_t> &allObserverValues,
                                                       std::vector<std::vector<size_t> >& observerValues)
{
    for(size_t i = 0; i < observers.size(); ++i)
    {
        observers[i]->fill_values_to_reduce(observerValues[i]);
        allObserverValues.insert(allObserverValues.end(), observerValues[i].begin(), observerValues[i].end());
    }
}
void ModificationNotifier::set_max_values_on_observers(const std::vector<size_t> &maxValues,
                                                       std::vector<std::vector<size_t> >& observerValues)
{
    std::vector<size_t>::const_iterator start = maxValues.begin();
    for(size_t i = 0; i < observers.size(); ++i)
    {
        observerValues[i].assign(start, start + observerValues[i].size());
        observers[i]->set_reduced_values(observerValues[i]);
        start += observerValues[i].size();
    }
}

}
}
