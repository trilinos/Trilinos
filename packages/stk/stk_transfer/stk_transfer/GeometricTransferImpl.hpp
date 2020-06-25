/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_STK_TRANSFER_STK_TRANSFER_GEOMETRICTRANSFERIMPL_HPP_
#define STK_STK_TRANSFER_STK_TRANSFER_GEOMETRICTRANSFERIMPL_HPP_

#include <algorithm>
#include <vector>
#include <stk_util/environment/Env.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>

namespace stk {
namespace transfer {

namespace impl {

template <class BoundingBoxType>
struct BoundingBoxCompare{

  bool operator()(const BoundingBoxType & a, const BoundingBoxType & b) const
  {
    return a.second.id() < b.second.id();
  }
};


template <class BoundingBoxB, class EntityProcB>
  struct compare {
    bool operator()(const BoundingBoxB &a, const EntityProcB  &b) const
    {
      return a.second < b;
    }

    bool operator()(const EntityProcB  &a, const BoundingBoxB &b) const
    {
      return a < b.second;
    }
  };

template <class ForwardIterator, class Compare>
bool local_is_sorted(ForwardIterator first, ForwardIterator last, Compare compare)
{
  if (first == last) return true;
  ForwardIterator next = first;
  while (++next!=last) {
      if ( compare(*next, *first) ) return false;
      ++first;
  }
  return true;
}


template <class INTERPOLATE>
void delete_range_points_found(
                               std::vector<typename INTERPOLATE::MeshB::BoundingBox>            &range_vector,
                               const typename INTERPOLATE::EntityProcRelationVec          &del)
{

  std::vector<typename INTERPOLATE::EntityProcB> range_entities_found;
  range_entities_found.reserve(del.size());
  for (auto && i : del) {
    range_entities_found.push_back(i.first);
  }
  {
    std::sort(range_entities_found.begin(), range_entities_found.end());
    const typename std::vector<typename INTERPOLATE::EntityProcB>::iterator it = std::unique(range_entities_found.begin(), range_entities_found.end());
    range_entities_found.resize(it-range_entities_found.begin());
  }

  std::vector<typename INTERPOLATE::MeshB::BoundingBox> difference(range_vector.size());
  {
    const auto it =
      std::set_difference(
        range_vector.        begin(), range_vector.        end(),
        range_entities_found.begin(), range_entities_found.end(),
        difference.begin(), compare<typename INTERPOLATE::MeshB::BoundingBox, typename INTERPOLATE::EntityProcB>());
    difference.resize(it-difference.begin());
  }
  swap(difference, range_vector);
}

template <typename T, typename=int>
struct CallOptionalPostCoarseSearchFilter
{
  static void call(typename T::EntityProcRelationVec &BtoA,
      const typename T::MeshA * mesha,
      const typename T::MeshB * meshb)
  {
    //Doesn't exist, so don't call
  }
};

template <typename T>
struct CallOptionalPostCoarseSearchFilter<T, decltype((void) T::post_coarse_search_filter, 0)>
{
  static void call(typename T::EntityProcRelationVec &BtoA,
      const typename T::MeshA * mesha,
      const typename T::MeshB * meshb)
  {
    ThrowRequire(mesha);
    ThrowRequire(meshb);
    T::post_coarse_search_filter(BtoA, *mesha, *meshb);
  }
};

template <typename T>
auto call_copy_entities(T & mesh, typename T::EntityProcVec & entities_to_copy, const std::string & name)
-> decltype( (void) mesh.copy_entities(entities_to_copy, name))
{
  mesh.copy_entities(entities_to_copy, name);
}

template <typename... T>
void call_copy_entities(T & ... t)
{
  //copy entities doesn't exist, so NO-OP
}

template <class INTERPOLATE>
void print_expansion_warnings(stk::ParallelMachine comm, int not_empty_count, size_t range_vector_size) {
  // sum and provide message to user
  size_t g_range_vector_size = 0;
  stk::all_reduce_max( comm, &range_vector_size, &g_range_vector_size, 1);
  sierra::Env::outputP0() << "GeometricTransfer<INTERPOLATE>::coarse_search(): Number of points not found: " << g_range_vector_size
      << " after expanding bounding boxes: " << not_empty_count << " time(s)" << std::endl;
  sierra::Env::outputP0() << "...will now expand the set of candidate bounding boxes and re-attempt the coarse search" << std::endl;
}

template <class INTERPOLATE>
void coarse_search_impl(typename INTERPOLATE::EntityProcRelationVec   &range_to_domain,
    stk::ParallelMachine comm,
    const typename INTERPOLATE::MeshA             * mesha,
    const typename INTERPOLATE::MeshB             * meshb,
    const stk::search::SearchMethod search_method,
    const double expansion_factor)
{
  using BoundingBoxA = typename INTERPOLATE::MeshA::BoundingBox;
  using BoundingBoxB = typename INTERPOLATE::MeshB::BoundingBox;
  std::vector<BoundingBoxB> range_vector;
  std::vector<BoundingBoxA> domain_vector;

  if (mesha)
    mesha->bounding_boxes(domain_vector);

  const bool domainVectorEmpty = stk::is_true_on_all_procs(comm, domain_vector.empty());
  if (domainVectorEmpty) {
    return;
  }

  if (meshb)
    meshb->bounding_boxes(range_vector);

  if( !local_is_sorted( domain_vector.begin(), domain_vector.end(), BoundingBoxCompare<BoundingBoxA>() ) )
    std::sort(domain_vector.begin(),domain_vector.end(),BoundingBoxCompare<BoundingBoxA>());
  if( !local_is_sorted( range_vector.begin(), range_vector.end(), BoundingBoxCompare<BoundingBoxB>() ) )
    std::sort(range_vector.begin(),range_vector.end(),BoundingBoxCompare<BoundingBoxB>());

  // check track of how many times the coarse search needs to exercise the expanstion_factor
  int not_empty_count = 0;

  bool range_vector_empty = stk::is_true_on_all_procs(comm, range_vector.empty());
  while (!range_vector_empty) { // Keep going until all range points are processed.
    // Slightly confusing: coarse_search documentation has domain->range
    // relations sorted by domain key.  We want range->domain type relations
    // sorted on range key. It might appear we have the arguments revered
    // in coarse_search call, but really, this is what we want.
    typename INTERPOLATE::EntityProcRelationVec rng_to_dom;
    search::coarse_search(range_vector, domain_vector, search_method, comm, rng_to_dom);

    // increment how many times we are within the while loop
    not_empty_count++;

    CallOptionalPostCoarseSearchFilter<INTERPOLATE>::call(rng_to_dom, mesha, meshb);
    range_to_domain.insert(range_to_domain.end(), rng_to_dom.begin(), rng_to_dom.end());

    impl::delete_range_points_found<INTERPOLATE>(range_vector, rng_to_dom);

    range_vector_empty = stk::is_true_on_all_procs(comm, range_vector.empty());

    if (expansion_factor <= 1.0) {
        if (!range_vector_empty) {
            size_t range_vector_size = range_vector.size();
            size_t g_range_vector_size = 0;
            stk::all_reduce_max( comm, &range_vector_size, &g_range_vector_size, 1);
            sierra::Env::outputP0() << "GeometricTransfer<INTERPOLATE>::coarse_search(): Number of points not found: "
                                    << g_range_vector_size
                                    << " in initial coarse search" << std::endl;
        }

        break;
    }

    if (!range_vector_empty) {
      for (BoundingBoxB& rangeBox : range_vector) {
        // If points were missed, increase search radius.
        search::scale_by(rangeBox.first, expansion_factor);
      }
      // If points were missed, increase search radius; extract the number of these points and tell the user
      for (BoundingBoxA& domainBox : domain_vector) {
        search::scale_by(domainBox.first, expansion_factor);
      }

      print_expansion_warnings<INTERPOLATE>(comm, not_empty_count, range_vector.size());
    }
  }
  sort (range_to_domain.begin(), range_to_domain.end());
}

template <class T>
typename T::EntityKeyMap::iterator  insert (typename T::EntityKeyMap &map,
                                      const typename T::EntityKeyMap::key_type &k,
                                      const typename T::EntityKeyA &a) {
  const typename T::EntityKeyMap::mapped_type m(a);
  const typename T::EntityKeyMap::value_type  v(k, m);
  const typename T::EntityKeyMap::iterator it = map.insert(v);
  return it;
}


template <class INTERPOLATE> void localize_entity_key_map(typename INTERPOLATE::MeshB & meshb,
    const typename INTERPOLATE::EntityProcRelationVec & global_range_to_domain,
    typename INTERPOLATE::EntityKeyMap & local_range_to_domain)
{

  ParallelMachine comm = meshb.comm();
  const unsigned my_rank = parallel_machine_rank(comm);

  local_range_to_domain.clear();
  for (auto && i : global_range_to_domain)
  {
    const unsigned range_owning_rank = i.first.proc();
    if (range_owning_rank == my_rank) insert<INTERPOLATE>(local_range_to_domain, i.first.id(), i.second.id());
  }
}
}

}
}




#endif /* STK_STK_TRANSFER_STK_TRANSFER_GEOMETRICTRANSFERIMPL_HPP_ */
