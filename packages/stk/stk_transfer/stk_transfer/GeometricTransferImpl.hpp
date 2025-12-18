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
#include <limits>
#include <vector>
#include <stk_util/parallel/OutputStreams.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include "stk_search/BoundingBox.hpp"
#include "stk_search/CoarseSearch.hpp"
#include "stk_util/parallel/Parallel.hpp"
#include "stk_util/util/ReportHandler.hpp"

namespace stk {
namespace transfer {

namespace impl {

template <class BoundingBoxType>
struct BoundingBoxCompare 
{
  bool operator()(const BoundingBoxType & a, const BoundingBoxType & b) const
  {
    return a.second.id() < b.second.id();
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

// requires domain boxes and domain to range to be sorted
template <class INTERPOLATE>
void get_remaining_domain_points(std::vector<typename INTERPOLATE::MeshB::BoundingBox> &domain_boxes,
    const typename INTERPOLATE::EntityProcRelationVec &this_pass_domain_to_range)
{
  std::vector<typename INTERPOLATE::MeshB::BoundingBox> domain_boxes_missing_a_range_candidate;

  auto domain_it = domain_boxes.begin();
  auto range_it = this_pass_domain_to_range.begin();

  while (domain_it != domain_boxes.end() && range_it != this_pass_domain_to_range.end()) {
    const auto domain_id_a = domain_it->second;
    const auto domain_id_b = range_it->first;

    if (domain_id_a < domain_id_b) {
      domain_boxes_missing_a_range_candidate.push_back(*domain_it);
      ++domain_it;
    } else if (domain_id_a > domain_id_b) {
      ++range_it;
    } else {
      ++domain_it;
      ++range_it;
    }
  }

  // remaining domain were all missed
  while (domain_it != domain_boxes.end()) {
      domain_boxes_missing_a_range_candidate.push_back(*domain_it);
      ++domain_it;
  }

  domain_boxes.swap(domain_boxes_missing_a_range_candidate);
}

template <typename T, typename=int>
struct OptionalPostCoarseSearchFilter
{
  static void call(typename T::EntityProcRelationVec & /*BtoA*/,
      const typename T::MeshA * /*mesha*/,
      const typename T::MeshB * /*meshb*/)
  {
    //Doesn't exist, so don't call
  }
};

template <typename T>
struct OptionalPostCoarseSearchFilter<T, decltype((void) T::post_coarse_search_filter, 0)>
{
  static void call(typename T::EntityProcRelationVec &BtoA,
      const typename T::MeshA * mesha,
      const typename T::MeshB * meshb)
  {
    STK_ThrowRequire(mesha);
    STK_ThrowRequire(meshb);
    T::post_coarse_search_filter(BtoA, *mesha, *meshb);
  }
};

template <typename INTERPOLATE, typename=int>
struct BoxExpansionHandler
{
  using MeshA = typename INTERPOLATE::MeshA;
  using MeshB = typename INTERPOLATE::MeshB;
  using BoundingBoxA = typename INTERPOLATE::MeshA::BoundingBox;
  using BoundingBoxB = typename INTERPOLATE::MeshB::BoundingBox;

  static void call(
    const MeshA* /*mesha*/,
    const MeshB* /*meshb*/,
    std::vector<BoundingBoxA>& range,  // from
    std::vector<BoundingBoxB>& domain, // to
    const double expansion_factor
  )
  {
    for(auto&& box : range) {
      search::scale_by(box.first, expansion_factor);
    }

    for(auto&& box : domain) {
      search::scale_by(box.first, expansion_factor);
    }
  }
};

template <typename INTERPOLATE>
struct BoxExpansionHandler<INTERPOLATE, decltype((void) INTERPOLATE::handle_box_expansions, 0)>
{
  using MeshA = typename INTERPOLATE::MeshA;
  using MeshB = typename INTERPOLATE::MeshB;
  using BoundingBoxA = typename INTERPOLATE::MeshA::BoundingBox;
  using BoundingBoxB = typename INTERPOLATE::MeshB::BoundingBox;

  static void call(
    const MeshA* mesha,
    const MeshB* meshb,
    std::vector<BoundingBoxA>& range,  // from
    std::vector<BoundingBoxB>& domain, // to
    const double expansion_factor
  )
  {
    INTERPOLATE::handle_box_expansions(mesha, meshb, range, domain, expansion_factor);
  }
};

template <typename INTERPOLATE, typename=int>
struct ExpansionWarningPrinter
{
  static void call(
    int num_coarse_search_passes,
    size_t global_domain_size
  )
  {
    stk::outputP0() << "GeometricTransfer<INTERPOLATE>::coarse_search(): Number of points not found: "
                    << global_domain_size << " after expanding bounding boxes " << num_coarse_search_passes
                    << " time(s)" << std::endl;
    stk::outputP0() << "...will now expand the set of candidate bounding boxes and re-attempt the coarse search"
                    << std::endl;
  }
};

template <typename INTERPOLATE>
struct ExpansionWarningPrinter<INTERPOLATE, decltype((void) INTERPOLATE::print_expansion_warnings, 0)>
{
  static void call(
    int num_coarse_search_passes,
    size_t global_domain_size
  )
  {
    INTERPOLATE::print_expansion_warnings(num_coarse_search_passes, global_domain_size);
  }
};


template <typename T>
auto call_copy_entities(T & mesh, typename T::EntityProcVec & entities_to_copy, const std::string & name)
-> decltype( (void) mesh.copy_entities(entities_to_copy, name))
{
  mesh.copy_entities(entities_to_copy, name);
}

template <typename... T>
void call_copy_entities(T & ... /*t*/)
{
  //copy entities doesn't exist, so NO-OP
}



template <class INTERPOLATE>
void coarse_search_impl(typename INTERPOLATE::EntityProcRelationVec &domain_to_range,
    stk::ParallelMachine comm,
    typename INTERPOLATE::MeshA *mesha,
    typename INTERPOLATE::MeshB *meshb,
    const stk::search::SearchMethod search_method,
    const double expansion_factor)
{
  using BoundingBoxA = typename INTERPOLATE::MeshA::BoundingBox;
  using BoundingBoxB = typename INTERPOLATE::MeshB::BoundingBox;

  using domain_range_pairs_t = typename INTERPOLATE::EntityProcRelationVec;

  domain_to_range.clear();

  std::vector<BoundingBoxB> domain;
  std::vector<BoundingBoxA> range;

  if (mesha) mesha->bounding_boxes(range);

  if (stk::is_true_on_all_procs(comm, range.empty())) {
    return;
  }

  if (meshb) meshb->bounding_boxes(domain);

  if (!local_is_sorted(range.begin(), range.end(), BoundingBoxCompare<BoundingBoxA>())) {
    std::sort(range.begin(), range.end(), BoundingBoxCompare<BoundingBoxA>());
  }

  if (!local_is_sorted(domain.begin(), domain.end(), BoundingBoxCompare<BoundingBoxB>())) {
    std::sort(domain.begin(), domain.end(), BoundingBoxCompare<BoundingBoxB>());
  }

  // coarse search continues until all domain points are mapped to range candidates
  // TODO where is the safety termination condition preventing infinite loop?
  int num_coarse_search_passes = 0;
  bool domain_empty = stk::is_true_on_all_procs(comm, domain.empty());
  while (!domain_empty) {
    domain_range_pairs_t this_pass_domain_to_range;

    search::coarse_search(domain, range, search_method, comm, this_pass_domain_to_range, true, true, true);

    num_coarse_search_passes++;
    OptionalPostCoarseSearchFilter<INTERPOLATE>::call(this_pass_domain_to_range, mesha, meshb);
    impl::get_remaining_domain_points<INTERPOLATE>(domain, this_pass_domain_to_range);

    domain_empty = stk::is_true_on_all_procs(comm, domain.empty());
    const bool terminate_on_first_pass = (domain_empty || expansion_factor <= 1.0) && num_coarse_search_passes == 1;

    size_t global_domain_size = 0;
    size_t domain_size = domain.size();

    if (stk::util::get_common_coupling_version() >= 16) {
      if (!domain_empty){
        stk::all_reduce_sum(comm, &domain_size, &global_domain_size, 1);
      }
    }

    if (terminate_on_first_pass) {
      domain_to_range.swap(this_pass_domain_to_range);

      if ( !domain_empty ){
        if (stk::util::get_common_coupling_version() < 16) {
          stk::all_reduce_max(comm, &domain_size, &global_domain_size, 1);
        }
        stk::outputP0() << "GeometricTransfer<INTERPOLATE>::coarse_search(): Number of points not found: "
                        << global_domain_size << " in initial coarse search" << std::endl;
      }

      break;
    }

    // otherwise extend current domain to range, expand domain / range boxes, print warnings and repeat coarse search
    domain_to_range.insert(domain_to_range.end(), this_pass_domain_to_range.begin(), this_pass_domain_to_range.end());

    if (!domain_empty) {
      BoxExpansionHandler<INTERPOLATE>::call(
        mesha, meshb, range, domain, expansion_factor
      );

      // preserve MPI pattern from old throw_error_on_coarse_expansion_limit
      if (stk::util::get_common_coupling_version() == 15) {

        std::string local_err_msg;
        std::ostringstream global_err_msg;
        stk::all_write_string(comm, global_err_msg, local_err_msg);

        int err = 0;
        int global_err = 0;
        stk::all_reduce_sum(comm, &err, &global_err, 1);

        STK_ThrowRequireMsg(!global_err, global_err_msg.str());
      }

      // preserve MPI pattern from old print_expansion_warnings
      if (stk::util::get_common_coupling_version() < 16) {
        stk::all_reduce_max(comm, &domain_size, &global_domain_size, 1);
      }

      ExpansionWarningPrinter<INTERPOLATE>::call(
        num_coarse_search_passes, global_domain_size
      );

      if (stk::util::get_common_coupling_version() >= 16) {
        // Custom expansion handlers may remove items from the domain
        domain_empty = stk::is_true_on_all_procs(comm, domain.empty());
      }
    }
  }
  sort(domain_to_range.begin(), domain_to_range.end());
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


template <typename T>
class HasNeedRepeatSearch
{
  private:
    using yes = std::array<int, 1>;
    using no  = std::array<int, 2>;

    template <typename C>
    static yes test(decltype(&C::need_repeat_search));

    template <typename C>
    static no test(...);

  public:
    static constexpr bool value = sizeof(test<T>(0)) == sizeof(yes);
};

template <typename Mesh>
constexpr bool has_need_repeat_search()
{
  return HasNeedRepeatSearch<Mesh>::value;
}

template <typename Mesh>
bool need_repeat_search(Mesh& mesh)
{
  if constexpr (has_need_repeat_search<Mesh>())
  {
    return mesh.need_repeat_search();
  } else
  {
    return false;
  }
}
}



}
}

#endif /* STK_STK_TRANSFER_STK_TRANSFER_GEOMETRICTRANSFERIMPL_HPP_ */
