#ifndef KOKKOS_GLOBAL_TO_LOCAL_IDS_HPP
#define KOKKOS_GLOBAL_TO_LOCAL_IDS_HPP

#include <Kokkos_Core.hpp>

#include <Kokkos_UnorderedMap.hpp>

#include <vector>
#include <algorithm>
#include <iomanip>

#include <impl/Kokkos_Timer.hpp>

// This test will simulate global ids

namespace G2L {

static const unsigned begin_id_size = 256u;
static const unsigned end_id_size = 1u << 25;
static const unsigned id_step = 2u;

//use to help generate global ids
union helper
{
  uint32_t word;
  uint8_t byte[4];
};


//generate a unique global id from the local id
template <typename Device>
struct generate_ids
{
  typedef Device device_type;
  typedef typename device_type::size_type size_type;
  typedef Kokkos::View<uint32_t*,device_type> local_id_view;

  local_id_view local_2_global;

  generate_ids( local_id_view & ids)
    : local_2_global(ids)
  {
    Kokkos::parallel_for(local_2_global.size(), *this);
  }


  KOKKOS_INLINE_FUNCTION
  void operator()(size_type i) const
  {

    helper x = {static_cast<uint32_t>(i)};

    // shuffle the bytes of i to create a unique, semi-random global_id
    x.word = ~x.word;

    uint8_t tmp = x.byte[3];
    x.byte[3] = x.byte[1];
    x.byte[1] = tmp;

    tmp = x.byte[2];
    x.byte[2] = x.byte[0];
    x.byte[0] = tmp;

    local_2_global[i] = x.word;
  }

};

// fill a map of global_id -> local_id
template <typename Device>
struct fill_map
{
  typedef Device device_type;
  typedef typename device_type::size_type size_type;
  typedef Kokkos::View<const uint32_t*,device_type, Kokkos::MemoryRandomAccess> local_id_view;
  typedef Kokkos::UnorderedMap<uint32_t,size_type,device_type> global_id_view;

  global_id_view global_2_local;
  local_id_view local_2_global;

  fill_map( global_id_view gIds, local_id_view lIds)
    : global_2_local(gIds) , local_2_global(lIds)
  {
    Kokkos::parallel_for(local_2_global.size(), *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(size_type i) const
  {
    global_2_local.insert( local_2_global[i], i);
  }

};

// check that the global id is found and that it maps to the local id
template <typename Device>
struct find_test
{
  typedef Device device_type;
  typedef typename device_type::size_type size_type;
  typedef Kokkos::View<const uint32_t*,device_type, Kokkos::MemoryRandomAccess> local_id_view;
  typedef Kokkos::UnorderedMap<const uint32_t, const size_type,device_type> global_id_view;

  global_id_view global_2_local;
  local_id_view local_2_global;

  typedef size_t value_type;

  find_test( global_id_view gIds, local_id_view lIds, value_type & num_errors)
    : global_2_local(gIds) , local_2_global(lIds)
  {
    Kokkos::parallel_reduce(local_2_global.size(), *this, num_errors);
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type & v) const
  { v = 0; }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type & dst, volatile value_type const & src) const
  { dst += src; }

  KOKKOS_INLINE_FUNCTION
  void operator()(size_type i, value_type & num_errors) const
  {
    uint32_t index = global_2_local.find( local_2_global[i] );

    if (  !global_2_local.valid_at(index)
        || global_2_local.key_at(index) != local_2_global[i]
        || global_2_local.value_at(index) != i)
      ++num_errors;
  }

};

// run test
template <typename Device>
size_t test_global_to_local_ids(unsigned num_ids, unsigned capacity, unsigned num_find_iterations)
{

  typedef Device device_type;
  typedef typename device_type::size_type size_type;

  typedef Kokkos::View<uint32_t*,device_type> local_id_view;
  typedef Kokkos::UnorderedMap<uint32_t,size_type,device_type> global_id_view;

  double elasped_time = 0;
  Kokkos::Impl::Timer timer;

  local_id_view local_2_global("local_ids", num_ids);
  global_id_view global_2_local(capacity);

  int shiftw = 15;

  //create
  elasped_time = timer.seconds();
  std::cout << std::setw(shiftw) <<  "allocate: " <<  elasped_time << std::endl;
  timer.reset();

  // generate unique ids
  {
    generate_ids<Device> gen(local_2_global);
  }

  // generate
  elasped_time = timer.seconds();
  std::cout << std::setw(shiftw) << "generate: " <<  elasped_time << std::endl;
  timer.reset();

  {
    fill_map<Device> fill(global_2_local, local_2_global);
  }

  // fill
  elasped_time = timer.seconds();
  std::cout << std::setw(shiftw) << "fill: " <<  elasped_time << std::endl;
  timer.reset();


  size_t num_errors = global_2_local.failed_insert();

  if (num_errors == 0u) {
    for (unsigned i=0; i<num_find_iterations; ++i)
    {
      find_test<Device> find(global_2_local, local_2_global,num_errors);
    }

    // find
    elasped_time = timer.seconds();
    std::cout << std::setw(shiftw) << "lookup: " <<  elasped_time << std::endl;
  }
  else {
    std::cout << "    !!! Fill Failed !!!" << std::endl;
  }

  return num_errors;
}

template <typename Device>
size_t run_test(unsigned num_ids, unsigned num_find_iterations)
{
  // expect to fail
  unsigned capacity = (num_ids*2u)/3u;
  std::cout << " 66% of needed capacity (should fail)" << std::endl;
  test_global_to_local_ids<Device>(num_ids, capacity, num_find_iterations);

  //should not fail
  std::cout << " 100% of needed capacity" << std::endl;
  capacity = num_ids;
  size_t num_errors = test_global_to_local_ids<Device>(num_ids, capacity, num_find_iterations);

  //should not fail
  std::cout << " 150% of needed capacity" << std::endl;
  capacity = (num_ids*3u)/2u;
  num_errors += test_global_to_local_ids<Device>(num_ids, capacity, num_find_iterations);

  return num_errors;
}


} // namespace G2L


#endif //KOKKOS_GLOBAL_TO_LOCAL_IDS_HPP

