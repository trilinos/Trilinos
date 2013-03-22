/*--------------------------------------------------------------------*/
/*    Copyright 2009, 2011 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>

#include <math.h>

#include <stk_mesh/base/FieldParallel.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_percept/PerceptMesh.hpp>

#include <stk_percept/Util.hpp>
#include <stk_percept/ExceptionWatch.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <stk_percept/fixtures/Fixture.hpp>
#include <stk_percept/fixtures/QuadFixture.hpp>
#include <stk_util/util/MallocUsed.h>


namespace stk
{
  namespace adapt
  {
    namespace regression_tests
    {


#define DO_MEMORY_ACCOUNTING 1

typedef uint64_t MemorySizeType;

struct MemoryInfo
{
  MemorySizeType m_malloc_used;
  MemorySizeType m_malloc_footprint;
  MemorySizeType m_malloc_max_footprint;
  static const MemorySizeType MB = 1024*1024;

  MemoryInfo() { get_memory_usage(); }

  void get_memory_usage() 
  {
#if DO_MEMORY_ACCOUNTING
#if !defined(SIERRA_PTMALLOC3_ALLOCATOR) && !defined(SIERRA_PTMALLOC2_ALLOCATOR)
          std::cout << "WARNING: ptmalloc2|3 not compiled in so malloc_used info unavailable.  Recompile with e.g. 'bake allocator=ptmalloc2 (or 3)'.  Printing zeros..." << std::endl;
#else

#endif
    m_malloc_used = malloc_used();
    m_malloc_footprint = malloc_footprint();
    m_malloc_max_footprint = malloc_max_footprint();
#else
    m_malloc_used = 0;
    m_malloc_footprint = 0;
    m_malloc_max_footprint = 0;
#endif
  }
  void set_state() { get_memory_usage(); }
  void get_increment() {
    MemoryInfo old_state = *this;
    get_memory_usage();
    m_malloc_used -= old_state.m_malloc_used;
    m_malloc_footprint -= old_state.m_malloc_footprint;
  }

};

inline double MegaByte(MemorySizeType x) { return  ((double)x/1024.0/1024.0); }

std::ostream& operator<<(std::ostream& os, const MemoryInfo& mem)
{
  char buf[1024];
  sprintf(buf, "\n%20s %20s %20s\n%20g %20g %20g\n", "used [MB]", "footprint [MB]", "max_footprint [MB]",
          MegaByte(mem.m_malloc_used), MegaByte(mem.m_malloc_footprint), MegaByte(mem.m_malloc_max_footprint) );
  os << buf;
  return os;
}


STKUNIT_UNIT_TEST(adapt, count_memory)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD ;

  const unsigned p_size = stk::parallel_machine_size( pm );
  if (p_size == 1)
    {
      const unsigned n = 20;
      //const unsigned nx = n , ny = n , nz = p_size*n ;
      const unsigned nx = n , ny = n;

      percept::QuadFixture<double, shards::Triangle<3> > fixture( pm , nx , ny, false);
      fixture.meta_data.commit();
      fixture.generate_mesh();

      percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);
      //eMesh.print_info("quad mesh",2);

      // see stk_samba/perf_test_tri_refine.cpp
      //const size_t num_new_tris = 2000*2000;
      const size_t num_new_tris = 20*20;
      const size_t num_nodes_per_tri = 3;
      const size_t num_new_nodes = num_new_tris*num_nodes_per_tri;
      MemoryInfo mem_delta_node;
      double time = -stk::percept::Util::cpu_time();

      std::vector<stk::mesh::Entity *> new_nodes, new_elements;

      eMesh.get_bulk_data()->modification_begin();
      eMesh.createEntities(eMesh.node_rank(), num_new_nodes, new_nodes);
      eMesh.get_bulk_data()->modification_end();

      mem_delta_node.get_increment();
      double mem_per_node = double(mem_delta_node.m_malloc_used)/double(num_new_nodes);
      std::cout << "\nstk_mesh count_memory mem_per_node = " << mem_per_node << "\n" << std::endl;

      MemoryInfo mem_delta_elem_0, mem_delta_elem_1;

      eMesh.get_bulk_data()->modification_begin();
      eMesh.createEntities(eMesh.element_rank(), num_new_tris, new_elements);
      eMesh.get_bulk_data()->modification_end();

      mem_delta_elem_0.get_increment();

      eMesh.get_bulk_data()->modification_begin();
      size_t i_node=0;
      for (size_t i=0; i<num_new_tris; ++i) {
        unsigned ordinal = 0;
        for (size_t j=0; j < num_nodes_per_tri; ++j) {
          eMesh.get_bulk_data()->declare_relation(*new_elements[i],*new_nodes[i_node],ordinal);
          ++ordinal;
          ++i_node;
        }
      }
      eMesh.get_bulk_data()->modification_end();

      mem_delta_elem_1.get_increment();
      double mem_per_elem_0 = double(mem_delta_elem_0.m_malloc_used)/double(num_new_tris);
      double mem_per_elem_1 = double(mem_delta_elem_1.m_malloc_used)/double(num_new_tris);

      time += stk::percept::Util::cpu_time();

      std::cout << "\nstk_mesh count_memory mem_per_elem (no connectivity) = " << mem_per_elem_0 << " with connectivity= " << mem_per_elem_1 << " cpu= " << time << std::endl;

    }
}


    }
  }
}
