/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/util/perf_util.hpp>

#include <stk_io/util/Gmesh_STKmesh_Fixture.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include <sstream>

STKUNIT_UNIT_TEST(many_parts, many_parts)
{
  // vector of mesh-dimensions holds the number of elements in each dimension.
  // Hard-wired to 3. This test can run with spatial-dimension less than 3,
  // (if generated-mesh can do that) but not greater than 3.
  //
  // Doesn't really matter what the mesh size is since we never commit
  std::vector<int> mesh_dims(3, 42);

  std::ostringstream oss;
  oss << mesh_dims[0] << "x" << mesh_dims[1] << "x" << mesh_dims[2];

  stk::io::util::Gmesh_STKmesh_Fixture fixture(MPI_COMM_WORLD, oss.str());
  stk::mesh::MetaData& meta = fixture.getMetaData();

  double start_time = stk::cpu_time();

  stk::mesh::Part& super1 = meta.declare_part("super1");
  stk::mesh::Part& super2 = meta.declare_part("super2");
  stk::mesh::Part& super3 = meta.declare_part("super3");
  stk::mesh::Part& super4 = meta.declare_part("super4");

  stk::mesh::Part& sub1 = meta.declare_part("sub1");
  stk::mesh::Part& sub2 = meta.declare_part("sub2");
  stk::mesh::Part& sub3 = meta.declare_part("sub3");
  stk::mesh::Part& sub4 = meta.declare_part("sub4");

  stk::mesh::PartVector parts;
  unsigned num_parts = 50000;
  for(unsigned i=0; i<num_parts; ++i) {
    std::ostringstream ossname;
    ossname << "part_"<<i;
    stk::mesh::Part& part = meta.declare_part(ossname.str());
    meta.declare_part_subset(super1, part);
    meta.declare_part_subset(super2, part);
    meta.declare_part_subset(super3, part);
    meta.declare_part_subset(part, sub1);
    meta.declare_part_subset(part, sub2);
    meta.declare_part_subset(part, sub3);
    parts.push_back(&part);
  }

  double part_create_time = stk::cpu_time() - start_time;
  start_time = stk::cpu_time();

  typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorField;
  VectorField& field = meta.declare_field<VectorField>(stk::topology::NODE_RANK, "field");
  for(size_t i=0; i<parts.size(); ++i) {
    const stk::mesh::Part& part = *parts[i];
    stk::mesh::put_field(field, part, 3);
  }

  double field_reg_time = stk::cpu_time() - start_time;
  start_time = stk::cpu_time();

  for(size_t i=0; i<parts.size(); ++i) {
    meta.declare_part_subset(super4, *parts[i]);
    meta.declare_part_subset(*parts[i], sub4);
  }

  double part_subset_time = stk::cpu_time() - start_time;

  double total_time    = part_create_time + field_reg_time + part_subset_time;

  static const int NUM_TIMERS = 4;
  const double timers[NUM_TIMERS] = {part_create_time, field_reg_time, part_subset_time, total_time};
  const char* timer_names[NUM_TIMERS] = {"Part create", "Field register", "Part subset", "Total time"};

  stk::print_timers_and_memory(&timer_names[0], &timers[0], NUM_TIMERS);

  stk::parallel_print_time_without_output_and_hwm(MPI_COMM_WORLD, total_time);
}
