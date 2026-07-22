#ifndef KRINO_KRINO_UNIT_TESTS_AKRI_UNIT_CREATEFACETEDSPHERE_HPP_
#define KRINO_KRINO_UNIT_TESTS_AKRI_UNIT_CREATEFACETEDSPHERE_HPP_

#include <stk_math/StkVector.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <vector>
#include <string>

namespace krino {

void fill_sphere_vertices_and_connectivity(const double radius,
  const double meshSize,
  std::vector<stk::math::Vector3d> & vertices,
  std::vector<std::array<unsigned,3>> & facetConnectivity);


void write_stl_for_sphere(const std::string &filename,
  const double radius,
  const double meshSize,
  const stk::ParallelMachine comm);

}

#endif /* KRINO_KRINO_UNIT_TESTS_AKRI_UNIT_CREATEFACETEDSPHERE_HPP_ */
