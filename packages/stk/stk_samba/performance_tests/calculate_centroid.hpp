#ifndef calculate_centroid_hpp
#define calculate_centroid_hpp

#include <vector>

namespace performance_tests {

/**
   num_nodes: number of nodes connected to the element
   elem_node_coords: array of length num_nodes*3, containing the
   coordinates for an element's nodes
   elem_centroid: array of length 3
*/
inline void
calculate_centroid_3d(
  size_t                num_nodes,
  const double *        elem_node_coords,
  double *              elem_centroid)
{
  //compute the element-centroid:
  for(size_t n = 0; n < num_nodes; ++n) {
    elem_centroid[0] += elem_node_coords[n*3 + 0];
    elem_centroid[1] += elem_node_coords[n*3 + 1];
    elem_centroid[2] += elem_node_coords[n*3 + 2];
  }
  elem_centroid[0] /= num_nodes;
  elem_centroid[1] /= num_nodes;
  elem_centroid[2] /= num_nodes;
}

inline void
calculate_centroid_hex_8(
  const double *        elem_node_coords,
  double *              elem_centroid)
{
  const size_t num_nodes = 8;
  //compute the element-centroid:
  for(size_t n = 0; n < num_nodes; ++n) {
    elem_centroid[0] += elem_node_coords[n*3 + 0];
    elem_centroid[1] += elem_node_coords[n*3 + 1];
    elem_centroid[2] += elem_node_coords[n*3 + 2];
  }
  elem_centroid[0] /= num_nodes;
  elem_centroid[1] /= num_nodes;
  elem_centroid[2] /= num_nodes;
}
} //performance_tests

#endif //calculate_centroid_hpp

