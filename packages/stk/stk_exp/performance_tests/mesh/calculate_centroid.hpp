#ifndef calculate_centroid_hpp
#define calculate_centroid_hpp

#include <vector>

namespace sierra {
namespace mesh {
namespace performance_tests {

/**
   num_nodes: number of nodes connected to the element
   elem_node_coords: array of length num_nodes*3, containing the
   coordinates for an element's nodes
   elem_centroid: array of length 3
*/
template<class T>
inline
void
calculate_centroid_3d(
  size_t                num_nodes,
  const T *        elem_node_coords,
  T *              elem_centroid)
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

/**
   num_elements: number of element
   num_nodes: number of nodes connected to the element
   elem_node_coords: array of length num_nodes*3, containing the
   coordinates for an element's nodes
   elem_centroid: array of length 3
*/
template<class T>
inline
void
calculate_centroid_3d(
  size_t                num_elements,
  size_t                num_nodes,
  const T *        elem_node_coords,
  T *              elem_centroid)
{
  //compute the element-centroid:
  for (size_t i = 0; i < num_elements; ++i) {
    for(size_t n = 0; n < num_nodes; ++n) {
      elem_centroid[i*3 + 0] += elem_node_coords[i*num_nodes*3 + n*3 + 0];
      elem_centroid[i*3 + 1] += elem_node_coords[i*num_nodes*3 + n*3 + 1];
      elem_centroid[i*3 + 2] += elem_node_coords[i*num_nodes*3 + n*3 + 2];
    }
    elem_centroid[i*3 + 0] /= num_nodes;
    elem_centroid[i*3 + 1] /= num_nodes;
    elem_centroid[i*3 + 2] /= num_nodes;
  }
}

}
}
}

#endif
