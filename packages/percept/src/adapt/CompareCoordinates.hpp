#ifndef adapt_CompareCoordinates_hpp
#define adapt_CompareCoordinates_hpp

#include <array>
#include <numeric>

inline bool float_less(double a, double b)
{
  return static_cast<float>(a) < static_cast<float>(b);
}

inline bool is_less_than_in_x_then_y_then_z(const double * A, const double * B)
{
    if (float_less(A[0], B[0])) return true;
    if (float_less(B[0], A[0])) return false;
    if (float_less(A[1], B[1])) return true;
    if (float_less(B[1], A[1])) return false;
    if (float_less(A[2], B[2])) return true;
    if (float_less(B[2], A[2])) return false;

    //  How do you want to handle a tie?
    return false;
}

template <size_t SIZE>
std::array<int,SIZE> get_rank_of_nodes_based_on_coordinates(std::array<double *,SIZE> & node_coords)
{
  // initialize original index locations
  std::array<size_t,SIZE> index;
  std::iota(index.begin(), index.end(), 0);

  std::stable_sort(index.begin(), index.end(),
       [&node_coords](size_t i1, size_t i2) {return is_less_than_in_x_then_y_then_z(node_coords[i1], node_coords[i2]);});

  // invert indices to get node rank by coordinates
  std::array<int,SIZE> rank;
  for (size_t i=0; i < SIZE; ++i)
    rank[index[i]] = i;
  return rank;
}
    
#endif
