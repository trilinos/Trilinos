#ifndef SEARCH_1D_H
#define SEARCH_1D_H

#include "discretization_1d.hpp"

namespace stk {
namespace middle_mesh {
namespace disc1d {
namespace impl {

class Search1D
{
  public:
    explicit Search1D(Discretization1DPtr disc)
      : m_disc(disc)
      , m_vertCoords(count_valid(m_disc->get_mesh()->get_vertices()))
    {
      setup_for_search();
    }

    mesh::MeshEntityPtr find_containing_element(double x);

    // given a point x that lies within el, compute its parametric coordinate
    double get_xi_coordinate(mesh::MeshEntityPtr el, double x);

  private:
    void setup_for_search();

    int clamp(int min, int val, int max);

    Discretization1DPtr m_disc;
    std::vector<double> m_vertCoords;
};

} // namespace impl

} // namespace disc1d
} // namespace middle_mesh
} // namespace stk
#endif
