#ifndef STENCIL_INTERPOLATOR_H
#define STENCIL_INTERPOLATOR_H

#include "simple_interpolator.hpp"
#include <vector>

namespace stk {
namespace middle_mesh {
namespace disc1d {
namespace impl {

class StencilInterpolator : public InterpolatorBase
{
  public:
    StencilInterpolator(Discretization1DPtr discSrc, Discretization1DPtr discDest, int stencilsize)
      : m_discSrc(discSrc)
      , m_discDest(discDest)
      , m_search(discSrc)
      , m_stencilsize(stencilsize)
    {}

    void interpolate(const std::vector<double>& srcFunc, std::vector<double>& destFunc) override;

    void get_interpolation_matrix(std::vector<std::vector<double>>& coeffs,
                                  std::vector<std::vector<int>>& indices) override;

  private:
    double interpolate_point(const std::vector<double>& srcFunc, double x);

    void get_interpolation_data(double x, std::vector<double>& coeffs, std::vector<int>& indices);

    std::vector<mesh::MeshEntityPtr> get_mesh_entities_in_stencil(double x);

    std::vector<double> get_node_coords(const std::vector<mesh::MeshEntityPtr>& elements);

    std::vector<int> get_node_dofs(const std::vector<mesh::MeshEntityPtr>& elements);

    std::vector<double> compute_coefficients(const std::vector<double>& nodeCoords, double x);

    Discretization1DPtr m_discSrc;
    Discretization1DPtr m_discDest;
    Search1D m_search;
    int m_stencilsize;
};
} // namespace impl

} // namespace disc1d
} // namespace middle_mesh
} // namespace stk
#endif
