#ifndef SIMPLE_INTERPOLATOR_H
#define SIMPLE_INTERPOLATOR_H

#include "discretization_1d.hpp"
#include "search_1d.hpp"

namespace stk {
namespace middle_mesh {
namespace disc1d {
namespace impl {

class InterpolatorBase
{
  public:
    virtual ~InterpolatorBase(){};

    virtual void interpolate(const std::vector<double>& srcFunc, std::vector<double>& destFunc) = 0;

    virtual void get_interpolation_matrix(std::vector<std::vector<double>>& coeffs,
                                          std::vector<std::vector<int>>& indices) = 0;
};

// use getInterpolationMatrix to do the interpolation
// Prefer the InterpolatorBase::interpolation() method, this method is here mostly
// for testing
void apply_interpolation_from_matrix(InterpolatorBase& interp, const std::vector<double>& srcFunc,
                                     std::vector<double>& destFunc);

class SimpleInterpolator : public InterpolatorBase
{
  public:
    SimpleInterpolator(Discretization1DPtr discSrc, Discretization1DPtr discDest)
      : m_discSrc(discSrc)
      , m_discDest(discDest)
      , m_search(discSrc)
    {}

    void interpolate(const std::vector<double>& srcFunc, std::vector<double>& destFunc) override;

    void get_interpolation_matrix(std::vector<std::vector<double>>& coeffs,
                                  std::vector<std::vector<int>>& indices) override;

  private:
    double interpolate_point(const std::vector<double>& srcFunc, double x);

    void get_interpolation_data(double x, std::vector<double>& coeffs, std::vector<int>& indices);

    Discretization1DPtr m_discSrc;
    Discretization1DPtr m_discDest;
    Search1D m_search;
};

} // namespace impl

} // namespace disc1d
} // namespace middle_mesh
} // namespace stk
#endif