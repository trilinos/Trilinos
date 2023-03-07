#ifndef CONSERVATIVE_TRANSFER_ADAPTOR_H
#define CONSERVATIVE_TRANSFER_ADAPTOR_H

#include "simple_interpolator.hpp"

namespace stk {
namespace middle_mesh {
namespace disc1d {
namespace impl {

struct CSRIndex
{
    int row;
    int colIdx;
};

class ConservativeTransferAdaptor : public InterpolatorBase
{
  public:
  public:
    ConservativeTransferAdaptor(Discretization1DPtr discSrc, Discretization1DPtr discDest,
                                std::shared_ptr<InterpolatorBase> interpolator)
      : m_discSrc(discSrc)
      , m_discDest(discDest)
      , m_interpolator(interpolator)
    {
      correct_coeffs();
    }

    void interpolate(const std::vector<double>& srcFunc, std::vector<double>& destFunc) override
    {
      // TODO: this is suboptimal because applyInterpolationFromMatrix uses getInterpolationMatrix,
      //       which makes an unnecessary copy
      apply_interpolation_from_matrix(*this, srcFunc, destFunc);
    }

    void get_interpolation_matrix(std::vector<std::vector<double>>& coeffs,
                                  std::vector<std::vector<int>>& indices) override
    {
      coeffs  = m_rowCoeffs;
      indices = m_rowIndices;
    }

  private:
    void correct_coeffs();

    void get_column_indices();

    // solves a quadratic program for the special case where G = I and A is a single row
    void solve_quadratic_program(const std::vector<double>& x0, const std::vector<double>& a, double b,
                                 std::vector<double>& x);

    double dot(const std::vector<double>& a, const std::vector<double>& b);

    Discretization1DPtr m_discSrc;
    Discretization1DPtr m_discDest;
    std::shared_ptr<InterpolatorBase> m_interpolator;
    std::vector<std::vector<double>> m_rowCoeffs;
    std::vector<std::vector<int>> m_rowIndices;
    std::vector<std::vector<CSRIndex>> m_columnIndices;

    std::vector<double> m_srcMass;
    std::vector<double> m_destMass;
};

} // namespace impl

} // namespace disc1d
} // namespace middle_mesh
} // namespace stk
#endif