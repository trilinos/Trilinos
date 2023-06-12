#include "stk_middle_mesh/conservative_transfer_adaptor.hpp"
#include "stk_middle_mesh/simple_interpolator.hpp"
#include "stk_middle_mesh/stencil_interpolator.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {

class SimpleInterpolatorTester : public ::testing::Test
{
  protected:
    template <typename T>
    void initialize(double xmin, double xmax, int numelSrc, int numelDest, int basisDegreeSrc, int basisDegreeDest,
                    int quadDegree, T func)
    {
      disc1d::impl::MeshSpec1D specSrc{numelSrc, xmin, xmax}, specDest{numelDest, xmin, xmax};

      mDiscSrc  = disc1d::impl::make_disc1_d(specSrc, basisDegreeSrc, quadDegree);
      mDiscDest = disc1d::impl::make_disc1_d(specDest, basisDegreeDest, quadDegree);

      mInterpolator = std::make_shared<disc1d::impl::SimpleInterpolator>(mDiscSrc, mDiscDest);

      mFSrc    = disc1d::impl::make_grid_function(mDiscSrc, func);
      mFDestEx = disc1d::impl::make_grid_function(mDiscDest, func);
    }

    template <typename T>
    void initialize(double xmin, double xmax, int numelSrc, int numelDest, int basisDegree, int quadDegree, T func)
    {
      initialize(xmin, xmax, numelSrc, numelDest, basisDegree, basisDegree, quadDegree, func);
    }

    disc1d::impl::Discretization1DPtr mDiscSrc;
    disc1d::impl::Discretization1DPtr mDiscDest;
    std::shared_ptr<disc1d::impl::SimpleInterpolator> mInterpolator;
    std::vector<double> mFSrc;
    std::vector<double> mFDestEx;
};

class StencilInterpolatorTester : public ::testing::Test
{
  protected:
    template <typename T>
    void initialize(double xmin, double xmax, int numelSrc, int numelDest, int basisDegreeSrc, int basisDegreeDest,
                    int quadDegree, int stencilsize, T func)
    {
      disc1d::impl::MeshSpec1D specSrc{numelSrc, xmin, xmax}, specDest{numelDest, xmin, xmax};

      mDiscSrc  = disc1d::impl::make_disc1_d(specSrc, basisDegreeSrc, quadDegree);
      mDiscDest = disc1d::impl::make_disc1_d(specDest, basisDegreeDest, quadDegree);

      mInterpolator = std::make_shared<disc1d::impl::StencilInterpolator>(mDiscSrc, mDiscDest, stencilsize);

      mFSrc    = disc1d::impl::make_grid_function(mDiscSrc, func);
      mFDestEx = disc1d::impl::make_grid_function(mDiscDest, func);
    }

    template <typename T>
    void initialize(double xmin, double xmax, int numelSrc, int numelDest, int basisDegree, int quadDegree,
                    int stencilsize, T func)
    {
      initialize(xmin, xmax, numelSrc, numelDest, basisDegree, basisDegree, quadDegree, stencilsize, func);
    }

    disc1d::impl::Discretization1DPtr mDiscSrc;
    disc1d::impl::Discretization1DPtr mDiscDest;
    std::shared_ptr<disc1d::impl::StencilInterpolator> mInterpolator;
    std::vector<double> mFSrc;
    std::vector<double> mFDestEx;
};

class ConvergenceRateCalculator
{
  public:
    void add_measurement(double h, double error)
    {
      m_sizes.push_back(h);
      m_errors.push_back(error);
    }

    void output_data()
    {
      std::cout << "mesh  h            error       ratio     slope" << std::endl;
      for (size_t i = 0; i < m_errors.size(); ++i)
      {
        // double ratio = i == 0 ? 0 : m_errors[i-1]/m_errors[i];
        // double slope = i == 0 ? 0 : (std::log(m_errors[i]) - std::log(m_errors[i-1]))/(std::log(m_sizes[i]) -
        // std::log(m_sizes[i-1]));
        std::printf("%d     %6.5e  %6.5e %9.5f %6.5f\n", int(i), m_sizes[i], m_errors[i], get_ratio(i), get_slope(i));
      }
    }

    double get_ratio(int mesh)
    {
      if (mesh == 0)
        return 0;
      else
        return m_errors[mesh - 1] / m_errors[mesh];
    }

    double get_slope(int mesh)
    {
      if (mesh == 0)
        return 0;
      else
        return (std::log(m_errors[mesh]) - std::log(m_errors[mesh - 1])) /
               (std::log(m_sizes[mesh]) - std::log(m_sizes[mesh - 1]));
    }

  private:
    std::vector<double> m_errors;
    std::vector<double> m_sizes;
};

double compute_error_norm(disc1d::impl::Discretization1DPtr disc, const std::vector<double>& fEx,
                          const std::vector<double>& f)
{
  std::vector<double> error(disc->get_num_dofs());
  for (int i = 0; i < disc->get_num_dofs(); ++i)
    error[i] = f[i] - fEx[i];

  auto mError = apply_mass_matrix(disc, error);

  double errorNorm = 0;
  for (int i = 0; i < disc->get_num_dofs(); ++i)
    errorNorm += error[i] * error[i];

  errorNorm = std::sqrt(errorNorm);

  return errorNorm;
}

} // namespace

TEST_F(SimpleInterpolatorTester, Polynomial)
{
  double xmin = 0, xmax = 2;
  int numelSrc = 5, numelDest = 11;
  int basisDegree = 2, quadDegree = 3;

  for (int polyDegree = 0; polyDegree <= 2; ++polyDegree)
  {
    auto f = [polyDegree](double x) { return std::pow(x, polyDegree); };
    initialize(xmin, xmax, numelSrc, numelDest, basisDegree, quadDegree, f);

    std::vector<double> fDest(mDiscDest->get_num_dofs(), 42);
    mInterpolator->interpolate(mFSrc, fDest);

    std::vector<double> fDestMatrix(mDiscDest->get_num_dofs(), 42);
    disc1d::impl::apply_interpolation_from_matrix(*mInterpolator, mFSrc, fDestMatrix);

    for (int i = 0; i < mDiscDest->get_num_dofs(); ++i)
    {
      EXPECT_NEAR(fDest[i], mFDestEx[i], 1e-13);
      EXPECT_NEAR(fDestMatrix[i], mFDestEx[i], 1e-13);
    }
  }
}

TEST_F(StencilInterpolatorTester, Polynomial)
{
  double xmin = 0, xmax = 2;
  int numelSrc = 5, numelDest = 11;
  int basisDegree = 1, quadDegree = 3;
  int stencilsize = 4;

  for (int polyDegree = 0; polyDegree <= 3; ++polyDegree)
  {
    auto f = [polyDegree](double x) { return std::pow(x, polyDegree); };
    initialize(xmin, xmax, numelSrc, numelDest, basisDegree, quadDegree, stencilsize, f);

    std::vector<double> fDest(mDiscDest->get_num_dofs(), 42);
    mInterpolator->interpolate(mFSrc, fDest);

    std::vector<double> fDestMatrix(mDiscDest->get_num_dofs(), 42);
    disc1d::impl::apply_interpolation_from_matrix(*mInterpolator, mFSrc, fDestMatrix);

    for (int i = 0; i < mDiscDest->get_num_dofs(); ++i)
    {
      EXPECT_NEAR(fDest[i], mFDestEx[i], 1e-13);
      EXPECT_NEAR(fDestMatrix[i], mFDestEx[i], 1e-13);
    }
  }
}

TEST_F(SimpleInterpolatorTester, ConvergenceRate)
{
  double xmin = 0, xmax = 2;
  int numelSrc = 5, numelDest = 11;
  int basisDegree = 2, quadDegree = 3;
  int nmeshes = 6;

  ConvergenceRateCalculator errorSlope, conservationSlope;
  for (int mesh = 0; mesh < nmeshes; ++mesh)
  {
    // std::cout << "mesh " << mesh << std::endl;
    double mean  = (xmax + xmin) / 2;
    double sigma = (xmax - xmin) / 4;
    double pi    = std::atan(1) * 4;

    auto f = [&](double x) {
      return x * std::exp(-(x - mean) * (x - mean) / (2 * sigma * sigma)) / (sigma * std::sqrt(2 * pi));
    };
    initialize(xmin, xmax, numelSrc, numelDest, basisDegree, quadDegree, f);

    std::vector<double> fDest(mDiscDest->get_num_dofs(), 42), error(mDiscDest->get_num_dofs());
    mInterpolator->interpolate(mFSrc, fDest);

    double errorNorm         = compute_error_norm(mDiscDest, mFDestEx, fDest);
    double integralSrc       = disc1d::impl::integrate_function(mDiscSrc, mFSrc);
    double integralDest      = disc1d::impl::integrate_function(mDiscDest, fDest);
    double conservationError = std::abs(integralSrc - integralDest);

    errorSlope.add_measurement(1.0 / numelSrc, errorNorm);
    conservationSlope.add_measurement(1.0 / numelSrc, conservationError);

    numelSrc *= 2;
    numelDest *= 2;
  }

  EXPECT_NEAR(errorSlope.get_slope(nmeshes - 1), 2.5, 0.1);
  EXPECT_NEAR(conservationSlope.get_slope(nmeshes - 1), 4.02, 0.1);

  std::cout << "\nsolution error:" << std::endl;
  errorSlope.output_data();

  std::cout << "\nconservation error:" << std::endl;
  conservationSlope.output_data();
}

TEST_F(StencilInterpolatorTester, ConvergenceRate)
{
  double xmin = 0, xmax = 2;
  int numelSrc = 5, numelDest = 11;
  int basisDegree = 2, quadDegree = 3;
  int stencilsize = 4;
  int nmeshes     = 6;

  ConvergenceRateCalculator errorSlope, conservationSlope;
  for (int mesh = 0; mesh < nmeshes; ++mesh)
  {
    // std::cout << "mesh " << mesh << std::endl;
    double mean  = (xmax + xmin) / 2;
    double sigma = (xmax - xmin) / 4;
    double pi    = std::atan(1) * 4;

    auto f = [&](double x) {
      return x * std::exp(-(x - mean) * (x - mean) / (2 * sigma * sigma)) / (sigma * std::sqrt(2 * pi));
    };
    initialize(xmin, xmax, numelSrc, numelDest, basisDegree, quadDegree, stencilsize, f);

    std::vector<double> fDest(mDiscDest->get_num_dofs(), 42), error(mDiscDest->get_num_dofs());
    mInterpolator->interpolate(mFSrc, fDest);

    double errorNorm         = compute_error_norm(mDiscDest, mFDestEx, fDest);
    double integralSrc       = disc1d::impl::integrate_function(mDiscSrc, mFSrc);
    double integralDest      = disc1d::impl::integrate_function(mDiscDest, fDest);
    double conservationError = std::abs(integralSrc - integralDest);

    errorSlope.add_measurement(1.0 / numelSrc, errorNorm);
    conservationSlope.add_measurement(1.0 / numelSrc, conservationError);

    numelSrc *= 2;
    numelDest *= 2;
  }

  EXPECT_NEAR(errorSlope.get_slope(nmeshes - 1), 4.5, 0.1);
  EXPECT_NEAR(conservationSlope.get_slope(nmeshes - 1), 3.9, 0.1);

  std::cout << "\nsolution error:" << std::endl;
  errorSlope.output_data();

  std::cout << "\nconservation error:" << std::endl;
  conservationSlope.output_data();
}

TEST_F(StencilInterpolatorTester, ConservativeConvergenceRate)
{
  double xmin = 0, xmax = 2;
  int numelSrc = 5, numelDest = 11;
  int basisDegreeSrc = 2, basisDegreeDest = 2;
  int quadDegree  = 5;
  int stencilsize = 4;
  int nmeshes     = 6;

  ConvergenceRateCalculator errorSlope, conservationSlope;
  for (int mesh = 0; mesh < nmeshes; ++mesh)
  {
    // std::cout << "mesh " << mesh << std::endl;
    double mean  = (xmax + xmin) / 2;
    double sigma = (xmax - xmin) / 4;
    double pi    = std::atan(1) * 4;

    auto f = [&](double x) {
      return x * std::exp(-(x - mean) * (x - mean) / (2 * sigma * sigma)) / (sigma * std::sqrt(2 * pi));
    };
    initialize(xmin, xmax, numelSrc, numelDest, basisDegreeSrc, basisDegreeDest, quadDegree, stencilsize, f);
    disc1d::impl::ConservativeTransferAdaptor interpolator(mDiscSrc, mDiscDest, mInterpolator);

    std::vector<double> fDest(mDiscDest->get_num_dofs(), 42), error(mDiscDest->get_num_dofs());
    interpolator.interpolate(mFSrc, fDest);

    double errorNorm         = compute_error_norm(mDiscDest, mFDestEx, fDest);
    double integralSrc       = disc1d::impl::integrate_function(mDiscSrc, mFSrc);
    double integralDest      = disc1d::impl::integrate_function(mDiscDest, fDest);
    double conservationError = std::abs(integralSrc - integralDest);

    errorSlope.add_measurement(1.0 / numelSrc, errorNorm);
    conservationSlope.add_measurement(1.0 / numelSrc, conservationError);

    EXPECT_NEAR(conservationError, 0, 1e-13);

    numelSrc *= 2;
    numelDest *= 2;
  }

  std::cout << "\nsolution error:" << std::endl;
  errorSlope.output_data();

  std::cout << "\nconservation error:" << std::endl;
  conservationSlope.output_data();
}

TEST_F(SimpleInterpolatorTester, Conservation)
{
  double xmin = 0, xmax = 2;
  int numelSrc = 5, numelDest = 11;
  int basisDegree = 1, quadDegree = 3;

  double mean  = (xmax + xmin) / 2;
  double sigma = (xmax - xmin) / 4;
  double pi    = std::atan(1) * 4;

  auto f = [&](double x) {
    return std::exp(-(x - mean) * (x - mean) / (2 * sigma * sigma)) / (sigma * std::sqrt(2 * pi));
  };
  initialize(xmin, xmax, numelSrc, numelDest, basisDegree, quadDegree, f);

  std::vector<double> fDest(mDiscDest->get_num_dofs(), 42);
  mInterpolator->interpolate(mFSrc, fDest);

  double integralSrc    = disc1d::impl::integrate_function(mDiscSrc, mFSrc);
  double integralDest   = disc1d::impl::integrate_function(mDiscDest, fDest);
  double integralDestEx = disc1d::impl::integrate_function(mDiscDest, mFDestEx);

  std::cout << "integral_src = " << integralSrc << ", integral_dest = " << integralDest
            << ", diff = " << std::abs(integralSrc - integralDest)
            << ", relative error = " << std::abs(integralSrc - integralDest) / integralSrc << std::endl;
  std::cout << "integral_dest_ex = " << integralDestEx << std::endl;
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
