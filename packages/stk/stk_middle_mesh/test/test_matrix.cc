#include "matrix.h"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

template <typename T>
using Matrix = utils::impl::Matrix<T>;

TEST(Matrix, Sizes)
{
  Matrix<double> a(2, 3);
  EXPECT_EQ(a.extent(0), 2);
  EXPECT_EQ(a.extent(1), 3);
  EXPECT_EQ(a.extent0(), 2);
  EXPECT_EQ(a.extent1(), 3);
}

TEST(Matrix, Indexing)
{
  Matrix<double> a(2, 3);
  a(0, 0) = 1;
  a(0, 1) = 2;
  a(0, 2) = 3;
  a(1, 0) = 4;
  a(1, 1) = 5;
  a(1, 2) = 6;

  EXPECT_EQ(a(0, 0), 1);
  EXPECT_EQ(a(0, 1), 2);
  EXPECT_EQ(a(0, 2), 3);
  EXPECT_EQ(a(1, 0), 4);
  EXPECT_EQ(a(1, 1), 5);
  EXPECT_EQ(a(1, 2), 6);

  double* p = a.data();
  EXPECT_NE(p, nullptr);

  const Matrix<double> b(2, 3);
  const double* p2 = b.data();
  EXPECT_NE(p2, nullptr);
}

TEST(Matrix, CBlas_dgemv)
{
  // regular
  {
    Matrix<double> a(2, 2, {1, 2, 3, 4});
    double x[2] = {1, 2};
    double y[2] = {0, 0};
    matvec(1, a, x, 0, y);

    EXPECT_FLOAT_EQ(y[0], 5);
    EXPECT_FLOAT_EQ(y[1], 11);
  }

  // transposed
  {
    Matrix<double> a(2, 2, {1, 2, 3, 4});
    double x[2] = {1, 2};
    double y[2] = {0, 0};
    matvec(1, a, x, 0, y, utils::impl::BlasTrans::Trans);

    EXPECT_FLOAT_EQ(y[0], 7);
    EXPECT_FLOAT_EQ(y[1], 10);
  }

  // rectangular
  {
    Matrix<double> a(2, 3, {1, 2, 3, 4, 5, 6});
    double x[3] = {1, 2, 3};
    double y[2] = {0, 0};
    matvec(1, a, x, 0, y);

    EXPECT_FLOAT_EQ(y[0], 14);
    EXPECT_FLOAT_EQ(y[1], 32);
  }

  {
    Matrix<double> a(3, 2, {1, 2, 3, 4, 5, 6});
    double x[2] = {1, 2};
    double y[3] = {
        0,
        0,
        0,
    };
    matvec(1, a, x, 0, y);

    EXPECT_FLOAT_EQ(y[0], 5);
    EXPECT_FLOAT_EQ(y[1], 11);
    EXPECT_FLOAT_EQ(y[2], 17);
  }
}

TEST(Matrix, CBlas_dgesv)
{
  // regular
  {
    Matrix<double> a(2, 2, {1, 2, 3, 4});
    Matrix<double> a0 = a;
    double b[2]       = {1, 2};
    double x[2]       = {b[0], b[1]};
    double b2[3]      = {0, 0};
    int ipiv[2]       = {0, 0};

    solve_linear_system(a, ipiv, x);
    matvec(1, a0, x, 0, b2);

    EXPECT_NEAR(b[0], b2[0], 1e-13);
    EXPECT_NEAR(b[1], b2[1], 1e-13);
  }
}

TEST(Matrix, lapack_qr)
{
  // regular
  {
    Matrix<double> a(2, 2, {1, 2, 3, 4});
    Matrix<double> a0 = a;
    double b[2]       = {1, 2};
    double x[2]       = {b[0], b[1]};
    double b2[3]      = {0, 0};

    Matrix<double> work(2, 2);
    std::vector<double> tau(2);
    utils::impl::compute_qr_factorization(a, work, tau.data());
    utils::impl::solve_qr_factorization(a, work, tau.data(), x);

    matvec(1, a0, x, 0, b2);

    EXPECT_NEAR(b[0], b2[0], 1e-13);
    EXPECT_NEAR(b[1], b2[1], 1e-13);
  }
}

TEST(Matrix, CBlas_least_squares)
{
  // overdetermined
  {
    Matrix<double> a(3, 2, {1, 2, 3, 4, 5, 6});
    Matrix<double> b(3, 1, {1, 2, 3});
    Matrix<double> work(a);
    solve_least_squares(a, b, work);

    EXPECT_NEAR(b(0, 0), 0.0, 1e-13);
    EXPECT_NEAR(b(1, 0), 0.5, 1e-13);
  }

  // underdetermined
  {
    Matrix<double> a(2, 3, {1, 2, 3, 4, 5, 6});
    Matrix<double> a0 = a;
    Matrix<double> b(3, 1, {1, 2, -666}); // we are solving a 2 x 3 linear system, but
                                          // B must be larged enough for both x and b
    double b0[2] = {b(0, 0), b(1, 0)};

    Matrix<double> work(a);
    solve_least_squares(a, b, work);

    double x[3] = {b(0, 0), b(1, 0), b(2, 0)};
    double b2[2];
    matvec(1, a0, x, 0, b2);

    EXPECT_NEAR(b2[0] - b0[0], 0, 1e-12);
    EXPECT_NEAR(b2[1] - b0[1], 0, 1e-12);
  }
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
