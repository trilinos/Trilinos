#include "gtest/gtest.h"

#include "complex_utils.h"
#include "mat2x2.h"

namespace stk {
namespace middle_mesh {
namespace impl {

using Complex = std::complex<double>;

template <typename T>
using Mat2x2 = utils::impl::Mat2x2<T>;

// func = f(A) -> T
template <typename T, typename Tfunc>
Mat2x2<T> compute_fd(Mat2x2<T>& a, Tfunc func, double eps = 1e-7)
{
  Mat2x2<T> b;
  T val0 = func(a);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
    {
      a(i, j) += eps;
      b(i, j) = (func(a) - val0) / eps;
      a(i, j) -= eps;
    }

  return b;
}

TEST(Mat2x2, Addition)
{
  // scalar
  {
    Mat2x2<double> a = {1, 2, 3, 4};
    auto b           = a + 2;
    EXPECT_FLOAT_EQ(b(0, 0), 3);
    EXPECT_FLOAT_EQ(b(0, 1), 4);
    EXPECT_FLOAT_EQ(b(1, 0), 5);
    EXPECT_FLOAT_EQ(b(1, 1), 6);

    auto c = 2 + a;
    EXPECT_FLOAT_EQ(c(0, 0), 3);
    EXPECT_FLOAT_EQ(c(0, 1), 4);
    EXPECT_FLOAT_EQ(c(1, 0), 5);
    EXPECT_FLOAT_EQ(c(1, 1), 6);
  }

  // matrix
  {
    // Mat2x2<double> A = {1, 2, 3, 4};
    // Mat2x2<double> B = {5, 6, 7, 8};

    Mat2x2<double> a, b;
    a(0, 0) = 1;
    a(0, 1) = 2;
    a(1, 0) = 3;
    a(1, 1) = 4;

    b(0, 0) = 5;
    b(0, 1) = 6;
    b(1, 0) = 7;
    b(1, 1) = 8;

    // std::cout << "A = " << A << std::endl;
    // std::cout << "B = " << B << std::endl;
    auto c = a + b;
    // std::cout << "C = " << C << std::endl;
    EXPECT_FLOAT_EQ(c(0, 0), 6);
    EXPECT_FLOAT_EQ(c(0, 1), 8);
    EXPECT_FLOAT_EQ(c(1, 0), 10);
    EXPECT_FLOAT_EQ(c(1, 1), 12);

    auto d = b + a;
    EXPECT_FLOAT_EQ(d(0, 0), 6);
    EXPECT_FLOAT_EQ(d(0, 1), 8);
    EXPECT_FLOAT_EQ(d(1, 0), 10);
    EXPECT_FLOAT_EQ(d(1, 1), 12);
  }
}

TEST(Mat2x2, Subtraction)
{
  // scalar
  {
    Mat2x2<double> a = {1, 2, 3, 4};
    auto b           = a - 2;
    EXPECT_FLOAT_EQ(b(0, 0), -1);
    EXPECT_FLOAT_EQ(b(0, 1), 0);
    EXPECT_FLOAT_EQ(b(1, 0), 1);
    EXPECT_FLOAT_EQ(b(1, 1), 2);

    auto c = 2 - a;
    EXPECT_FLOAT_EQ(c(0, 0), 1);
    EXPECT_FLOAT_EQ(c(0, 1), 0);
    EXPECT_FLOAT_EQ(c(1, 0), -1);
    EXPECT_FLOAT_EQ(c(1, 1), -2);
  }

  // matrix
  {
    Mat2x2<double> a = {1, 2, 3, 4};
    Mat2x2<double> b = {6, 8, 10, 12};
    auto c           = a - b;
    EXPECT_FLOAT_EQ(c(0, 0), -5);
    EXPECT_FLOAT_EQ(c(0, 1), -6);
    EXPECT_FLOAT_EQ(c(1, 0), -7);
    EXPECT_FLOAT_EQ(c(1, 1), -8);

    auto d = -a;
    EXPECT_FLOAT_EQ(d(0, 0), -1);
    EXPECT_FLOAT_EQ(d(0, 1), -2);
    EXPECT_FLOAT_EQ(d(1, 0), -3);
    EXPECT_FLOAT_EQ(d(1, 1), -4);
  }
}

TEST(Mat2x2, Multiplication)
{
  // scalar
  {
    Mat2x2<double> a = {1, 2, 3, 4};
    auto b           = a * 2;
    EXPECT_FLOAT_EQ(b(0, 0), 2);
    EXPECT_FLOAT_EQ(b(0, 1), 4);
    EXPECT_FLOAT_EQ(b(1, 0), 6);
    EXPECT_FLOAT_EQ(b(1, 1), 8);

    auto c = 2 * a;
    EXPECT_FLOAT_EQ(c(0, 0), 2);
    EXPECT_FLOAT_EQ(c(0, 1), 4);
    EXPECT_FLOAT_EQ(c(1, 0), 6);
    EXPECT_FLOAT_EQ(c(1, 1), 8);
  }

  // matrix
  {
    Mat2x2<double> a = {1, 2, 3, 4};
    Mat2x2<double> b = {5, 6, 7, 8};
    auto c           = a * b;
    EXPECT_FLOAT_EQ(c(0, 0), 19);
    EXPECT_FLOAT_EQ(c(0, 1), 22);
    EXPECT_FLOAT_EQ(c(1, 0), 43);
    EXPECT_FLOAT_EQ(c(1, 1), 50);
  }
}

TEST(Mat2x2, Division)
{
  // scalar
  {
    Mat2x2<double> a = {1, 2, 3, 4};
    auto b           = a / 2;
    EXPECT_FLOAT_EQ(b(0, 0), 0.5);
    EXPECT_FLOAT_EQ(b(0, 1), 1);
    EXPECT_FLOAT_EQ(b(1, 0), 1.5);
    EXPECT_FLOAT_EQ(b(1, 1), 2);
  }
}

TEST(Mat2x2, Inverse)
{
  {
    Mat2x2<double> a = {1, 2, 3, 4};

    inverse2x2(a);

    EXPECT_FLOAT_EQ(a(0, 0), -2.0);
    EXPECT_FLOAT_EQ(a(0, 1), 1.0);
    EXPECT_FLOAT_EQ(a(1, 0), 1.5);
    EXPECT_FLOAT_EQ(a(1, 1), -0.5);
  }

  {
    Mat2x2<double> a = {1, 2, 3, 4};
    double x[2], b[2] = {1.0, 2.0};
    matsolve2x2(a, x, b);

    EXPECT_FLOAT_EQ(x[0], 0.0);
    EXPECT_FLOAT_EQ(x[1], 0.5);
  }

  {
    Mat2x2<double> a = {0, 2, 3, 4};
    double x[2], b[2] = {1.0, 2.0};
    matsolve2x2(a, x, b);

    EXPECT_FLOAT_EQ(x[0], 0.0);
    EXPECT_FLOAT_EQ(x[1], 0.5);
  }
}

namespace {
void test_matsolve_fd(Mat2x2<double> a, double b[2])
{
  Mat2x2<double> aDot = {0, 0, 0, 0};
  double x[2], x2[2], x3[2], xDot[2], bDot[2] = {0, 0}, derivFd[2];

  // test derivative wrt A
  double eps = 1e-7;
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
    {
      matsolve2x2(a, x, b);

      a(i, j) += eps;
      matsolve2x2(a, x2, b);
      a(i, j) -= eps;

      derivFd[0] = (x2[0] - x[0]) / eps;
      derivFd[1] = (x2[1] - x[1]) / eps;

      aDot(i, j) = 1;
      matsolve2x2_dot(a, aDot, x3, xDot, b, bDot);
      aDot(i, j) = 0;

      EXPECT_NEAR(xDot[0], derivFd[0], 1e-6);
      EXPECT_NEAR(xDot[1], derivFd[1], 1e-6);
      EXPECT_NEAR(x[0], x3[0], 1e-13);
      EXPECT_NEAR(x[1], x3[1], 1e-13);
    }

  // test derivative wrt b
  for (int i = 0; i < 2; ++i)
  {
    matsolve2x2(a, x, b);

    b[i] += eps;
    matsolve2x2(a, x2, b);
    b[i] -= eps;

    derivFd[0] = (x2[0] - x[0]) / eps;
    derivFd[1] = (x2[1] - x[1]) / eps;

    bDot[i] = 1;
    matsolve2x2_dot(a, aDot, x, xDot, b, bDot);
    bDot[i] = 0;

    EXPECT_NEAR(xDot[0], derivFd[0], 1e-6);
    EXPECT_NEAR(xDot[1], derivFd[1], 1e-6);
  }
}
} // namespace

TEST(Mat2x2, solve_dot)
{
  {
    Mat2x2<double> a = {1, 2, 3, 4};
    double b[2]      = {1.0, 2.0};
    test_matsolve_fd(a, b);
  }

  {
    Mat2x2<double> a = {3, 4, 1, 2};
    double b[2]      = {2.0, 1.0};
    test_matsolve_fd(a, b);
  }
}

TEST(Mat2x2, MatVec)
{
  Mat2x2<double> a = {1, 2, 3, 4};
  double x[2]      = {1, 2};
  double b[2];

  matvec2x2(a, x, b);

  EXPECT_FLOAT_EQ(b[0], 5.0);
  EXPECT_FLOAT_EQ(b[1], 11.0);
}

TEST(Mat2x2, Norm_F)
{
  Mat2x2<double> a = {1, 2, 3, 4};

  EXPECT_FLOAT_EQ(norm_f(a), std::sqrt(30.0));

  Mat2x2<double> aBar({0, 0, 0, 0});
  double nBar = 2;
  utils::impl::norm_f_rev(a, aBar, nBar);
  Mat2x2<double> aBarFd = compute_fd(a, utils::impl::norm_f<double>);

  double tol = 1e-6;
  EXPECT_NEAR(aBar(0, 0), nBar * aBarFd(0, 0), tol);
  EXPECT_NEAR(aBar(0, 1), nBar * aBarFd(0, 1), tol);
  EXPECT_NEAR(aBar(1, 0), nBar * aBarFd(1, 0), tol);
  EXPECT_NEAR(aBar(1, 1), nBar * aBarFd(1, 1), tol);
}

TEST(Mat2x2, Det)
{
  Mat2x2<double> a = {1, 2, 3, 4};

  EXPECT_FLOAT_EQ(det2x2(a), -2);

  Mat2x2<double> aBar({0, 0, 0, 0});
  double dBar = 2;
  det2x2_rev(a, aBar, dBar);
  Mat2x2<double> aBarFd = compute_fd(a, utils::impl::det2x2<double>);

  double tol = 1e-6;
  EXPECT_NEAR(aBar(0, 0), dBar * aBarFd(0, 0), tol);
  EXPECT_NEAR(aBar(0, 1), dBar * aBarFd(0, 1), tol);
  EXPECT_NEAR(aBar(1, 0), dBar * aBarFd(1, 0), tol);
  EXPECT_NEAR(aBar(1, 1), dBar * aBarFd(1, 1), tol);
}

TEST(Mat2x2, matmat_dot)
{
  double eps = 1e-20;
  Complex pert(0, eps);

  Mat2x2<Complex> a = {Complex(1, 0), Complex(2, 0), Complex(3, 0), Complex(4, 0)};
  Mat2x2<Complex> b = {Complex(5, 0), Complex(6, 0), Complex(7, 0), Complex(8, 0)};

  Mat2x2<std::array<Complex, 2>> aDot = {{Complex(9, 0), Complex(10, 0)},
                                         {Complex(11, 0), Complex(12, 0)},
                                         {Complex(13, 0), Complex(14, 0)},
                                         {Complex(15, 0), Complex(16, 0)}};

  auto cDot = matmat_dot(aDot, b);

  // complex step
  // auto C = A*B;
  for (int i = 0; i < 2; ++i)
  {
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        a(j, k) += pert * aDot(j, k)[i];

    auto c2 = a * b;

    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
      {
        auto val = c2(i, j).imag() / eps;
        EXPECT_FLOAT_EQ(val, cDot(i, j)[i].real());
        a(j, k) -= pert * aDot(j, k)[i];
      }
  }
}

TEST(Mat2x2, norm_f_dot)
{
  double eps = 1e-20;
  Complex pert(0, eps);

  Mat2x2<Complex> a                   = {Complex(1, 0), Complex(2, 0), Complex(3, 0), Complex(4, 0)};
  Mat2x2<std::array<Complex, 2>> aDot = {{Complex(9, 0), Complex(10, 0)},
                                         {Complex(11, 0), Complex(12, 0)},
                                         {Complex(13, 0), Complex(14, 0)},
                                         {Complex(15, 0), Complex(16, 0)}};

  auto valDot = utils::impl::norm_f_dot(a, aDot);

  // complex step
  // auto val0 = utils::impl::norm_f(A);
  for (int i = 0; i < 2; ++i)
  {
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        a(j, k) += pert * aDot(j, k)[i];

    double val = utils::impl::norm_f(a).imag() / eps;
    EXPECT_FLOAT_EQ(val, valDot[i].real());

    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        a(j, k) -= pert * aDot(j, k)[i];
  }
}

TEST(Mat2x2, det2x2_dot)
{
  double eps = 1e-20;
  Complex pert(0, eps);

  Mat2x2<Complex> a                   = {Complex(1, 0), Complex(2, 0), Complex(3, 0), Complex(4, 0)};
  Mat2x2<std::array<Complex, 2>> aDot = {{Complex(9, 0), Complex(10, 0)},
                                         {Complex(11, 0), Complex(12, 0)},
                                         {Complex(13, 0), Complex(14, 0)},
                                         {Complex(15, 0), Complex(16, 0)}};

  auto valDot = det2x2_dot(a, aDot);

  // complex step
  auto val0 = utils::impl::norm_f(a);
  for (int i = 0; i < 2; ++i)
  {
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        a(j, k) += pert * aDot(j, k)[i];

    double val = (det2x2(a) - val0).imag() / eps;
    EXPECT_FLOAT_EQ(val, valDot[i].real());

    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        a(j, k) -= pert * aDot(j, k)[i];
  }
}

TEST(Mat2x2, det2x2_rev_dot)
{
  double eps = 1e-20;
  Complex pert(0, eps);

  Mat2x2<Complex> a                   = {Complex(1, 0), Complex(2, 0), Complex(3, 0), Complex(4, 0)};
  Mat2x2<std::array<Complex, 2>> aDot = {{Complex(9, 0), Complex(10, 0)},
                                         {Complex(11, 0), Complex(12, 0)},
                                         {Complex(13, 0), Complex(14, 0)},
                                         {Complex(15, 0), Complex(16, 0)}};
  Complex dBar(5, 0);
  std::array<Complex, 2> dBarDot = {Complex(6, 0), Complex(7, 0)};

  Mat2x2<std::array<Complex, 2>> aBarDot;
  det2x2_rev_dot(a, aDot, dBar, dBarDot, aBarDot);

  for (int i = 0; i < 2; ++i)
  {
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        a(j, k) += pert * aDot(j, k)[i];

    dBar += pert * dBarDot[i];

    Mat2x2<Complex> aBar;
    det2x2_rev(a, aBar, dBar);

    dBar -= pert * dBarDot[i];
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
      {
        a(j, k) -= pert * aDot(j, k)[i];
        EXPECT_FLOAT_EQ(aBar(j, k).imag() / eps, aBarDot(j, k)[i].real());
      }
  }
}

TEST(Mat2x2, norm_f_rev_dot)
{
  double eps = 1e-20;
  Complex pert(0, eps);

  Mat2x2<Complex> a                   = {Complex(1, 0), Complex(2, 0), Complex(3, 0), Complex(4, 0)};
  Mat2x2<std::array<Complex, 2>> aDot = {{Complex(9, 0), Complex(10, 0)},
                                         {Complex(11, 0), Complex(12, 0)},
                                         {Complex(13, 0), Complex(14, 0)},
                                         {Complex(15, 0), Complex(16, 0)}};
  Complex dBar(5, 0);
  std::array<Complex, 2> dBarDot = {Complex(6, 0), Complex(7, 0)};

  Mat2x2<std::array<Complex, 2>> aBarDot;
  utils::impl::norm_f_rev_dot(a, aDot, dBar, dBarDot, aBarDot);

  for (int i = 0; i < 2; ++i)
  {
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        a(j, k) += pert * aDot(j, k)[i];

    dBar += pert * dBarDot[i];

    Mat2x2<Complex> aBar;
    norm_f_rev(a, aBar, dBar);

    dBar -= pert * dBarDot[i];
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
      {
        a(j, k) -= pert * aDot(j, k)[i];
        EXPECT_FLOAT_EQ(aBar(j, k).imag() / eps, aBarDot(j, k)[i].real());
      }
  }
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
