#ifndef MAT_2X2_H
#define MAT_2X2_H

#include <array>
#include <cmath>
#include <complex> // needed for std::sqrt(std::complex)
#include <initializer_list>

#include <iostream> //TODO: DEBUGGING

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

#if defined(__GNUC__) && !defined(__llvm__) && !defined(__INTEL_COMPILER)
#define FORCE_NOINLINE_GCC __attribute__((noinline))
#else
#define FORCE_NOINLINE_GCC
#endif

template <typename T>
class Mat2x2
{
  public:
    // Mat2x2() : m_data{{0, 0, 0, 0}} {}
    Mat2x2()
    {
      for (int i = 0; i < 4; ++i)
        m_data[i] = T{};
    }

    Mat2x2(std::initializer_list<T> list)
    {
      // std::array doesn't support initializer lists
      auto it = list.begin();
      for (int i = 0; i < 4; ++i)
      {
        m_data[i] = *it;
        it++;
      }
    }

    template <typename T2>
    Mat2x2(const Mat2x2<T2>& other) :
      m_data{other(0, 0), other(0, 1), other(1, 0), other(1, 1)}
    {}

    template <typename T2>
    Mat2x2<T> operator=(const Mat2x2<T2>& other)
    {
      for (int i=0; i < 2; ++i)
        for (int j=0; j < 2; ++j)
          (*this)(i, j) = other(i, j);
    }

    using value_type = T;

    T& operator()(const int i, const int j)
    {
#ifdef NDEBUG
      return m_data[2 * i + j];
#else
      return m_data.at(2 * i + j);
#endif
    }

    const T& operator()(const int i, const int j) const
    {
#ifdef NDEBUG
      return m_data[2 * i + j];
#else
      return m_data.at(2 * i + j);
#endif
    }

    static int extent(const int /*dim*/) { return 2; }

    template <typename T2>
    FORCE_NOINLINE_GCC Mat2x2<T>& operator+=(const T2& val)
    {
      for (int i = 0; i < 4; ++i)
        m_data[i] += val;

      return *this;
    }

    template <typename T2>
    FORCE_NOINLINE_GCC Mat2x2<T>& operator+=(const Mat2x2<T2>& a)
    {
      for (int i = 0; i < 2; ++i)
        for (int j=0; j < 2; ++j)
          (*this)(i, j) += a(i, j);

      return *this;
    }

    template <typename T2>
    FORCE_NOINLINE_GCC Mat2x2<T>& operator-=(const T2& val)
    {
      for (int i = 0; i < 4; ++i)
        m_data[i] -= val;

      return *this;
    }

    template <typename T2>
    FORCE_NOINLINE_GCC Mat2x2<T>& operator-=(const Mat2x2<T2>& a)
    {
      for (int i=0; i < 2; ++i)
        for (int j=0; j < 2; ++j)
          (*this)(i, j) -= a(i, j);

      return *this;
    }

    template <typename T2>
    FORCE_NOINLINE_GCC Mat2x2<T>& operator*=(const T2& val)
    {
      for (int i = 0; i < 4; ++i)
        m_data[i] *= val;

      return *this;
    }

    template <typename T2>
    FORCE_NOINLINE_GCC Mat2x2<T>& operator/=(const T2& val)
    {
      for (int i = 0; i < 4; ++i)
        m_data[i] /= val;

      return *this;
    }

  private:
    std::array<T, 4> m_data;
};

template <typename T, typename T2>
Mat2x2<std::common_type_t<T, T2>> operator+(const Mat2x2<T>& a, const T2& val)
{
  Mat2x2<std::common_type_t<T, T2>> c = a;
  c += val;
  return c;
}

template <typename T, typename T2>
Mat2x2<std::common_type_t<T, T2>> operator+(const T2& val, const Mat2x2<T>& a)
{
  return a + val;
}

template <typename T, typename T2>
Mat2x2<std::common_type_t<T, T2>> operator-(const Mat2x2<T>& a, const T2& val)
{
  Mat2x2<std::common_type_t<T, T2>> c = a;
  c -= val;
  return c;
}

template <typename T, typename T2>
Mat2x2<std::common_type_t<T, T2>> operator-(const T2& val, const Mat2x2<T>& a)
{
  Mat2x2<std::common_type_t<T, T2>> c = -a;
  return c + val;
}

template <typename T>
Mat2x2<T> operator-(const Mat2x2<T>& a)
{
  return a * -1;
}

template <typename T, typename T2>
Mat2x2<std::common_type_t<T, T2>> operator*(const Mat2x2<T>& a, const T2& val)
{
  Mat2x2<std::common_type_t<T, T2>> c = a;
  c *= val;
  return c;
}

template <typename T, typename T2>
Mat2x2<std::common_type_t<T, T2>> operator*(const T2& val, const Mat2x2<T>& a)
{
  return a * val;
}

template <typename T, typename T2>
Mat2x2<std::common_type_t<T, T2>> operator*(const Mat2x2<T>& a, const Mat2x2<T2>& b)
{
  Mat2x2<std::common_type_t<T, T2>> c;
  c(0, 0) = a(0, 0) * b(0, 0) + a(0, 1) * b(1, 0);
  c(0, 1) = a(0, 0) * b(0, 1) + a(0, 1) * b(1, 1);
  c(1, 0) = a(1, 0) * b(0, 0) + a(1, 1) * b(1, 0);
  c(1, 1) = a(1, 0) * b(0, 1) + a(1, 1) * b(1, 1);
  return c;
}

template <typename T, typename T2>
Mat2x2<std::common_type_t<T, T2>> operator/(const Mat2x2<T>& a, const T2& val)
{
  Mat2x2<std::common_type_t<T, T2>> c = a;
  c /= val;
  return c;
}

template <typename T, typename T2>
Mat2x2<std::common_type_t<T, T2>> operator+(const Mat2x2<T>& a, const Mat2x2<T2>& b)
{
  Mat2x2<std::common_type_t<T, T2>> c = a;
  c += b;
  return c;
}

template <typename T, typename T2>
Mat2x2<std::common_type_t<T, T2>> operator-(const Mat2x2<T>& a, const Mat2x2<T2>& b)
{
  Mat2x2<std::common_type_t<T, T2>> c = a;
  c -= b;
  return c;
}

template <typename T>
void inverse2x2(Mat2x2<T>& a)
{
  auto a11 = a(0, 0);
  auto a12 = a(0, 1);
  auto a21 = a(1, 0);
  auto a22 = a(1, 1);

  auto detA = det2x2(a);

  a(0, 0) = a22 / detA;
  a(0, 1) = -a12 / detA;
  a(1, 0) = -a21 / detA;
  a(1, 1) = a11 / detA;
}

// solve Ax = b
template <typename T, typename T2>
void matsolve2x2(const Mat2x2<T>& a, std::common_type_t<T, T2> x[2], const T2 b[2])
{
  using Tcommon = std::common_type_t<T, T2>;
  if (std::abs(a(0, 0)) < std::abs(a(1, 0)))
  {
    T fac = a(0, 0) / a(1, 0);
    // T A00 = A(0, 0) - fac*A(1, 0); // = 0
    T a01 = a(0, 1) - fac * a(1, 1);
    Tcommon b0  = b[0] - fac * b[1];

    x[1] = b0 / a01;
    x[0] = (b[1] - a(1, 1) * x[1]) / a(1, 0);
  } else
  {
    // Gaussian elimination
    T fac = a(1, 0) / a(0, 0);
    // T A10 = A(1, 0) - fac*A(0, 0);  // = 0
    T a11 = a(1, 1) - fac * a(0, 1);
    Tcommon b1  = b[1] - fac * b[0];

    x[1] = b1 / a11;
    x[0] = (b[0] - a(0, 1) * x[1]) / a(0, 0);
  }
}

template <typename T, typename T2, typename T3, typename T4, typename T5, typename T6>
void matsolve2x2_dot(const Mat2x2<T>& a, const Mat2x2<T2>& aDot, T3 x[2], T4 xDot[2], const T5 b[2], const T6 bDot[2])
{
  using Tcommon = std::common_type_t<T, T2>;
  if (std::abs(a(0, 0)) < std::abs(a(1, 0)))
  {
    Tcommon fac    = a(0, 0) / a(1, 0);
    Tcommon facDot = aDot(0, 0) / a(1, 0) - aDot(1, 0) * a(0, 0) / (a(1, 0) * a(1, 0));

    // T A00 = A(0, 0) - fac*A(1, 0); // = 0
    Tcommon a01    = a(0, 1) - fac * a(1, 1);
    Tcommon a01Dot = aDot(0, 1) - fac * aDot(1, 1) - facDot * a(1, 1);

    Tcommon b0    = b[0] - fac * b[1];
    Tcommon b0Dot = bDot[0] - fac * bDot[1] - facDot * b[1];

    x[1]    = b0 / a01;
    xDot[1] = b0Dot / a01 - b0 * a01Dot / (a01 * a01);

    x[0]    = (b[1] - a(1, 1) * x[1]) / a(1, 0);
    xDot[0] = (bDot[1] - aDot(1, 1) * x[1] - a(1, 1) * xDot[1]) / a(1, 0) -
              aDot(1, 0) * (b[1] - a(1, 1) * x[1]) / (a(1, 0) * a(1, 0));
  } else
  {
    // Gaussian elimination
    Tcommon fac    = a(1, 0) / a(0, 0);
    Tcommon facDot = aDot(1, 0) / a(0, 0) - aDot(0, 0) * a(1, 0) / (a(0, 0) * a(0, 0));

    // T A10 = A(1, 0) - fac*A(0, 0);  // = 0
    Tcommon a11    = a(1, 1) - fac * a(0, 1);
    Tcommon a11Dot = aDot(1, 1) - fac * aDot(0, 1) - facDot * (a(0, 1));

    Tcommon b1    = b[1] - fac * b[0];
    Tcommon b1Dot = bDot[1] - fac * bDot[0] - facDot * b[0];

    x[1]    = b1 / a11;
    xDot[1] = b1Dot / a11 - a11Dot * b1 / (a11 * a11);

    x[0]    = (b[0] - a(0, 1) * x[1]) / a(0, 0);
    xDot[0] = (bDot[0] - aDot(0, 1) * x[1] - a(0, 1) * xDot[1]) / a(0, 0) -
              aDot(0, 0) * (b[0] - a(0, 1) * x[1]) / (a(0, 0) * a(0, 0));
  }
}

template <typename T, typename T2>
void matvec2x2(const Mat2x2<T>& a, const T2 x[2], std::common_type_t<T, T2> b[2])
{
  b[0] = a(0, 0) * x[0] + a(0, 1) * x[1];
  b[1] = a(1, 0) * x[0] + a(1, 1) * x[1];
}

template <typename T>
T det2x2(const Mat2x2<T>& a)
{
  return a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0);
}

template <typename T, typename T2, typename T3>
void det2x2_rev(const Mat2x2<T>& a, Mat2x2<T2>& aBar, const T3 dBar)
{
  aBar(0, 0) = dBar * a(1, 1);
  aBar(0, 1) = -dBar * a(1, 0);
  aBar(1, 0) = -dBar * a(0, 1);
  aBar(1, 1) = dBar * a(0, 0);
}

// Frobenious norm
template <typename T>
T norm_f(const Mat2x2<T>& a)
{
  T val = a(0, 0) * a(0, 0) + a(0, 1) * a(0, 1) + a(1, 0) * a(1, 0) + a(1, 1) * a(1, 1);
  return std::sqrt(val);
}

template <typename T, typename T2, typename T3>
void norm_f_rev(const Mat2x2<T>& a, Mat2x2<T2>& aBar, const T3 dBar)
{
  T val = a(0, 0) * a(0, 0) + a(0, 1) * a(0, 1) + a(1, 0) * a(1, 0) + a(1, 1) * a(1, 1);

  auto valBar = 0.5 * dBar / std::sqrt(val);

  aBar(0, 0) = 2.0 * a(0, 0) * valBar;
  aBar(0, 1) = 2.0 * a(0, 1) * valBar;
  aBar(1, 0) = 2.0 * a(1, 0) * valBar;
  aBar(1, 1) = 2.0 * a(1, 1) * valBar;
}

template <typename T>
void transpose(Mat2x2<T>& a)
{
  T tmp   = a(0, 1);
  a(0, 1) = a(1, 0);
  a(1, 0) = tmp;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Mat2x2<T>& a)
{
  os << "Mat2x2:" << std::endl;
  ;
  os << "[" << a(0, 0) << ", " << a(0, 1) << "]" << std::endl;
  os << "[" << a(1, 0) << ", " << a(1, 1) << "]";

  return os;
}

// forward mode differentiation of mat-mat, where only A has derivatives
template <typename T, typename T2, std::size_t N>
Mat2x2<std::array<std::common_type_t<T, T2>, N>> matmat_dot( // const Mat2x2<T>& A,
    const Mat2x2<std::array<T, N>>& aDot, const Mat2x2<T2>& b)
{
  Mat2x2<std::array<std::common_type_t<T, T2>, N>> cDot;
  for (unsigned int i = 0; i < N; ++i)
  {
    cDot(0, 0)[i] = aDot(0, 0)[i] * b(0, 0) + aDot(0, 1)[i] * b(1, 0);
    cDot(0, 1)[i] = aDot(0, 0)[i] * b(0, 1) + aDot(0, 1)[i] * b(1, 1);
    cDot(1, 0)[i] = aDot(1, 0)[i] * b(0, 0) + aDot(1, 1)[i] * b(1, 0);
    cDot(1, 1)[i] = aDot(1, 0)[i] * b(0, 1) + aDot(1, 1)[i] * b(1, 1);
  }

  return cDot;
}

// forward mode differentiation of Frobenious norm
template <typename T, typename T2, std::size_t N>
std::array<std::common_type_t<T, T2>, N> norm_f_dot(const Mat2x2<T>& a, const Mat2x2<std::array<T2, N>>& aDot)
{
  using Tcommon = std::common_type_t<T, T2>;
  T val = a(0, 0) * a(0, 0) + a(0, 1) * a(0, 1) + a(1, 0) * a(1, 0) + a(1, 1) * a(1, 1);

  T tmp = 1.0 / (2.0 * std::sqrt(val));

  std::array<Tcommon, N> valDot;
  for (unsigned int i = 0; i < N; ++i)
  {
    Tcommon v1Dot = 2.0 * a(0, 0) * aDot(0, 0)[i] + 2.0 * a(0, 1) * aDot(0, 1)[i] + 2.0 * a(1, 0) * aDot(1, 0)[i] +
              2.0 * a(1, 1) * aDot(1, 1)[i];
    valDot[i] = tmp * v1Dot;
  }

  return valDot;
}

// forward mode differentiation of determinant
template <typename T, typename T2, std::size_t N>
std::array<std::common_type_t<T, T2>, N> det2x2_dot(const Mat2x2<T>& a, const Mat2x2<std::array<T2, N>>& aDot)
{
  // T val = A(0, 0)*A(1, 1) - A(0, 1)*A(1, 0);
  std::array<std::common_type_t<T, T2>, N> valDot;
  for (unsigned int i = 0; i < N; ++i)
    valDot[i] = a(1, 1) * aDot(0, 0)[i] + a(0, 0) * aDot(1, 1)[i] + -a(1, 0) * aDot(0, 1)[i] - a(0, 1) * aDot(1, 0)[i];

  return valDot;
}

// foward mode over reverse mode
// A_bar_dot is overwritten with the result
template <typename T, typename T2, typename T3, typename T4, typename T5, std::size_t N>
void det2x2_rev_dot(const Mat2x2<T>& a, Mat2x2<std::array<T2, N>>& aDot, const T3 dBar, const std::array<T4, N>& dBarDot,
                    Mat2x2<std::array<T5, N>>& aBarDot)
{
  /*
  A_bar(0, 0) =  d_bar*A(1, 1);
  A_bar(0, 1) = -d_bar*A(1, 0);
  A_bar(1, 0) = -d_bar*A(0, 1);
  A_bar(1, 1) =  d_bar*A(0, 0);
  */

  for (unsigned int i = 0; i < N; ++i)
  {
    aBarDot(0, 0)[i] = dBarDot[i] * a(1, 1) + dBar * aDot(1, 1)[i];
    aBarDot(0, 1)[i] = -(dBarDot[i] * a(1, 0) + dBar * aDot(1, 0)[i]);
    aBarDot(1, 0)[i] = -(dBarDot[i] * a(0, 1) + dBar * aDot(0, 1)[i]);
    aBarDot(1, 1)[i] = dBarDot[i] * a(0, 0) + dBar * aDot(0, 0)[i];
  }
}

// forward mode over reverse mode
// A_bar_dot is updated (not overwritten) with the result
template <typename T, typename T2, typename T3, typename T4, typename T5, std::size_t N>
void norm_f_rev_dot(const Mat2x2<T>& a, const Mat2x2<std::array<T2, N>>& aDot, const T3 dBar,
                    const std::array<T4, N>& dBarDot, Mat2x2<std::array<T5, N>>& aBarDot)
{
  T val = a(0, 0) * a(0, 0) + a(0, 1) * a(0, 1) + a(1, 0) * a(1, 0) + a(1, 1) * a(1, 1);

  std::array<std::common_type_t<T, T2>, N> valDot;
  for (unsigned int i = 0; i < N; ++i)
    valDot[i] = 2.0 * a(0, 0) * aDot(0, 0)[i] + 2.0 * a(0, 1) * aDot(0, 1)[i] + 2.0 * a(1, 0) * aDot(1, 0)[i] +
                2.0 * a(1, 1) * aDot(1, 1)[i];

  auto valBar = 0.5 * dBar / std::sqrt(val);
  std::array<std::common_type_t<T, T2, T4>, N> valBarDot;
  for (unsigned int i = 0; i < N; ++i)
    valBarDot[i] = 0.5 * dBarDot[i] / std::sqrt(val) - 0.25 * dBar * valDot[i] / std::pow(val, 1.5);
  /*
    A_bar(0, 0) = 2*A(0, 0)*val_bar;
    A_bar(0, 1) = 2*A(0, 1)*val_bar;
    A_bar(1, 0) = 2*A(1, 0)*val_bar;
    A_bar(1, 1) = 2*A(1, 1)*val_bar;
  */
  for (unsigned int i = 0; i < N; ++i)
  {
    aBarDot(0, 0)[i] += 2.0 * aDot(0, 0)[i] * valBar + 2.0 * a(0, 0) * valBarDot[i];
    aBarDot(0, 1)[i] += 2.0 * aDot(0, 1)[i] * valBar + 2.0 * a(0, 1) * valBarDot[i];
    aBarDot(1, 0)[i] += 2.0 * aDot(1, 0)[i] * valBar + 2.0 * a(1, 0) * valBarDot[i];
    aBarDot(1, 1)[i] += 2.0 * aDot(1, 1)[i] * valBar + 2.0 * a(1, 1) * valBarDot[i];
  }
}

} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
#endif
