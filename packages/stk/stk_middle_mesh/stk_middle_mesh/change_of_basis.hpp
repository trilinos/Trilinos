#ifndef CHANGE_OF_BASIS
#define CHANGE_OF_BASIS

#include "projection.hpp"
#include <array>
#include <cassert>

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

class ChangeOfBasis
{
  public:
    ChangeOfBasis()
      : m_basis{Point(), Point(), Point()}
    {}
    // basis vectors
    explicit ChangeOfBasis(const std::array<Point, 3>& basis)
      : m_basis(basis)
    {
#ifndef NDEBUG
      // std::cout << "dot1 = " << std::abs(dot(basis[0], basis[1])) << std::endl;
      assert(std::abs(dot(m_basis[0], m_basis[1])) < 1e-11);
      assert(std::abs(dot(m_basis[0], m_basis[2])) < 1e-11);
      assert(std::abs(dot(m_basis[1], m_basis[2])) < 1e-11);
#endif
    }

    // given a point in xyz space, returns the same point on the
    // basis space
    // (more generally, the point can be in whatever coordinate system the
    //  basis vectors are in)
    template <typename T>
    PointT<T> project_forward(const PointT<T>& pt)
    {
      Point ptPrime;
      ptPrime.x = dot(m_basis[0], pt);
      ptPrime.y = dot(m_basis[1], pt);
      ptPrime.z = dot(m_basis[2], pt);

      return ptPrime;
    }

    template <typename T>
    PointT<T> project_back(const PointT<T>& ptPrime)
    {
      Point pt;
      pt.x += ptPrime.x * m_basis[0].x;
      pt.y += ptPrime.x * m_basis[0].y;
      pt.z += ptPrime.x * m_basis[0].z;

      pt.x += ptPrime.y * m_basis[1].x;
      pt.y += ptPrime.y * m_basis[1].y;
      pt.z += ptPrime.y * m_basis[1].z;

      pt.x += ptPrime.z * m_basis[2].x;
      pt.y += ptPrime.z * m_basis[2].y;
      pt.z += ptPrime.z * m_basis[2].z;

      return pt;
    }

  private:
    std::array<Point, 3> m_basis;

    friend std::ostream& operator<<(std::ostream& os, const ChangeOfBasis& proj);
};

std::array<Point, 3> compute_basis(Point norm);

inline std::ostream& operator<<(std::ostream& os, const ChangeOfBasis& proj)
{
  for (int i = 0; i < 3; ++i)
    os << "basis vector " << i << " = " << proj.m_basis[i] << std::endl;

  return os;
}

} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
#endif
