#include <cassert>
#include <stdexcept>
#include <vector>

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

// TODO: untested
template <typename T>
class BandedMatrix
{
  public:
    BandedMatrix(int m, int bandwidth)
      : m_data(m * bandwidth)
      , m_bandwidth(m)
      , m_halfBandwidth(m / 2)
      , m_m(m)
    {
      if (bandwidth % 2 != 1)
        throw std::runtime_error("bandwidth must be odd number");
    }

    T& operator()(int i, int j) { return m_data[get_idx(i, j)]; }

    const T& operator()(int i, int j) const { return m_data[get_idx(i, j)]; }

    int extent(const int dim) const { return m_m; }

    void fill(const T& val)
    {
      for (auto& v : m_data)
        v = val;
    }

    int extent0() const { return m_m; }

    int extent1() const { return m_m; }

  private:
    int get_idx(int i, int j)
    {
      assert(std::abs(i - j) <= m_half_bandwidth);
      assert(i >= 0 && i <= m_m);
      assert(j >= 0 && j <= m_m);

      int rowStart    = m_bandwidth * i;
      int idxOfCenter = rowStart + m_halfBandwidth;
      return idxOfCenter - (i - j);
    }

    std::vector<T> m_data;
    int m_bandwidth;
    int m_halfBandwidth;
    int m_m;
};

} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
