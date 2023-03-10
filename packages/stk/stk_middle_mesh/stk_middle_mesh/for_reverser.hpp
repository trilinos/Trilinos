#ifndef FOR_REVERSER_H
#define FOR_REVERSER_H

#include <cstddef>

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

// provides a utility for runtime reversing for loops
// std::vector<T> vals{1, 2, 3}
// for (std::size_t i=0; i < vals.size(); ++i)  should be written as:
//
// ForReverser r(vals.size());
// for (std::size_t i=r.init(), r.condition(), r.increment(i))
class ForReverser
{
  public:
    explicit ForReverser(const std::size_t len, const bool reversed = false)
      : m_len(len)
      , m_reversed(reversed)
      , m_increment(reversed ? -1 : 1)
    {}

    std::size_t init()
    {
      if (m_reversed)
        return m_len - 1;
      else
        return 0;
    }

    bool condition(const std::size_t idx)
    {
      if (m_reversed)
        return idx < m_len; // std::size_to rolls over to a large positive value
      else
        return idx < m_len;
    }

    void increment(std::size_t& idx) { idx += m_increment; }

  private:
    const std::size_t m_len;
    const bool m_reversed;
    const int m_increment;
};

} // namespace impl
} // namespace utils
} // namespace middle_mesh
} // namespace stk
#endif
