#include <stk_percept/function/internal/Dimensions.hpp>

namespace stk_classic
{
  namespace percept
  {

    std::ostream &operator<<(std::ostream& out, const Dimensions& dim)
    {
      //out << "Dimensions: ";
      out << "[ ";
      for (unsigned int i = 0; i < dim.size(); i++)
        {
          out << " " << dim[i];
        }
      out << "] ";
      return out;
    }

  }
}
