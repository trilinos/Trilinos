
#include <stk_percept/function/internal/GenericFunction.hpp>

namespace stk
{
  namespace percept
  {
    std::ostream &operator<<(std::ostream& out,  GenericFunction& func)
    {
      out <<  "GenericFunction:: domain dims: " << func.getDomainDimensions() << " codomain dims: " << func.getCodomainDimensions();
      return out;
    }
  }
}
