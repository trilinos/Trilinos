
#include "Tifpack_Factory.hpp"

namespace Tifpack {

bool supportsUnsymmetric(const std::string& prec_type)
{
  bool result = false;
  if (prec_type == "POINT_RELAXATION") result = true;
  else if (prec_type == "ILUT") result = true;
  else
    throw std::runtime_error("Tifpack::supportsUnsymmetric ERROR, unrecognized prec_type");

  return result;
}

}//namespace Tifpack

