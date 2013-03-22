
#include "Ifpack2_Factory.hpp"

namespace Ifpack2 {

bool supportsUnsymmetric(const std::string& prec_type)
{
  bool result = false;
  if (prec_type == "RELAXATION" ||
      prec_type == "CHEBYSHEV"  ||
      prec_type == "DIAGONAL"   ||
      prec_type == "RILUK"      ||
      prec_type == "ILUT"       ||
      prec_type == "SCHWARZ"    ||
      prec_type == "KRYLOV")
  {
    result = true;
  }
  else {
    throw std::runtime_error("Ifpack2::supportsUnsymmetric ERROR, unrecognized prec_type");
  }

  return result;
}

}//namespace Ifpack2

