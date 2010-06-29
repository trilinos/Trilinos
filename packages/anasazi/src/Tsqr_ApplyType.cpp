#include <Tsqr_ApplyType.hpp>
#include <stdexcept>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  const ApplyType ApplyType::NoTranspose = ApplyType ("N");
  const ApplyType ApplyType::Transpose = ApplyType ("T");
  const ApplyType ApplyType::ConjugateTranspose = ApplyType ("H");

  bool 
  ApplyType::decide_transposed (const std::string& op) const 
  {
    if (op[0] == 'N' || op[0] == 'n')
      return false;
    else if (op[0] == 'T' || op[0] == 't' || op[0] == 'H' || op[0] == 'h')
      return true;
    else
      throw std::invalid_argument ("Invalid \"op\" argument \"" + op + "\"");
  }

  ApplyType::ApplyType_
  ApplyType::decide_apply_type (const std::string& op) const 
  {
    if (op[0] == 'T' || op[0] == 't')
      return Transpose_;
    else if (op[0] == 'N' || op[0] == 'n')
      return NoTranspose_;
    else if (op[0] == 'H' || op[0] == 'h')
      return ConjugateTranspose_;
    else
      throw std::invalid_argument ("Invalid \"op\" argument \"" + op + "\"");
  }
}

