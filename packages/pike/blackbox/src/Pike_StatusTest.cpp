#include "Pike_StatusTest.hpp"

namespace pike {

  std::ostream& operator<<(std::ostream& os, const pike::StatusTest& st)
  {
    st.describe(os,st.getVerbLevel());
    return os;
  }

  std::string statusToString(const pike::SolveStatus& s)
  {
    std::string value;

    if (s == pike::UNCHECKED)
      value = "??...........";
    else if (s == pike::UNCONVERGED)
      value = "**...........";
    else if (s == pike::CONVERGED)
      value = "CONVERGED....";
    else if (s == pike::FAILED)
      value = "FAILED.......";
    
    return value;
  }
}
