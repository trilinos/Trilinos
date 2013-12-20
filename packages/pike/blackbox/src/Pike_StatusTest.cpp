#include "Pike_StatusTest.hpp"

namespace pike {

  std::ostream& operator<<(std::ostream& os, const pike::StatusTest& st)
  {
    st.describe(os,st.getVerbLevel());
    return os;
  }

}
