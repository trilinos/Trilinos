#include <string>
#include "Teuchos_Comm.hpp"

namespace buildDemo {
  int src_file(std::string& infile, Teuchos::Comm<int>& comm);
}
