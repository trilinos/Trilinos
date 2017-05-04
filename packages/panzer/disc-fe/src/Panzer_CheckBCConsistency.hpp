#ifndef PANZER_CHECK_BC_CONSISTENCY
#define PANZER_CHECK_BC_CONSISTENCY

#include <vector>
#include <string>

namespace panzer {

  class BC;

  void checkBCConsistency(const std::vector<std::string>& element_block_names,
                          const std::vector<std::string>& sideset_names,
                          const std::vector<panzer::BC>& bcs);

}

#endif
