#ifndef PACKAGES_STK_STK_TOOLS_STK_TOOLS_BLOCK_EXTRACTOR_EXTRACT_BLOCKS_HPP_
#define PACKAGES_STK_STK_TOOLS_STK_TOOLS_BLOCK_EXTRACTOR_EXTRACT_BLOCKS_HPP_

#include <vector>
#include <string>

namespace stk { namespace mesh { class BulkData; } }

namespace stk {
namespace tools {

void extract_blocks_from_file(const std::string &inFile,
                              const std::string &outFile,
                              const std::vector<std::string> &blockNames,
                              MPI_Comm comm);
void extract_blocks(stk::mesh::BulkData &oldBulk, stk::mesh::BulkData &newBulk, const std::vector<std::string> &blockNames);

}
}

#endif
