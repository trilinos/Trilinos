#ifndef MESH_INPUT_H
#define MESH_INPUT_H

#include <string>
#include <utility>
#include <vector>

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

struct MeshInput
{
    using NamePair = std::pair<std::string, std::string>;
    explicit MeshInput(const std::string& fnameIn_ = "", const std::string& fnameOut_ = "",
                       const std::string& fnameOut2_            = "",
                       const std::vector<NamePair>& interfaces_ = std::vector<NamePair>())
      : fnameIn(fnameIn_)
      , fnameOut(fnameOut_)
      , fnameOut2(fnameOut2_)
      , interfaces(interfaces_)
    {}

    std::string fnameIn;              // name of input meshfile
    std::string fnameOut;             // name to copy input mesh to (with updated
                                      // coordinates
    std::string fnameOut2;            // name of meshfile to create with nonconformal
                                      // interface mesh
    std::vector<NamePair> interfaces; // pair of names of sidesets that
                                      // interface
                                      // meshes will be created for
};

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
