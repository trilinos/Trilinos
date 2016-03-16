
#ifndef FACECREATOR_HPP_
#define FACECREATOR_HPP_

#include <vector>
#include "stk_mesh/base/Types.hpp"

namespace stk {
namespace mesh {

class BulkData;
class ElemElemGraph;
class SideSetEntry;
class sharing_info;

class FaceCreator {
public:
    FaceCreator(stk::mesh::BulkData& bulkData, stk::mesh::ElemElemGraph& elemElemGraph);
    ~FaceCreator();

    void create_side_entities_given_sideset(const std::vector<SideSetEntry> &skinnedSideSet, const stk::mesh::PartVector& skinParts);
private:
    void fill_side_ordinals(size_t element_side_index, const std::vector<SideSetEntry> &skinnedSideSet, std::vector<int>& ordinals);
    std::vector<int> get_side_ordinals_of_element(size_t element_side_index, const std::vector<SideSetEntry> &skinnedSideSet);
    size_t create_face_entities_per_element(size_t element_side_index, const std::vector<SideSetEntry> &skinnedSideSet, const stk::mesh::PartVector& skinParts, std::vector<stk::mesh::sharing_info> &sharedModified);
    stk::mesh::BulkData& m_bulkData;
    stk::mesh::ElemElemGraph& m_eeGraph;
};

}
}


#endif /* FACECREATOR_HPP_ */
