
#ifndef FACECREATOR_HPP_
#define FACECREATOR_HPP_

#include <vector>

namespace stk {
namespace mesh {

class BulkData;
class ElemElemGraph;

class FaceCreator {
public:
    FaceCreator(stk::mesh::BulkData& bulkData, stk::mesh::ElemElemGraph& elemElemGraph);
    ~FaceCreator();

    void fill_side_ordinals(size_t element_side_index, const std::vector<SideSetEntry> &skinnedSideSet, std::vector<int>& ordinals);
    std::vector<int> get_side_ordinals_of_element(size_t element_side_index, const std::vector<SideSetEntry> &skinnedSideSet);
    size_t create_face_entities_per_element(size_t element_side_index, const std::vector<SideSetEntry> &skinnedSideSet, const stk::mesh::PartVector& skinParts, std::vector<stk::mesh::sharing_info> &sharedModified);
    void create_side_entities_given_sideset(const std::vector<SideSetEntry> &skinnedSideSet, const stk::mesh::PartVector& skinParts, std::vector<stk::mesh::sharing_info>& sharedModified);
private:
    stk::mesh::BulkData& m_bulkData;
    stk::mesh::ElemElemGraph& m_eeGraph;
};

}
}


#endif /* FACECREATOR_HPP_ */
