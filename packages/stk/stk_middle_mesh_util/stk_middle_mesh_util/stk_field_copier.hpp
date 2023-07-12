#ifndef STK_MIDDLE_MESH_UTILS_STK_FIELD_COPIER
#define STK_MIDDLE_MESH_UTILS_STK_FIELD_COPIER

#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_middle_mesh/field.hpp"
#include "stk_middle_mesh/mesh.hpp"
#include "create_stk_mesh.hpp"

namespace stk {
namespace middle_mesh {
namespace stk_interface {

class StkFieldCopier
{
  public:
    using VertIdType = StkMeshCreator::VertIdType;

    StkFieldCopier(std::shared_ptr<stk::mesh::BulkData> bulkDataPtr, stk::mesh::Part* part,
                   std::shared_ptr<mesh::Mesh> middleMesh,
                   stk::mesh::Field<VertIdType>* stkVertField);

    mesh::FieldPtr<double> create_middle_mesh_field(const stk::mesh::Field<double>& stkField);

    stk::mesh::Field<double>* create_stk_field(mesh::FieldPtr<double> middleMeshField, const std::string& name);

    void copy(const stk::mesh::Field<double>& stkField, mesh::FieldPtr<double>& middleMeshFieldPtr);

    void copy(const mesh::FieldPtr<double> middleMeshFieldPtr, stk::mesh::Field<double>& stkField);


  private:
    void check_field_shapes(const stk::mesh::Field<double>& stkField, const mesh::FieldPtr<double> meshField);

    std::pair<mesh::FieldShape, int> get_field_shape_and_num_components(const stk::mesh::Field<double>& stkField);

    std::pair<int, int> get_field_shape(const stk::mesh::Field<double>& stkField);

    std::shared_ptr<::stk::mesh::BulkData> m_bulkDataPtr;
    stk::mesh::Part* m_part;
    std::shared_ptr<mesh::Mesh> m_middleMesh;
    stk::mesh::Field<VertIdType>* m_stkVertField;

};

}
}
}

#endif
