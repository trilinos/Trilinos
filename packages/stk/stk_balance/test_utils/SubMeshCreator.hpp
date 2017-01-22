#ifndef SubMeshCreator_hpp
#define SubMeshCreator_hpp

namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Selector; } }

namespace stk {
namespace balance {
namespace internal {

void createNewSubMesh(stk::mesh::MetaData &oldMeta,
                      stk::mesh::BulkData &oldBulk,
                      stk::mesh::Selector subMeshSelector,
                      stk::mesh::MetaData &newMeta,
                      stk::mesh::BulkData &newBulk);

}
}
}
#endif
