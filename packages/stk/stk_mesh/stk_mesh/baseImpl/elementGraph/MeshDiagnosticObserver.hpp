#ifndef STK_MESH_DIAGNOSTIC_OBSERVER_HPP
#define STK_MESH_DIAGNOSTIC_OBSERVER_HPP

#include <stk_mesh/base/ModificationObserver.hpp>

namespace stk { namespace mesh {

class BulkData;

class MeshDiagnosticObserver : public stk::mesh::ModificationObserver
{
public:
    MeshDiagnosticObserver(stk::mesh::BulkData& bulkData) : mBulkData(bulkData){}
    virtual void finished_modification_end_notification();

private:
    stk::mesh::BulkData &mBulkData;
};

} }

#endif
