#ifndef STK_MESH_DIAGNOSTIC_OBSERVER_HPP
#define STK_MESH_DIAGNOSTIC_OBSERVER_HPP

#include <stk_mesh/base/ModificationObserver.hpp>
#include <unordered_set>

namespace stk { namespace mesh {

class BulkData;

class MeshDiagnosticObserver : public stk::mesh::ModificationObserver
{
public:
    MeshDiagnosticObserver(stk::mesh::BulkData& bulkData)
      : m_initialized(false),
        m_bulkData(bulkData) {}

    virtual void finished_modification_end_notification();
    void gather_new_errors(std::ofstream & out, const std::vector<std::string> & errorList);

private:
    bool m_initialized;
    stk::mesh::BulkData &m_bulkData;
    std::unordered_set<std::string> m_errorDatabase;
};

} }

#endif
