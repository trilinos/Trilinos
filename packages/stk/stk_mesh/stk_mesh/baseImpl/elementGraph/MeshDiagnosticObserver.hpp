#ifndef STK_MESH_DIAGNOSTIC_OBSERVER_HPP
#define STK_MESH_DIAGNOSTIC_OBSERVER_HPP

#include <stk_mesh/base/ModificationObserver.hpp>
#include <unordered_set>

namespace stk { namespace mesh {

class BulkData;

enum MeshDiagnosticFlag
{
    NO_RULE   = 0     ,
    RULE_1    = 1L<<1 ,
    RULE_2    = 1L<<2 ,
    RULE_3    = 1L<<3 ,
    ALL_RULES = RULE_1 | RULE_2 | RULE_3,
    SOLO_SIDES = 1L<<4
};

class MeshDiagnosticObserver : public stk::mesh::ModificationObserver
{
public:
    MeshDiagnosticObserver(stk::mesh::BulkData& bulkData)
      : m_initialized(false),
        m_throwOnError(false),
        m_bulkData(bulkData),
        m_diagnosticFlag(NO_RULE) {}

    MeshDiagnosticObserver(stk::mesh::BulkData& bulkData, stk::mesh::MeshDiagnosticFlag diagnosticFlag)
      : m_initialized(false),
        m_throwOnError(false),
        m_bulkData(bulkData),
        m_diagnosticFlag(diagnosticFlag) {}

    virtual void finished_modification_end_notification();
    void enable_rule(stk::mesh::MeshDiagnosticFlag flag);
    unsigned get_number_of_errors() const {return m_errorDatabase.size();}
    void throw_if_errors_exist() {m_throwOnError = true;}
private:
    bool m_initialized;
    bool m_throwOnError;
    stk::mesh::BulkData &m_bulkData;
    std::unordered_set<std::string> m_errorDatabase;
    stk::mesh::MeshDiagnosticFlag m_diagnosticFlag;

    void gather_new_errors(std::ofstream & out, const std::vector<std::string> & errorList);
};

} }

#endif
