#include "MeshDiagnosticObserver.hpp"
#include <map>
#include <utility>
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MeshDiagnostics.hpp"
#include <cstdio>

namespace stk { namespace mesh {

void MeshDiagnosticObserver::gather_new_errors(std::ofstream & out, const std::vector<std::string> & errorList)
{
    for (const std::string & error : errorList) {
        std::pair<std::unordered_set<std::string>::iterator, bool> result = m_errorDatabase.insert(error);
        if (result.second) {
            out << error;

            if(m_throwOnError)
                STK_ThrowRequireMsg( false,
                                 "Mesh diagnostic rule failure: " + error );
        }
    }
}

void MeshDiagnosticObserver::finished_modification_end_notification()
{
    if(m_diagnosticFlag)
    {
        int my_proc_id = m_bulkData.parallel_rank();
        const std::string outFileName = "mesh_diagnostics_failures_" + std::to_string(my_proc_id) + ".txt";
        if (!m_initialized) {
            std::remove(outFileName.c_str());
            m_initialized = true;
        }
        std::ofstream out(outFileName, std::ios::app);
        out << "Reason for mesh modification: " << m_bulkData.get_last_modification_description() << std::endl;

        if(m_diagnosticFlag & RULE_1)
        {
            stk::mesh::SplitCoincidentInfo splitCoincidentElements = stk::mesh::get_split_coincident_elements(m_bulkData);
            std::vector<std::string> splitCoincidentErrors = stk::mesh::get_messages_for_split_coincident_elements(m_bulkData, splitCoincidentElements);
            gather_new_errors(out, splitCoincidentErrors);
        }

        if(m_diagnosticFlag & RULE_2)
        {
            std::vector<stk::mesh::EntityKeyProc> badKeyProcs = stk::mesh::get_non_unique_key_procs(m_bulkData);
            std::vector<std::string> nonUniqueKeyErrors = stk::mesh::get_non_unique_key_messages(m_bulkData, badKeyProcs);
            gather_new_errors(out, nonUniqueKeyErrors);
        }

        if(m_diagnosticFlag & RULE_3)
        {
            std::vector<stk::mesh::Entity> orphanedSides = stk::mesh::get_orphaned_sides_with_attached_element_on_different_proc(m_bulkData);
            std::vector<std::string> orphanedSideErrors = stk::mesh::get_messages_for_orphaned_owned_sides(m_bulkData, orphanedSides);
            gather_new_errors(out, orphanedSideErrors);
        }

        if(m_diagnosticFlag & SOLO_SIDES)
        {
            std::vector<stk::mesh::Entity> orphanedSides = stk::mesh::get_solo_sides_without_element_on_different_proc(m_bulkData);
            std::vector<std::string> orphanedSideErrors = stk::mesh::get_messages_for_solo_sides(m_bulkData, orphanedSides);
            gather_new_errors(out, orphanedSideErrors);
        }

        out.close();
    }

//    throw_if_any_proc_has_false(m_bulkData.parallel(), badKeyProcs.empty() && splitCoincidentElements.empty() && orphanedSides.empty());
}

void MeshDiagnosticObserver::enable_rule(stk::mesh::MeshDiagnosticFlag flag)
{
    m_diagnosticFlag = static_cast<stk::mesh::MeshDiagnosticFlag> (m_diagnosticFlag | flag);
}

}}
