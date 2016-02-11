#include <stddef.h>                     // for size_t, nullptr
#include <string>                       // for string
#include <stk_util/parallel/ParallelReduce.hpp>
#include <map>
#include <string>
#include "../../base/BulkData.hpp"
#include "../../base/MetaData.hpp"
#include "ElemElemGraph.hpp"
#include "MeshDiagnosticObserver.hpp"
#include "../EquivalentEntityBlocks.hpp"

namespace stk { namespace mesh {

class ElemGraphForDiagnostics : public stk::mesh::ElemElemGraph
{
public:
    ElemGraphForDiagnostics(stk::mesh::BulkData& bulkData, const stk::mesh::Selector &selector, const stk::mesh::Selector *air = nullptr) :
        stk::mesh::ElemElemGraph(bulkData, selector, air) { }

    stk::mesh::impl::ParallelGraphInfo& get_parallel_info() { return m_parallelInfoForGraphEdges.get_parallel_graph_info(); }
    stk::mesh::Entity get_entity(stk::mesh::impl::LocalId local_id) const { return m_local_id_to_element_entity[local_id]; }
    const stk::mesh::impl::SparseGraph& get_coincident_graph() const { return m_coincidentGraph; }
};

std::map<stk::mesh::EntityId, std::pair<stk::mesh::EntityId, int> > get_split_coincident_elements(stk::mesh::BulkData& bulkData)
{
    stk::mesh::Selector sel = bulkData.mesh_meta_data().locally_owned_part();
    ElemGraphForDiagnostics graph(bulkData, sel);
    const stk::mesh::impl::SparseGraph& coingraph = graph.get_coincident_graph();

    std::map<stk::mesh::EntityId, std::pair<stk::mesh::EntityId, int> > badElements;

    for(const stk::mesh::impl::SparseGraph::value_type& extractedEdgesForElem : coingraph)
    {
        //const stk::mesh::impl::LocalId possibleMCE = extractedEdgesForElem.first;
        const std::vector<stk::mesh::GraphEdge>& coincidentEdgesForElem = extractedEdgesForElem.second;
        for(const stk::mesh::GraphEdge& edge : coincidentEdgesForElem)
        {
            if(edge.elem2 < 0)
            {
                stk::mesh::Entity entity = graph.get_entity(edge.elem1);
                stk::mesh::EntityId id = bulkData.identifier(entity);
                stk::mesh::impl::ParallelGraphInfo& par_info = graph.get_parallel_info();
                stk::mesh::impl::ParallelGraphInfo::iterator iter = par_info.find(edge);
                ThrowRequireMsg(iter!=par_info.end(), "Program error. Contact sierra-help@sandia.gov for support.");
                badElements[id] = std::make_pair(-edge.elem2, iter->second.m_other_proc);
            }
        }
    }
    return badElements;
}

void write_mesh_diagnostics(const stk::mesh::BulkData& bulkData, const std::map<stk::mesh::EntityId, std::pair<stk::mesh::EntityId, int> > & splitCoincidentElements)
{
    int my_proc_id = bulkData.parallel_rank();
    std::ofstream out("mesh_diagnostics_failures_" + std::to_string(my_proc_id) + ".txt");
    for(const auto& item : splitCoincidentElements) {
        stk::mesh::Entity element = bulkData.get_entity(stk::topology::ELEM_RANK,item.first);
        const stk::mesh::PartVector& elementParts = bulkData.bucket(element).supersets();
        std::string blockNames;
        blockNames = "{";
        for (const stk::mesh::Part* part : elementParts) {
            if (stk::mesh::impl::is_element_block(*part)) {
                blockNames += " " + part->name();
            }
        }
        blockNames += " }";
        out << "[" << my_proc_id << "] Element " << item.first << " (" << bulkData.bucket(element).topology() << ") in blocks " << blockNames << " is coincident with element " << item.second.first << " on processor " << item.second.second << std::endl;
    }
    out.close();
}

void print_and_throw_if_elements_are_split(const stk::mesh::BulkData& bulkData, const std::map<stk::mesh::EntityId, std::pair<stk::mesh::EntityId, int> > &splitCoincidentElements)
{
    bool is_all_ok_locally = splitCoincidentElements.empty();
    bool is_all_ok_globally = stk::is_true_on_all_procs(bulkData.parallel(), is_all_ok_locally);
    if(!is_all_ok_locally)
        write_mesh_diagnostics(bulkData, splitCoincidentElements);
    ThrowRequireMsg(is_all_ok_globally, "Mesh diagnostics failed.");
}

} }
