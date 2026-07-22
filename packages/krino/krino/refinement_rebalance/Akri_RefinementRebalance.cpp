#include <Akri_RefinementRebalance.hpp>

#include <Akri_Refinement.hpp>
#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>

namespace krino {

class RefinementRebalance : public stk::balance::GraphCreationSettings
{
public:
    RefinementRebalance(const krino::Refinement & refinement)
        : mRefinement(refinement)
    {
      setVertexWeightMethod(stk::balance::VertexWeightMethod::FIELD);
    }

    ~RefinementRebalance() = default;

    void modifyDecomposition(stk::balance::DecompositionChangeList & decomp_changes) const override;
    bool shouldPrintMetrics() const override { return true; }
    virtual double getFieldVertexWeight(const stk::mesh::BulkData &/*bulkData*/, stk::mesh::Entity entity, int /*criteria_index*/) const override
    {
      return 1.0 * mRefinement.rebalance_element_count_incorporating_parallel_owner_constraints(entity);
    }

private:
    void remove_decomp_changes_for_entities_with_parallel_ownership_constraints(stk::balance::DecompositionChangeList & decomp_changes) const;
    void add_parallel_ownership_constraints_into_decomp_changes(stk::balance::DecompositionChangeList & decomp_changes) const;
    const krino::Refinement & mRefinement;
};

void RefinementRebalance::remove_decomp_changes_for_entities_with_parallel_ownership_constraints(stk::balance::DecompositionChangeList & decompChanges) const
{
    for(auto && change : decompChanges.get_all_partition_changes())
    {
        stk::mesh::Entity entity = change.first;
        if(mRefinement.has_parallel_owner_rebalance_constraint(entity))
            decompChanges.delete_entity(entity);
    }
}

void RefinementRebalance::add_parallel_ownership_constraints_into_decomp_changes(stk::balance::DecompositionChangeList & decompChanges) const
{
    stk::mesh::EntityVector adaptChildren;
    for(auto && change : decompChanges.get_all_partition_changes())
    {
        stk::mesh::Entity entity = change.first;
        const auto dest = change.second;

        mRefinement.fill_child_elements_that_must_stay_on_same_proc_as_parent(entity, adaptChildren);
        for(auto && child : adaptChildren)
            decompChanges.set_entity_destination(child, dest);
    }
}

void RefinementRebalance::modifyDecomposition(stk::balance::DecompositionChangeList & decompChanges) const
{
    remove_decomp_changes_for_entities_with_parallel_ownership_constraints(decompChanges);
    add_parallel_ownership_constraints_into_decomp_changes(decompChanges);
}

bool rebalance_refined_mesh(const Refinement & refinement, stk::mesh::BulkData & mesh)
{
    const RefinementRebalance balancer(refinement);
    return stk::balance::balanceStkMesh(balancer, mesh);
}

}
