#include <stk_mesh/base/GetEntities.hpp>

#include <stk_unit_test_utils/unittestMeshUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>

#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>
#include <stk_balance/internal/GeometricVertices.hpp>
#include <stk_balance/internal/ZoltanGeometricAdapter.hpp>
#include <stk_balance/internal/StkGeometricMethodViaZoltan.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

namespace
{

class GeometricBalanceSettingsTester : public stk::balance::GraphCreationSettings
{
public:
    GeometricBalanceSettingsTester(const std::string& decompMethod)
    : method(decompMethod) { }
    virtual ~GeometricBalanceSettingsTester() = default;

    virtual std::string getDecompMethod() const { return method; }

private:
    const std::string& method;
};

class ZoltanGeometricMethods : public stk::unit_test_util::MeshFixture
{
protected:

    std::vector<double> get_weights_per_proc(const std::vector<unsigned> &processorOntoWhichEntityBelongs,
                                             const stk::balance::internal::GeometricVertices& vertexInfo,
                                             const stk::mesh::EntityVector &entitiesToBalance)
    {
        std::vector<double> weightsPerProc(get_bulk().parallel_size(),0);
        const std::vector<double> &vertWeights = vertexInfo.get_vertex_weights();

        for(size_t i=0;i<entitiesToBalance.size();++i)
            weightsPerProc[processorOntoWhichEntityBelongs[i]] += vertWeights[i];
        return weightsPerProc;
    }

    std::vector<double> get_summed_weights_per_proc(const std::vector<unsigned> &processorOntoWhichEntityBelongs, const stk::balance::internal::GeometricVertices& vertexInfo, const stk::mesh::EntityVector &entitiesToBalance)
    {
        std::vector<double> weightsPerProc = get_weights_per_proc(processorOntoWhichEntityBelongs, vertexInfo, entitiesToBalance);
        std::vector<double> sumWeightsPerProc(get_bulk().parallel_size(),0);
        stk::all_reduce_sum(get_bulk().parallel(), weightsPerProc.data(), sumWeightsPerProc.data(), get_bulk().parallel_size());
        return sumWeightsPerProc;
    }

    double  get_total_weights_per_proc(const std::vector<double> &sumWeightsPerProc)
    {
        double sum = 0;
        for(size_t i=0;i<sumWeightsPerProc.size();++i)
            sum += sumWeightsPerProc[i];
        return sum;
    }

    void test_decomposition(const std::vector<unsigned> &processorOntoWhichEntityBelongs, const stk::balance::internal::GeometricVertices& vertexInfo, const stk::mesh::EntityVector &entitiesToBalance)
    {
        std::vector<double> sumWeightsPerProc = get_summed_weights_per_proc(processorOntoWhichEntityBelongs, vertexInfo, entitiesToBalance);
        double totalWeightAcrossProcs = get_total_weights_per_proc(sumWeightsPerProc);
        double averagePerProc = totalWeightAcrossProcs/get_bulk().parallel_size();
        EXPECT_EQ(sumWeightsPerProc[get_bulk().parallel_rank()], averagePerProc);
    }

    void run_decomp_with_method(const std::string& method, stk::mesh::Selector selector, stk::mesh::EntityRank primary_rank)
    {
        stk::mesh::EntityVector entitiesToBalance = stk::balance::internal::get_entities_to_balance(selector, primary_rank, get_bulk());
        stk::balance::GraphCreationSettings balanceSettings;
        stk::balance::internal::GeometricVertices vertexInfo(balanceSettings, get_bulk(), entitiesToBalance, {selector});
        for(size_t i=0;i<entitiesToBalance.size();++i)
        {
            vertexInfo.set_vertex_weight(i, get_bulk().identifier(entitiesToBalance[i]));
        }
        GeometricBalanceSettingsTester geometricBalanceSettings(method);
        std::vector<unsigned> processorOntoWhichEntityBelongs = stk::balance::get_decomposition(vertexInfo, geometricBalanceSettings, get_bulk().parallel_size(), get_bulk().parallel());
        test_decomposition(processorOntoWhichEntityBelongs, vertexInfo, entitiesToBalance);
    }

    void run_geometric_decomp_tests_for_elements_1_3_4()
    {
        stk::mesh::Selector selector = *part;
        const std::vector<std::string> geometric_decomposition_methods = { "rcb", "rib", "block", "multijagged" };
        for(const std::string geometric_method : geometric_decomposition_methods)
          run_decomp_with_method(geometric_method, selector, stk::topology::ELEMENT_RANK);
    }

    void create_decomposition_part()
    {
        part = &get_meta().declare_part("balance_this");
    }

    void move_element_1_3_4_into_decomposition_part()
    {
        std::vector<unsigned> ids = { 1, 3, 4 };
        get_bulk().modification_begin();
        change_parts_on_entities_with_these_ids(ids);
        get_bulk().modification_end();
    }

    void change_parts_on_entities_with_these_ids(const std::vector<unsigned>& ids)
    {
        for(unsigned id : ids)
        {
            stk::mesh::Entity element = get_bulk().get_entity(stk::topology::ELEM_RANK, id);
            if(get_bulk().is_valid(element) && get_bulk().bucket(element).owned())
                get_bulk().change_entity_parts(element, stk::mesh::ConstPartVector{part});
        }
    }

    void set_up_1x1x6_mesh_with_elements_1_3_4_in_decomposition_part()
    {
        setup_mesh("generated:1x1x6", stk::mesh::BulkData::AUTO_AURA);
        create_decomposition_part();
        move_element_1_3_4_into_decomposition_part();
    }

private:
    stk::mesh::Part* part = nullptr;
};


TEST_F(ZoltanGeometricMethods, balanceOnlySelectedElements)
{
    if(stk::parallel_machine_size(get_comm())==2)
    {
        set_up_1x1x6_mesh_with_elements_1_3_4_in_decomposition_part();
        run_geometric_decomp_tests_for_elements_1_3_4();
    }
}

}
