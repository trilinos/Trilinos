#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include "stk_mesh/base/Comm.hpp"
#include "stk_mesh/base/GetEntities.hpp"

namespace
{
//BEGINRcbSettings
class RcbSettings : public stk::balance::BalanceSettings
{
public:
    RcbSettings() {}
    virtual ~RcbSettings() {}

    virtual bool isIncrementalRebalance() const { return false; }
    virtual std::string getDecompMethod() const { return std::string("rcb"); }
    virtual std::string getCoordinateFieldName() const { return std::string("coordinates"); }
    virtual bool shouldPrintMetrics() const { return true; }
};
//ENDRcbSettings

class StkBalanceHowTo : public stk::unit_test_util::MeshFixture
{};

bool is_mesh_balanced(const stk::mesh::BulkData& bulk)
{
    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(bulk, counts);

    size_t numElementBalanced = counts[stk::topology::ELEM_RANK]/bulk.parallel_size();
    size_t numLocalElements = stk::mesh::count_selected_entities(bulk.mesh_meta_data().locally_owned_part(), bulk.buckets(stk::topology::ELEM_RANK));
    return numElementBalanced == numLocalElements;
}

bool is_mesh_balanced_wrt_weight(const stk::mesh::BulkData& bulk, stk::mesh::Field<double>& weightField)
{
    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(bulk, counts);

    size_t numElements = counts[stk::topology::ELEM_RANK];
    double totalWeight = (numElements*(numElements+1))/2;
    double balancedWeight = totalWeight/bulk.parallel_size();

    double weightThisProc = 0;
    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK, elements);
    for(stk::mesh::Entity element : elements)
    {
        if(bulk.bucket(element).owned())
        {
            double *data = stk::mesh::field_data(weightField, element);
            weightThisProc += *data;
        }
    }

    double tolerance = 20;
    return fabs(weightThisProc - balancedWeight)  < tolerance;
}

//BEGINBalanceTest1
TEST_F(StkBalanceHowTo, UseRebalanceWithGeometricMethods)
{
    if(stk::parallel_machine_size(get_comm()) == 2)
    {
        setup_mesh("generated:4x4x4|sideset:xX", stk::mesh::BulkData::NO_AUTO_AURA);

        RcbSettings balanceSettings;
        stk::balance::balanceStkMesh(balanceSettings, get_bulk());

        EXPECT_TRUE(is_mesh_balanced(get_bulk()));
    }
}
//ENDBalanceTest1


//BEGINParmetisSettings
class ParmetisSettings : public stk::balance::GraphCreationSettings
{
public:
    virtual std::string getDecompMethod() const { return "parmetis"; }

    size_t getNumNodesRequiredForConnection(stk::topology element1Topology, stk::topology element2Topology) const
    {
        const int noConnection = 1000;
        const int s = noConnection;
        const static int connectionTable[7][7] = {
            {1, 1, 1, 1, 1, 1, s}, // 0 dim
            {1, 1, 1, 1, 1, 1, s}, // 1 dim
            {1, 1, 2, 3, 2, 3, s}, // 2 dim linear
            {1, 1, 3, 3, 3, 3, s}, // 3 dim linear
            {1, 1, 2, 3, 3, 4, s}, // 2 dim higher-order
            {1, 1, 3, 3, 4, 4, s}, // 3 dim higher-order
            {s, s, s, s, s, s, s}  // super element
        };

        int element1Index = getConnectionTableIndex(element1Topology);
        int element2Index = getConnectionTableIndex(element2Topology);

        return connectionTable[element1Index][element2Index];
    }

    virtual double getGraphEdgeWeight(stk::topology element1Topology, stk::topology element2Topology) const
    {
        const double noConnection = 0;
        const double s = noConnection;
        const double largeWeight = 1000;
        const double L = largeWeight;
        const double twoDimWeight = 5;
        const double q = twoDimWeight;
        const double defaultWeight = 1.0;
        const double D = defaultWeight;
        const static double weightTable[7][7] = {
            {L, L, L, L, L, L, s}, // 0 dim
            {L, L, L, L, L, L, s}, // 1 dim
            {L, L, q, q, q, q, s}, // 2 dim linear
            {L, L, q, D, q, D, s}, // 3 dim linear
            {L, L, q, q, q, q, s}, // 2 dim higher-order
            {L, L, q, D, q, D, s}, // 3 dim higher-order
            {s, s, s, s, s, s, s}  // super element
        };

        int element1Index = getConnectionTableIndex(element1Topology);
        int element2Index = getConnectionTableIndex(element2Topology);

        return weightTable[element1Index][element2Index];
    }

    using BalanceSettings::getGraphVertexWeight;

    virtual int getGraphVertexWeight(stk::topology type) const
    {
        switch(type)
        {
            case stk::topology::PARTICLE:
            case stk::topology::LINE_2:
            case stk::topology::BEAM_2:
                return 1;
            case stk::topology::SHELL_TRIANGLE_3:
                return 3;
            case stk::topology::SHELL_TRIANGLE_6:
                return 6;
            case stk::topology::SHELL_QUADRILATERAL_4:
                return 6;
            case stk::topology::SHELL_QUADRILATERAL_8:
                return 12;
            case stk::topology::HEXAHEDRON_8:
                return 3;
            case stk::topology::HEXAHEDRON_20:
                return 12;
            case stk::topology::TETRAHEDRON_4:
                return 1;
            case stk::topology::TETRAHEDRON_10:
                return 3;
            case stk::topology::WEDGE_6:
                return 2;
            case stk::topology::WEDGE_15:
                return 12;
            default:
                if ( type.is_superelement( ))
                {
                    return 10;
                }
                throw("Invalid Element Type In WeightsOfElement");
        }
    }
};
//ENDParmetisSettings

//BEGINBalanceTest2
TEST_F(StkBalanceHowTo, UseRebalanceWithParmetis)
{
    if(stk::parallel_machine_size(get_comm()) == 2)
    {
        setup_mesh("generated:4x4x4|sideset:xX", stk::mesh::BulkData::NO_AUTO_AURA);

        ParmetisSettings balanceSettings;
        stk::balance::balanceStkMesh(balanceSettings, get_bulk());

        EXPECT_TRUE(is_mesh_balanced(get_bulk()));
    }
}
//ENDBalanceTest2

//BEGINParmeticSearchSettings
class ParmetisWithSearchSettings : public ParmetisSettings
{
    using ParmetisSettings::getToleranceForFaceSearch;
    virtual bool includeSearchResultsInGraph() const { return true; }
    virtual double getToleranceForFaceSearch() const { return 0.0001; }
    virtual double getVertexWeightMultiplierForVertexInSearch() const { return 6.0; }
    virtual double getGraphEdgeWeightForSearch() const { return 1000; }
};
//ENDParmeticSearchSettings


//BEGINBalanceTest3
TEST_F(StkBalanceHowTo, UseRebalanceWithParmetisAugmentedWithSearch)
{
    if(stk::parallel_machine_size(get_comm()) == 2)
    {
        setup_mesh("generated:4x4x4|sideset:xX", stk::mesh::BulkData::NO_AUTO_AURA);

        ParmetisWithSearchSettings balanceSettings;
        stk::balance::balanceStkMesh(balanceSettings, get_bulk());

        EXPECT_TRUE(is_mesh_balanced(get_bulk()));
    }
}
//ENDBalanceTest3

//BEGINRcbFieldSettings
class FieldVertexWeightSettings : public stk::balance::GraphCreationSettings
{
public:
    FieldVertexWeightSettings(const stk::balance::DoubleFieldType &weightField,
                              const double defaultWeight = 0.0)
  {
    setVertexWeightMethod(stk::balance::VertexWeightMethod::FIELD);
    setVertexWeightFieldName(weightField.name());
    setDefaultFieldWeight(defaultWeight);
  }

    virtual ~FieldVertexWeightSettings() = default;

    virtual double getGraphEdgeWeight(stk::topology element1Topology, stk::topology element2Topology) const { return 1.0; }

    virtual int getGraphVertexWeight(stk::topology type) const { return 1; }
    virtual double getImbalanceTolerance() const { return 1.0001; }
    virtual std::string getDecompMethod() const { return "rcb"; }

protected:
    FieldVertexWeightSettings() = delete;
    FieldVertexWeightSettings(const FieldVertexWeightSettings&) = delete;
    FieldVertexWeightSettings& operator=(const FieldVertexWeightSettings&) = delete;
};

//ENDRcbFieldSettings

void set_vertex_weights(const stk::mesh::BulkData& bulk, stk::mesh::Selector selector, stk::mesh::Field<double>& weightField)
{
    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK, selector, elements);
    for(stk::mesh::Entity element : elements)
    {
        if(bulk.bucket(element).owned())
        {
            double *data = stk::mesh::field_data(weightField, element);
            *data = static_cast<double>(bulk.identifier(element));
        }
    }
}

//BEGINBalanceTest4
TEST_F(StkBalanceHowTo, UseRebalanceWithFieldSpecifiedVertexWeights)
{
    if(stk::parallel_machine_size(get_comm()) == 2)
    {
        setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        stk::mesh::Field<double> &weightField = get_meta().declare_field<double>(stk::topology::ELEM_RANK, "vertex_weights");
        stk::mesh::put_field_on_mesh(weightField, get_meta().universal_part(), nullptr);
        stk::io::fill_mesh("generated:4x4x4|sideset:xX", get_bulk());
        set_vertex_weights(get_bulk(), get_meta().locally_owned_part(), weightField);

        FieldVertexWeightSettings balanceSettings(weightField);
        stk::balance::balanceStkMesh(balanceSettings, get_bulk());

        EXPECT_TRUE(is_mesh_balanced_wrt_weight(get_bulk(), weightField));
    }
}
//ENDBalanceTest4

TEST_F(StkBalanceHowTo, DISABLED_UseRebalanceWithFieldSpecifiedVertexWeightsOnLocallyOwnedPart)
{
    if(stk::parallel_machine_size(get_comm()) == 2)
    {
        setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        stk::mesh::Field<double> &weightField = get_meta().declare_field<double>(stk::topology::ELEM_RANK, "vertex_weights");
        stk::mesh::put_field_on_mesh(weightField, get_meta().locally_owned_part(), nullptr);
        stk::io::fill_mesh("generated:4x4x4|sideset:xX", get_bulk());
        set_vertex_weights(get_bulk(), get_meta().locally_owned_part(), weightField);

        FieldVertexWeightSettings balanceSettings(weightField);
        stk::balance::balanceStkMesh(balanceSettings, get_bulk());

        EXPECT_TRUE(is_mesh_balanced_wrt_weight(get_bulk(), weightField));
    }
}


//BEGINMultiCriteriaSelectorSettings
class MultipleCriteriaSelectorSettings : public ParmetisSettings
{
public:
    MultipleCriteriaSelectorSettings() { }
    virtual ~MultipleCriteriaSelectorSettings() = default;

    virtual bool isMultiCriteriaRebalance() const { return true;}

protected:
    MultipleCriteriaSelectorSettings(const MultipleCriteriaSelectorSettings&) = delete;
    MultipleCriteriaSelectorSettings& operator=(const MultipleCriteriaSelectorSettings&) = delete;
};
//ENDMultiCriteriaSelectorSettings

void put_elements_in_different_parts(stk::mesh::BulkData &bulk, stk::mesh::Part &part1, stk::mesh::Part &part2)
{
    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK, bulk.mesh_meta_data().locally_owned_part(), elements);
    bulk.modification_begin();
    for(stk::mesh::Entity element : elements)
    {
        stk::mesh::EntityId id = bulk.identifier(element);
        if(id%2 == 0)
            bulk.change_entity_parts(element, stk::mesh::ConstPartVector{&part1});
        else
            bulk.change_entity_parts(element, stk::mesh::ConstPartVector{&part2});
    }
    bulk.modification_end();
}

void verify_mesh_balanced_wrt_selectors(const stk::mesh::BulkData& bulk, const std::vector<stk::mesh::Selector> &selectors)
{
    std::vector<size_t> counts;
    for(const stk::mesh::Selector & sel : selectors)
    {
        stk::mesh::EntityVector elements;
        size_t num_elements = stk::mesh::count_selected_entities(sel, bulk.buckets(stk::topology::ELEM_RANK));
        counts.clear();
        stk::mesh::comm_mesh_counts(bulk, counts, &sel);

        size_t goldNumElements = counts[stk::topology::ELEM_RANK]/bulk.parallel_size();

        EXPECT_EQ(goldNumElements, num_elements);
    }
}

//BEGINBalanceTest5
TEST_F(StkBalanceHowTo, UseRebalanceWithMultipleCriteriaWithSelectors)
{
    if(stk::parallel_machine_size(get_comm()) == 2)
    {
        setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        stk::mesh::Part &part1 = get_meta().declare_part("madeup_part_1", stk::topology::ELEM_RANK);
        stk::mesh::Part &part2 = get_meta().declare_part("part_2", stk::topology::ELEM_RANK);
        stk::io::fill_mesh("generated:4x4x4|sideset:xX", get_bulk());

        put_elements_in_different_parts(get_bulk(), part1, part2);

        std::vector<stk::mesh::Selector> selectors = { part1, part2 };

        MultipleCriteriaSelectorSettings balanceSettings;
        stk::balance::balanceStkMesh(balanceSettings, get_bulk(), selectors);

        verify_mesh_balanced_wrt_selectors(get_bulk(), selectors);
    }
}
//ENDBalanceTest5

//BEGINMultiCriteriaFieldSettings
class MultipleCriteriaFieldSettings : public ParmetisSettings
{
public:
    MultipleCriteriaFieldSettings(const std::vector<stk::mesh::Field<double>*> critFields,
                                  const double default_weight = 0.0)
    {
      setNumCriteria(critFields.size());
      setVertexWeightMethod(stk::balance::VertexWeightMethod::FIELD);
      for (unsigned i = 0; i < critFields.size(); ++i) {
        setVertexWeightFieldName(critFields[i]->name(), i);
      }
      setDefaultFieldWeight(default_weight);
    }
    virtual ~MultipleCriteriaFieldSettings() override = default;

    virtual bool isMultiCriteriaRebalance() const override { return true;}

protected:
    MultipleCriteriaFieldSettings() = delete;
    MultipleCriteriaFieldSettings(const MultipleCriteriaFieldSettings&) = delete;
    MultipleCriteriaFieldSettings& operator=(const MultipleCriteriaFieldSettings&) = delete;
};
//ENDMultiCriteriaFieldSettings

void verify_mesh_balanced_wrt_fields(const stk::mesh::BulkData& bulk, const std::vector<stk::mesh::Field<double>*> &critFields)
{
    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK, bulk.mesh_meta_data().locally_owned_part(), elements);

    std::vector<double> sums(critFields.size(),0);

    for(size_t i=0;i<critFields.size();++i)
    {
        for(stk::mesh::Entity element : elements)
        {
            double *data = stk::mesh::field_data(*critFields[i], element);
            sums[i] += *data;
        }
    }

    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(bulk, counts);

    double gold_value = static_cast<double>(counts[stk::topology::ELEM_RANK]/bulk.parallel_size());
    double gold_value_per_field = static_cast<double>(gold_value/critFields.size());

    for(size_t i=0;i<sums.size();++i)
        EXPECT_EQ(gold_value_per_field, sums[i]);
}

void set_vertex_weights_checkerboard(stk::mesh::BulkData& bulk, stk::mesh::Selector selector, stk::mesh::Field<double> &weightField1, stk::mesh::Field<double> &weightField2)
{
    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK, selector, elements);
    for(stk::mesh::Entity element : elements)
    {
        double *data1 = stk::mesh::field_data(weightField1, element);
        double *data2 = stk::mesh::field_data(weightField2, element);
        stk::mesh::EntityId id = bulk.identifier(element);
        if(id%2==0)
        {
            *data1 = 1.0;
            *data2 = 0.0;
        }
        else
        {
            *data1 = 0.0;
            *data2 = 1.0;
        }
    }
}

//BEGINBalanceTest6
TEST_F(StkBalanceHowTo, UseRebalanceWithMultipleCriteriaWithFields)
{
    if(stk::parallel_machine_size(get_comm()) == 2)
    {
        setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        stk::mesh::Field<double> &weightField1 = get_meta().declare_field<double>(stk::topology::ELEM_RANK, "vertex_weights1");
        stk::mesh::put_field_on_mesh(weightField1, get_meta().universal_part(), nullptr);

        stk::mesh::Field<double> &weightField2 = get_meta().declare_field<double>(stk::topology::ELEM_RANK, "vertex_weights2");
        stk::mesh::put_field_on_mesh(weightField2, get_meta().universal_part(), nullptr);

        stk::io::fill_mesh("generated:4x4x4|sideset:xX", get_bulk());

        set_vertex_weights_checkerboard(get_bulk(), get_meta().locally_owned_part(), weightField1, weightField2);

        std::vector<stk::mesh::Field<double>*> critFields = { &weightField1, &weightField2 };
        MultipleCriteriaFieldSettings balanceSettings(critFields);
        stk::balance::balanceStkMesh(balanceSettings, get_bulk());

        verify_mesh_balanced_wrt_fields(get_bulk(), critFields);
    }
}
//ENDBalanceTest6
}
