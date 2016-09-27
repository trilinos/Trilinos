#ifndef STK_BALANCE_UTILS
#define STK_BALANCE_UTILS

#include <mpi.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_topology/topology.hpp>
#include "stk_mesh/base/Field.hpp"  // for field_data

namespace stk
{
namespace balance
{
//rcb, multijagged, rib, hsfc, patoh, phg, metis, parmetis, parma, scotch, ptscotch, block, cyclic, random, zoltan, nd

typedef std::vector<int> ElementDecomposition;
typedef stk::mesh::Field<double> DoubleFieldType;

class DecompositionChangeList
{
public:
    DecompositionChangeList(stk::mesh::BulkData &stkMeshBulkData, const stk::mesh::EntityProcVec& decomposition);
    ~DecompositionChangeList();

    bool has_entity(stk::mesh::Entity entity);
    int  get_entity_destination(stk::mesh::Entity entity);
    void set_entity_destination(stk::mesh::Entity entity, const int destination);
    void delete_entity(stk::mesh::Entity entity);
    stk::mesh::EntityProcVec get_all_partition_changes();
    stk::mesh::BulkData & get_bulk();
    size_t get_num_global_entity_migrations() const;
    size_t get_max_global_entity_migrations() const;

    class Impl;
    Impl *pImpl = nullptr;
};


class BalanceSettings
{
public:
    BalanceSettings():
    fifteen(15),
    onePlus(1.05){}
    virtual ~BalanceSettings() {}

    enum GraphOption
    {
        LOADBALANCE = 0,
        COLORING
    };

    virtual size_t getNumNodesRequiredForConnection(stk::topology element1Topology, stk::topology element2Topology) const = 0;
    virtual double getGraphEdgeWeight(stk::topology element1Topology, stk::topology element2Topology) const = 0;
    virtual int getGraphVertexWeight(stk::topology type) const = 0;
    virtual double getGraphVertexWeight(stk::mesh::Entity entity, int criteria_index = 0) const = 0;
    virtual GraphOption getGraphOption() const = 0;


    // Graph (parmetis) based options only
    virtual bool includeSearchResultsInGraph() const { return false; }
    virtual double getToleranceForFaceSearch() const { return 0.0; }
    virtual double getToleranceForParticleSearch() const { return 0.0; }
    virtual double getGraphEdgeWeightForSearch() const { return 1.0; }
    virtual bool getEdgesForParticlesUsingSearch() const { return false; }
    virtual double getVertexWeightMultiplierForVertexInSearch() const { return fifteen; }

    virtual bool isIncrementalRebalance() const { return false; }
    virtual bool isMultiCriteriaRebalance() const { return false; }

    virtual bool userSpecifiedVertexWeights() const { return false; }
    virtual bool fieldSpecifiedVertexWeights() const { return false; }
    virtual std::vector<double> getVertexWeights() const { return std::vector<double>(); }

    virtual double getImbalanceTolerance() const { return onePlus; }

    virtual void setDecompMethod(const std::string& method) {}
    virtual std::string getDecompMethod() const { return std::string("parmetis"); }

    virtual std::string getCoordinateFieldName() const { return std::string("coordinates"); }

    virtual bool shouldPrintMetrics() const { return false; }

    virtual int getNumCriteria() const { return 1; }

    // Given an element/proc pair, can modify decomposition before elements are migrated
    virtual void modifyDecomposition(DecompositionChangeList & decomp) const {};

    // Need to implement getting particle radius from field
    virtual double getParticleRadius(stk::mesh::Entity particle) const { return 0.5; }

    // experimental
    virtual bool setVertexWeightsBasedOnNumberAdjacencies() const { return false; }
private:
    unsigned fifteen;
    double onePlus;
};

class ColoringSettings : public BalanceSettings
{
public:
    ColoringSettings() {}

    virtual size_t getNumNodesRequiredForConnection(stk::topology element1Topology, stk::topology element2Topology) const
    {
        return 1;
    }
    virtual double getGraphEdgeWeight(stk::topology element1Topology, stk::topology element2Topology) const
    {
        return 1.0;
    }
    virtual int getGraphVertexWeight(stk::topology type) const
    {
        return 1;
    }
    virtual GraphOption getGraphOption() const
    {
        return BalanceSettings::LOADBALANCE;
    }

    virtual double getGraphVertexWeight(stk::mesh::Entity entity, int criteria_index = 0) const
    {
        return 1.0;
    }
};

class OurColoringSettings : public ColoringSettings
{
public:
    OurColoringSettings() {}

    virtual GraphOption getGraphOption() const
    {
        return BalanceSettings::COLORING;
    }
};

class GraphCreationSettings : public BalanceSettings
{
public:
    GraphCreationSettings(): mToleranceForFaceSearch(0.0001),
                             mToleranceForParticleSearch( 0.0001),
                             edgeWeightForSearch (1000),
                             method(std::string("parmetis")),
                             vertexWeightMultiplierForVertexInSearch(6.0)
    {}

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

    virtual double getGraphEdgeWeightForSearch() const
    {
        return edgeWeightForSearch;
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

    virtual double getGraphVertexWeight(stk::mesh::Entity entity, int criteria_index = 0) const { return 1.0; }

    virtual int getGraphVertexWeight(stk::topology type) const
    {
        switch(type)
        {
            case stk::topology::PARTICLE:
            case stk::topology::LINE_2:
            case stk::topology::BEAM_2:
                return 1;
                break;
            case stk::topology::SHELL_TRIANGLE_3:
                return 3;
                break;
            case stk::topology::SHELL_TRIANGLE_6:
                return 6;
                break;
            case stk::topology::SHELL_QUADRILATERAL_4:
                return 6;
                break;
            case stk::topology::SHELL_QUADRILATERAL_8:
                return 12;
                break;
            case stk::topology::HEXAHEDRON_8:
                return 3;
                break;
            case stk::topology::HEXAHEDRON_20:
                return 12;
                break;
            case stk::topology::TETRAHEDRON_4:
                return 1;
                break;
            case stk::topology::TETRAHEDRON_10:
                return 3;
                break;
            case stk::topology::WEDGE_6:
                return 2;
                break;
            case stk::topology::WEDGE_15:
                return 12;
                break;
            default:
                if ( type.is_superelement( ))
                {
                    return 10;
                }
                throw("Invalid Element Type In WeightsOfElement");
                break;
        }
        return 0;
    }

    virtual GraphOption getGraphOption() const
    {
        return BalanceSettings::LOADBALANCE;
    }

    virtual bool includeSearchResultsInGraph() const { return true; }
    virtual double getToleranceForParticleSearch() const { return mToleranceForParticleSearch; }
    virtual double getToleranceForFaceSearch() const { return mToleranceForFaceSearch; }
    virtual bool getEdgesForParticlesUsingSearch() const { return true; }
    virtual double getVertexWeightMultiplierForVertexInSearch() const { return vertexWeightMultiplierForVertexInSearch; }
    virtual std::string getDecompMethod() const { return method; }

    virtual void setDecompMethod(const std::string& input_method) { method = input_method;}
    virtual void setToleranceForFaceSearch(double tol) { mToleranceForFaceSearch = tol; }
    virtual void setToleranceForParticleSearch(double tol) { mToleranceForParticleSearch = tol; }

private:
    int getConnectionTableIndex(stk::topology elementTopology) const
    {
        int tableIndex = -1;
        switch(elementTopology)
        {
            case stk::topology::PARTICLE:
                tableIndex = 0;
                break;
            case stk::topology::LINE_2:
            case stk::topology::LINE_2_1D:
            case stk::topology::LINE_3_1D:
            case stk::topology::BEAM_2:
            case stk::topology::BEAM_3:
            case stk::topology::SHELL_LINE_2:
            case stk::topology::SHELL_LINE_3:
                tableIndex = 1;
                break;
            case stk::topology::TRI_3_2D:
            case stk::topology::TRI_4_2D:
            case stk::topology::QUAD_4_2D:
            case stk::topology::SHELL_TRI_3:
            case stk::topology::SHELL_TRI_4:
            case stk::topology::SHELL_QUAD_4:
                tableIndex = 2;
                break;
            case stk::topology::TET_4:
            case stk::topology::PYRAMID_5:
            case stk::topology::WEDGE_6:
            case stk::topology::HEX_8:
                tableIndex = 3;
                break;
            case stk::topology::TRI_6_2D:
            case stk::topology::QUAD_8_2D:
            case stk::topology::QUAD_9_2D:
            case stk::topology::SHELL_TRI_6:
            case stk::topology::SHELL_QUAD_8:
            case stk::topology::SHELL_QUAD_9:
                tableIndex = 4;
                break;
            case stk::topology::TET_8:
            case stk::topology::TET_10:
            case stk::topology::TET_11:
            case stk::topology::PYRAMID_13:
            case stk::topology::PYRAMID_14:
            case stk::topology::WEDGE_15:
            case stk::topology::WEDGE_18:
            case stk::topology::HEX_20:
            case stk::topology::HEX_27:
                tableIndex = 5;
                break;
            default:
                if(elementTopology.is_superelement())
                {
                    tableIndex = 6;
                }
                else
                {
                    std::cerr << "Topology is " << elementTopology << std::endl;
                    throw("Invalid Element Type in GetDimOfElement");
                }
                break;
        };
        return tableIndex;
    }
    double mToleranceForFaceSearch;
    double mToleranceForParticleSearch;
    double edgeWeightForSearch;
    std::string method;
    double vertexWeightMultiplierForVertexInSearch;

};

class GraphCreationSettingsWithCustomTolerances : public GraphCreationSettings
{
public:
    GraphCreationSettingsWithCustomTolerances() : mToleranceForFaceSearch(0.1), mToleranceForParticleSearch(1.0) { }

    virtual double getToleranceForFaceSearch() const { return mToleranceForFaceSearch; }
    void setToleranceForFaceSearch(double tol) { mToleranceForFaceSearch = tol; }

    virtual double getToleranceForParticleSearch() const { return mToleranceForParticleSearch; }
    void setToleranceForParticleSearch(double tol)
    {
        mToleranceForParticleSearch = tol;
    }
    virtual bool getEdgesForParticlesUsingSearch() const { return true; }
    virtual bool setVertexWeightsBasedOnNumberAdjacencies() const { return true; }

private:
    double mToleranceForFaceSearch;
    double mToleranceForParticleSearch;
};

class BasicZoltan2Settings : public GraphCreationSettings
{
public:
    virtual bool includeSearchResultsInGraph() const { return false; }

    virtual double getToleranceForFaceSearch() const { return 0.0005 ; }
    virtual double getToleranceForParticleSearch() const { return 0.3; }
    virtual double getGraphEdgeWeightForSearch() const { return 100.0; }
    virtual bool getEdgesForParticlesUsingSearch() const { return true; }
    virtual double getVertexWeightMultiplierForVertexInSearch() const { return 6.0; }
    //virtual double getImbalanceTolerance() const { return 1.05; }
    virtual std::string getDecompMethod() const { return std::string("rcb"); }
};

class UserSpecifiedVertexWeightsSetting : public GraphCreationSettings
{
public:
    virtual double getGraphEdgeWeight(stk::topology element1Topology, stk::topology element2Topology) const { return 1.0; }
    virtual bool userSpecifiedVertexWeights() const { return true; }
    virtual bool includeSearchResultsInGraph() const { return false; }
    void setVertexWeights(const std::vector<double>& weights) { vertex_weights = weights; }
    virtual std::vector<double> getVertexWeights() const { return vertex_weights; }
    virtual int getGraphVertexWeight(stk::topology type) const { return 1; }
    //virtual double getImbalanceTolerance() const { return 1.05; }
    virtual void setDecompMethod(const std::string& input_method) { method = input_method;}
    virtual std::string getDecompMethod() const { return method; }
    void setCoordinateFieldName(const std::string& field_name) { m_field_name = field_name; }
    virtual std::string getCoordinateFieldName() const { return m_field_name; }

private:
    std::vector<double> vertex_weights;
    std::string method = std::string("parmetis");
    std::string m_field_name = std::string("coordinates");
};

class GraphCreationSettingsForZoltan2 : public GraphCreationSettingsWithCustomTolerances
{
public:
    virtual bool setVertexWeightsBasedOnNumberAdjacencies() const { return false; }
};

class FieldVertexWeightSettings : public GraphCreationSettings
{
public:
    FieldVertexWeightSettings(stk::mesh::BulkData &stkMeshBulkData,
                              const DoubleFieldType &weightField,
                              const double defaultWeight = 0.0)
      : m_stkMeshBulkData(stkMeshBulkData),
        m_weightField(weightField),
        m_defaultWeight(defaultWeight) { }
    virtual ~FieldVertexWeightSettings() = default;

    virtual double getGraphEdgeWeight(stk::topology element1Topology, stk::topology element2Topology) const { return 1.0; }
    virtual bool userSpecifiedVertexWeights() const { return false; }
    virtual bool fieldSpecifiedVertexWeights() const { return true; }
    virtual bool includeSearchResultsInGraph() const { return false; }
    virtual int getGraphVertexWeight(stk::topology type) const { return 1; }
    virtual double getImbalanceTolerance() const { return 1.05; }
    virtual void setDecompMethod(const std::string& input_method) { method = input_method;}
    virtual std::string getDecompMethod() const { return method; }

    virtual double getGraphVertexWeight(stk::mesh::Entity entity, int criteria_index = 0) const
    {
        const double *weight = stk::mesh::field_data(m_weightField, entity);
        if(weight) return *weight;

        return m_defaultWeight;
    }

protected:
    FieldVertexWeightSettings() = default;
    FieldVertexWeightSettings(const FieldVertexWeightSettings&) = delete;
    FieldVertexWeightSettings& operator=(const FieldVertexWeightSettings&) = delete;

    const stk::mesh::BulkData & m_stkMeshBulkData;
    const DoubleFieldType &m_weightField;
    const double m_defaultWeight;
    std::string method = std::string("parmetis");
};

class GraphEdge
{
public:
    GraphEdge(const stk::mesh::Entity element1, const stk::mesh::EntityId element2, int vertex2ProcOwner, double edgeWeight, bool isEdgeFromSearchArg = false) :
        mVertex1(element1), mVertex2(element2), mVertex2OwningProc(vertex2ProcOwner), mWeight(edgeWeight), mIsEdgeFromSearch(isEdgeFromSearchArg)
    {}

    ~GraphEdge() {}

    stk::mesh::Entity vertex1() const { return mVertex1; }
    stk::mesh::EntityId vertex2() const { return mVertex2; }
    int vertex2OwningProc() const { return mVertex2OwningProc; }
    double weight() const { return mWeight; }
    bool isEdgeFromSearch() const { return mIsEdgeFromSearch; }

private:
    stk::mesh::Entity mVertex1;
    stk::mesh::EntityId mVertex2;
    int mVertex2OwningProc;
    double mWeight;
    bool mIsEdgeFromSearch;
};

inline bool operator<(const GraphEdge &a, const GraphEdge &b)
{
    bool aLessB = (a.vertex1().m_value < b.vertex1().m_value);
    if(a.vertex1().m_value == b.vertex1().m_value)
    {
        aLessB = (a.vertex2() < b.vertex2());
        if(a.vertex2() == b.vertex2())
        {
            aLessB = (a.weight() < b.weight());
        }
    }
    return aLessB;
}

inline bool operator==(const GraphEdge &a, const GraphEdge &b)
{
    return (a.vertex1().m_value == b.vertex1().m_value) && (a.vertex2() == b.vertex2());
}

}
}
#endif
