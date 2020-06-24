// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 // 
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 // 
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 // 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#ifndef STK_BALANCE_UTILS
#define STK_BALANCE_UTILS

#include "mpi.h"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_topology/topology.hpp>
#include "stk_mesh/base/Field.hpp"  // for field_data
#include "stk_mesh/base/FieldBase.hpp"
#include "stk_balance/setup/DefaultSettings.hpp"
#include <memory>

namespace stk
{
namespace balance
{
//rcb, multijagged, rib, hsfc, patoh, phg, metis, parmetis, parma, scotch, ptscotch, block, cyclic, random, zoltan, nd

class FaceSearchTolerance;

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
    BalanceSettings();
    virtual ~BalanceSettings() {}

    enum GraphOption
    {
        LOAD_BALANCE = 0,
        COLOR_MESH,
        COLOR_MESH_BY_TOPOLOGY,
        COLOR_MESH_AND_OUTPUT_COLOR_FIELDS,
    };

    virtual size_t getNumNodesRequiredForConnection(stk::topology element1Topology, stk::topology element2Topology) const;
    virtual double getGraphEdgeWeight(stk::topology element1Topology, stk::topology element2Topology) const;
    virtual int getGraphVertexWeight(stk::topology type) const;
    virtual double getGraphVertexWeight(stk::mesh::Entity entity, int criteria_index = 0) const ;
    virtual GraphOption getGraphOption() const;

    // Graph (parmetis) based options only
    virtual bool includeSearchResultsInGraph() const;
    virtual void setIncludeSearchResultsInGraph(bool doContactSearch);

    virtual double getToleranceForFaceSearch(const stk::mesh::BulkData & mesh,
                                             const stk::mesh::FieldBase & coordField,
                                             const stk::mesh::Entity * faceNodes,
                                             const unsigned numFaceNodes) const;
    virtual void setToleranceFunctionForFaceSearch(std::shared_ptr<stk::balance::FaceSearchTolerance> faceSearchTolerance);
    virtual bool isConstantFaceSearchTolerance() const;

    virtual double getToleranceForParticleSearch() const;
    virtual double getAbsoluteToleranceForParticleSearch(stk::mesh::Entity particle) const;
    virtual double getGraphEdgeWeightForSearch() const;
    virtual bool getEdgesForParticlesUsingSearch() const;
    virtual double getVertexWeightMultiplierForVertexInSearch() const;

    virtual bool isIncrementalRebalance() const;
    virtual bool isMultiCriteriaRebalance() const;

    virtual bool areVertexWeightsProvidedInAVector() const;
    virtual std::vector<double> getVertexWeightsViaVector() const;
    virtual bool areVertexWeightsProvidedViaFields() const;

    virtual double getImbalanceTolerance() const;

    virtual void setDecompMethod(const std::string& method) ;
    virtual std::string getDecompMethod() const ;

    virtual std::string getInitialDecompMethod() const ;

    virtual std::string getCoordinateFieldName() const ;

    virtual bool shouldPrintMetrics() const;

    virtual int getNumCriteria() const;

    // Given an element/proc pair, can modify decomposition before elements are migrated
    virtual void modifyDecomposition(DecompositionChangeList & decomp) const ;

    // Need to implement getting particle radius from field
    virtual double getParticleRadius(stk::mesh::Entity particle) const ;

    // experimental
    virtual bool setVertexWeightsBasedOnNumberAdjacencies() const;

    virtual bool allowModificationOfVertexWeightsForSmallMeshes() const;

    virtual bool shouldFixMechanisms() const;
    virtual bool shouldFixSpiders() const;
    virtual std::string getSpiderBeamConnectivityCountFieldName() const;
    virtual std::string getSpiderVolumeConnectivityCountFieldName() const;
    virtual const stk::mesh::Field<int> * getSpiderBeamConnectivityCountField(const stk::mesh::BulkData & stkMeshBulkData) const;
    virtual const stk::mesh::Field<int> * getSpiderVolumeConnectivityCountField(const stk::mesh::BulkData & stkMeshBulkData) const;
    virtual bool useLocalIds() const;

    virtual bool useNodeBalancer() const;
    virtual double getNodeBalancerTargetLoadBalance() const;
    virtual unsigned getNodeBalancerMaxIterations() const;

    virtual void set_input_filename(const std::string& filename);
    virtual std::string get_input_filename() const;

    virtual void set_output_filename(const std::string& filename);
    virtual std::string get_output_filename() const;

    virtual void setShouldFixSpiders(bool fixSpiders) { }
    virtual void setEdgeWeightForSearch(double w) { }
    virtual void setVertexWeightMultiplierForVertexInSearch(double w) { }
    virtual void setToleranceForFaceSearch(double tol) { }

protected:
    const DefaultSettings defaults;

private:
    std::string m_initialDecompMethod;
    std::string m_inputFilename;
    std::string m_outputFilename;
};

class BasicGeometricSettings : public BalanceSettings
{
public:
    virtual std::string getDecompMethod() const { return "rcb"; }
};

class GraphCreationSettings : public BalanceSettings
{
public:
    GraphCreationSettings()
      : mToleranceForFaceSearch(defaults.faceSearchAbsTol),
        mToleranceForParticleSearch(defaults.particleSearchTol),
        edgeWeightForSearch(defaults.faceSearchEdgeWeight),
        method(defaults.decompMethod),
        vertexWeightMultiplierForVertexInSearch(defaults.faceSearchVertexMultiplier),
        m_UseConstantToleranceForFaceSearch(true),
        m_shouldFixSpiders(false),
        m_spiderBeamConnectivityCountField(nullptr),
        m_spiderVolumeConnectivityCountField(nullptr),
        m_includeSearchResultInGraph(defaults.useContactSearch),
        m_useNodeBalancer(false),
        m_nodeBalancerTargetLoadBalance(1.0),
        m_nodeBalancerMaxIterations(5)
    {}

    GraphCreationSettings(double faceSearchTol, double particleSearchTol, double edgeWeightSearch, const std::string& decompMethod, double multiplierVWSearch)
      : GraphCreationSettings()
    {
      mToleranceForFaceSearch = faceSearchTol;
      mToleranceForParticleSearch = particleSearchTol;
      edgeWeightForSearch = edgeWeightSearch;
      method = decompMethod;
      vertexWeightMultiplierForVertexInSearch = multiplierVWSearch;
    }

    virtual ~GraphCreationSettings() = default;

    size_t getNumNodesRequiredForConnection(stk::topology element1Topology, stk::topology element2Topology) const override;

    virtual double getGraphEdgeWeightForSearch() const override;

    virtual double getGraphEdgeWeight(stk::topology element1Topology, stk::topology element2Topology) const override;

    virtual double getGraphVertexWeight(stk::mesh::Entity entity, int criteria_index = 0) const override;

    virtual int getGraphVertexWeight(stk::topology type) const override;

    virtual GraphOption getGraphOption() const override;
    virtual bool includeSearchResultsInGraph() const override;
    virtual void setIncludeSearchResultsInGraph(bool doContactSearch) override;
    virtual double getToleranceForParticleSearch() const override;

    virtual double getToleranceForFaceSearch(const stk::mesh::BulkData & mesh,
                                             const stk::mesh::FieldBase & coordField,
                                             const stk::mesh::Entity * faceNodes,
                                             const unsigned numFaceNodes) const override;
    virtual void setToleranceFunctionForFaceSearch(std::shared_ptr<stk::balance::FaceSearchTolerance> faceSearchTolerance) override;
    virtual bool isConstantFaceSearchTolerance() const override;

    virtual bool getEdgesForParticlesUsingSearch() const override;
    virtual double getVertexWeightMultiplierForVertexInSearch() const override;
    virtual std::string getDecompMethod() const override;

    virtual void setDecompMethod(const std::string& input_method) override;
    virtual void setToleranceForFaceSearch(double tol);
    virtual void setToleranceForParticleSearch(double tol);
    virtual void setEdgeWeightForSearch(double w);
    virtual void setVertexWeightMultiplierForVertexInSearch(double w);
    virtual void setShouldFixSpiders(bool fixSpiders);

    virtual bool shouldFixMechanisms() const override;
    virtual bool shouldFixSpiders() const override;
    virtual const stk::mesh::Field<int> * getSpiderBeamConnectivityCountField(const stk::mesh::BulkData & stkMeshBulkData) const override;
    virtual const stk::mesh::Field<int> * getSpiderVolumeConnectivityCountField(const stk::mesh::BulkData & stkMeshBulkData) const override;

    virtual void setUseNodeBalancer(bool useBalancer);
    virtual void setNodeBalancerTargetLoadBalance(double targetLoadBalance);
    virtual void setNodeBalancerMaxIterations(unsigned maxIterations);

    virtual bool useNodeBalancer() const override;
    virtual double getNodeBalancerTargetLoadBalance() const override;
    virtual unsigned getNodeBalancerMaxIterations() const override;

protected:
    int getConnectionTableIndex(stk::topology elementTopology) const;
    int getEdgeWeightTableIndex(stk::topology elementTopology) const;
    double mToleranceForFaceSearch;
    double mToleranceForParticleSearch;
    double edgeWeightForSearch;
    std::string method;
    double vertexWeightMultiplierForVertexInSearch;
    bool m_UseConstantToleranceForFaceSearch;
    bool m_shouldFixSpiders;
    mutable const stk::mesh::Field<int> * m_spiderBeamConnectivityCountField;
    mutable const stk::mesh::Field<int> * m_spiderVolumeConnectivityCountField;
    std::shared_ptr<stk::balance::FaceSearchTolerance> m_faceSearchToleranceFunction;
    bool m_includeSearchResultInGraph;
    bool m_useNodeBalancer;
    double m_nodeBalancerTargetLoadBalance;
    unsigned m_nodeBalancerMaxIterations;
};

class GraphCreationSettingsWithCustomTolerances : public GraphCreationSettings
{
public:
    GraphCreationSettingsWithCustomTolerances()
      : GraphCreationSettings()
    {
        mToleranceForFaceSearch = 0.1;
        mToleranceForParticleSearch = 1.0;
    }

    virtual bool getEdgesForParticlesUsingSearch() const { return true; }
    virtual bool setVertexWeightsBasedOnNumberAdjacencies() const { return true; }
};

class StkBalanceSettings : public GraphCreationSettings
{
public:
    StkBalanceSettings()
      : GraphCreationSettings()
    {}

    ~StkBalanceSettings() = default;

    std::string getCoordinateFieldName() const override {
      return "balance_coordinates";
    }
};

class BasicZoltan2Settings : public GraphCreationSettings
{
public:
    BasicZoltan2Settings()
      : GraphCreationSettings(0.0005, 0.3, 100.0, "rcb", 6.0) 
      {
        m_includeSearchResultInGraph = false;
      }
    virtual bool getEdgesForParticlesUsingSearch() const { return true; }
};

class UserSpecifiedVertexWeightsSetting : public GraphCreationSettings
{
public:
    UserSpecifiedVertexWeightsSetting()
    {
      method = "parmetis";
      m_includeSearchResultInGraph = false;
    }
    virtual double getGraphEdgeWeight(stk::topology element1Topology, stk::topology element2Topology) const { return 1.0; }
    virtual bool areVertexWeightsProvidedInAVector() const { return true; }
    void setVertexWeights(const std::vector<double>& weights) { vertex_weights = weights; }
    virtual std::vector<double> getVertexWeightsViaVector() const { return vertex_weights; }
    virtual int getGraphVertexWeight(stk::topology type) const { return 1; }
    virtual double getGraphVertexWeight(stk::mesh::Entity entity, int criteria_index) const { return 1.0; }
    //virtual double getImbalanceTolerance() const { return 1.05; }
    virtual void setDecompMethod(const std::string& input_method) { method = input_method;}
    virtual std::string getDecompMethod() const { return method; }
    void setCoordinateFieldName(const std::string& field_name) { m_field_name = field_name; }
    virtual std::string getCoordinateFieldName() const { return m_field_name; }
    virtual bool shouldFixMechanisms() const { return false; }

private:
    std::vector<double> vertex_weights;
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
        m_defaultWeight(defaultWeight) 
        { 
            method = "parmetis";
            m_includeSearchResultInGraph = false;
        }
    virtual ~FieldVertexWeightSettings() = default;

    virtual double getGraphEdgeWeight(stk::topology element1Topology, stk::topology element2Topology) const { return 1.0; }
    virtual bool areVertexWeightsProvidedInAVector() const { return false; }
    virtual bool areVertexWeightsProvidedViaFields() const { return true; }
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
    FieldVertexWeightSettings() = delete;
    FieldVertexWeightSettings(const FieldVertexWeightSettings&) = delete;
    FieldVertexWeightSettings& operator=(const FieldVertexWeightSettings&) = delete;

    const stk::mesh::BulkData & m_stkMeshBulkData;
    const DoubleFieldType &m_weightField;
    const double m_defaultWeight;
};

class MultipleCriteriaSettings : public GraphCreationSettings
{
public:
    MultipleCriteriaSettings(const std::vector<const stk::mesh::Field<double>*> critFields,
                             const double default_weight = 0.0)
      : m_critFields(critFields), m_defaultWeight(default_weight)
    { 
        m_includeSearchResultInGraph = false;
    }

    MultipleCriteriaSettings(double faceSearchTol, double particleSearchTol, double edgeWeightSearch,
                             const std::string& decompMethod, double multiplierVWSearch,
                             const std::vector<const stk::mesh::Field<double>*> critFields,
                             bool includeSearchResults, const double default_weight = 0.0)
      : GraphCreationSettings(faceSearchTol, particleSearchTol, edgeWeightSearch, decompMethod, multiplierVWSearch),
        m_critFields(critFields), m_defaultWeight(default_weight)
    { 
        m_includeSearchResultInGraph = includeSearchResults;
    }

    virtual ~MultipleCriteriaSettings() = default;

    virtual double getGraphEdgeWeight(stk::topology element1Topology,
      stk::topology element2Topology) const override { return 1.0; }
    virtual bool areVertexWeightsProvidedViaFields() const override { return true; }
    virtual int getGraphVertexWeight(stk::topology type) const override { return 1; }
    virtual double getImbalanceTolerance() const override { return 1.05; }
    virtual int getNumCriteria() const override { return m_critFields.size(); }
    virtual bool isMultiCriteriaRebalance() const override { return true;}
    virtual bool isIncrementalRebalance() const override { return true; }


    virtual double getGraphVertexWeight(stk::mesh::Entity entity, int criteria_index) const override
    {
        ThrowRequireWithSierraHelpMsg(criteria_index>=0 && static_cast<size_t>(criteria_index)<m_critFields.size());
        const double *weight = stk::mesh::field_data(*m_critFields[criteria_index], entity);
        if(weight != nullptr)
        {
            ThrowRequireWithSierraHelpMsg(*weight >= 0);
            return *weight;
        }
        else
        {
            return m_defaultWeight;
        }
    }

protected:
    MultipleCriteriaSettings() = delete;
    MultipleCriteriaSettings(const MultipleCriteriaSettings&) = delete;
    MultipleCriteriaSettings& operator=(const MultipleCriteriaSettings&) = delete;
    const std::vector<const stk::mesh::Field<double>*> m_critFields;
    const double m_defaultWeight;
};

class BasicColoringSettings : public BalanceSettings
{
public:
    BasicColoringSettings() {}

    virtual GraphOption getGraphOption() const
    {
        return BalanceSettings::COLOR_MESH;
    }
};

class BasicColoringByTopologySettings : public BalanceSettings
{
public:
    BasicColoringByTopologySettings() {}

    virtual GraphOption getGraphOption() const
    {
        return BalanceSettings::COLOR_MESH_BY_TOPOLOGY;
    }
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

const std::string& get_coloring_part_base_name();
stk::mesh::Part* get_coloring_part(const stk::mesh::BulkData& bulk, const stk::mesh::Entity& entity);

stk::mesh::PartVector get_root_topology_parts_for_rank(const stk::mesh::BulkData& bulk, stk::mesh::EntityRank rank);

}
}
#endif
