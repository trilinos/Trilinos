// Copyright (c) 2013, Sandia Corporation.
 // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 // the U.S. Government retains certain rights in this software.
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
 //     * Neither the name of Sandia Corporation nor the names of its
 //       contributors may be used to endorse or promote products derived
 //       from this software without specific prior written permission.
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
    BalanceSettings() {}
    virtual ~BalanceSettings() {}

    enum GraphOption
    {
        LOADBALANCE = 0,
        COLORING
    };

    virtual size_t getNumNodesRequiredForConnection(stk::topology element1Topology, stk::topology element2Topology) const;
    virtual double getGraphEdgeWeight(stk::topology element1Topology, stk::topology element2Topology) const;
    virtual int getGraphVertexWeight(stk::topology type) const;
    virtual double getGraphVertexWeight(stk::mesh::Entity entity, int criteria_index = 0) const ;
    virtual GraphOption getGraphOption() const;

    // Graph (parmetis) based options only
    virtual bool includeSearchResultsInGraph() const;

    virtual double getToleranceForFaceSearch(const stk::mesh::BulkData & mesh,
                                             const stk::mesh::FieldBase & coordField,
                                             const stk::mesh::Entity * faceNodes,
                                             const unsigned numFaceNodes) const;
    virtual void setToleranceFunctionForFaceSearch(std::shared_ptr<stk::balance::FaceSearchTolerance> faceSearchTolerance);

    virtual double getToleranceForParticleSearch() const;
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
    virtual std::string getSpiderConnectivityCountFieldName() const;
    virtual const stk::mesh::Field<int> * getSpiderConnectivityCountField(const stk::mesh::BulkData & stkMeshBulkData) const;
};

class BasicGeomtricSettings : public BalanceSettings
{
    virtual std::string getDecompMethod() const { return "rcb"; }
};

class GraphCreationSettings : public BalanceSettings
{
public:
    GraphCreationSettings(): GraphCreationSettings(0.0001, 3, 15, "parmetis", 5.0)
    {}

    GraphCreationSettings(double faceSearchTol, double particleSearchTol, double edgeWeightSearch, const std::string& decompMethod, double multiplierVWSearch)
                           : mToleranceForFaceSearch(faceSearchTol),
                             mToleranceForParticleSearch(particleSearchTol),
                             edgeWeightForSearch (edgeWeightSearch),
                             method(decompMethod),
                             vertexWeightMultiplierForVertexInSearch(multiplierVWSearch),
                             m_UseConstantToleranceForFaceSearch(true),
                             m_shouldFixSpiders(false),
                             m_spiderConnectivityCountField(nullptr)
    {}

    virtual ~GraphCreationSettings() {}

    size_t getNumNodesRequiredForConnection(stk::topology element1Topology, stk::topology element2Topology) const;

    virtual double getGraphEdgeWeightForSearch() const;

    virtual double getGraphEdgeWeight(stk::topology element1Topology, stk::topology element2Topology) const;

    virtual double getGraphVertexWeight(stk::mesh::Entity entity, int criteria_index = 0) const;

    virtual int getGraphVertexWeight(stk::topology type) const;

    virtual GraphOption getGraphOption() const;
    virtual bool includeSearchResultsInGraph() const ;
    virtual double getToleranceForParticleSearch() const ;

    virtual double getToleranceForFaceSearch(const stk::mesh::BulkData & mesh,
                                             const stk::mesh::FieldBase & coordField,
                                             const stk::mesh::Entity * faceNodes,
                                             const unsigned numFaceNodes) const;
    virtual void setToleranceFunctionForFaceSearch(std::shared_ptr<stk::balance::FaceSearchTolerance> faceSearchTolerance);

    virtual bool getEdgesForParticlesUsingSearch() const ;
    virtual double getVertexWeightMultiplierForVertexInSearch() const ;
    virtual std::string getDecompMethod() const ;

    virtual void setDecompMethod(const std::string& input_method);
    virtual void setToleranceForFaceSearch(double tol);
    virtual void setToleranceForParticleSearch(double tol) ;
    virtual void setEdgeWeightForSearch(double w) ;
    virtual void setVertexWeightMultiplierForVertexInSearch(double w) ;
    virtual void setShouldFixSpiders(bool fixSpiders);

    virtual bool shouldFixMechanisms() const override;
    virtual bool shouldFixSpiders() const override;
    virtual const stk::mesh::Field<int> * getSpiderConnectivityCountField(const stk::mesh::BulkData & stkMeshBulkData) const override;

protected:
    int getConnectionTableIndex(stk::topology elementTopology) const;
    double mToleranceForFaceSearch;
    double mToleranceForParticleSearch;
    double edgeWeightForSearch;
    std::string method;
    double vertexWeightMultiplierForVertexInSearch;
    bool m_UseConstantToleranceForFaceSearch;
    bool m_shouldFixSpiders;
    mutable const stk::mesh::Field<int> * m_spiderConnectivityCountField;
    std::shared_ptr<stk::balance::FaceSearchTolerance> m_faceSearchToleranceFunction;
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

class BasicZoltan2Settings : public GraphCreationSettings
{
public:
    BasicZoltan2Settings()
      : GraphCreationSettings(0.0005, 0.3, 100.0, "rcb", 6.0) {}
    virtual bool includeSearchResultsInGraph() const { return false; }
    virtual bool getEdgesForParticlesUsingSearch() const { return true; }
};

class UserSpecifiedVertexWeightsSetting : public GraphCreationSettings
{
public:
    virtual double getGraphEdgeWeight(stk::topology element1Topology, stk::topology element2Topology) const { return 1.0; }
    virtual bool areVertexWeightsProvidedInAVector() const { return true; }
    virtual bool includeSearchResultsInGraph() const { return false; }
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
    virtual bool areVertexWeightsProvidedInAVector() const { return false; }
    virtual bool areVertexWeightsProvidedViaFields() const { return true; }
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

class MultipleCriteriaSettings : public stk::balance::GraphCreationSettings
{
public:
    MultipleCriteriaSettings(const std::vector<const stk::mesh::Field<double>*> critFields,
                             const double default_weight = 0.0)
      : m_critFields(critFields), m_defaultWeight(default_weight)
    { }

    MultipleCriteriaSettings(double faceSearchTol, double particleSearchTol, double edgeWeightSearch,
                             const std::string& decompMethod, double multiplierVWSearch,
                             const std::vector<const stk::mesh::Field<double>*> critFields,
                             bool includeSearchResults, const double default_weight = 0.0)
      : GraphCreationSettings(faceSearchTol, particleSearchTol, edgeWeightSearch, decompMethod, multiplierVWSearch),
        m_critFields(critFields), m_includeSearchResults(includeSearchResults), m_defaultWeight(default_weight)
    { }

    virtual ~MultipleCriteriaSettings() = default;

    virtual double getGraphEdgeWeight(stk::topology element1Topology, stk::topology element2Topology) const { return 1.0; }
    virtual bool areVertexWeightsProvidedViaFields() const { return true; }
    virtual bool includeSearchResultsInGraph() const { return m_includeSearchResults; }
    virtual int getGraphVertexWeight(stk::topology type) const { return 1; }
    virtual double getImbalanceTolerance() const { return 1.05; }
    virtual int getNumCriteria() const { return m_critFields.size(); }
    virtual bool isMultiCriteriaRebalance() const { return true;}
    virtual bool isIncrementalRebalance() const override { return true; }


    virtual double getGraphVertexWeight(stk::mesh::Entity entity, int criteria_index) const
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
    MultipleCriteriaSettings() = default;
    MultipleCriteriaSettings(const MultipleCriteriaSettings&) = delete;
    MultipleCriteriaSettings& operator=(const MultipleCriteriaSettings&) = delete;

    const std::vector<const stk::mesh::Field<double>*> m_critFields;
    bool m_includeSearchResults = false;
    const double m_defaultWeight;
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
