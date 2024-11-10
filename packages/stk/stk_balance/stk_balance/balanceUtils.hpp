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
#include <stk_mesh/base/Entity.hpp>
#include <stk_topology/topology.hpp>
#include "stk_mesh/base/Field.hpp"  // for field_data
#include "stk_mesh/base/FieldBase.hpp"
#include "stk_balance/setup/DefaultSettings.hpp"
#include <memory>

namespace stk
{
namespace mesh { class BulkData; }

namespace balance
{

class FaceSearchTolerance;

using ElementDecomposition = std::vector<int>;
using DoubleFieldType =  stk::mesh::Field<double>;
using BlockWeightMultipliers =  std::map<std::string, double>;

class DecompositionChangeList
{
public:
  explicit DecompositionChangeList(stk::mesh::BulkData & bulk);
  DecompositionChangeList(stk::mesh::BulkData & bulk, const stk::mesh::EntityProcVec & decomposition);
  ~DecompositionChangeList() = default;

  bool has_entity(stk::mesh::Entity entity) const;
  int  get_entity_destination(stk::mesh::Entity entity) const;
  void set_entity_destination(stk::mesh::Entity entity, const int destination);
  void delete_entity(stk::mesh::Entity entity);
  stk::mesh::EntityProcVec get_all_partition_changes();
  stk::mesh::BulkData & get_bulk();
  size_t get_num_global_entity_migrations() const;
  size_t get_max_global_entity_migrations() const;
  stk::mesh::EntityProcVec get_decomposition() const;
  stk::mesh::EntityProcVec get_decomposition_with_full_closure() const;
  auto begin() { return m_dataMap.begin(); }
  auto end() { return m_dataMap.end(); }

private:
  void fill_change_list_from_raw_decomposition(const stk::mesh::EntityProcVec& decomposition);
  void add_downward_closure_for_entity(const std::pair<const stk::mesh::Entity, int> & entityProc,
                                       stk::mesh::EntityProcVec & finalDecomposition) const;

  stk::mesh::BulkData * m_bulk;
  std::map<stk::mesh::Entity, int> m_dataMap;
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

  virtual void set_num_input_processors(unsigned numInputProcs) { m_numInputProcessors = numInputProcs; }
  virtual unsigned get_num_input_processors() const { return m_numInputProcessors; }

  virtual void set_num_output_processors(unsigned numOutputProcs) { m_numOutputProcessors = numOutputProcs; }
  virtual unsigned get_num_output_processors() const { return m_numOutputProcessors; }

  virtual void set_is_rebalancing(bool isRebalancing) { m_isRebalancing = isRebalancing; }
  virtual bool get_is_rebalancing() const { return m_isRebalancing; }

  virtual size_t getNumNodesRequiredForConnection(stk::topology element1Topology, stk::topology element2Topology) const;
  virtual double getGraphEdgeWeight(stk::topology element1Topology, stk::topology element2Topology) const;
  virtual int getGraphVertexWeight(stk::topology type) const;

  virtual double getFieldVertexWeight(const stk::mesh::BulkData &bulkData, stk::mesh::Entity entity, int criteria_index = 0) const;
  virtual GraphOption getGraphOption() const;
  virtual double getGraphEdgeWeightMultiplier() const;

  virtual void setVertexWeightMethod(VertexWeightMethod method);
  virtual VertexWeightMethod getVertexWeightMethod() const;

  virtual bool shouldFixCoincidentElements() const;
  virtual void setShouldFixCoincidentElements(bool fixCoincidentElements);

  // Graph based options only
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

  virtual void setVertexWeightBlockMultiplier(const std::string & blockName, double multiplier);
  virtual const BlockWeightMultipliers & getVertexWeightBlockMultipliers() const;

  virtual bool isIncrementalRebalance() const;
  virtual bool isMultiCriteriaRebalance() const;

  virtual void setVertexWeightFieldName(std::string field_name, unsigned criteria_index = 0);
  virtual std::string getVertexWeightFieldName(unsigned criteria_index = 0) const;
  virtual const stk::mesh::Field<double> * getVertexWeightField(const stk::mesh::BulkData & stkMeshBulkData, unsigned criteria_index = 0) const;
  virtual void setDefaultFieldWeight(double defaultFieldWeight);

  virtual double getImbalanceTolerance() const;

  virtual void setDecompMethod(const std::string& method) ;
  virtual std::string getDecompMethod() const ;

  virtual std::string getInitialDecompMethod() const ;

  virtual std::string getCoordinateFieldName() const ;

  virtual void setShouldPrintDiagnostics(bool shouldPrint);

  virtual bool shouldPrintDiagnostics() const;
  virtual bool shouldPrintMetrics() const;

  virtual int getNumCriteria() const { return m_numCriteria; }
  virtual void setNumCriteria(int num_criteria);

  // Given an element/proc pair, can modify decomposition before elements are migrated
  virtual void modifyDecomposition(DecompositionChangeList & decomp) const ;

  // Need to implement getting particle radius from field
  virtual double getParticleRadius(stk::mesh::Entity particle) const ;

  // experimental
  virtual bool setVertexWeightsBasedOnNumberAdjacencies() const;

  virtual bool allowModificationOfVertexWeightsForSmallMeshes() const;

  virtual bool shouldFixMechanisms() const;
  virtual bool shouldFixSpiders() const;
  virtual std::string getSpiderPartName() const;
  virtual std::string getSpiderVolumeConnectivityCountFieldName() const;
  virtual std::string getOutputSubdomainFieldName() const;
  virtual std::string getDiagnosticElementWeightFieldName() const;
  virtual std::string getVertexConnectivityWeightFieldName() const;
  virtual stk::mesh::Part * getSpiderPart(const stk::mesh::BulkData & stkMeshBulkData) const;
  virtual const stk::mesh::Field<int> * getSpiderVolumeConnectivityCountField(const stk::mesh::BulkData & stkMeshBulkData) const;
  virtual const stk::mesh::Field<int> * getOutputSubdomainField(const stk::mesh::BulkData & stkMeshBulkData) const;
  virtual const stk::mesh::Field<double> * getDiagnosticElementWeightField(const stk::mesh::BulkData & stkMeshBulkData) const;
  virtual const stk::mesh::Field<double> * getVertexConnectivityWeightField(const stk::mesh::BulkData & stkMeshBulkData) const;
  virtual bool usingColoring() const;

  virtual bool useNodeBalancer() const;
  virtual double getNodeBalancerTargetLoadBalance() const;
  virtual unsigned getNodeBalancerMaxIterations() const;

  virtual void set_input_filename(const std::string& filename);
  virtual std::string get_input_filename() const;

  virtual void set_output_filename(const std::string& filename);
  virtual std::string get_output_filename() const;

  virtual void set_log_filename(const std::string& filename);
  virtual std::string get_log_filename() const;

  virtual void setShouldFixMechanisms(bool fixMechanisms) { }
  virtual void setShouldFixSpiders(bool fixSpiders) { }
  virtual void setEdgeWeightForSearch(double w) { }
  virtual void setVertexWeightMultiplierForVertexInSearch(double w) { }
  virtual void setToleranceForFaceSearch(double tol) { }
  virtual void setGraphEdgeWeightMultiplier(double multiplier);

  void set_use_nested_decomp(bool useNestedDecomp) { m_useNestedDecomp = useNestedDecomp; }
  bool get_use_nested_decomp() const { return m_useNestedDecomp; }

private:
  unsigned m_numInputProcessors;
  unsigned m_numOutputProcessors;
  bool m_isRebalancing;
  bool m_shouldFixCoincidentElements;
  std::string m_initialDecompMethod;
  std::string m_inputFilename;
  std::string m_outputFilename;
  std::string m_logFilename;
  int m_numCriteria;
  double m_defaultFieldWeight;
  std::vector<std::string> m_vertexWeightFieldNames;
  BlockWeightMultipliers m_vertexWeightBlockMultipliers;
  bool m_useNestedDecomp;
  bool m_shouldPrintDiagnostics;
  mutable const stk::mesh::Field<double> * m_diagnosticElementWeightsField;
  mutable const stk::mesh::Field<double> * m_vertexConnectivityWeightField;
  mutable std::vector<const stk::mesh::Field<double> *> m_vertexWeightFields;
  VertexWeightMethod m_vertexWeightMethod;
  double m_graphEdgeWeightMultiplier;
};

class BasicGeometricSettings : public BalanceSettings
{
public:
  virtual std::string getDecompMethod() const { return "rcb"; }
};

class GraphCreationSettings : public BalanceSettings
{
public:
  GraphCreationSettings();
  GraphCreationSettings(double faceSearchTol, double particleSearchTol, double edgeWeightSearch,
                        const std::string& decompMethod, double multiplierVWSearch);

  virtual ~GraphCreationSettings() = default;

  size_t getNumNodesRequiredForConnection(stk::topology element1Topology, stk::topology element2Topology) const override;

  virtual double getGraphEdgeWeightForSearch() const override;

  virtual double getGraphEdgeWeight(stk::topology element1Topology, stk::topology element2Topology) const override;

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
  virtual void setToleranceForFaceSearch(double tol) override;
  virtual void setToleranceForParticleSearch(double tol);
  virtual void setEdgeWeightForSearch(double w) override;
  virtual void setVertexWeightMultiplierForVertexInSearch(double w) override;
  virtual void setShouldFixSpiders(bool fixSpiders) override;
  virtual void setShouldFixMechanisms(bool fixMechanisms) override;

  virtual bool shouldFixMechanisms() const override;
  virtual bool shouldFixSpiders() const override;
  virtual stk::mesh::Part * getSpiderPart(const stk::mesh::BulkData & stkMeshBulkData) const override;
  virtual const stk::mesh::Field<int> * getSpiderVolumeConnectivityCountField(const stk::mesh::BulkData & stkMeshBulkData) const override;
  virtual const stk::mesh::Field<int> * getOutputSubdomainField(const stk::mesh::BulkData & stkMeshBulkData) const override;

  virtual void setUseNodeBalancer(bool useBalancer);
  virtual void setNodeBalancerTargetLoadBalance(double targetLoadBalance);
  virtual void setNodeBalancerMaxIterations(unsigned maxIterations);

  virtual bool useNodeBalancer() const override;
  virtual double getNodeBalancerTargetLoadBalance() const override;
  virtual unsigned getNodeBalancerMaxIterations() const override;

protected:
  int getConnectionTableIndex(stk::topology elementTopology) const;
  int getEdgeWeightTableIndex(stk::topology elementTopology) const;

  std::string m_method;
  double m_ToleranceForFaceSearch;
  double m_ToleranceForParticleSearch;
  double m_vertexWeightMultiplierForVertexInSearch;
  double m_edgeWeightForSearch;
  bool m_UseConstantToleranceForFaceSearch;
  bool m_shouldFixSpiders;
  bool m_shouldFixMechanisms;
  mutable stk::mesh::Part * m_spiderPart;
  mutable const stk::mesh::Field<int> * m_spiderVolumeConnectivityCountField;
  mutable const stk::mesh::Field<int> * m_outputSubdomainField;
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
    m_ToleranceForFaceSearch = 0.1;
    m_ToleranceForParticleSearch = 1.0;
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
    : m_stkMeshBulkData(stkMeshBulkData)
  {
    m_method = "parmetis";
    m_includeSearchResultInGraph = false;
    setVertexWeightMethod(VertexWeightMethod::FIELD);
    setVertexWeightFieldName(weightField.name());
    setDefaultFieldWeight(defaultWeight);
  }
  virtual ~FieldVertexWeightSettings() = default;

  virtual double getGraphEdgeWeight(stk::topology element1Topology, stk::topology element2Topology) const { return 1.0; }
  virtual int getGraphVertexWeight(stk::topology type) const { return 1; }
  virtual double getImbalanceTolerance() const { return 1.05; }
  virtual void setDecompMethod(const std::string& input_method) { m_method = input_method;}
  virtual std::string getDecompMethod() const { return m_method; }

protected:
  FieldVertexWeightSettings() = delete;
  FieldVertexWeightSettings(const FieldVertexWeightSettings&) = delete;
  FieldVertexWeightSettings& operator=(const FieldVertexWeightSettings&) = delete;

  const stk::mesh::BulkData & m_stkMeshBulkData;
};

class MultipleCriteriaSettings : public GraphCreationSettings
{
public:
  MultipleCriteriaSettings(const std::vector<const stk::mesh::Field<double>*> critFields,
                           const double default_weight = 0.0)
  {
    setIncludeSearchResultsInGraph(false);
    setNumCriteria(critFields.size());
    setVertexWeightMethod(VertexWeightMethod::FIELD);
    for (unsigned i = 0; i < critFields.size(); ++i) {
      setVertexWeightFieldName(critFields[i]->name(), i);
    }
    setDefaultFieldWeight(default_weight);
  }

  MultipleCriteriaSettings(double faceSearchTol, double particleSearchTol, double edgeWeightSearch,
                           const std::string& decompMethod, double multiplierVWSearch,
                           const std::vector<const stk::mesh::Field<double>*> critFields,
                           bool includeSearchResults, const double default_weight = 0.0)
    : GraphCreationSettings(faceSearchTol, particleSearchTol, edgeWeightSearch, decompMethod, multiplierVWSearch)
  {
    setIncludeSearchResultsInGraph(includeSearchResults);
    setNumCriteria(critFields.size());
    setVertexWeightMethod(VertexWeightMethod::FIELD);
    for (unsigned i = 0; i < critFields.size(); ++i) {
      setVertexWeightFieldName(critFields[i]->name(), i);
    }
    setDefaultFieldWeight(default_weight);
  }

  virtual ~MultipleCriteriaSettings() override = default;

  virtual double getGraphEdgeWeight(stk::topology element1Topology,
                                    stk::topology element2Topology) const override { return 1.0; }
  virtual int getGraphVertexWeight(stk::topology type) const override { return 1; }
  virtual double getImbalanceTolerance() const override { return 1.05; }
  virtual bool isMultiCriteriaRebalance() const override { return true;}
  virtual bool isIncrementalRebalance() const override { return true; }

protected:
  MultipleCriteriaSettings() = delete;
  MultipleCriteriaSettings(const MultipleCriteriaSettings&) = delete;
  MultipleCriteriaSettings& operator=(const MultipleCriteriaSettings&) = delete;
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

class M2NBalanceSettings : public GraphCreationSettings
{
public:
  M2NBalanceSettings()
    : GraphCreationSettings(),
      m_numOutputProcs(0),
      m_useNestedDecomp(false)
  {}

  M2NBalanceSettings(const std::string & inputFileName,
                     unsigned numOutputProcs,
                     bool useNestedDecomp = false)
    : GraphCreationSettings(),
      m_numOutputProcs(numOutputProcs),
      m_useNestedDecomp(useNestedDecomp)
  {
    set_input_filename(inputFileName);
  }

  ~M2NBalanceSettings() = default;

  void set_num_output_processors(unsigned numOutputProcs) { m_numOutputProcs = numOutputProcs; }
  unsigned get_num_output_processors() const { return m_numOutputProcs; }

  void set_use_nested_decomp(bool useNestedDecomp) { m_useNestedDecomp = useNestedDecomp; }
  bool get_use_nested_decomp() const { return m_useNestedDecomp; }

protected:
  unsigned m_numOutputProcs;
  bool m_useNestedDecomp;
};

class GraphEdge
{
public:
  GraphEdge() = default;

  GraphEdge(const stk::mesh::Entity & element1, const stk::mesh::EntityId & element2Id,
            int vertex2ProcOwner, double edgeWeight, bool isEdgeFromSearchArg = false)
    : m_vertex1(element1),
      m_vertex2Id(element2Id),
      m_vertex2OwningProc(vertex2ProcOwner),
      m_weight(edgeWeight),
      m_isEdgeFromSearch(isEdgeFromSearchArg)
  {}

  ~GraphEdge() = default;

  stk::mesh::Entity vertex1() const { return m_vertex1; }
  stk::mesh::EntityId vertex2_id() const { return m_vertex2Id; }
  int vertex2_owning_proc() const { return m_vertex2OwningProc; }
  double weight() const { return m_weight; }
  bool is_edge_from_search() const { return m_isEdgeFromSearch; }

private:
  stk::mesh::Entity m_vertex1;
  stk::mesh::EntityId m_vertex2Id;
  int m_vertex2OwningProc;
  double m_weight;
  bool m_isEdgeFromSearch;
};

inline bool operator<(const GraphEdge &a, const GraphEdge &b)
{
  bool aLessB = (a.vertex1().m_value < b.vertex1().m_value);
  if(a.vertex1().m_value == b.vertex1().m_value)
  {
    aLessB = (a.vertex2_id() < b.vertex2_id());
    if(a.vertex2_id() == b.vertex2_id())
    {
      aLessB = (a.weight() < b.weight());
    }
  }
  return aLessB;
}

inline bool operator==(const GraphEdge &a, const GraphEdge &b)
{
  return (a.vertex1().m_value == b.vertex1().m_value) && (a.vertex2_id() == b.vertex2_id());
}

const std::string& get_coloring_part_base_name();
stk::mesh::Part* get_coloring_part(const stk::mesh::BulkData& bulk, const stk::mesh::Entity& entity);

stk::mesh::PartVector get_root_topology_parts_for_rank(const stk::mesh::BulkData& bulk, stk::mesh::EntityRank rank);

}
}
#endif
