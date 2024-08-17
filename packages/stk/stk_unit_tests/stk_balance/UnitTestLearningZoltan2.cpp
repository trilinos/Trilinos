#include <stk_mesh/base/GetEntities.hpp>

#include <stk_unit_test_utils/unittestMeshUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>

#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>

namespace
{

// First adjacency examples are: element to vertex, element to node
// Second adjacency examples are: element to element

stk::mesh::EntityRank get_mapped_topology(Zoltan2::MeshEntityType etype)
{
  stk::mesh::EntityRank entityRank = stk::topology::INVALID_RANK;
  if(etype==Zoltan2::MESH_VERTEX)
  {
    entityRank = stk::topology::NODE_RANK;
  }
  else if(etype==Zoltan2::MESH_EDGE)
  {
    entityRank = stk::topology::EDGE_RANK;
  }
  else if(etype==Zoltan2::MESH_FACE)
  {
    entityRank = stk::topology::FACE_RANK;
  }
  else if(etype==Zoltan2::MESH_REGION)
  {
    entityRank = stk::topology::ELEM_RANK;
  }
  return entityRank;
}

std::string get_string_for_rank(stk::mesh::EntityRank rank)
{
  std::string rank_name = "unknown";
  if(rank==stk::topology::NODE_RANK)
    rank_name = "vertex";
  else if(rank==stk::topology::EDGE_RANK)
    rank_name = "edge";
  else if(rank==stk::topology::FACE_RANK)
    rank_name = "face";
  else if(rank==stk::topology::ELEM_RANK)
    rank_name = "region";
  return rank_name;
}

Zoltan2::MeshEntityType map_rank_to_mesh_entity_type(stk::mesh::EntityRank rank)
{
  std::string rank_name = "unknown";
  if(rank==stk::topology::NODE_RANK)
    return Zoltan2::MESH_VERTEX;
  else if(rank==stk::topology::EDGE_RANK)
    return Zoltan2::MESH_EDGE;
  else if(rank==stk::topology::FACE_RANK)
    return Zoltan2::MESH_FACE;

  return Zoltan2::MESH_REGION;
}

typedef Zoltan2::BasicUserTypes<double, BalanceLocalNumber, BalanceGlobalNumber> learning_data_t;

class LearningZoltan2Adapter : public Zoltan2::MeshAdapter<learning_data_t>
{
  struct mygraph
  {
    std::vector<double> mVertexCoordinates;
    std::vector<BalanceGlobalNumber> mVertexIds;
    std::vector<double> mVertexWeights;
    unsigned mSpatialDim = 0;
    unsigned mNumFieldCriteria = 1;
    std::vector<double> mEdgeWeights;
    std::vector<BalanceLocalNumber> mOffsets;
    std::vector<BalanceGlobalNumber> mAdjacency;
    size_t mNumGlobalElements = 0;
  };

public:

  typedef Zoltan2::MeshAdapter<learning_data_t> base_adapter_t;

  LearningZoltan2Adapter(stk::mesh::BulkData& bulkData,
                         stk::mesh::EntityRank primary_rank,
                         stk::mesh::EntityRank secondary_rank)
    : m_bulk_data(bulkData),
      m_primary_rank(primary_rank),
      m_secondary_rank(secondary_rank),
      m_primary_entities(),
      m_secondary_entities()
  {
    fillPrimaryAndSecondaryEntities();
    fillVertexIds();
    setEntityTypes(get_string_for_rank(m_primary_rank), get_string_for_rank(m_secondary_rank), get_string_for_rank(m_secondary_rank));
    fillAdjacency();
    fillCoordinates();
    fillVertexWeights();
  }

  virtual ~LearningZoltan2Adapter() { }

  virtual size_t getLocalNumOf(Zoltan2::MeshEntityType etype) const
  {
    size_t num_of = 0;
    if(get_mapped_topology(etype) == m_primary_rank)
      num_of = m_primary_entities.size();
    else if(get_mapped_topology(etype) == m_secondary_rank)
      num_of = m_secondary_entities.size();
    else
      STK_ThrowRequireMsg(false, "getLocalNumOf couldn't return valid answer.");
    // std::cerr<<"getLocalNumOf: "<<num_of<< " for rank: " << get_mapped_topology(etype) << std::endl;
    return num_of;
  }

  virtual void getIDsViewOf(Zoltan2::MeshEntityType etype, BalanceGlobalNumber const *&Ids) const
  {
    Ids = nullptr;
    if(get_mapped_topology(etype) == m_primary_rank)
      Ids = m_graph.mVertexIds.data();
    else if(get_mapped_topology(etype) == m_secondary_rank)
      Ids = m_secondary_ids.data();
  }

  virtual int getDimension() const
  {
    return m_bulk_data.mesh_meta_data().spatial_dimension();
  }

  virtual void getCoordinatesViewOf(Zoltan2::MeshEntityType etype, const scalar_t *&coords, int &stride, int coordDim) const
  {
    coords = NULL;
    stride = getDimension();
    STK_ThrowRequireMsg(get_mapped_topology(etype)== m_primary_rank, "Error!");
    if(get_mapped_topology(etype)== m_primary_rank)
    {
      if(!m_graph.mVertexCoordinates.empty())
      {
        coords = &m_graph.mVertexCoordinates[coordDim];
      }
    }
  }

  virtual int getNumWeightsPerOf(Zoltan2::MeshEntityType etype) const
  {
    STK_ThrowRequireMsg(get_mapped_topology(etype)== m_primary_rank, "Error!");
    size_t numWeightsPerVertex = 0;
    if(get_mapped_topology(etype)==m_primary_rank)
    {
      numWeightsPerVertex = 1;
    }
    return numWeightsPerVertex;
  }

  virtual void getWeightsViewOf(Zoltan2::MeshEntityType etype, const scalar_t *&weights, int &stride, int idx = 0) const
  {
    STK_ThrowRequireMsg(get_mapped_topology(etype)== m_primary_rank, "Error!");
    weights = NULL;
    stride = 1;
    if(get_mapped_topology(etype)==m_primary_rank)
    {
      if(!m_graph.mVertexWeights.empty())
      {
        weights = m_graph.mVertexWeights.data();
      }
    }
  }

  // Fill in 2nd adjs if needed

  virtual bool avail2ndAdjs(Zoltan2::MeshEntityType sourcetarget, Zoltan2::MeshEntityType through) const
  {
    bool avail2ndAdj = false;
    //if(get_mapped_topology(source)==m_primary_rank && get_mapped_topology(target)==m_secondary_rank)
    //    avail2ndAdj = true;
    return avail2ndAdj;
  }

  virtual size_t getLocalNum2ndAdjs(Zoltan2::MeshEntityType sourcetarget, Zoltan2::MeshEntityType through) const
  {
    // size of the adjacency list (for 2nd adjacency)
    return 0;
  }

  virtual void get2ndAdjsView(Zoltan2::MeshEntityType sourcetarget, Zoltan2::MeshEntityType through, const BalanceLocalNumber *&offsets, const BalanceGlobalNumber *&adjacencyIds) const
  {
    //ThrowRequireMsg(get_mapped_topology(sourcetarget)==m_primary_rank && get_mapped_topology(through)==m_secondary_rank, "Error!");
    STK_ThrowRequireMsg(false, "Error!");
  }

  virtual int getNumWeightsPer2ndAdj(Zoltan2::MeshEntityType sourcetarget, Zoltan2::MeshEntityType through) const
  {
    //ThrowRequireMsg(get_mapped_topology(sourcetarget)==m_primary_rank && get_mapped_topology(through)==m_secondary_rank, "Error!");
    STK_ThrowRequireMsg(false, "Error!");
    return 1;
  }

  virtual void get2ndAdjWeightsView(Zoltan2::MeshEntityType sourcetarget, Zoltan2::MeshEntityType through, const scalar_t *&weights, int &stride, int idx) const
  {
    STK_ThrowRequireMsg(false, "Error!");
    // return edge weights (one per edge)
  }

  // first adjacency API below

  virtual bool availAdjs(Zoltan2::MeshEntityType source, Zoltan2::MeshEntityType target) const
  {
    bool availAdj = false;
    if(get_mapped_topology(source)==m_primary_rank && get_mapped_topology(target)==m_secondary_rank)
      availAdj = true;
    //std::cerr<<"availAdj: " << availAdj << ", get_mapped_topology(source)=" << get_mapped_topology(source) << ", get_mapped_topology(target)=" << get_mapped_topology(target) << std::endl;
    return availAdj;
  }

  virtual size_t getLocalNumAdjs(Zoltan2::MeshEntityType source, Zoltan2::MeshEntityType target) const
  {
    size_t numAdjs = 0;
    STK_ThrowRequireMsg(get_mapped_topology(source)==m_primary_rank && get_mapped_topology(target)==m_secondary_rank, "Error!");
    if(get_mapped_topology(source)==m_primary_rank && get_mapped_topology(target)==m_secondary_rank)
    {
      numAdjs = m_graph.mAdjacency.size();
    }
    //std::cerr<<"getLocalNumAdjs: "<<numAdjs<<std::endl;
    return numAdjs;
  }

  virtual void getAdjsView(Zoltan2::MeshEntityType source, Zoltan2::MeshEntityType target, const BalanceLocalNumber *&offsets, const BalanceGlobalNumber *& adjacencyIds) const
  {
    STK_ThrowRequireMsg(get_mapped_topology(source)==m_primary_rank && get_mapped_topology(target)==m_secondary_rank, "Error!");
    offsets = NULL;
    adjacencyIds = NULL;
    if(get_mapped_topology(source)==m_primary_rank && get_mapped_topology(target)==m_secondary_rank)
    {
      if(!m_graph.mAdjacency.empty())
      {
        adjacencyIds = m_graph.mAdjacency.data();
      }
      if(!m_graph.mOffsets.empty())
      {
        offsets = m_graph.mOffsets.data();
      }
    }
  }

  virtual bool useDegreeAsWeightOf(Zoltan2::MeshEntityType etype, int idx) const
  {
    return false;
  }

private:

  void fillPrimaryAndSecondaryEntities()
  {
    stk::mesh::get_selected_entities(m_bulk_data.mesh_meta_data().locally_owned_part(), m_bulk_data.buckets(m_primary_rank), m_primary_entities);
    stk::mesh::get_selected_entities(m_bulk_data.mesh_meta_data().locally_owned_part(), m_bulk_data.buckets(m_secondary_rank), m_secondary_entities);
  }

  void fillVertexWeights()
  {
    m_graph.mVertexWeights.resize(m_primary_entities.size(), 1.0);
    for(size_t i=0;i<m_primary_entities.size();++i)
      m_graph.mVertexWeights[i] = m_bulk_data.identifier(m_primary_entities[i]);
  }

  void fillVertexIds()
  {
    m_graph.mVertexIds.resize(m_primary_entities.size());
    for(size_t i=0;i<m_primary_entities.size();++i)
      m_graph.mVertexIds[i] = m_bulk_data.identifier(m_primary_entities[i]);
    m_secondary_ids.resize(m_secondary_entities.size());
    for(size_t i=0;i<m_secondary_entities.size();++i)
      m_secondary_ids[i] = m_bulk_data.identifier(m_secondary_entities[i]);
  }

  void fillAdjacency()
  {
    unsigned counter = 0;
    m_graph.mOffsets.push_back(0);

    for(stk::mesh::Entity p_entity : m_primary_entities )
    {
      unsigned num_adj_entities = m_bulk_data.num_connectivity(p_entity, m_secondary_rank);
      const stk::mesh::Entity* adj_entities = m_bulk_data.begin(p_entity, m_secondary_rank);
      for(unsigned j=0;j<num_adj_entities;++j)
      {
        BalanceGlobalNumber id = m_bulk_data.identifier(adj_entities[j]);
        m_graph.mAdjacency.push_back(id);
      }

      counter++;
      m_graph.mOffsets.push_back(m_graph.mOffsets[counter-1] + num_adj_entities);
    }

    //        for(size_t i=0;i<m_primary_entities.size();++i)
    //        {
    //            std::cerr << "Entity with index: " << i << " is connected to: ";
    //            for(BalanceLocalNumber j=m_graph.mOffsets[i];j<m_graph.mOffsets[i+1];j++)
    //            {
    //                std::cerr << m_graph.mAdjacency[j] << "\t";
    //            }
    //            std::cerr << std::endl;
    //        }
  }


  void fillCoordinates()
  {
    m_graph.mVertexCoordinates.resize(m_primary_entities.size()*getDimension(), 0);
    const stk::mesh::FieldBase * coord = m_bulk_data.mesh_meta_data().get_field(stk::topology::NODE_RANK,"coordinates");

    for(size_t i=0;i<m_primary_entities.size();++i)
      stk::balance::internal::fillEntityCentroid(m_bulk_data, coord, m_primary_entities[i], &m_graph.mVertexCoordinates[getDimension()*i]);

    //        size_t num_entities = m_primary_entities.size();
    //        for(size_t i=0;i<num_entities;++i)
    //        {
    //            std::cerr << m_bulk_data.entity_key(m_primary_entities[i]) << " has coordinates: ";
    //            std::cerr << m_graph.mVertexCoordinates[3*i] << "\t" << m_graph.mVertexCoordinates[3*i+1] << "\t" << m_graph.mVertexCoordinates[3*i+2] << std::endl;
    //        }
  }

  stk::mesh::BulkData &m_bulk_data;
  stk::mesh::EntityRank m_primary_rank;
  stk::mesh::EntityRank m_secondary_rank;
  stk::mesh::EntityVector m_primary_entities;
  stk::mesh::EntityVector m_secondary_entities;
  std::vector<BalanceGlobalNumber> m_secondary_ids;
  mygraph m_graph;

public:

  const stk::mesh::EntityVector& get_entities() const
  {
    return m_primary_entities;
  }
};

////////////////////////////////////////////////////////////////////////////////////////////

class UsingZoltan2 : public stk::unit_test_util::MeshFixture
{
protected:
  void run_decomp_with_method(const std::string& method, int nparts, stk::mesh::EntityRank primary_rank, stk::mesh::EntityRank secondary_rank)
  {
    stk::io::write_mesh("junk.exo", get_bulk());
    LearningZoltan2Adapter adapter(get_bulk(), primary_rank, secondary_rank);
    std::vector<int> elem2proc(adapter.getLocalNumOf(map_rank_to_mesh_entity_type(primary_rank)), get_bulk().parallel_rank());
    use_zoltan2_with_adapter(adapter, method, nparts, elem2proc);

    const stk::mesh::EntityVector entities = adapter.get_entities();

    std::ostringstream os;
    os << "For processor: " << get_bulk().parallel_rank() << " using " << method << ":\n";
    for(size_t i=0;i<entities.size();++i)
      os << get_bulk().entity_key(entities[i]) << " goes to processor " << elem2proc[i] << std::endl;
    std::cerr << os.str();
  }

  void use_zoltan2_with_adapter(LearningZoltan2Adapter& adapter, const std::string &method, int nparts, std::vector<int> &elem2proc)
  {
    Teuchos::ParameterList params("test params");
    params.set("debug_level", "basic_status");

    double imbalance_allowed = 1.1;
    params.set("imbalance_tolerance", imbalance_allowed);
    params.set("num_global_parts", nparts);
    params.set("algorithm", method);

    Teuchos::ParameterList &zparams = params.sublist("zoltan_parameters", false);
    zparams.set("LB_METHOD", "PHG");
    zparams.set("LB_APPROACH", "PARTITION");

    Zoltan2::PartitioningProblem<LearningZoltan2Adapter> problem(&adapter, &params, get_bulk().parallel());
    problem.solve();

    const StkMeshZoltanAdapter::part_t *processorOntoWhichEntityBelongs = problem.getSolution().getPartListView();
    for(size_t i=0;i<elem2proc.size();++i)
      elem2proc[i] = processorOntoWhichEntityBelongs[i];
  }
};

const std::vector<std::string> get_decomp_list()
{
  const std::vector<std::string> decomp_options = { "rcb", "block", "multijagged", "parmetis", "zoltan" };
  return decomp_options;
}

TEST_F(UsingZoltan2, testWithAuraElem2Node)
{
  if(stk::parallel_machine_size(get_comm())<=2)
  {
    int numParts = 2;
    setup_mesh("generated:1x1x6", stk::mesh::BulkData::AUTO_AURA);
    const std::vector<std::string> decomp_options = get_decomp_list();
    for(size_t i=0;i<decomp_options.size();++i)
      run_decomp_with_method(decomp_options[i], numParts, stk::topology::ELEMENT_RANK, stk::topology::NODE_RANK);
  }
}


TEST_F(UsingZoltan2, testWithAuraNode2Elem)
{
  if(stk::parallel_machine_size(get_comm())<=2)
  {
    int numParts = 2;
    setup_mesh("generated:1x1x6", stk::mesh::BulkData::AUTO_AURA);
    const std::vector<std::string> decomp_options = get_decomp_list();
    for(size_t i=0;i<decomp_options.size();++i)
      run_decomp_with_method(decomp_options[i], numParts, stk::topology::NODE_RANK, stk::topology::ELEMENT_RANK);
  }
}

}
