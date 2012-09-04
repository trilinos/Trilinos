/**
   1. class/struct for the data type stored on a sub-dim entity (value)
   2. key (subDimEntity) gives value& 
        2a. getDefaultValueFromKeyIfKeyNotPresent
        2b. value * isKeyInMap(key) [returns 0 if not present]
        2c. insert (key,value) pair, with overwrite
        2d. for a given key, reset the data (value type) to "empty"

        is_value_null, is_value_cleared, is_value_null_or_cleared
        insert(key,value)
        getValueWithDefaultIfNotPresent(key)

        Paradigm: avoid determining if it's null - just ask for is_value_null_or_cleared, if it is,
           then getValueWithDefaultIfNotPresent(key), else

        Or, always do is_value_null, then insertAndReturnDefaultValue(), else getValue().  This is needed
           for performance since we don't want to always construct a default value, as the operator[] does.

        


   3. we treat key,value as value-types, but we need to be aware of performance
        and pass/return references as needed
   4. notify delete a node: find subDimEntity containing the node and remove from the value
   5. notify delete an element: if an element is about to be deleted from stk_mesh, notify
        the DB so the appropriate owning element flags can be reset
   6.
 */

#ifndef stk_adapt_NodeRegistry_hpp
#define stk_adapt_NodeRegistry_hpp

#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <limits>
#include <cmath>
#include <utility>
#include <math.h>
#include <map>
#include <set>
#include <vector>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

#include <stk_mesh/base/EntityKey.hpp>

#include <stk_percept/stk_mesh.hpp>

#include <stk_util/environment/CPUTime.hpp>

#include <stk_percept/NoMallocArray.hpp>
#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/Util.hpp>

#include <stk_percept/PerceptBoostArray.hpp>

#include <boost/tuple/tuple_io.hpp>
#include <boost/tuple/tuple_comparison.hpp>

/// define only one of these to be 1
/// current best setting is NODE_REGISTRY_MAP_TYPE_BOOST = 1
#define NODE_REGISTRY_MAP_TYPE_BOOST 1
#define NODE_REGISTRY_MAP_TYPE_TR1 0
#define NODE_REGISTRY_MAP_TYPE_STD 0
#define NODE_REGISTRY_MAP_TYPE_GOOGLE 0

// next two are deprecated
#define NODE_REGISTRY_MAP_TYPE_PERCEPT_HASHTABLE 0
#define NODE_REGISTRY_MAP_TYPE_TEUCHOS_HASHTABLE 0

/// current best is STK_ADAPT_NODEREGISTRY_USE_ENTITY_REPO 1, DO_REHASH 1
#define STK_ADAPT_NODEREGISTRY_USE_ENTITY_REPO 0
#define STK_ADAPT_NODEREGISTRY_DO_REHASH 1

#if NODE_REGISTRY_MAP_TYPE_BOOST
#include <boost/unordered_map.hpp>
#endif

#if NODE_REGISTRY_MAP_TYPE_STD
#include <map>
#endif

#if NODE_REGISTRY_MAP_TYPE_TR1
#include <tr1/unordered_map>
#endif


#include <stk_adapt/SubDimCell.hpp>

#if NODE_REGISTRY_MAP_TYPE_GOOGLE
#include <google/sparse_hash_map>
#include <google/dense_hash_map>
#endif

// set to maximum number of procs
#define PSEUDO_ELEMENT_MAGIC_NUMBER 1000000
// set to skip over the FAMILY_TREE_RANK (which is 1 + element rank)
#define PSEUDO_ELEMENT_RANK_SHIFT 2


namespace stk {
  namespace adapt {

    using namespace stk::percept;
    using std::vector;
    using std::map;
    using std::set;

    // FIXME
    static bool s_allow_empty_sub_dims = true; // for uniform refine, this should be false for testing

    typedef std::pair<stk::mesh::EntityRank, unsigned> NeededEntityType;

    // using tuple here instead of pair to allow for future expansion

    // pair of node id and the owning element for a node on a sub-dimensional entity (like a face or edge)
    enum SubDimCellDataEnum {
      SDC_DATA_GLOBAL_NODE_IDS,
      SDC_DATA_OWNING_ELEMENT_KEY
    };

    typedef  stk::mesh::Entity* NodeIdsOnSubDimEntityTypeQuantum;

    /// data on a sub-dim entity (global node ids on the entity, the owning element's id)
    struct NodeIdsOnSubDimEntityType : public std::vector<NodeIdsOnSubDimEntityTypeQuantum>
    {
      typedef std::vector<NodeIdsOnSubDimEntityTypeQuantum> base_type;
      typedef std::vector<stk::mesh::EntityId> entity_id_vector_type;
      entity_id_vector_type m_entity_id_vector;
      unsigned m_mark;

      NodeIdsOnSubDimEntityType(unsigned sz=1, NodeIdsOnSubDimEntityTypeQuantum allValues=0) : base_type(sz,allValues),
                                                                                               m_entity_id_vector(sz,0u),
                                                                                               m_mark(0u)  {}

      void resize(size_t sz)
      {
        m_entity_id_vector.resize(sz);
        base_type::resize(sz);
      }

      void pack(CommBuffer& buff)
      {
        buff.pack< unsigned > ( this->size() );
        m_entity_id_vector.resize( this->size() );
        for (unsigned ii = 0; ii < this->size(); ii++)
          {
            if (!(*this)[ii])
              throw std::logic_error("logic err in NodeIdsOnSubDimEntityType::pack");
            stk::mesh::EntityId id = ((*this)[ii])->identifier();
            m_entity_id_vector[ii] = id;
            buff.pack<stk::mesh::EntityId>( id );
          }
      }
      void unpack(PerceptMesh& eMesh, CommBuffer& buff)
      {
        unsigned sz;
        buff.unpack< unsigned > ( sz );
        this->resize( sz );
        m_entity_id_vector.resize( sz );
        for (unsigned ii = 0; ii < this->size(); ii++)
          {
            stk::mesh::EntityId id=0;
            buff.unpack<stk::mesh::EntityId>( id );
            m_entity_id_vector[ii] = id;
          }
      }
    };

    inline std::ostream &operator<<(std::ostream& out, const boost::array<stk::mesh::EntityId, 1>& arr)
    {
      out << arr[0];
      return out;
    }


    typedef boost::tuple<NodeIdsOnSubDimEntityType, stk::mesh::EntityKey> SubDimCellData;

#define SDS_ENTITY_TYPE_ID 0
#if SDS_ENTITY_TYPE_ID
    typedef stk::mesh::EntityId SDSEntityType;

    struct CompareSDSEntityType {
      bool operator() (SDSEntityType i, SDSEntityType j) { return (i < j) ; }
    };

#else
    typedef stk::mesh::Entity * SDSEntityType;

    struct CompareSDSEntityType {
      bool operator() (SDSEntityType i, SDSEntityType j) { return (i->identifier() < j->identifier());}
    };

    template<>
    inline int SubDimCell<SDSEntityType, 4, SubDimCellCompare<SDSEntityType> >::hashCode()
#if 1
    {
      typedef stk::percept::NoMallocArray<SDSEntityType,4> base_type;

      std::size_t sum = 0;

      for ( base_type::iterator i = this->begin(); i != this->end(); i++)
        {
          sum += static_cast<std::size_t>((*i)->identifier());
        }
      return sum;
    }
#endif

#endif

    typedef SubDimCell<SDSEntityType> SubDimCell_SDSEntityType;

    inline std::ostream& operator<<(std::ostream& out,  SubDimCellData& val)
    {
      out << "SDC:: node ids= " << val.get<SDC_DATA_GLOBAL_NODE_IDS>()
          << " owning element rank= " << stk::mesh::entity_rank(val.get<SDC_DATA_OWNING_ELEMENT_KEY>())
          << " owning element id= " << stk::mesh::entity_id(val.get<SDC_DATA_OWNING_ELEMENT_KEY>());
      return out;
    }

  }
}

namespace stk {
  namespace adapt {

    typedef stk::mesh::Entity *EntityPtr;

    /// map of the node ids on a sub-dim entity to the data on the sub-dim entity

#if NODE_REGISTRY_MAP_TYPE_BOOST
#  ifdef STK_HAVE_TBB

    typedef tbb::scalable_allocator<std::pair<SubDimCell_SDSEntityType const, SubDimCellData> > RegistryAllocator;
    typedef boost::unordered_map<SubDimCell_SDSEntityType, SubDimCellData, my_fast_hash<SDSEntityType, 4>, my_fast_equal_to<SDSEntityType, 4>, RegistryAllocator > SubDimCellToDataMap;

#  else

    typedef boost::unordered_map<SubDimCell_SDSEntityType, SubDimCellData, my_fast_hash<SDSEntityType, 4>, my_fast_equal_to<SDSEntityType, 4> > SubDimCellToDataMap;
    typedef boost::unordered_map<stk::mesh::EntityId, EntityPtr > EntityRepo;

    typedef boost::tuple<const stk::mesh::Entity *, stk::mesh::EntityRank, unsigned> ElementSideTuple;

    struct my_tuple_hash : public std::unary_function< ElementSideTuple, std::size_t>
    {
      inline std::size_t
      operator()(const ElementSideTuple& x) const
      {
        return std::size_t(x.get<0>())+std::size_t(x.get<1>())+std::size_t(x.get<2>());
      }
    };

    typedef boost::unordered_map<ElementSideTuple, SubDimCellData, my_tuple_hash > ElementSideMap;

#  endif
#endif

#if NODE_REGISTRY_MAP_TYPE_GOOGLE

    //         typedef google::sparse_hash_map<unsigned, unsigned *> google_sparse_map_type;
    //         google_sparse_map_type google_sparse_map1(init_capacity);
    //         google_sparse_map1.set_empty_key(2*N);

#  ifdef STK_HAVE_TBB

    typedef tbb::scalable_allocator<std::pair<SubDimCell_SDSEntityType const, SubDimCellData> > RegistryAllocator;
    typedef google::sparse_hash_map<SubDimCell_SDSEntityType, SubDimCellData, my_hash<stk::mesh::EntityId,4>, my_equal_to<stk::mesh::EntityId,4>, RegistryAllocator > SubDimCellToDataMap;
    typedef google::sparse_hash_map<stk::mesh::EntityId, EntityPtr > EntityRepo;

#  else

    typedef google::sparse_hash_map<SubDimCell_SDSEntityType, SubDimCellData, my_hash<stk::mesh::EntityId,4>, my_equal_to<stk::mesh::EntityId,4> > SubDimCellToDataMap;
    typedef boost::unordered_map<stk::mesh::EntityId, EntityPtr > EntityRepo;

#  endif
#endif

#if NODE_REGISTRY_MAP_TYPE_TR1
    typedef tr1::unordered_map<SubDimCell_SDSEntityType, SubDimCellData, my_hash<stk::mesh::EntityId,4>, my_equal_to<stk::mesh::EntityId,4> > SubDimCellToDataMap;
    typedef tr1::unordered_map<stk::mesh::EntityId, EntityPtr > EntityRepo;
#endif

#if NODE_REGISTRY_MAP_TYPE_STD

#  if STK_ADAPT_SUBDIMCELL_USES_STL_SET
    typedef map<SubDimCell_SDSEntityType, SubDimCellData> SubDimCellToDataMap;
#  endif

#  if STK_ADAPT_SUBDIMCELL_USES_STL_VECTOR
    typedef map<SubDimCell_SDSEntityType, SubDimCellData, SubDimCell_compare<stk::mesh::EntityId> > SubDimCellToDataMap;
#  endif

#  if STK_ADAPT_SUBDIMCELL_USES_NO_MALLOC_ARRAY
    typedef map<SubDimCell_SDSEntityType, SubDimCellData, SubDimCell_compare<stk::mesh::EntityId> > SubDimCellToDataMap;
#  endif

    typedef map<stk::mesh::EntityId, EntityPtr > EntityRepo;
#endif

    // Rank of sub-dim cells needing new nodes, which sub-dim entity, one non-owning element identifier, nodeId_elementOwnderId.first
    // FIXME - consider using bitfields to pack first two entries into a short - does this save anything on buffer size?
    enum CommDataTypeEnum {
      CDT_NEEDED_ENTITY_RANK,
      CDT_SUB_DIM_ENTITY_ORDINAL,
      CDT_NON_OWNING_ELEMENT_KEY
      //,NEW_NODE_IDS
    };
    enum
      {
        NEW_NODE_IDS
      };
    typedef boost::tuple<stk::mesh::EntityRank, unsigned, stk::mesh::EntityKey> CommDataType;

    enum {
      MODE_SERIAL,
      MODE_BUFFER_SIZING,
      MODE_BUFFERS_ALLOCD,
      MODE_SEND_DONE
      // after unpack, reset to MODE_SERIAL, i.e. it's cyclical
    };

    enum NodeRegistryState {
      NRS_NONE,
      NRS_START_REGISTER_NODE,
      NRS_END_REGISTER_NODE,
      NRS_START_CHECK_FOR_REMOTE,
      NRS_END_CHECK_FOR_REMOTE,
      NRS_START_PARALLEL_OPS,
      NRS_START_GET_FROM_REMOTE,
      NRS_END_GET_FROM_REMOTE,
      NRS_END_PARALLEL_OPS
    };


    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    class NodeRegistry
    {
      static const unsigned NR_MARK_NONE = 1u;
      static const unsigned NR_MARK = 2u;

    public:
      // FIXME use unordered_set
      typedef std::set<stk::mesh::Entity *> SetOfEntities;

      //========================================================================================================================
      // high-level interface
      //NodeRegistry(percept::PerceptMesh& eMesh) : m_eMesh(eMesh), m_comm_all(eMesh.get_bulk_data()->parallel()), m_gee_cnt(0), m_gen_cnt(0),
      //m_entity_repo(stk::mesh::stk::percept::EntityRankEnd)

      NodeRegistry(percept::PerceptMesh& eMesh, bool useCustomGhosting = false) : m_eMesh(eMesh), 
                                                  //m_comm_all(0),
                                                  m_comm_all( new stk::CommAll(eMesh.get_bulk_data()->parallel()) ),
                                                  //m_comm_all(eMesh.get_bulk_data()->parallel()),
                                                  // why does this cause failures? 
                                                  //m_cell_2_data_map(eMesh.get_number_elements()*8u),
                                                  m_useCustomGhosting(useCustomGhosting),
                                                  m_gee_cnt(0), m_gen_cnt(0),
                                                  m_entity_repo(stk::percept::EntityRankEnd),
                                                  m_debug(false),
                                                  m_state(NRS_NONE)
      {
#if !PERCEPT_USE_PSEUDO_ELEMENTS
        m_useCustomGhosting = true;
#endif
        //m_comm_all( new stk::CommAll(eMesh.get_bulk_data()->parallel()) ),

#if NODE_REGISTRY_MAP_TYPE_GOOGLE
        //SubDimCell_SDSEntityType empty_key;
        //empty_key.insert( std::numeric_limits<stk::mesh::EntityId>::max() );
        //m_cell_2_data_map.set_empty_key(empty_key);

        SubDimCell_SDSEntityType deleted_key;
        deleted_key.insert( std::numeric_limits<stk::mesh::EntityId>::max() - 1u );
        m_cell_2_data_map.set_deleted_key(deleted_key);
#endif
      }

      ~NodeRegistry() {
        if (m_comm_all)
          delete m_comm_all;
      }

      void init_comm_all()
      {
        //std::cout << "tmp &m_eMesh = " << &m_eMesh << std::endl;
        if (m_comm_all)
          delete m_comm_all;
        m_comm_all = new stk::CommAll(m_eMesh.get_bulk_data()->parallel());
      }

      void init_entity_repo()
      {
        for (unsigned i = 0; i < stk::percept::EntityRankEnd; i++) m_entity_repo[i].clear();
      }
      
      void clear_dangling_nodes(SetOfEntities* nodes_to_be_deleted)
      {
        const bool debug = false;
        if (debug) std::cout <<  "tmp srk NodeRegistry::clear_dangling_nodes start" << std::endl;
        double cpu_0 = stk::cpu_time();

        SubDimCellToDataMap::iterator iter;
        SubDimCellToDataMap& map = getMap();

        std::vector<SubDimCell_SDSEntityType> to_erase;
        int num_delete=0;

        for (iter = map.begin(); iter != map.end(); ++iter)
          {
            SubDimCellData& nodeId_elementOwnderId = (*iter).second;
            NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
            VERIFY_OP_ON(nodeIds_onSE.size(), ==, nodeIds_onSE.m_entity_id_vector.size(), "NodeRegistry::clear_dangling_nodes id vector/size mismatch");
            unsigned nnodes = nodeIds_onSE.size();
            NodeIdsOnSubDimEntityType node_to_keep(0);
            //std::vector<stk::mesh::Entity *> node_to_keep;
            std::vector<stk::mesh::EntityId> node_id_to_keep(0);
            for (unsigned inode=0; inode < nnodes; inode++)
              {
                if (!nodeIds_onSE[inode]) continue;
                stk::mesh::EntityId id = nodeIds_onSE.m_entity_id_vector[inode];
                stk::mesh::EntityId id_check = nodeIds_onSE[inode]->identifier();
                VERIFY_OP_ON(id_check, ==, id, "NodeRegistry::clear_dangling_nodes id");

                //if (  stk::mesh::EntityLogDeleted == nodeIds_onSE[inode]->log_query() )
                if (nodes_to_be_deleted && nodes_to_be_deleted->find(nodeIds_onSE[inode]) != nodes_to_be_deleted->end())
                  {
                    ++num_delete;
                  }
                else if (!nodes_to_be_deleted && stk::mesh::EntityLogDeleted == nodeIds_onSE[inode]->log_query() )
                  {
                    ++num_delete;
                  }
                else
                  {
                    node_to_keep.push_back(nodeIds_onSE[inode]);
                    node_id_to_keep.push_back(id);
                  }
              }
            nodeIds_onSE = node_to_keep;
            nodeIds_onSE.m_entity_id_vector = node_id_to_keep;
            if (nodeIds_onSE.size() != nodeIds_onSE.m_entity_id_vector.size())
              {
                std::cout << "NodeRegistry::clear_dangling_nodes id vector/size mismatch 1 size= " << nodeIds_onSE.size() << " id.size= " << nodeIds_onSE.m_entity_id_vector.size() << std::endl;
              }
            VERIFY_OP_ON(nodeIds_onSE.size(), ==, nodeIds_onSE.m_entity_id_vector.size(), "NodeRegistry::clear_dangling_nodes id vector/size mismatch 1");

            if (nodeIds_onSE.size() == 0)
              to_erase.push_back(iter->first);

          }
        if (debug) std::cout << "tmp srk NodeRegistry::clear_dangling_nodes num_delete= " << num_delete <<  std::endl;
        if (to_erase.size())
          {
            if (debug) std::cout << "tmp srk NodeRegistry::clear_dangling_nodes nodeIds_onSE.size() != node_to_keep.size()), to_erase= " << to_erase.size() <<  std::endl;
            for (unsigned i=0; i < to_erase.size(); i++)
              {
                map.erase(to_erase[i]);
              }
          }

        // check
        if (1)
          {
            for (iter = map.begin(); iter != map.end(); ++iter)
              {
                SubDimCellData& nodeId_elementOwnderId = (*iter).second;
                NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
                VERIFY_OP_ON(nodeIds_onSE.size(), ==, nodeIds_onSE.m_entity_id_vector.size(), "NodeRegistry::clear_dangling_nodes id vector/size mismatch after erase");
              }
          }
        double cpu_1 = stk::cpu_time();
        if (debug) std::cout <<  "tmp srk NodeRegistry::clear_dangling_nodes end, time= " << (cpu_1-cpu_0) << std::endl;
      }

      void initialize()
      {
        //std::cout << "tmp &m_eMesh = " << &m_eMesh << std::endl;
        //delete m_comm_all;
        //m_comm_all = new stk::CommAll(m_eMesh.get_bulk_data()->parallel());
        m_cell_2_data_map.clear();
        init_entity_repo();
      }

      void //NodeRegistry::
      beginRegistration()
      {
        m_nodes_to_ghost.resize(0);
        m_pseudo_entities.clear();
        m_state = NRS_START_REGISTER_NODE;
        if (m_debug)
          std::cout << "P[" << m_eMesh.get_rank() << "] tmp NodeRegistry::beginRegistration" << std::endl;
      }

      void //NodeRegistry::
      endRegistration()
      {
        if (m_debug)
          std::cout << "P[" << m_eMesh.get_rank() << "] tmp NodeRegistry::endRegistration start" << std::endl;

        //putInESMap();

        removeUnmarkedSubDimEntities();

        m_eMesh.get_bulk_data()->modification_begin();
        this->createNewNodesInParallel();
        m_nodes_to_ghost.resize(0);

#if STK_ADAPT_NODEREGISTRY_DO_REHASH
        m_cell_2_data_map.rehash(m_cell_2_data_map.size());
#endif
        if (m_debug)
          std::cout << "P[" << m_eMesh.get_rank() << "] tmp NodeRegistry::endRegistration end" << std::endl;

        m_state = NRS_END_REGISTER_NODE;

      }

      void //NodeRegistry::
      beginLocalMeshMods()
      {
      }

      void //NodeRegistry::
      endLocalMeshMods()
      {
      }

      void //NodeRegistry::
      beginCheckForRemote()
      {
        m_state = NRS_START_CHECK_FOR_REMOTE;
        if (m_debug)
          std::cout << "P[" << m_eMesh.get_rank() << "] tmp NodeRegistry::beginCheckForRemote " << std::endl;
      }

      void //NodeRegistry::
      endCheckForRemote()
      {
        if (m_debug)
          std::cout << "P[" << m_eMesh.get_rank() << "] tmp NodeRegistry::endCheckForRemote start " << std::endl;
        stk::ParallelMachine pm = m_eMesh.get_bulk_data()->parallel();
        int failed = 0;
        stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );

        this->allocateBuffers();

#if STK_ADAPT_NODEREGISTRY_DO_REHASH
        m_cell_2_data_map.rehash(m_cell_2_data_map.size());
#endif

        if (m_debug)
          std::cout << "P[" << m_eMesh.get_rank() << "] tmp NodeRegistry::endCheckForRemote end " << std::endl;

        m_state = NRS_END_CHECK_FOR_REMOTE;

      }

      void //NodeRegistry::
      beginGetFromRemote()
      {
        m_state = NRS_START_GET_FROM_REMOTE;
        if (m_debug)
          std::cout << "P[" << m_eMesh.get_rank() << "] tmp NodeRegistry::beginGetFromRemote  " << std::endl;

      }
      void //NodeRegistry::
      endGetFromRemote()
      {
        if (m_debug)
          std::cout << "P[" << m_eMesh.get_rank() << "] tmp NodeRegistry::endGetFromRemote start " << std::endl;
        stk::ParallelMachine pm = m_eMesh.get_bulk_data()->parallel();
        int failed = 0;
        stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );

        this->communicate();

        failed = 0;
        stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );

        if (m_useCustomGhosting)
        {
          //std::cout << "m_useCustomGhosting= " << m_useCustomGhosting << std::endl;
          stk::mesh::Ghosting & ghosting = m_eMesh.get_bulk_data()->create_ghosting( std::string("new_nodes") );

          vector<stk::mesh::Entity*> receive;

          ghosting.receive_list( receive );
          //if (receive.size()) std::cout << "NodeRegistry::endGetFromRemote receive.size() = " << receive.size() << std::endl;

          m_eMesh.get_bulk_data()->change_ghosting( ghosting, m_nodes_to_ghost, receive);

        }

        failed = 0;
        stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );

#if PERCEPT_USE_PSEUDO_ELEMENTS
        if (!m_useCustomGhosting) setAllReceivedNodeData();
#endif

        m_eMesh.get_bulk_data()->modification_end();

        if (m_useCustomGhosting) setAllReceivedNodeData();

        if (0 && !m_useCustomGhosting)
          {
            m_eMesh.get_bulk_data()->modification_begin();
            removePseudoEntities();
            m_eMesh.get_bulk_data()->modification_end();
          }

        if (m_debug)
          std::cout << "P[" << m_eMesh.get_rank() << "] tmp NodeRegistry::endGetFromRemote end " << std::endl;

        m_state = NRS_END_GET_FROM_REMOTE;
      }

      void removePseudoEntities()
      {
        
        for (SetOfEntities::iterator it = m_pseudo_entities.begin(); it != m_pseudo_entities.end(); ++it)
          {
            stk::mesh::Entity *pseudo_elem = *it;
            bool did_destroy = m_eMesh.get_bulk_data()->destroy_entity(pseudo_elem);
            VERIFY_OP_ON(did_destroy, ==, true, "NodeRegistry::removePseudoEntities couldn't destroy");
          }
      }

      void setAllReceivedNodeData()
      {
        EXCEPTWATCH;
        //m_eMesh.get_bulk_data()->modification_begin();
        SubDimCellToDataMap::iterator iter;
        stk::mesh::PartVector empty_parts;

        for (iter = m_cell_2_data_map.begin(); iter != m_cell_2_data_map.end(); ++iter)
          {
            //const SubDimCell_SDSEntityType& subDimEntity = (*iter).first;
            SubDimCellData& nodeId_elementOwnderId = (*iter).second;

            NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();

            if (nodeIds_onSE.m_entity_id_vector.size() != nodeIds_onSE.size())
              {
                std::cout << " nodeIds_onSE.m_entity_id_vector.size() = " << nodeIds_onSE.m_entity_id_vector.size()
                          << " nodeIds_onSE.size() = " <<  nodeIds_onSE.size() << std::endl;

                throw std::logic_error("NodeRegistry:: setAllReceivedNodeData logic err #0");
              }

            for (unsigned ii = 0; ii < nodeIds_onSE.size(); ii++)
              {
                if (!nodeIds_onSE.m_entity_id_vector[ii])
                  {
                    nodeIds_onSE[ii] = 0;
                  }
                else
                  {
                    stk::mesh::Entity *node = get_entity_node_I(*m_eMesh.get_bulk_data(), stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE.m_entity_id_vector[ii]);  // FIXME
                    if (!node)
                      {
                        if (m_useCustomGhosting)
                          {
                            throw std::logic_error("NodeRegistry:: setAllReceivedNodeData logic err #3");
                          }
                        else
                          {
                            //std::cout << "tmp P[" << m_eMesh.get_rank() << "] NodeRegistry::setAllReceivedNodeData id= " << nodeIds_onSE.m_entity_id_vector[ii] << std::endl;
                            node = & m_eMesh.get_bulk_data()->declare_entity(m_eMesh.node_rank(), nodeIds_onSE.m_entity_id_vector[ii], empty_parts);
#if PERCEPT_USE_PSEUDO_ELEMENTS
                            stk::mesh::Entity *elem = & m_eMesh.get_bulk_data()->declare_entity(m_eMesh.element_rank()+PSEUDO_ELEMENT_RANK_SHIFT, 
                                                                                              nodeIds_onSE.m_entity_id_vector[ii]*PSEUDO_ELEMENT_MAGIC_NUMBER+m_eMesh.get_rank(), 
                                                                                              empty_parts);
                            m_pseudo_entities.insert(elem);
                            m_eMesh.get_bulk_data()->declare_relation(*elem, *node, 0);
#endif
                            if (!node) throw std::logic_error("NodeRegistry:: setAllReceivedNodeData logic err #3.1");
                          }
                      }
                    nodeIds_onSE[ii] = node;
                  }
              }
          }
        //m_eMesh.get_bulk_data()->modification_end();
      }

      /// when a sub-dim entity is visited during node registration but is flagged as not being marked, and thus not requiring 
      ///   any new nodes, we flag it with NR_MARK_NONE, then remove it here
      void removeUnmarkedSubDimEntities()
      {
        EXCEPTWATCH;
        SubDimCellToDataMap::iterator iter;

        for (iter = m_cell_2_data_map.begin(); iter != m_cell_2_data_map.end(); ++iter)
          {
            //const SubDimCell_SDSEntityType& subDimEntity = (*iter).first;
            SubDimCellData& nodeId_elementOwnderId = (*iter).second;

            NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
            if (nodeIds_onSE.size())
              {
                unsigned mark = nodeIds_onSE.m_mark;
                unsigned is_marked = mark & NR_MARK;
                unsigned is_not_marked = mark & NR_MARK_NONE;
#if 0
                if (!(is_marked || is_not_marked)) 
                  {
                    std::cout << "is_not_marked = " << is_not_marked << " is_marked= " << is_marked << std::endl;
                    throw std::logic_error("removeUnmarkedSubDimEntities::err1");
                  }
#endif
                //nodeIds_onSE.m_mark = 0u;
                if (!is_marked && is_not_marked)
                  {
                    //std::cout << "tmp SRK FOUND NR_MARK " << NR_MARK << " " << NR_MARK_NONE << std::endl;
                    nodeIds_onSE.resize(0);
                  }
              }
          }
      }

      bool is_empty( const stk::mesh::Entity& element, stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd);

      /// Register the need for a new node on the sub-dimensional entity @param subDimEntity on element @param element.
      /// If the element is a ghost element, the entity is still registered: the locality/ownership of the new entity
      /// can be determined by the locality of the element (ghost or not).
      bool registerNeedNewNode(const stk::mesh::Entity& element, NeededEntityType& needed_entity_rank, unsigned iSubDimOrd, bool needNodes)
      {
        static SubDimCell_SDSEntityType subDimEntity;
        getSubDimEntity(subDimEntity, element, needed_entity_rank.first, iSubDimOrd);

        static SubDimCellData new_SubDimCellData;
        static SubDimCellData empty_SubDimCellData;

        SubDimCellData* nodeId_elementOwnderId_ptr = getFromMapPtr(subDimEntity);
        SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
        bool is_empty = nodeId_elementOwnderId_ptr == 0;
        bool is_not_empty_but_data_cleared = (!is_empty && nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>().size() == 0);

        // if empty or if my id is the smallest, make this element the owner
        bool should_put_in =
          (element.identifier()  < stk::mesh::entity_id(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>()))
          || (element.entity_rank() > stk::mesh::entity_rank(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>()));

#define DEBUG_NR_UNREF 0
        if (DEBUG_NR_UNREF)
          {
            std::cout << "registerNeedNewNode:: is_empty= " << is_empty << " should_put_in= " << should_put_in << " needed_entity_rank= " 
                      << needed_entity_rank.first << " subDimEntity= ";
                          
            for (unsigned k=0; k < subDimEntity.size(); k++)
              {
                std::cout << " " << subDimEntity[k]->identifier() << " ";
              }
#if 0
            for (unsigned k=0; k < subDimEntity.size(); k++)
              {
                std::cout << " " << (stk::mesh::EntityId) subDimEntity[k] << " ";
              }
#endif
            std::cout << std::endl;
          }

        if (!is_empty)
          {
            unsigned& mark = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>().m_mark;
            if (needNodes)
              mark |= NR_MARK;
            else
              mark |= NR_MARK_NONE;
          }


        /// once it's in, the assertion should be:
        ///   owning_elementId < non_owning_elementId && owning_elementRank >= non_owning_elementRank
        ///
        if (is_empty || is_not_empty_but_data_cleared || should_put_in)
          {
            // new SubDimCellData SDC_DATA_OWNING_ELEMENT_KEY
            // CHECK

            unsigned numNewNodes = needed_entity_rank.second;

            NodeIdsOnSubDimEntityType nid_new(numNewNodes, 0u);
            if (needNodes)
              nid_new.m_mark |= NR_MARK;
            else
              nid_new.m_mark |= NR_MARK_NONE;

            SubDimCellData data(nid_new, stk::mesh::EntityKey(element.entity_rank(), element.identifier()) );
            putInMap(subDimEntity,  data);

            if (0 && DEBUG_NR_UNREF)
              {
                std::cout << "registerNeedNewNode:: is_empty= " << is_empty << " should_put_in= " << should_put_in << " subDimEntity= ";
                for (unsigned k=0; k < subDimEntity.size(); k++)
                  {
                    std::cout << " " << subDimEntity[k]->identifier() << " ";
                  }
#if 0
                for (unsigned k=0; k < subDimEntity.size(); k++)
                  {
                    std::cout << " " << (stk::mesh::EntityId) subDimEntity[k] << " ";
                  }
#endif
                std::cout << std::endl;
              }

            // new_SubDimCellData.get<0>().init(needed_entity_rank.second, 0u);
            // new_SubDimCellData.get<1>() = newEntityKey;
            //m_cell_2_data_map.insert(std::pair<SubDimCell_SDSEntityType, SubDimCellData>(subDimEntity, new_SubDimCellData) );
            //m_cell_2_data_map.insert(std::pair<SubDimCell_SDSEntityType, SubDimCellData>(subDimEntity, data) );

            return true;
          }
        return false;
      }

      /// check the newly registered node from the registry, which does one of three things, depending on what mode we are in:
      ///   1. counts buffer in prep for sending (just does a pack)
      ///   2. packs the buffer (after buffers are alloc'd)
      ///   3. returns the new node after all communications are done
      bool checkForRemote(const stk::mesh::Entity& element, NeededEntityType& needed_entity_rank, unsigned iSubDimOrd, bool needNodes_notUsed)
      {
        EXCEPTWATCH;
        static SubDimCellData empty_SubDimCellData;
        static CommDataType buffer_entry;

        bool isGhost = m_eMesh.isGhostElement(element);

        if (!isGhost) return true;

        static SubDimCell_SDSEntityType subDimEntity;
        getSubDimEntity(subDimEntity, element, needed_entity_rank.first, iSubDimOrd);

        stk::CommAll& comm_all = *m_comm_all;
        unsigned proc_size = comm_all.parallel_size();
        unsigned proc_rank = comm_all.parallel_rank();
        unsigned owner_proc_rank = element.owner_rank();


        SubDimCellData* nodeId_elementOwnderId_ptr = getFromMapPtr(subDimEntity);
        SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
        bool is_empty = nodeId_elementOwnderId_ptr == 0;

        if (is_empty)
          {
            std::cout << "element= " << element
                      << " needed_entity_rank= " << needed_entity_rank.first<< " " << needed_entity_rank.second << std::endl;
            std::cout << "subDimEntity= " << subDimEntity << std::endl;
            std::cout << "nodeId_elementOwnderId= " << nodeId_elementOwnderId << std::endl;
            std::cout << "empty_SubDimCellData= " << empty_SubDimCellData << std::endl;
            throw std::logic_error("NodeRegistry::checkForRemote no data (is_empty=true) - logic error.");
            return false;
          }
        else
          {

            unsigned owning_elementId = stk::mesh::entity_id(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>());
            unsigned owning_elementRank = stk::mesh::entity_rank(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>());

            NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
            unsigned nidsz = nodeIds_onSE.size();

            // error check
            bool isNotOK = (element.identifier() < owning_elementId) && (element.entity_rank() > owning_elementRank) ;

            if ( isNotOK )
              {
                std::cout << "P[" << proc_rank << "] elem id = " << element.identifier()
                          << " nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>() = "
                          << owning_elementId
                          << std::endl;
                throw std::logic_error("NodeRegistry::checkForRemote logic: owning element info is wrong");
              }

            unsigned erank = m_eMesh.element_rank();
            erank = owning_elementRank;
            //VERIFY_OP(erank, <=, owning_elementRank , "erank...");
            stk::mesh::Entity * owning_element = get_entity_element(*m_eMesh.get_bulk_data(), erank, owning_elementId);

            if (!owning_element)
              throw std::logic_error("NodeRegistry::checkForRemote logic: owning_element is null");

            bool owning_element_is_ghost = m_eMesh.isGhostElement(*owning_element);

            // if this element is a ghost, and the owning element of the node is not a ghost, send info
            //   to ghost element's owner proc
            if (!owning_element_is_ghost && isGhost)
              {
                buffer_entry = CommDataType(
                                            needed_entity_rank.first,
                                            iSubDimOrd,
                                            element.key()
                                            );

                if (nodeIds_onSE.m_entity_id_vector.size() != nodeIds_onSE.size())
                  {
                    throw std::logic_error("NodeRegistry::checkForRemote logic err #0.1");
                  }

                for (unsigned iid = 0; iid < nidsz; iid++)
                  {
                    if (nodeIds_onSE[iid] == 0)
                      {
                        throw std::logic_error("logic: hmmm #5.0");
                      }

                    if (!nodeIds_onSE.m_entity_id_vector[iid])
                      {
                        throw std::logic_error("NodeRegistry::checkForRemote logic err #0.2");
                      }

                    //stk::mesh::Entity * new_node = get_entity_node_Ia(*m_eMesh.get_bulk_data(), Node, nodeIds_onSE, iid);
                    stk::mesh::Entity * new_node = nodeIds_onSE[iid];

                    if (0)
                      {
                        stk::mesh::Entity * new_node_1 = get_entity_node_I(*m_eMesh.get_bulk_data(), stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE.m_entity_id_vector[iid]);
                        if (new_node != new_node_1)
                          {
                            throw std::logic_error("NodeRegistry::checkForRemote logic err #0.3");
                          }
                      }

                    if (!new_node)
                      {
                        throw std::logic_error("NodeRegistry::checkForRemote logic: new_node is null");
                      }

                    m_nodes_to_ghost.push_back( stk::mesh::EntityProc(new_node, owner_proc_rank) );
                  }

                if (isParallelRun(proc_size))
                  {
                    //std::cout << "P[" << proc_rank << "] : pack " << buffer_entry << " owner_proc_rank= " << owner_proc_rank << std::endl;
                    m_comm_all->send_buffer( owner_proc_rank ).pack< CommDataType > (buffer_entry);
                    NodeIdsOnSubDimEntityType& nids = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
                    nids.pack(m_comm_all->send_buffer( owner_proc_rank ));
                  }
                else
                  {
                    // FIXME createNodeAndConnect(buffer_entry);
                  }
              }
          }
        return true; // FIXME
      }

      bool getFromRemote(const stk::mesh::Entity& element, NeededEntityType& needed_entity_rank, unsigned iSubDimOrd, bool needNodes_notUsed)
      {
        return checkForRemote(element, needed_entity_rank, iSubDimOrd, needNodes_notUsed);
      }


#define LOCAL_INLINE inline

      LOCAL_INLINE stk::mesh::Entity* get_entity_using_find(stk::mesh::EntityRank& rank, const stk::mesh::EntityId& id) const
      {
        const EntityRepo::const_iterator i = m_entity_repo[rank].find( id );
        return i != m_entity_repo[rank].end() ? i->second : NULL ;
      }

      LOCAL_INLINE stk::mesh::Entity *get_entity(stk::mesh::BulkData& bulk, stk::mesh::EntityRank rank, stk::mesh::EntityId id)
      {

#if STK_ADAPT_NODEREGISTRY_USE_ENTITY_REPO
        stk::mesh::Entity* entity = get_entity_using_find(rank, id);
        if (entity)
          return entity;
        else
          {
            entity = bulk.get_entity(rank, id);
            m_entity_repo[rank][id] = entity;
            return entity;
          }
#else
        return bulk.get_entity(rank, id);
#endif
      }

      LOCAL_INLINE stk::mesh::Entity *get_entity_I(stk::mesh::BulkData& bulk, stk::mesh::EntityRank rank, stk::mesh::EntityId id)
      {

#if STK_ADAPT_NODEREGISTRY_USE_ENTITY_REPO
        stk::mesh::Entity* entity = get_entity_using_find(rank, id);
        if (entity)
          return entity;
        else
          {
            entity = bulk.get_entity(rank, id);
            m_entity_repo[rank][id] = entity;
            return entity;
          }
#else
        return bulk.get_entity(rank, id);
#endif
      }

      LOCAL_INLINE stk::mesh::Entity *get_entity_Ia(stk::mesh::BulkData& bulk, stk::mesh::EntityRank rank, stk::mesh::EntityId id)
      {

#if STK_ADAPT_NODEREGISTRY_USE_ENTITY_REPO
        stk::mesh::Entity* entity = get_entity_using_find(rank, id);
        if (entity)
          return entity;
        else
          {
            entity = bulk.get_entity(rank, id);
            m_entity_repo[rank][id] = entity;
            return entity;
          }
#else
        return bulk.get_entity(rank, id);
#endif
      }

      LOCAL_INLINE stk::mesh::Entity *get_entity_Ib(stk::mesh::BulkData& bulk, stk::mesh::EntityRank rank, stk::mesh::EntityId id)
      {
#if STK_ADAPT_NODEREGISTRY_USE_ENTITY_REPO
        stk::mesh::Entity* entity = get_entity_using_find(rank, id);
        if (entity)
          return entity;
        else
          {
            entity = bulk.get_entity(rank, id);
            m_entity_repo[rank][id] = entity;
            return entity;
          }
#else
        return bulk.get_entity(rank, id);
#endif
      }


      LOCAL_INLINE stk::mesh::Entity *get_entity_element(stk::mesh::BulkData& bulk, stk::mesh::EntityRank rank, stk::mesh::EntityId id)
      {
        //m_gee_cnt++;
        return get_entity(bulk, rank, id);
      }

      LOCAL_INLINE stk::mesh::Entity *get_entity_node_I(stk::mesh::BulkData& bulk, stk::mesh::EntityRank rank, stk::mesh::EntityId id)
      {
        //m_gen_cnt++;
        return get_entity_I(bulk, rank, id);
      }

      LOCAL_INLINE stk::mesh::Entity *get_entity_node_Ia(stk::mesh::BulkData& bulk, stk::mesh::EntityRank rank, NodeIdsOnSubDimEntityType& nodeIds_onSE, unsigned index)
      {
        //m_gen_cnt++;
        //stk::mesh::Entity * entity = get_entity(bulk, rank, nodeIds_onSE.m_entity_id_vector[index]);
        stk::mesh::Entity * entity = nodeIds_onSE[index];
        return entity;
      }

      LOCAL_INLINE stk::mesh::Entity *get_entity_node_Ib(stk::mesh::BulkData& bulk, stk::mesh::EntityRank rank, stk::mesh::EntityId id)
      {
        //m_gen_cnt++;
        return get_entity_Ib(bulk, rank, id);
      }

    NodeIdsOnSubDimEntityType* getNewNodesOnSubDimEntity(const stk::mesh::Entity& element,  stk::mesh::EntityRank& needed_entity_rank, unsigned iSubDimOrd);

      /// makes coordinates of this new node be the centroid of its sub entity
      void makeCentroidCoords(const stk::mesh::Entity& element,  stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd)
      {
        makeCentroidField(element, needed_entity_rank, iSubDimOrd, m_eMesh.get_coordinates_field());
      }

      void makeCentroidField(const stk::mesh::Entity& element,  stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd, stk::mesh::FieldBase *field)
      {
        //EXCEPTWATCH;

        int spatialDim = m_eMesh.get_spatial_dim();
        stk::mesh::EntityRank field_rank = stk::mesh::fem::FEMMetaData::NODE_RANK;
        {
          unsigned nfr = field->restrictions().size();
          //if (print_info) std::cout << "P[" << p_rank << "] info>    number of field restrictions= " << nfr << std::endl;
          for (unsigned ifr = 0; ifr < nfr; ifr++)
            {
              const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
              //mesh::Part& frpart = metaData.get_part(fr.ordinal());
              field_rank = fr.entity_rank();
              spatialDim = fr.dimension() ;
            }
        }

        if (field_rank != stk::mesh::fem::FEMMetaData::NODE_RANK)
          {
            return;
          }

        unsigned *null_u = 0;
        static SubDimCell_SDSEntityType subDimEntity;
        //subDimEntity.clear();
        getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
        static SubDimCellData empty_SubDimCellData;

        SubDimCellData* nodeId_elementOwnderId_ptr = getFromMapPtr(subDimEntity);
        SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
        bool is_empty = nodeId_elementOwnderId_ptr == 0;

        if (s_allow_empty_sub_dims && is_empty)
          {
            return;
          }
        if (is_empty)
          {
            const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);
            shards::CellTopology cell_topo(cell_topo_data);

            std::cout << "NodeRegistry::makeCentroidField: no node found, cell_topo = " << cell_topo.getName()
                      << "\n subDimEntity= " << subDimEntity
                      << "\n element= " << element
                      << "\n element.entity_rank() = " << element.entity_rank()
                      << "\n needed_entity_rank= " << needed_entity_rank
                      << "\n iSubDimOrd= " << iSubDimOrd << std::endl;
            throw std::runtime_error("makeCentroidField: no node found");
          }
        NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
        if (nodeIds_onSE.size() != 1)
          throw std::runtime_error("makeCentroidField not ready for multiple nodes");
        //stk::mesh::Entity * c_node = m_eMesh.get_bulk_data()->get_entity(Node, nodeIds_onSE[0]);
        //stk::mesh::Entity * c_node = get_entity_node(*m_eMesh.get_bulk_data(), stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE[0]);
        //stk::mesh::Entity * c_node = get_entity_node_Ia(*m_eMesh.get_bulk_data(), stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE, 0u);
        stk::mesh::Entity * c_node = nodeIds_onSE[0];


        if (!c_node)
          {
            throw std::runtime_error("makeCentroidField: bad node found 0");
          }

        //std::vector<double> c_p(spatialDim, 0.0);
        double c_p[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        bool doPrint = false;

        if (needed_entity_rank == m_eMesh.element_rank())
          {
            const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);
            unsigned npts = elem_nodes.size();
            //if (npts == 2) doPrint=true;
            double dnpts = elem_nodes.size();
            for (unsigned ipts = 0; ipts < npts; ipts++)
              {
                //stk::mesh::EntityId nodeId = elem_nodes[ipts];
                stk::mesh::Entity * node = elem_nodes[ipts].entity(); //m_eMesh.get_bulk_data()->get_entity(Node, nodeId);
                if (!node)
                  {
                    throw std::runtime_error("makeCentroidField: bad node found 1");
                  }
                //double * const coord = stk::mesh::field_data( *field , *node );
                double *  coord = m_eMesh.field_data(field, *node, null_u);

                if (doPrint && coord)
                  {
                    const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);
                    shards::CellTopology cell_topo(cell_topo_data);

                    std::cout << "tmp NodeRegistry::makeCentroidField cell_topo = " << cell_topo.getName() << " ipts= " << ipts
                              << " coord= " << coord[0] << " " << coord[1] << " " << coord[2] << std::endl;
                  }

                if (coord)
                  {
                    for (int isp = 0; isp < spatialDim; isp++)
                      {
                        c_p[isp] += coord[isp]/dnpts;
                      }
                  }
              }

          }
        else
          {
            double dnpts = subDimEntity.size();
            for (SubDimCell_SDSEntityType::iterator ids = subDimEntity.begin(); ids != subDimEntity.end(); ids++)
              {
                SDSEntityType nodeId = *ids;

                //stk::mesh::Entity * node = m_eMesh.get_bulk_data()->get_entity(Node, nodeId);
                //!!stk::mesh::Entity * node = get_entity_node_II(*m_eMesh.get_bulk_data(),Node, nodeId);
                stk::mesh::Entity * node = nodeId;
                if (!node)
                  {
                    throw std::runtime_error("makeCentroidField: bad node found 2");
                  }
                //double * const coord = stk::mesh::field_data( *field, *node );
                double *  coord = m_eMesh.field_data(field, *node, null_u);
                if (coord)
                  {
                    for (int isp = 0; isp < spatialDim; isp++)
                      {
                        c_p[isp] += coord[isp]/dnpts;
                      }
                  }
              }
          }
        //double * c_coord = stk::mesh::field_data( *field , *c_node );
        double *  c_coord = m_eMesh.field_data(field, *c_node, null_u);
        if (c_coord)
          {
            //std::string coord_str;
            for (int isp = 0; isp < spatialDim; isp++)
              {
                c_coord[isp] = c_p[isp];
                //if (doPrint)
                //  coord_str += toString(c_coord[isp])+ " ";
              }
            if (doPrint)
              {
                const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);
                shards::CellTopology cell_topo(cell_topo_data);

                std::cout << "tmp NodeRegistry::makeCentroidField cell_topo = " << cell_topo.getName()
                      << "\n subDimEntity= " << subDimEntity
                      << "\n element= " << element
                      << "\n element.entity_rank() = " << element.entity_rank()
                      << "\n needed_entity_rank= " << needed_entity_rank
                      << "\n iSubDimOrd= " << iSubDimOrd << std::endl;
                //std::cout << "P[" << m_eMesh.get_rank() << "] needed_entity_rank= " << needed_entity_rank << " coord= " << coord_str << std::endl;
              }
          }

      }

      static void normalize_spacing(unsigned nsz, unsigned nsp, double spc[8][3], double den_xyz[3])
      {
        //double m_min_spacing_factor = 0.5;
        double m_min_spacing_factor = 0.0;
        for (unsigned isp = 0; isp < nsp; isp++)
          {
            double den = 0.0;
            unsigned ipts=0;
            for (ipts=0; ipts < nsz; ipts++)
              {
                spc[ipts][isp] = 1.0/spc[ipts][isp]; 
                den += spc[ipts][isp];
              }
            for (ipts=0; ipts < nsz; ipts++)
              {
                spc[ipts][isp] /= den;
              }
            // now it's a fraction [0,1], check if it's too big
            for (ipts=0; ipts < nsz; ipts++)
              {
                if ( spc[ipts][isp] > 1.0 - m_min_spacing_factor)
                  {
                    spc[ipts][isp] = 1.0 - m_min_spacing_factor;
                    for (unsigned jpts=0; jpts < nsz; jpts++)
                      {
                        if (ipts != jpts) 
                          spc[ipts][isp] = m_min_spacing_factor/((double)(nsz-1));
                      }

                    break;
                  }
              }         
            // now renormalize it
            den = 0.0;
            for (ipts=0; ipts < nsz; ipts++)
              {
                den += spc[ipts][isp];
              }
            for (ipts=0; ipts < nsz; ipts++)
              {
                spc[ipts][isp] /= den;
              }

          }
      }

      /// makes coordinates of this new node be the centroid of its sub entity - this version does it for all new nodes
      void makeCentroid(stk::mesh::FieldBase *field)
      {
        EXCEPTWATCH;
        //unsigned *null_u = 0;
        stk::mesh::FieldBase *spacing_field    = m_eMesh.get_field("ref_spacing_field");

        int spatialDim = m_eMesh.get_spatial_dim();
        stk::mesh::EntityRank field_rank = stk::mesh::fem::FEMMetaData::NODE_RANK;
        {
          EXCEPTWATCH;
          unsigned nfr = field->restrictions().size();
          //if (print_info) std::cout << "P[" << p_rank << "] info>    number of field restrictions= " << nfr << std::endl;
          for (unsigned ifr = 0; ifr < nfr; ifr++)
            {
              const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
              //mesh::Part& frpart = metaData.get_part(fr.ordinal());
              field_rank = fr.entity_rank();
              spatialDim = fr.dimension() ;
            }
        }
        if (field_rank != stk::mesh::fem::FEMMetaData::NODE_RANK)
          {
            if (field_rank == m_eMesh.element_rank())
              {
                
              }
            return;
          }

        SubDimCellToDataMap::iterator iter;

        for (iter = m_cell_2_data_map.begin(); iter != m_cell_2_data_map.end(); ++iter)
          {
            const SubDimCell_SDSEntityType& subDimEntity = (*iter).first;
            SubDimCellData& nodeId_elementOwnderId = (*iter).second;

            NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
            static const SubDimCellData empty_SubDimCellData;

            bool is_empty = (nodeId_elementOwnderId == empty_SubDimCellData);

            if (s_allow_empty_sub_dims && is_empty)
              {
                return;
              }

            if (is_empty)
              {
                throw std::runtime_error("makeCentroid(field) empty cell found");
              }

            if (nodeIds_onSE.size() != 1)
              {
                continue;
              }

            stk::mesh::EntityRank needed_entity_rank = stk::mesh::fem::FEMMetaData::NODE_RANK;
            // SPECIAL CASE
            // SPECIAL CASE
            // SPECIAL CASE
            if (subDimEntity.size() == 1)
              {
                needed_entity_rank = m_eMesh.element_rank();
              }

            if (nodeIds_onSE[0] == 0)
              {
                continue;
              }

            stk::mesh::Entity * c_node = nodeIds_onSE[0];

            if (!c_node)
              {
                throw std::runtime_error("makeCentroid(field): bad node found 0.0");
              }

            double c_p[] = {0.0, 0.0, 0.0};
            bool doPrint = false;
            std::vector<stk::mesh::Entity *> nodes(8,(stk::mesh::Entity *)0);
            unsigned nsz = 0;

            if (needed_entity_rank == m_eMesh.element_rank())
              {
                EXCEPTWATCH;
                stk::mesh::Entity *element_p = 0;
                {
                  SDSEntityType elementId = *subDimEntity.begin();
                  //!!element_p = get_entity_element(*m_eMesh.get_bulk_data(), m_eMesh.element_rank(), elementId);
                  element_p = elementId;
                  if (!element_p)
                    {
                      throw std::runtime_error("makeCentroid(field): bad elem found 2");
                    }
                }

                stk::mesh::Entity& element = *element_p;
                bool element_is_ghost = m_eMesh.isGhostElement(element);
                if (element_is_ghost)
                  {
                    //std::cout << "tmp found ghost" << std::endl;
                  }
                else
                  {
                    const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);
                    unsigned npts = elem_nodes.size();
                    nsz = npts;
                    nodes.resize(nsz, (stk::mesh::Entity *)0);
                    //if (npts == 2) doPrint=true;
                    //double dnpts = elem_nodes.size();
                    for (unsigned ipts = 0; ipts < npts; ipts++)
                      {
                        stk::mesh::Entity * node = elem_nodes[ipts].entity();
                        if (!node)
                          {
                            throw std::runtime_error("makeCentroid(field): bad node found 1.0");
                          }
                        nodes[ipts] = node;
                      }
                  }
              }
            else
              {
                unsigned ipts=0;
                nsz = subDimEntity.size();
                nodes.resize(nsz, (stk::mesh::Entity *)0);

                for (SubDimCell_SDSEntityType::const_iterator ids = subDimEntity.begin(); ids != subDimEntity.end(); ++ids, ++ipts)
                  {
                    SDSEntityType nodeId = *ids;
                    stk::mesh::Entity * node = nodeId;
                    nodes[ipts]=node;
                  }
              }

            {
              //if ( (spacing_field && (spacing_field != field) && subDimEntity.size() == 2))
              bool do_spacing=false;
              if (do_spacing && (spacing_field && (spacing_field != field) ) )
                {
                  EXCEPTWATCH;
                  unsigned ipts=0;

                  double * coord[8] = {0,0,0,0,0,0,0,0};
                  double * spacing[8] = {0,0,0,0,0,0,0,0};

                  for (ipts=0; ipts < nsz; ipts++)
                    {
                      coord[ipts] = m_eMesh.field_data_inlined(field, *nodes[ipts]);
                      spacing[ipts] = m_eMesh.field_data_inlined(spacing_field, *nodes[ipts]);
                    }

                  double spc[8][3] = {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
                  double den = 0.0;
                  double den_xyz[3] = {0,0,0};
                  if (nsz == 2) 
                    {
                      den = 0.0;
                      for (ipts=0; ipts < nsz; ipts++)
                        {
                          for (int isp = 0; isp < spatialDim; isp++)
                            {
                              spc[ipts][0] += (coord[1][isp] - coord[0][isp])*spacing[ipts][isp];
                            }                            
                        }
                      for (ipts=0; ipts < nsz; ipts++)
                        {
                          for (int isp = 1; isp < spatialDim; isp++)
                            {
                              spc[ipts][isp] = spc[ipts][0];
                            }                            
                        }
                    }
                  else
                    {
                      for (ipts=0; ipts < nsz; ipts++)
                        {
                          for (int isp = 0; isp < spatialDim; isp++)
                            {
                              spc[ipts][isp] = spacing[ipts][isp];
                            }
                        }
                    }
                  normalize_spacing(nsz, spatialDim, spc, den_xyz);
                  if (nsz==2 && (coord[1][0] < 1.e-3 && coord[0][0] < 1.e-3))
                    for (ipts=0; ipts < nsz; ipts++)
                      for (int isp = 0; isp < spatialDim; isp++)
                        {
                          std::cout << "y = " << coord[ipts][1] << " new spc[" << ipts << "]= " << spc[ipts][1] << std::endl;
                        }


                  for (ipts=0; ipts < nsz; ipts++)
                    {
                      for (int isp = 0; isp < spatialDim; isp++)
                        {
                          c_p[isp] += coord[ipts][isp]*spc[ipts][isp];
                        }
                    }

                }
              else
                {
                  EXCEPTWATCH;
                  double dnpts = nsz;
                  unsigned ipts=0;
                  for (ipts=0; ipts < nsz; ipts++)
                    {
                      stk::mesh::Entity * node = nodes[ipts];
                      if (!node)
                        {
                          throw std::runtime_error("makeCentroid(field): bad node found 2.0");
                        }
                      //double *  coord = m_eMesh.field_data(field, *node, null_u);
                      double *  coord = m_eMesh.field_data_inlined(field, *node);

                      if (doPrint && coord)
                        {
                          //const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);
                          //CellTopology cell_topo(cell_topo_data);

                          std::cout << "tmp NodeRegistry::makeCentroid(field) npts= " << subDimEntity.size() << " ipts= " << ipts
                                    << " coord= " << coord[0] << " " << coord[1] << " " << coord[2] << std::endl;
                        }

                      if (coord)
                        {
                          for (int isp = 0; isp < spatialDim; isp++)
                            {
                              c_p[isp] += coord[isp]/dnpts;
                            }
                        }
                    }
                }
            }

            // set coords
            {
              EXCEPTWATCH;

              //double *  c_coord = m_eMesh.field_data(field, *c_node, null_u);
              double *  c_coord = m_eMesh.field_data_inlined(field, *c_node);

              if (c_coord)
                {
                  for (int isp = 0; isp < spatialDim; isp++)
                    {
                      c_coord[isp] = c_p[isp];
                    }

                  if (doPrint)
                    std::cout << "tmp NodeRegistry::makeCentroid(field) c_coord= " << c_coord[0] << " " << c_coord[1] << " " << c_coord[2] << std::endl;


                }
            }
          }
      } // makeCentroid(stk::mesh::FieldBase *)

      /// do interpolation for all fields
      void interpolateFields(const stk::mesh::Entity& element,  stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd)
      {
        const stk::mesh::FieldVector & fields = m_eMesh.get_fem_meta_data()->get_fields();
        unsigned nfields = fields.size();
        //std::cout << "P[" << p_rank << "] info>    Number of fields = " << fields.size() << std::endl;
        for (unsigned ifld = 0; ifld < nfields; ifld++)
          {
            stk::mesh::FieldBase *field = fields[ifld];
            //std::cout << "P[" << m_eMesh.get_rank() << "] field = " << field->name() << std::endl;
            makeCentroidField(element, needed_entity_rank, iSubDimOrd, field);
          }
      }

      /// do interpolation for all fields
      void interpolateFields()
      {
        const stk::mesh::FieldVector & fields = m_eMesh.get_fem_meta_data()->get_fields();
        unsigned nfields = fields.size();
        //std::cout << "P[" << p_rank << "] info>    Number of fields = " << fields.size() << std::endl;
        for (unsigned ifld = 0; ifld < nfields; ifld++)
          {
            stk::mesh::FieldBase *field = fields[ifld];
            //std::cout << "P[" << m_eMesh.get_rank() << "] field = " << field->name() << std::endl;
            makeCentroid(field);
          }
      }


      /// check for adding new nodes to existing parts based on sub-entity part ownership

      void addToExistingParts(const stk::mesh::Entity& element,  stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd)
      {
        const std::vector< stk::mesh::Part * > & parts = m_eMesh.get_fem_meta_data()->get_parts();

        unsigned nparts = parts.size();

        //CHECK
        static SubDimCell_SDSEntityType subDimEntity;
        //subDimEntity.clear();
        getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
        static  SubDimCellData empty_SubDimCellData;
        SubDimCellData* nodeId_elementOwnderId_ptr = getFromMapPtr(subDimEntity);
        SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
        bool is_empty = nodeId_elementOwnderId_ptr == 0;

        if (s_allow_empty_sub_dims && is_empty)
          {
            return;
          }

        if (is_empty)
          {
            throw std::runtime_error("addToExistingParts: no node found");
          }
        NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
        unsigned nidsz = nodeIds_onSE.size();

        for (unsigned i_nid = 0; i_nid < nidsz; i_nid++)
          {
            //stk::mesh::Entity * c_node = get_entity_node(*m_eMesh.get_bulk_data(), stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE[i_nid]);
            //stk::mesh::Entity * c_node = get_entity_node_Ia(*m_eMesh.get_bulk_data(), stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE, i_nid);
            stk::mesh::Entity * c_node = nodeIds_onSE[i_nid];

            if (!c_node)
              {
                std::cout << "addToExistingParts: " <<  nodeIds_onSE[i_nid] << " i_nid= " << i_nid << " nidsz= " << nidsz
                          << " needed_entity_rank= " << needed_entity_rank << " iSubDimOrd= " << iSubDimOrd << std::endl;
                throw std::runtime_error("addToExistingParts: bad node found 0.1");
              }

            for (unsigned ipart=0; ipart < nparts; ipart++)
              {
                stk::mesh::Part& part = *parts[ipart];
                stk::mesh::Selector selector(part);

                //std::cout << "P[" << m_eMesh.get_rank() << "] NodeRegistry::addToExistingParts Part[" << ipart << "]= " << part.name() << std::endl;
                //std::string part_name = part.name();

                // FIXME - is there a better way to determine if a part is one of the "standard" parts?
                if (stk::mesh::is_auto_declared_part(part))  //part_name[0] == '{')  //  is_auto_declared_part
                  continue;

                //std::cout << "P[" << p_rank << "] info>     Part[" << ipart << "]= " << part.name()
                //              << " topology = " << (topology?shards::CellTopology(topology).getName():"null")
                //              << std::endl;


                bool found = true;
                if (needed_entity_rank == m_eMesh.element_rank())
                  {
                    const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);
                    unsigned npts = elem_nodes.size();
                    for (unsigned ipts = 0; ipts < npts; ipts++)
                      {
                        stk::mesh::Entity * node = elem_nodes[ipts].entity();
                        if (!node)
                          {
                            throw std::runtime_error("addToExistingParts: bad node found 1.1");
                          }
                        if (!selector(*node))
                          {
                            found = false;
                            break;
                          }

                      }
                  }
                else
                  {
                    //double dnpts = subDimEntity.size();
                    for (SubDimCell_SDSEntityType::iterator ids = subDimEntity.begin(); ids != subDimEntity.end(); ++ids)
                      {
                        SDSEntityType nodeId = *ids;
                        //!!Entity * node = get_entity_node_II(*m_eMesh.get_bulk_data(),Node, nodeId);
                        stk::mesh::Entity * node = nodeId;
                        if (!node)
                          {
                            throw std::runtime_error("addToExistingParts: bad node found 2.1");
                          }
                        if (!selector(*node))
                          {
                            found = false;
                            break;
                          }
                      }
                  }
                if (found)
                  {
                    // add to part
                    std::vector<stk::mesh::Part*> add_parts(1, &part);
                    std::vector<stk::mesh::Part*> remove_parts;
                    const CellTopologyData *const topology = stk::percept::PerceptMesh::get_cell_topology(part);
                    const unsigned part_rank = part.primary_entity_rank();

                    //if (!topology)
                    if (part_rank == stk::mesh::fem::FEMMetaData::NODE_RANK)
                      {
                        m_eMesh.get_bulk_data()->change_entity_parts( *c_node, add_parts, remove_parts );
                        if (0)
                          {
                            std::cout << "P[" << m_eMesh.get_rank() << "] adding node " << c_node->identifier() << " to   Part[" << ipart << "]= " << part.name()
                                      << " topology = " << (topology ? shards::CellTopology(topology).getName() : "null")
                                      << std::endl;
                          }
                      }

                  }
              }
          }
      }

      /// Check for adding new nodes to existing parts based on sub-entity part ownership.
      /// This version does it in bulk and thus avoids repeats on shared sub-dim entities.

      void addToExistingPartsNew()
      {
        static std::vector<stk::mesh::Part*> add_parts(1, static_cast<stk::mesh::Part*>(0));
        static std::vector<stk::mesh::Part*> remove_parts;

        //std::cout << "tmp addToExistingPartsNew... " << std::endl;
        const std::vector< stk::mesh::Part * > & parts = m_eMesh.get_fem_meta_data()->get_parts();

        unsigned nparts = parts.size();
        for (unsigned ipart=0; ipart < nparts; ipart++)
          {
            stk::mesh::Part& part = *parts[ipart];

            //std::cout << "P[" << m_eMesh.get_rank() << "] NodeRegistry::addToExistingParts Part[" << ipart << "]= " << part.name() << std::endl;
            //std::string part_name = part.name();

            // FIXME - is there a better way to determine if a part is one of the "standard" parts?
            if (stk::mesh::is_auto_declared_part(part)) //part_name[0] == '{')  // is_auto_declared_part
              continue;

            const CellTopologyData *const topology = stk::percept::PerceptMesh::get_cell_topology(part);
            const unsigned part_rank = part.primary_entity_rank();

            if (part_rank == stk::mesh::fem::FEMMetaData::NODE_RANK)
              {
                stk::mesh::Selector selector(part);

                //std::cout << "P[" << m_eMesh.get_rank() << "] NodeRegistry::addToExistingPartsNew rank=Node = Part[" << ipart << "]= " << part.name() << std::endl;
                add_parts[0] = &part;

                SubDimCellToDataMap::iterator iter;
                //std::cout << "tmp m_cell_2_data_map.size() = " << m_cell_2_data_map.size() << std::endl;
                for (iter = m_cell_2_data_map.begin(); iter != m_cell_2_data_map.end(); ++iter)
                  {
                    const SubDimCell_SDSEntityType& subDimEntity = (*iter).first;
                    SubDimCellData& nodeId_elementOwnderId = (*iter).second;

                    NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
                    //std::cout << "P[" << p_rank << "] info>     Part[" << ipart << "]= " << part.name()
                    //              << " topology = " << (topology?shards::CellTopology(topology).getName():"null")
                    //              << std::endl;

                    bool found = true;
                    stk::mesh::EntityRank needed_entity_rank = stk::mesh::fem::FEMMetaData::NODE_RANK;
                    //
                    // SPECIAL CASE
                    // SPECIAL CASE
                    // SPECIAL CASE
                    //
                    if( subDimEntity.size() == 1)
                      {
                        needed_entity_rank = m_eMesh.element_rank();
                      }

                    if (needed_entity_rank == m_eMesh.element_rank())
                      {
                        stk::mesh::Entity *element_p = 0;
                        {
                          SDSEntityType elementId = *subDimEntity.begin();
                          element_p = elementId;
                          if (!element_p)
                            {
                              throw std::runtime_error("addToExistingPartsNew: bad elem found 2");
                            }
                        }

                        stk::mesh::Entity& element = *element_p;

                        const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);
                        unsigned npts = elem_nodes.size();
                        for (unsigned ipts = 0; ipts < npts; ipts++)
                          {
                            stk::mesh::Entity * node = elem_nodes[ipts].entity();
                            if (!selector(*node))
                              {
                                found = false;
                                break;
                              }
                          }
                      }
                    else
                      {
                        for (SubDimCell_SDSEntityType::const_iterator ids = subDimEntity.begin(); ids != subDimEntity.end(); ++ids)
                          {
                            SDSEntityType nodeId = *ids;
                            stk::mesh::Entity * node = nodeId;
                            if (!selector(*node))
                              {
                                found = false;
                                break;
                              }
                          }
                      }
                    if (found)
                      {
                        // add to part

                        unsigned nidsz = nodeIds_onSE.size();

                        for (unsigned i_nid = 0; i_nid < nidsz; i_nid++)
                          {
                            stk::mesh::Entity * c_node = nodeIds_onSE[i_nid];

                            if (!c_node)
                              {
                                // note, this is ok - a null node here can come from a ghost element
                                if (1)
                                {
                                  continue;
                                }
                                else
                                {
                                  std::cout << "addToExistingPartsNew: " <<  nodeIds_onSE[i_nid] << " i_nid= " << i_nid << " nidsz= " << nidsz
                                            << std::endl;
                                  throw std::runtime_error("addToExistingParts: bad node found 0.3");
                                }
                              }
  
                            // only try to add to part if I am the owner
                            if (c_node->owner_rank() == m_eMesh.get_parallel_rank())
                              m_eMesh.get_bulk_data()->change_entity_parts( *c_node, add_parts, remove_parts );

                            if (0)
                              {
                                std::cout << "P[" << m_eMesh.get_rank() << "] adding node " << c_node->identifier() << " to   Part[" << ipart << "]= " << part.name()
                                          << " topology = " << (topology ? shards::CellTopology(topology).getName() : "null")
                                          << std::endl;
                              }
                          }

                      }
                  }
              }
          }
        //std::cout << "tmp addToExistingPartsNew...done " << std::endl;
      }

      SubDimCellData& getNewNodeAndOwningElement(SubDimCell_SDSEntityType& subDimEntity)
      {
        return m_cell_2_data_map[subDimEntity];
      }


#define NODE_REGISTRY_MAP_ACCESSORS_INLINED 1

      SubDimCellData * getFromMapPtr(const SubDimCell_SDSEntityType& subDimEntity) const
#if NODE_REGISTRY_MAP_ACCESSORS_INLINED
        {
          const SubDimCellToDataMap::const_iterator i = m_cell_2_data_map.find( subDimEntity );
          return i != m_cell_2_data_map.end() ? (const_cast<SubDimCellData *>(&(i->second))) : 0 ;
        }
#else
      ;
#endif

      SubDimCellData& getFromMap(const SubDimCell_SDSEntityType& subDimEntity) const
#if NODE_REGISTRY_MAP_ACCESSORS_INLINED
        {
          const SubDimCellToDataMap::const_iterator i = m_cell_2_data_map.find( subDimEntity );
          return * (const_cast<SubDimCellData *>(&(i->second)));
        }
#else
      ;
#endif
      void putInMap(SubDimCell_SDSEntityType& subDimEntity, SubDimCellData& data)
#if NODE_REGISTRY_MAP_ACCESSORS_INLINED
      {
        m_cell_2_data_map[subDimEntity] = data;
      }
#else
      ;
#endif

      /// @param needNodes should be true in general; it's used by registerNeedNewNode to generate actual data or not on the subDimEntity
      ///   For local refinement, subDimEntity's needs are not always known uniquely by the pair {elementId, iSubDimOrd}; for example, in
      ///   an element-based marking scheme, the shared face between two elements may be viewed differently.  So, we need the ability to
      ///   override the default behavior of always creating new nodes on the subDimEntity, but still allow the entity to be created in 
      ///   the NodeRegistry databse.

      typedef bool (NodeRegistry::*ElementFunctionPrototype)( const stk::mesh::Entity& element, NeededEntityType& needed_entity_rank, unsigned iSubDimOrd, bool needNodes);

      /// this is a helper method that loops over all sub-dimensional entities whose rank matches on of those in @param needed_entity_ranks
      ///    and registers that sub-dimensional entity as needing a new node.
      /// @param isGhost should be true if this element is a ghost, in which case this will call the appropriate method to set up for
      //     communications

      void //NodeRegistry::
      doForAllSubEntities(ElementFunctionPrototype function, const stk::mesh::Entity& element, vector<NeededEntityType>& needed_entity_ranks)
      {
        const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);

        shards::CellTopology cell_topo(cell_topo_data);
        const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);

        for (unsigned ineed_ent=0; ineed_ent < needed_entity_ranks.size(); ineed_ent++)
          {
            unsigned numSubDimNeededEntities = 0;
            stk::mesh::EntityRank needed_entity_rank = needed_entity_ranks[ineed_ent].first;

            if (needed_entity_rank == m_eMesh.edge_rank())
              {
                numSubDimNeededEntities = cell_topo_data->edge_count;
              }
            else if (needed_entity_rank == m_eMesh.face_rank())
              {
                numSubDimNeededEntities = cell_topo_data->side_count;
              }
            else if (needed_entity_rank == m_eMesh.element_rank())
              {
                numSubDimNeededEntities = 1;
              }

            for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
              {
                /// note: at this level of granularity we can do single edge refinement, hanging nodes, etc.
                //SubDimCell_SDSEntityType subDimEntity;
                //getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
                (this ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, true);

              } // iSubDimOrd
          } // ineed_ent
      }

      void //NodeRegistry::
      noInline_getSubDimEntity(SubDimCell_SDSEntityType& subDimEntity, const stk::mesh::Entity& element, stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd);

      /// fill
      ///    @param subDimEntity with the stk::mesh::EntityId's of
      ///    the ordinal @param iSubDimOrd sub-dimensional entity of
      ///    @param element of rank
      ///    @param needed_entity_rank
      ///
      void //NodeRegistry::
      getSubDimEntity(SubDimCell_SDSEntityType& subDimEntity, const stk::mesh::Entity& element, stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd)
      {
        subDimEntity.clear();
        // in the case of elements, we don't share any nodes so we just make a map of element id to node
        if (needed_entity_rank == m_eMesh.element_rank())
          {
            subDimEntity.insert( const_cast<stk::mesh::Entity*>(&element) );
            //!!subDimEntity.insert(element.identifier());
            return;
          }

        const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);

        //CellTopology cell_topo(cell_topo_data);
        const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);

        const unsigned *  inodes = 0;
        unsigned nSubDimNodes = 0;
        static const unsigned edge_nodes_2[2] = {0,1};
        static const unsigned face_nodes_3[3] = {0,1,2};
        static const unsigned face_nodes_4[4] = {0,1,2,3};

        // special case for faces in 3D
        if (needed_entity_rank == m_eMesh.face_rank() && needed_entity_rank == element.entity_rank())
          {
            nSubDimNodes = cell_topo_data->vertex_count;

            // note, some cells have sides with both 3 and 4 nodes (pyramid, prism)
            if (nSubDimNodes ==3 )
              inodes = face_nodes_3;
            else
              inodes = face_nodes_4;

          }
        // special case for edges in 2D
        else if (needed_entity_rank == m_eMesh.edge_rank() && needed_entity_rank == element.entity_rank())
          {
            nSubDimNodes = cell_topo_data->vertex_count;

            if (nSubDimNodes == 2 )
              {
                inodes = edge_nodes_2;
              }
            else
              {
                throw std::runtime_error("NodeRegistry bad for edges");
              }
          }
        else if (needed_entity_rank == m_eMesh.edge_rank())
          {
            inodes = cell_topo_data->edge[iSubDimOrd].node;
            nSubDimNodes = 2;
          }
        else if (needed_entity_rank == m_eMesh.face_rank())
          {
            nSubDimNodes = cell_topo_data->side[iSubDimOrd].topology->vertex_count;
            // note, some cells have sides with both 3 and 4 nodes (pyramid, prism)
            inodes = cell_topo_data->side[iSubDimOrd].node;
          }

        //subDimEntity.reserve(nSubDimNodes);
        for (unsigned jnode = 0; jnode < nSubDimNodes; jnode++)
          {
            subDimEntity.insert( elem_nodes[inodes[jnode]].entity() );
            //!!subDimEntity.insert(elem_nodes[inodes[jnode]].entity()->identifier());
          }

      }


      // FIXME
      unsigned total_size() {
        //throw std::runtime_error("not ready");
        //return m_cell_2_data_map.size();
        unsigned sz=0;

        for (SubDimCellToDataMap::iterator cell_iter = m_cell_2_data_map.begin(); cell_iter != m_cell_2_data_map.end(); ++cell_iter)
          {
            SubDimCellData& data = (*cell_iter).second;
            NodeIdsOnSubDimEntityType& nodeIds_onSE = data.get<SDC_DATA_GLOBAL_NODE_IDS>();

            sz += nodeIds_onSE.size();
          }
        return sz;
      }

      unsigned local_size()
      {
        unsigned sz=0;
        for (SubDimCellToDataMap::iterator cell_iter = m_cell_2_data_map.begin(); cell_iter != m_cell_2_data_map.end(); ++cell_iter)
          {
            SubDimCellData& data = (*cell_iter).second;

            stk::mesh::EntityId owning_elementId = stk::mesh::entity_id(data.get<SDC_DATA_OWNING_ELEMENT_KEY>());

            NodeIdsOnSubDimEntityType& nodeIds_onSE = data.get<SDC_DATA_GLOBAL_NODE_IDS>();
            if (nodeIds_onSE.size())
              {
                //!
                unsigned erank = m_eMesh.element_rank();
                erank = stk::mesh::entity_rank(data.get<SDC_DATA_OWNING_ELEMENT_KEY>());
                stk::mesh::Entity * owning_element = get_entity_element(*m_eMesh.get_bulk_data(), erank, owning_elementId);
                //!

                if (!owning_element)
                  {
                    std::cout << "tmp owning_element = null, owning_elementId= " << owning_elementId
                      //<< " nodeIds_onSE= " << nodeIds_onSE
                              << std::endl;
                    throw std::logic_error("logic: hmmm #5.2");
                  }
                if (!m_eMesh.isGhostElement(*owning_element))
                  {
                    //sz += 1;
                    sz += nodeIds_onSE.size();
                  }
              }
          }
        return sz;
      }

      //========================================================================================================================
      // low-level interface
      // FIXME
      bool isParallelRun(unsigned size) { return true; }

      void checkDB(std::string msg="")
      {
        if (0)
          {
            unsigned sz=0;
            for (SubDimCellToDataMap::iterator cell_iter = m_cell_2_data_map.begin(); cell_iter != m_cell_2_data_map.end(); ++cell_iter)
              {
                SubDimCellData& data = (*cell_iter).second;
                stk::mesh::EntityId owning_elementId = stk::mesh::entity_id(data.get<SDC_DATA_OWNING_ELEMENT_KEY>());
                NodeIdsOnSubDimEntityType& nodeIds_onSE = data.get<SDC_DATA_GLOBAL_NODE_IDS>();

                stk::mesh::Entity * owning_element = m_eMesh.get_bulk_data()->get_entity(m_eMesh.element_rank(), owning_elementId);
                if (!owning_element)
                  throw std::logic_error("logic: hmmm #5.3");
                bool isGhost = m_eMesh.isGhostElement(*owning_element);
                if (!m_eMesh.isGhostElement(*owning_element))
                  {
                    ++sz;
                  }
                if (!isGhost)
                  std::cout << "P[" << m_eMesh.get_rank() << "] owning_elementId = "  << owning_elementId << " isGhostElement = " << isGhost
                            << " nodeId = " << nodeIds_onSE << std::endl;
              }
          }
        if (1)
          {
            std::cout << "NodeRegistry::checkDB start msg= " << msg << std::endl;
            for (SubDimCellToDataMap::iterator cell_iter = m_cell_2_data_map.begin(); cell_iter != m_cell_2_data_map.end(); ++cell_iter)
              {
                //const SubDimCell_SDSEntityType& subDimEntity = (*iter).first;
                SubDimCellData&            subDimCellData      = (*cell_iter).second;
                //stk::mesh::EntityId        owning_elementId    = stk::mesh::entity_id(subDimCellData.get<SDC_DATA_OWNING_ELEMENT_KEY>());
                //stk::mesh::EntityRank      owning_element_rank = stk::mesh::entity_rank(subDimCellData.get<SDC_DATA_OWNING_ELEMENT_KEY>());
                NodeIdsOnSubDimEntityType& nodeIds_onSE        = subDimCellData.get<SDC_DATA_GLOBAL_NODE_IDS>();

                for (unsigned i=0; i < nodeIds_onSE.size(); i++)
                  {
                    stk::mesh::Entity *node = nodeIds_onSE[i];
                    stk::mesh::EntityId nodeId = nodeIds_onSE.m_entity_id_vector[i];
                    VERIFY_OP_ON(node, !=, 0, "checkDB #11.1");
                    VERIFY_OP_ON(nodeId, !=, 0, "checkDB #11.1.1");
                    VERIFY_OP_ON(node->identifier(), ==, nodeId, "checkDB #11.2");
                    stk::mesh::Entity *node_0 = m_eMesh.get_bulk_data()->get_entity(0, nodeId);
                    
                    VERIFY_OP_ON(node, ==, node_0, "checkDB #11.3");
                    VERIFY_OP_ON(node_0->identifier(), ==, nodeId, "checkDB #11.4");
                  }
              }
            std::cout << "NodeRegistry::checkDB end msg= " << msg << std::endl;
          }

        if (0)
          {
            std::cout << "NodeRegistry::checkDB start 1 msg= " << msg << std::endl;
            for (SubDimCellToDataMap::iterator cell_iter = m_cell_2_data_map.begin(); cell_iter != m_cell_2_data_map.end(); ++cell_iter)
              {
                //const SubDimCell_SDSEntityType& subDimEntity = (*iter).first;
                SubDimCellData&            subDimCellData      = (*cell_iter).second;
                stk::mesh::EntityId        owning_elementId    = stk::mesh::entity_id(subDimCellData.get<SDC_DATA_OWNING_ELEMENT_KEY>());
                stk::mesh::EntityRank      owning_element_rank = stk::mesh::entity_rank(subDimCellData.get<SDC_DATA_OWNING_ELEMENT_KEY>());
                NodeIdsOnSubDimEntityType& nodeIds_onSE        = subDimCellData.get<SDC_DATA_GLOBAL_NODE_IDS>();

                if (owning_elementId == 0) continue;
                stk::mesh::Entity * owning_element = m_eMesh.get_bulk_data()->get_entity(owning_element_rank, owning_elementId);
                if (!owning_element)
                  throw std::logic_error("logic: checkDB hmmm #11.0");
                //bool isGhost = m_eMesh.isGhostElement(*owning_element);
                if (!m_eMesh.isGhostElement(*owning_element))
                  {
                    for (unsigned i=0; i < nodeIds_onSE.size(); i++)
                      {
                        stk::mesh::Entity *node = nodeIds_onSE[i];
                        stk::mesh::EntityId nodeId = nodeIds_onSE.m_entity_id_vector[i];
                        VERIFY_OP_ON(node, !=, 0, "checkDB #11.1");
                        VERIFY_OP_ON(nodeId, !=, 0, "checkDB #11.1.1");
                        VERIFY_OP_ON(node->identifier(), ==, nodeId, "checkDB #11.2");
                        stk::mesh::Entity *node_0 = m_eMesh.get_bulk_data()->get_entity(0, nodeId);
                    
                        VERIFY_OP_ON(node, ==, node_0, "checkDB #11.3");
                        VERIFY_OP_ON(node_0->identifier(), ==, nodeId, "checkDB #11.4");
                        unsigned owner_rank = node->owner_rank();
                        VERIFY_OP_ON(owner_rank, ==, owning_element->owner_rank(), "checkDB #11.6");
                      }
                  }
              }
            std::cout << "NodeRegistry::checkDB end 1 msg= " << msg << std::endl;

          }

      }

      /// allocate the send/recv buffers for all-to-all communication
      bool allocateBuffers()
      {
        stk::CommAll& comm_all = *m_comm_all;
        unsigned proc_size = comm_all.parallel_size();
        unsigned proc_rank = comm_all.parallel_rank();

        // FIXME - add some error checking

#if 0
        if (!isParallelRun(proc_size))
          {
            return false;
          }
#endif

        bool local = true; // FIXME
        unsigned num_msg_bounds = proc_size < 4 ? proc_size : proc_size/4 ;
        bool global = comm_all.allocate_buffers(num_msg_bounds , false, local );
        if ( not global )
          {
            std::cout << "P[" << proc_rank << "] : not global" << std::endl;
            return false;
          }
        return true;
      }

      void communicate()
      {
        stk::CommAll& comm_all = *m_comm_all;
        //unsigned proc_size = comm_all.parallel_size();
        //unsigned proc_rank = comm_all.parallel_rank();

#if 0
        for (unsigned i_proc_rank = 0; i_proc_rank < proc_size; i_proc_rank++)
          {
            std::cout << "P[" << proc_rank << "] : i_proc_rank = " << i_proc_rank << " send buf size =  " 
                      <<   m_comm_all->send_buffer( i_proc_rank ).size() << " num in buf= " 
                      <<   m_comm_all->send_buffer( i_proc_rank ).size() / sizeof(CommDataType) <<  std::endl;
          }
#endif
        comm_all.communicate();

        stk::ParallelMachine pm = m_eMesh.get_bulk_data()->parallel();
        int failed = 0;
        stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );

        unpack();

      }

      void
      unpack()
      {
        stk::CommAll& comm_all = *m_comm_all;

        int failed = 0;
        std::string msg;

        stk::ParallelMachine pm = m_eMesh.get_bulk_data()->parallel();
        unsigned proc_size = m_eMesh.get_bulk_data()->parallel_size();
        unsigned proc_rank = comm_all.parallel_rank();

        vector<stk::mesh::EntityProc> nodes_to_ghost;

        CommDataType buffer_entry;
        NodeIdsOnSubDimEntityType nodeIds_onSE;
        try
          {
            for(unsigned from_proc = 0; from_proc < proc_size; ++from_proc )
              {
                stk::CommBuffer & recv_buffer = comm_all.recv_buffer( from_proc );

                //unsigned num_in_buffer = recv_buffer.size() / sizeof(CommDataType);
                //std::cout << "for proc= " << from_proc << " recv_buffer.size()= " << recv_buffer.size() << " num_in_buffer = " << num_in_buffer << std::endl;

                while ( recv_buffer.remaining() )
                  {
                    //
                    // this->unpack( this->container(), p, recv_buffer );
                    //
                    // Rank of sub-dim cells needing new nodes, which sub-dim entity, one non-owning element identifier, nodeId_elementOwnderId.first
                    // typedef boost::tuple::tuple<stk::mesh::EntityRank, unsigned, stk::mesh::EntityId, stk::mesh::EntityId> CommDataType;


                    recv_buffer.unpack< CommDataType >( buffer_entry );
                    nodeIds_onSE.unpack(m_eMesh, recv_buffer);
                    //recv_buffer.unpack< NodeIdsOnSubDimEntityType > (nodeIds_onSE);


                    //std::cout << "P[" << proc_rank << "] unpack for buffer from proc= " << from_proc << " " << buffer_entry << std::endl;
                    createNodeAndConnect(buffer_entry, nodeIds_onSE, from_proc, nodes_to_ghost);
                  }
              }
          }
        catch ( std::exception &x )
          {
            failed = 1;
            msg = std::string("unpack error: ")+x.what();
          }

        stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );
        if ( failed )
          {
            throw std::runtime_error( msg+" from unpack error, rank = "+toString(proc_rank) );
          }

        if (nodes_to_ghost.size())
          {
            stk::mesh::Ghosting & ghosting = m_eMesh.get_bulk_data()->create_ghosting( std::string("new_nodes") );

            vector<stk::mesh::Entity*> receive;
            ghosting.receive_list( receive );
            //if (receive.size()) std::cout << "NodeRegistry::endGetFromRemote receive.size() = " << receive.size() << std::endl;
            m_eMesh.get_bulk_data()->change_ghosting( ghosting, nodes_to_ghost, receive);
          }

      }// unpack


      /// after registering all needed nodes, this method is used to request new nodes on this processor
      void createNewNodesInParallel()
      {
        static int num_times_called = 0;
        ++num_times_called;

        stk::mesh::Part* new_nodes_part = m_eMesh.get_non_const_part("refine_new_nodes_part");
        if (new_nodes_part)
          {
            // first remove any existing nodes
            //unsigned proc_rank = m_eMesh.get_rank();
            //std::cout << "P[" << proc_rank << "] remove existing nodes... " <<  std::endl;
            std::vector<stk::mesh::Part*> remove_parts(1, new_nodes_part);
            std::vector<stk::mesh::Part*> add_parts;
            std::vector<stk::mesh::Entity *> node_vec;

            stk::mesh::Selector removePartSelector(*new_nodes_part & m_eMesh.get_fem_meta_data()->locally_owned_part() );
            const std::vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.node_rank() );
            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                stk::mesh::Bucket & bucket = **k ;
                if (removePartSelector(bucket))
                  {
                    const unsigned num_entity_in_bucket = bucket.size();
                    for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                      {
                        stk::mesh::Entity& node = bucket[ientity];
                        node_vec.push_back(&node);
                      }
                  }
              }
            for (unsigned ii=0; ii < node_vec.size(); ii++)
              {
                m_eMesh.get_bulk_data()->change_entity_parts( *node_vec[ii], add_parts, remove_parts );
              }

            //std::cout << "P[" << proc_rank << "] remove existing nodes...done " <<  std::endl;
          }

        unsigned num_nodes_needed = local_size();

        // FIXME
        // assert( bulk data is in modifiable mode)
        // create new entities on this proc
        vector<stk::mesh::Entity *> new_nodes;

        if (m_useCustomGhosting)
          m_eMesh.createEntities( stk::mesh::fem::FEMMetaData::NODE_RANK, num_nodes_needed, new_nodes);

        std::vector<stk::mesh::EntityId> ids(num_nodes_needed);

        if (!m_useCustomGhosting)
          {
#define NR_GEN_OWN_IDS 0
#if NR_GEN_OWN_IDS
            new_nodes.resize(num_nodes_needed);
#else
            m_eMesh.createEntities( stk::mesh::fem::FEMMetaData::NODE_RANK, num_nodes_needed, new_nodes);
#endif

            // bogus, but just for testing - we delete the just-created entities and re-declare them without ownership info
            //   so that modification_end takes care of assigning ownership
            stk::mesh::PartVector empty_parts;
            for (unsigned i=0; i < num_nodes_needed; i++)
              {
#if NR_GEN_OWN_IDS
                ids[i] = (num_times_called*100000) + i + (num_times_called*100000)*1000*m_eMesh.get_parallel_rank();
#else
                ids[i] = new_nodes[i]->identifier();
                bool did_destroy = m_eMesh.get_bulk_data()->destroy_entity(new_nodes[i]);
                VERIFY_OP_ON(did_destroy, ==, true, "createNewNodesInParallel couldn't destroy");
#endif
                new_nodes[i] = & m_eMesh.get_bulk_data()->declare_entity(m_eMesh.node_rank(), ids[i], empty_parts);
#if PERCEPT_USE_PSEUDO_ELEMENTS
                unsigned proc_rank = m_eMesh.get_rank();
                stk::mesh::Entity *elem = & m_eMesh.get_bulk_data()->declare_entity(m_eMesh.element_rank()+PSEUDO_ELEMENT_RANK_SHIFT, ids[i]*PSEUDO_ELEMENT_MAGIC_NUMBER+proc_rank, empty_parts);
                m_pseudo_entities.insert(elem);
                m_eMesh.get_bulk_data()->declare_relation(*elem, *new_nodes[i], 0);
#endif
              }
          }
        
        if (new_nodes_part)
          {
            //unsigned proc_rank = m_eMesh.get_rank();
            //std::cout << "P[" << proc_rank << "] add new nodes... " <<  std::endl;
            stk::mesh::Selector selector(m_eMesh.get_fem_meta_data()->locally_owned_part() );
            std::vector<stk::mesh::Part*> add_parts(1, new_nodes_part);
            std::vector<stk::mesh::Part*> remove_parts;
            for (unsigned ind = 0; ind < new_nodes.size(); ind++)
              {
                  
                m_eMesh.get_bulk_data()->change_entity_parts( *new_nodes[ind], add_parts, remove_parts );
              }
            //std::cout << "P[" << proc_rank << "] add new nodes...done " <<  std::endl;
          }

        // set map values to new node id's
        unsigned inode=0;

        for (SubDimCellToDataMap::iterator cell_iter = m_cell_2_data_map.begin(); cell_iter != m_cell_2_data_map.end(); ++cell_iter)
          {
            SubDimCellData& data = (*cell_iter).second;
            NodeIdsOnSubDimEntityType& nodeIds_onSE = data.get<SDC_DATA_GLOBAL_NODE_IDS>();
            if (!nodeIds_onSE.size())
              continue;

            stk::mesh::EntityId owning_elementId = stk::mesh::entity_id(data.get<SDC_DATA_OWNING_ELEMENT_KEY>());

            if (!owning_elementId)
              {
                throw std::logic_error("logic: hmmm #5.4.0");
              }

            //!
            unsigned erank = m_eMesh.element_rank();
            erank = stk::mesh::entity_rank(data.get<SDC_DATA_OWNING_ELEMENT_KEY>());
            stk::mesh::Entity * owning_element = get_entity_element(*m_eMesh.get_bulk_data(), erank, owning_elementId);
            //!

            if (!owning_element)
              {
                throw std::logic_error("logic: hmmm #5.4");
              }
            if (!m_eMesh.isGhostElement(*owning_element))
              {
                if (nodeIds_onSE.m_entity_id_vector.size() != nodeIds_onSE.size())
                  {
                    throw std::logic_error("NodeRegistry:: createNewNodesInParallel logic err #0.0");
                  }

                for (unsigned ii = 0; ii < nodeIds_onSE.size(); ii++)
                  {
                    //nodeIds_onSE[ii] = new_nodes[inode]->identifier();

                    VERIFY_OP(inode, < , num_nodes_needed, "UniformRefiner::doBreak() too many nodes");
                    if ( DEBUG_NR_UNREF)
                      {
                        std::cout << "tmp createNewNodesInParallel: old node id= " << (nodeIds_onSE[ii] ? toString(nodeIds_onSE[ii]->identifier()) : std::string("null")) << std::endl;
                        std::cout << "tmp createNewNodesInParallel: new node=";
                        m_eMesh.print_entity(std::cout, *new_nodes[inode]);
                      }

                    // if already exists from a previous iteration/call to doBreak, don't reset it and just use the old node
                    if (nodeIds_onSE[ii])
                      {
                        if (DEBUG_NR_UNREF)
                          {
                            std::cout << "tmp createNewNodesInParallel: old node id is no-null, re-using it= " << (nodeIds_onSE[ii] ? toString(nodeIds_onSE[ii]->identifier()) : std::string("null")) << std::endl;
                            std::cout << "tmp createNewNodesInParallel: new node=";
                            m_eMesh.print_entity(std::cout, *new_nodes[inode]);
                          }
                      }
                    else
                      {
                        nodeIds_onSE[ii] = new_nodes[inode];
                        nodeIds_onSE.m_entity_id_vector[ii] = new_nodes[inode]->identifier();
                      }

                    //nodeIds_onSE.m_entity_vector[ii] = new_nodes[inode];
                    inode++;
                  }
                //data.get<SDC_DATA_GLOBAL_NODE_IDS>()[0] = new_nodes[inode]->identifier();
              }
          }
      }


      /// unpacks the incoming information in @param buffer_entry and adds that information to my local node registry
      /// (i.e. the map of sub-dimensional entity to global node id is updated)
      void
      createNodeAndConnect(CommDataType& buffer_entry, NodeIdsOnSubDimEntityType& nodeIds_onSE, unsigned from_proc, vector<stk::mesh::EntityProc>& nodes_to_ghost)
      {
        //stk::mesh::EntityId&                  non_owning_elementId = buffer_entry.get<CDT_NON_OWNING_ELEMENT_KEY>();

        stk::mesh::EntityRank& needed_entity_rank                    = buffer_entry.get<CDT_NEEDED_ENTITY_RANK>();
        unsigned    iSubDimOrd                            = buffer_entry.get<CDT_SUB_DIM_ENTITY_ORDINAL>();
        stk::mesh::EntityKey&  non_owning_elementKey                 = buffer_entry.get<CDT_NON_OWNING_ELEMENT_KEY>();
        stk::mesh::EntityId    non_owning_elementId                  = stk::mesh::entity_id(non_owning_elementKey);
        stk::mesh::EntityRank  non_owning_elementRank                = stk::mesh::entity_rank(non_owning_elementKey);

        // create a new relation here?  no, we are going to delete this element, so we just register that the new node is attached to
        //stk::mesh::Entity * element = m_eMesh.get_bulk_data()->get_entity(m_eMesh.element_rank(), non_owning_elementId);

        //!
        unsigned erank = m_eMesh.element_rank();
        erank = non_owning_elementRank;
        stk::mesh::Entity * element = get_entity_element(*m_eMesh.get_bulk_data(), erank, non_owning_elementId);
        //!

        for (unsigned iid = 0; iid < nodeIds_onSE.size(); iid++)
          {
            //stk::mesh::Entity * node = get_entity_node_I(*m_eMesh.get_bulk_data(),stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE[iid]);
            //nodeIds_onSE.m_entity_vector[iid] = node;
            //stk::mesh::Entity * node = get_entity_node_Ia(*m_eMesh.get_bulk_data(),stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE, iid);
            stk::mesh::Entity * node = nodeIds_onSE[iid];

            // has to be null, right?
            if (node)
              {
                throw std::logic_error("logic: node should be null in createNodeAndConnect");
              }
          }
        if (!element)
          {
            throw std::logic_error("logic: element shouldn't be null in createNodeAndConnect");
          }

        static SubDimCell_SDSEntityType subDimEntity;
        getSubDimEntity(subDimEntity, *element, needed_entity_rank, iSubDimOrd);
        SubDimCellData& subDimCellData = getNewNodeAndOwningElement(subDimEntity);
        // assert it is empty?

        subDimCellData.get<SDC_DATA_GLOBAL_NODE_IDS>() = nodeIds_onSE;

#ifndef NDEBUG
        stk::mesh::EntityId owning_element_id = stk::mesh::entity_id(subDimCellData.get<SDC_DATA_OWNING_ELEMENT_KEY>());
        stk::mesh::EntityRank owning_element_rank = stk::mesh::entity_rank(subDimCellData.get<SDC_DATA_OWNING_ELEMENT_KEY>());
        //VERIFY_OP(owning_element_id, !=, non_owning_elementId, "createNodeAndConnect:: bad elem ids");
        //VERIFY_OP(owning_element_id, < , non_owning_elementId, "createNodeAndConnect:: bad elem ids 2");

        // owning_elementId < non_owning_elementId && owning_elementRank >= non_owning_elementRank
        if (owning_element_rank >= non_owning_elementRank)
          {
            VERIFY_OP(owning_element_id, !=, non_owning_elementId, "createNodeAndConnect:: bad elem ids");
            VERIFY_OP(owning_element_id, < , non_owning_elementId, "createNodeAndConnect:: bad elem ids 2");
          }
#endif

      }

    public:
      SubDimCellToDataMap& getMap() { return  m_cell_2_data_map; }
      PerceptMesh& getMesh() { return m_eMesh; }
      bool getUseCustomGhosting() { return m_useCustomGhosting; }

      // remove any sub-dim entities from the map that have a node in deleted_nodes
      void cleanDeletedNodes(std::set<stk::mesh::Entity *>& deleted_nodes, bool debug=false)
      {
        std::set<stk::mesh::Entity *> deleted_nodes_copy = deleted_nodes;

        if (DEBUG_NR_UNREF)
          std::cout << "tmp cleanDeletedNodes deleted_nodes size: " << deleted_nodes_copy.size() << std::endl;

        SubDimCellToDataMap::iterator iter;
        std::vector<SubDimCellToDataMap::iterator> to_delete;

        SubDimCellToDataMap& map = m_cell_2_data_map;
        if (DEBUG_NR_UNREF)
          std::cout << "tmp cleanDeletedNodes map size: " << map.size() << std::endl;

        for (iter = map.begin(); iter != map.end(); ++iter)
          {
            //const SubDimCell_SDSEntityType& subDimEntity = (*iter).first;
            SubDimCellData& nodeId_elementOwnderId = (*iter).second;
            
            NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();

            unsigned jj = 0; 
            bool found = false;
            for (unsigned ii = 0; ii < nodeIds_onSE.size(); ii++)
              {
                if (deleted_nodes.find(nodeIds_onSE[ii]) != deleted_nodes.end())
                  {
                    found = true;
                    jj = ii;
                    deleted_nodes_copy.erase(nodeIds_onSE[ii]);
                    break;
                  }
              }
            if (found)
              {
                if (DEBUG_NR_UNREF)
                  {
                    std::cout << "tmp cleanDeletedNodes:: removing node id= " << nodeIds_onSE[jj]->identifier() << std::endl;
                    std::cout << "Node: ";
                    m_eMesh.print_entity(std::cout, *nodeIds_onSE[jj]);
                  }
                if (!debug)
                  {
                    //unsigned owning_elementId = stk::mesh::entity_id(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>());
                    //unsigned owning_elementRank = stk::mesh::entity_rank(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>());
                    //nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>() = stk::mesh::EntityKey(0u, 0u);
                    //nodeIds_onSE.resize(0);
                    to_delete.push_back(iter);
                  }
              }
          }

        //std::cout << "tmp cleanDeletedNodes to_delete.size()= " << to_delete.size() << " map.size()= " << map.size() << std::endl;
        for (unsigned itd=0; itd < to_delete.size(); itd++)
          {
            map.erase(to_delete[itd]);
          }

        if (DEBUG_NR_UNREF && deleted_nodes_copy.size())
          {
            std::cout << "tmp cleanDeletedNodes some deleted nodes not found, size()=: " << deleted_nodes_copy.size() << " nodes= " << std::endl;
            std::set<stk::mesh::Entity *>::iterator it;
            for (it = deleted_nodes_copy.begin(); it != deleted_nodes_copy.end(); ++it)
              {
                stk::mesh::Entity *node = *it;
                std::cout << "Node: ";
                m_eMesh.print_entity(std::cout, *node);
              }
            
          }
      }

      // further cleanup of the NodeRegistry db - some elements get deleted on some procs but the ghost elements
      //   are still in the db - the easiest way to detect this is as used here: during modification_begin(), 
      //   cleanup of cached transactions are performed, then we find these by seeing if our element id's are
      //   no longer in the stk_mesh db.

      void clear_element_owner_data_phase_2()
      {
        m_eMesh.get_bulk_data()->modification_begin();

        SubDimCellToDataMap::iterator iter;

        SubDimCellToDataMap& map = m_cell_2_data_map;

        for (iter = map.begin(); iter != map.end(); ++iter)
          {
            //const SubDimCell_SDSEntityType& subDimEntity = (*iter).first;
            SubDimCellData& nodeId_elementOwnderId = (*iter).second;
            
            unsigned owning_elementId = stk::mesh::entity_id(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>());
            unsigned owning_elementRank = stk::mesh::entity_rank(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>());

            if (owning_elementId)
              {
                //stk::mesh::Entity * owning_element = get_entity_element(*m_eMesh.get_bulk_data(), owning_elementRank, owning_elementId);
                stk::mesh::Entity * owning_element = m_eMesh.get_bulk_data()->get_entity(owning_elementRank, owning_elementId);

#if 0
                if (owning_element != owning_element_1) 
                  {
                    
                    std::cout << "P[" << m_eMesh.get_rank() << "] NR::clear_2 error # 1= " << owning_element << " " << owning_element_1 <<  std::endl;
                    throw std::logic_error("NodeRegistry:: clear_element_owner_data_phase_2 error # 1");
                  }
#endif

                if (!owning_element)
                  {
                    NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
                    nodeIds_onSE.resize(0);

                    nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>() = stk::mesh::EntityKey(owning_elementRank, 0u);

                  }
              }
          }
        m_eMesh.get_bulk_data()->modification_end();

      }

      // remove/zero any data that points to a deleted element
      // called by Refiner with "children_to_be_removed_with_ghosts"
      void clear_element_owner_data( std::set<stk::mesh::Entity *>& elems_to_be_deleted)
      {
        SubDimCellToDataMap::iterator iter;

        SubDimCellToDataMap& map = m_cell_2_data_map;

        for (iter = map.begin(); iter != map.end(); ++iter)
          {
            //const SubDimCell_SDSEntityType& subDimEntity = (*iter).first;
            SubDimCellData& nodeId_elementOwnderId = (*iter).second;
            
            unsigned owning_elementId = stk::mesh::entity_id(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>());
            unsigned owning_elementRank = stk::mesh::entity_rank(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>());

            if (owning_elementId)
              {
                stk::mesh::Entity * owning_element = get_entity_element(*m_eMesh.get_bulk_data(), owning_elementRank, owning_elementId);
                stk::mesh::Entity * owning_element_1 = m_eMesh.get_bulk_data()->get_entity(owning_elementRank, owning_elementId);

                if (owning_element != owning_element_1) 
                  {
                    std::cout << "P[" << m_eMesh.get_rank() << "] NR::clear_1 error # 1= " << owning_element << " " << owning_element_1 <<  std::endl;
                    throw std::logic_error("NodeRegistry:: clear_element_owner_data_phase_1 error # 1");
                  }

                if (owning_element)
                  {

                    bool isGhost = m_eMesh.isGhostElement(*owning_element);
                        
                    bool in_deleted_list = elems_to_be_deleted.find(owning_element) != elems_to_be_deleted.end();

                    if (in_deleted_list)
                      {

                        if (0)
                          std::cout << "clear_element_owner_data: owning_elementId = " << owning_elementId 
                                    << " isGhost= " << isGhost << std::endl;

                        // FIXME
                        NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
                        nodeIds_onSE.resize(0);
                        // FIXME
                        nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>() = stk::mesh::EntityKey(owning_elementRank, 0u);
                      }
                  }
              }
          }
      }

      void dumpDB(std::string msg="")
      {
        if (!DEBUG_NR_UNREF) return;
        SubDimCellToDataMap::iterator iter;
        SubDimCellToDataMap& map = m_cell_2_data_map;
        std::cout << msg << " tmp dumpDB map size: " << map.size() << std::endl;

        for (iter = map.begin(); iter != map.end(); ++iter)
          {
            const SubDimCell_SDSEntityType& subDimEntity = (*iter).first;
            SubDimCellData& nodeId_elementOwnderId = (*iter).second;
            
            NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();

            for (unsigned ii = 0; ii < nodeIds_onSE.size(); ii++)
              {
                std::cout << "tmp ddb:: node id= " << nodeIds_onSE[ii]->identifier() << std::endl;
                std::cout << "subDimEntity= ";
                for (unsigned k=0; k < subDimEntity.size(); k++)
                  {
                    std::cout << " " << subDimEntity[k]->identifier() << " ";
                  }
                std::cout << "Node: ";
                m_eMesh.print_entity(std::cout, *nodeIds_onSE[ii]);
              }
          }
      }

      // estimate of memory used by this object
      unsigned get_memory_usage()
      {
        SubDimCellToDataMap::iterator iter;
        SubDimCellToDataMap& map = m_cell_2_data_map;

        unsigned mem=0;

        for (iter = map.begin(); iter != map.end(); ++iter)
          {
            const SubDimCell_SDSEntityType& subDimEntity = (*iter).first;
            SubDimCellData& nodeId_elementOwnderId = (*iter).second;
            
            mem += sizeof(SDSEntityType)*subDimEntity.size();
            mem += sizeof(stk::mesh::EntityKey);
            NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();

            unsigned mem1 = (sizeof(NodeIdsOnSubDimEntityTypeQuantum)+
                             sizeof(stk::mesh::EntityId))*nodeIds_onSE.size() +sizeof(unsigned);

            mem += mem1;
          }
        return mem;
      }


    private:
      percept::PerceptMesh& m_eMesh;
      stk::CommAll * m_comm_all;
      SubDimCellToDataMap m_cell_2_data_map;

      vector<stk::mesh::EntityProc> m_nodes_to_ghost;
      SetOfEntities m_pseudo_entities;
      bool m_useCustomGhosting;

    public:
      int m_gee_cnt;
      int m_gen_cnt;
      //const CellTopologyData * const m_cell_topo_data;
      std::vector<EntityRepo> m_entity_repo;

      bool m_debug;

      NodeRegistryState m_state;
    };

  }
}
#endif
