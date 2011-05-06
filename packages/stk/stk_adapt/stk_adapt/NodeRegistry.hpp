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

#include <stk_percept/NoMallocArray.hpp>
#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/Teuchos_My_Hashtable.hpp>
#include <stk_percept/Hashtable.hpp>

#include <boost/array.hpp>

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
#define STK_ADAPT_NODEREGISTRY_USE_ENTITY_REPO 1
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


namespace stk {
  namespace adapt {

    using namespace stk::percept;
    using std::vector;
    using std::map;
    using std::set;

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

      NodeIdsOnSubDimEntityType(unsigned sz=1, NodeIdsOnSubDimEntityTypeQuantum allValues=0) : base_type(sz,allValues),
                                                                                               m_entity_id_vector(sz,0u) {}
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

    public:
      //========================================================================================================================
      // high-level interface
      //NodeRegistry(percept::PerceptMesh& eMesh) : m_eMesh(eMesh), m_comm_all(eMesh.getBulkData()->parallel()), m_gee_cnt(0), m_gen_cnt(0),
      //m_entity_repo(stk::mesh::stk::percept::EntityRankEnd)

      NodeRegistry(percept::PerceptMesh& eMesh) : m_eMesh(eMesh), m_comm_all(eMesh.getBulkData()->parallel()),
                                                  // why does this cause failures? 
                                                  //m_cell_2_data_map(eMesh.getNumberElements()*8u),
                                                  m_gee_cnt(0), m_gen_cnt(0),
                                                  m_entity_repo(stk::percept::EntityRankEnd),
                                                  m_debug(false),
                                                  m_state(NRS_NONE)
      {
#if NODE_REGISTRY_MAP_TYPE_GOOGLE
        //SubDimCell_SDSEntityType empty_key;
        //empty_key.insert( std::numeric_limits<stk::mesh::EntityId>::max() );
        //m_cell_2_data_map.set_empty_key(empty_key);

        SubDimCell_SDSEntityType deleted_key;
        deleted_key.insert( std::numeric_limits<stk::mesh::EntityId>::max() - 1u ); 
        m_cell_2_data_map.set_deleted_key(deleted_key);
#endif
      }


      void initialize() 
      {
        m_cell_2_data_map.clear();
        for (unsigned i = 0; i < stk::percept::EntityRankEnd; i++) m_entity_repo[i].clear();
      }

      void //NodeRegistry::
      beginRegistration()
      {
        m_state = NRS_START_REGISTER_NODE;
        if (m_debug)
          std::cout << "P[" << m_eMesh.getRank() << "] tmp NodeRegistry::beginRegistration" << std::endl;
      }

      void //NodeRegistry::
      endRegistration()
      {
        if (m_debug)
          std::cout << "P[" << m_eMesh.getRank() << "] tmp NodeRegistry::endRegistration start" << std::endl;

        //putInESMap();

        m_eMesh.getBulkData()->modification_begin();  
        this->createNewNodesInParallel(); 
        m_nodes_to_ghost.resize(0);

#if STK_ADAPT_NODEREGISTRY_DO_REHASH
        m_cell_2_data_map.rehash(m_cell_2_data_map.size());
#endif
        if (m_debug)
          std::cout << "P[" << m_eMesh.getRank() << "] tmp NodeRegistry::endRegistration end" << std::endl;

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
          std::cout << "P[" << m_eMesh.getRank() << "] tmp NodeRegistry::beginCheckForRemote " << std::endl;
      }

      void //NodeRegistry::
      endCheckForRemote()
      {
        if (m_debug)
          std::cout << "P[" << m_eMesh.getRank() << "] tmp NodeRegistry::endCheckForRemote start " << std::endl;
        stk::ParallelMachine pm = m_eMesh.getBulkData()->parallel();
        int failed = 0;
        stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );

        this->allocateBuffers();  

#if STK_ADAPT_NODEREGISTRY_DO_REHASH
        m_cell_2_data_map.rehash(m_cell_2_data_map.size());
#endif

        if (m_debug)
          std::cout << "P[" << m_eMesh.getRank() << "] tmp NodeRegistry::endCheckForRemote end " << std::endl;

        m_state = NRS_END_CHECK_FOR_REMOTE;

      }

      void //NodeRegistry::
      beginGetFromRemote()
      {
        m_state = NRS_START_GET_FROM_REMOTE;
        if (m_debug)
          std::cout << "P[" << m_eMesh.getRank() << "] tmp NodeRegistry::beginGetFromRemote  " << std::endl;

      }
      void //NodeRegistry::
      endGetFromRemote()
      {
        if (m_debug)
          std::cout << "P[" << m_eMesh.getRank() << "] tmp NodeRegistry::endGetFromRemote start " << std::endl;
        stk::ParallelMachine pm = m_eMesh.getBulkData()->parallel();
        int failed = 0;
        stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );

        this->communicate();  

        failed = 0;
        stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );

        {
          stk::mesh::Ghosting & ghosting = m_eMesh.getBulkData()->create_ghosting( std::string("new_nodes") );
          
          vector<stk::mesh::Entity*> receive;

          ghosting.receive_list( receive );

          m_eMesh.getBulkData()->change_ghosting( ghosting, m_nodes_to_ghost, receive);

        }

        failed = 0;
        stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );

        m_eMesh.getBulkData()->modification_end();  

        setAllReceivedNodeData();

        if (m_debug)
          std::cout << "P[" << m_eMesh.getRank() << "] tmp NodeRegistry::endGetFromRemote end " << std::endl;

        m_state = NRS_END_GET_FROM_REMOTE;
      }

      void setAllReceivedNodeData()
      {
        EXCEPTWATCH;
        SubDimCellToDataMap::iterator iter;

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
                    stk::mesh::Entity *node = get_entity_node_I(*m_eMesh.getBulkData(), stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE.m_entity_id_vector[ii]);  // FIXME
                    if (!node)
                      {
                        throw std::logic_error("NodeRegistry:: setAllReceivedNodeData logic err #3");
                      }
                    nodeIds_onSE[ii] = node;
                  }
              }
          }
      }

      void putInESMap(const stk::mesh::Entity& element, stk::mesh::EntityRank& needed_entity_rank, unsigned iSubDimOrd, SubDimCellData& data) 
      {
        ElementSideTuple est(&element, needed_entity_rank, iSubDimOrd);
        m_element_side_map[est] = data;
      }

      /// Register the need for a new node on the sub-dimensional entity @param subDimEntity on element @param element.
      /// If the element is a ghost element, the entity is still registered: the locality/ownership of the new entity
      /// can be determined by the locality of the element (ghost or not).
      bool registerNeedNewNode(const stk::mesh::Entity& element, NeededEntityType& needed_entity_rank, unsigned iSubDimOrd) 
      {
        static SubDimCell_SDSEntityType subDimEntity;
        getSubDimEntity(subDimEntity, element, needed_entity_rank.first, iSubDimOrd);
        
        static SubDimCellData new_SubDimCellData;
        static SubDimCellData empty_SubDimCellData;

        SubDimCellData* nodeId_elementOwnderId_ptr = getFromMapPtr(subDimEntity);
        SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
        bool is_empty = nodeId_elementOwnderId_ptr == 0;

        // if empty or if my id is the smallest, make this element the owner
        bool should_put_in = 
          (element.identifier()  < stk::mesh::entity_id(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>()))
          || (element.entity_rank() > stk::mesh::entity_rank(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>()));

        /// once it's in, the assertion should be:
        ///   owning_elementId < non_owning_elementId && owning_elementRank >= non_owning_elementRank
        ///
        if (is_empty || should_put_in)
          {
            // new SubDimCellData SDC_DATA_OWNING_ELEMENT_KEY
            // CHECK
            SubDimCellData data(NodeIdsOnSubDimEntityType(needed_entity_rank.second, 0u), stk::mesh::EntityKey(element.entity_rank(), element.identifier()) );
            putInMap(subDimEntity,  data);
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
      bool checkForRemote(const stk::mesh::Entity& element, NeededEntityType& needed_entity_rank, unsigned iSubDimOrd)
      {
        EXCEPTWATCH;
        static SubDimCellData empty_SubDimCellData;
        static CommDataType buffer_entry;

        bool isGhost = m_eMesh.isGhostElement(element);

        if (!isGhost) return true;

        static SubDimCell_SDSEntityType subDimEntity;
        getSubDimEntity(subDimEntity, element, needed_entity_rank.first, iSubDimOrd);

        stk::CommAll& comm_all = m_comm_all;
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
            stk::mesh::Entity * owning_element = get_entity_element(*m_eMesh.getBulkData(), erank, owning_elementId);

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

                unsigned nidsz = nodeIds_onSE.size();
                for (unsigned iid = 0; iid < nidsz; iid++)
                  {
                    if (nodeIds_onSE[iid] == 0)
                      {
                        throw std::logic_error("logic: hmmm #5.0");
                      }

                    if (nodeIds_onSE.m_entity_id_vector.size() != nodeIds_onSE.size())
                      {
                        throw std::logic_error("NodeRegistry::checkForRemote logic err #0.1");
                      }
                    if (!nodeIds_onSE.m_entity_id_vector[iid])
                      {
                        throw std::logic_error("NodeRegistry::checkForRemote logic err #0.2");
                      }

                    //stk::mesh::Entity * new_node = get_entity_node_Ia(*m_eMesh.getBulkData(), Node, nodeIds_onSE, iid);
                    stk::mesh::Entity * new_node = nodeIds_onSE[iid];

                    if (0)
                      {
                        stk::mesh::Entity * new_node_1 = get_entity_node_I(*m_eMesh.getBulkData(), stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE.m_entity_id_vector[iid]);
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
                    m_comm_all.send_buffer( owner_proc_rank ).pack< CommDataType > (buffer_entry);
                    NodeIdsOnSubDimEntityType& nids = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
                    nids.pack(m_comm_all.send_buffer( owner_proc_rank ));
                  }
                else
                  {
                    // FIXME createNodeAndConnect(buffer_entry);
                  }
              }
          }
        return true; // FIXME
      }

      bool getFromRemote(const stk::mesh::Entity& element, NeededEntityType& needed_entity_rank, unsigned iSubDimOrd)
      {
        return checkForRemote(element, needed_entity_rank, iSubDimOrd);
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

      NodeIdsOnSubDimEntityType& getNewNodesOnSubDimEntity(const stk::mesh::Entity& element,  stk::mesh::EntityRank& needed_entity_rank, unsigned iSubDimOrd)
      {
        EXCEPTWATCH;
        static SubDimCell_SDSEntityType subDimEntity;
        getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
        static  SubDimCellData empty_SubDimCellData;

        SubDimCellData* nodeId_elementOwnderId_ptr = getFromMapPtr(subDimEntity);
        SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
        bool is_empty = nodeId_elementOwnderId_ptr == 0;

        if (is_empty)
          {
            const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);
            CellTopology cell_topo(cell_topo_data);

            std::cout << "NodeRegistry::getNewNodesOnSubDimEntity: no node found, cell_topo = " << cell_topo.getName()
                      << "\n subDimEntity= " << subDimEntity 
                      << "\n element= " << element 
                      << "\n element.entity_rank() = " << element.entity_rank()
                      << "\n needed_entity_rank= " << needed_entity_rank
                      << "\n iSubDimOrd= " << iSubDimOrd << std::endl;
            throw std::runtime_error("NodeRegistry::getNewNodesOnSubDimEntity: no node found");

            //return 0;
          }
        NodeIdsOnSubDimEntityType& nodeId = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
        return nodeId;
      }

      /// makes coordinates of this new node be the centroid of its sub entity
      void makeCentroidCoords(const stk::mesh::Entity& element,  stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd)
      {
        makeCentroidField(element, needed_entity_rank, iSubDimOrd, m_eMesh.getCoordinatesField());
      }

      void makeCentroidField(const stk::mesh::Entity& element,  stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd, stk::mesh::FieldBase *field)
      {
        //EXCEPTWATCH;

        int spatialDim = m_eMesh.getSpatialDim();
        stk::mesh::EntityRank field_rank = stk::mesh::fem::FEMMetaData::NODE_RANK;
        {
          unsigned nfr = field->restrictions().size();
          //if (printInfo) std::cout << "P[" << p_rank << "] info>    number of field restrictions= " << nfr << std::endl;
          for (unsigned ifr = 0; ifr < nfr; ifr++)
            {
              const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
              //mesh::Part& frpart = metaData.get_part(fr.ordinal());
              field_rank = fr.rank();
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

        if (is_empty)
          {
            const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);
            CellTopology cell_topo(cell_topo_data);

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
        //stk::mesh::Entity * c_node = m_eMesh.getBulkData()->get_entity(Node, nodeIds_onSE[0]);
        //stk::mesh::Entity * c_node = get_entity_node(*m_eMesh.getBulkData(), stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE[0]);
        //stk::mesh::Entity * c_node = get_entity_node_Ia(*m_eMesh.getBulkData(), stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE, 0u);
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
                stk::mesh::Entity * node = elem_nodes[ipts].entity(); //m_eMesh.getBulkData()->get_entity(Node, nodeId);
                if (!node)
                  {
                    throw std::runtime_error("makeCentroidField: bad node found 1");
                  }
                //double * const coord = stk::mesh::field_data( *field , *node );
                double *  coord = m_eMesh.field_data(field, *node, null_u);

                if (doPrint && coord)
                  {
                    const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);
                    CellTopology cell_topo(cell_topo_data);

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
                
                //stk::mesh::Entity * node = m_eMesh.getBulkData()->get_entity(Node, nodeId);
                //!!stk::mesh::Entity * node = get_entity_node_II(*m_eMesh.getBulkData(),Node, nodeId);
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
                CellTopology cell_topo(cell_topo_data);

                std::cout << "tmp NodeRegistry::makeCentroidField cell_topo = " << cell_topo.getName()
                      << "\n subDimEntity= " << subDimEntity 
                      << "\n element= " << element 
                      << "\n element.entity_rank() = " << element.entity_rank()
                      << "\n needed_entity_rank= " << needed_entity_rank
                      << "\n iSubDimOrd= " << iSubDimOrd << std::endl;
                //std::cout << "P[" << m_eMesh.getRank() << "] needed_entity_rank= " << needed_entity_rank << " coord= " << coord_str << std::endl;
              }
          }
        
      }

      /// makes coordinates of this new node be the centroid of its sub entity - this version does it for all new nodes
      void makeCentroid(stk::mesh::FieldBase *field)
      {
        EXCEPTWATCH;
        //unsigned *null_u = 0;

        int spatialDim = m_eMesh.getSpatialDim();
        stk::mesh::EntityRank field_rank = stk::mesh::fem::FEMMetaData::NODE_RANK;
        {
          EXCEPTWATCH;
          unsigned nfr = field->restrictions().size();
          //if (printInfo) std::cout << "P[" << p_rank << "] info>    number of field restrictions= " << nfr << std::endl;
          for (unsigned ifr = 0; ifr < nfr; ifr++)
            {
              const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
              //mesh::Part& frpart = metaData.get_part(fr.ordinal());
              field_rank = fr.rank();
              spatialDim = fr.dimension() ;
            }
        }
        if (field_rank != stk::mesh::fem::FEMMetaData::NODE_RANK)
          {
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

            double c_p[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            bool doPrint = false;

            if (needed_entity_rank == m_eMesh.element_rank())
              {
                EXCEPTWATCH;
                stk::mesh::Entity *element_p = 0; 
                {
                  SDSEntityType elementId = *subDimEntity.begin();
                  //!!element_p = get_entity_element(*m_eMesh.getBulkData(), m_eMesh.element_rank(), elementId);
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
                    //if (npts == 2) doPrint=true;
                    double dnpts = elem_nodes.size();
                    for (unsigned ipts = 0; ipts < npts; ipts++)
                      {
                        stk::mesh::Entity * node = elem_nodes[ipts].entity(); 
                        if (!node)
                          {
                            throw std::runtime_error("makeCentroid(field): bad node found 1.0");
                          }

                        //double *  coord = m_eMesh.field_data(field, *node, null_u);
                        double *  coord = m_eMesh.field_data_inlined(field, *node);

                        if (doPrint && coord)
                          {
                            const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);
                            CellTopology cell_topo(cell_topo_data);

                            std::cout << "tmp NodeRegistry::makeCentroid(field) cell_topo = " << cell_topo.getName() << " ipts= " << ipts
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
            else
              {
                EXCEPTWATCH;
                double dnpts = subDimEntity.size();

                for (SubDimCell_SDSEntityType::const_iterator ids = subDimEntity.begin(); ids != subDimEntity.end(); ++ids)
                  {
                    SDSEntityType nodeId = *ids;

                    //!!stk::mesh::Entity * node = get_entity_node_II(*m_eMesh.getBulkData(), mesh::Node, nodeId);
                    stk::mesh::Entity * node = nodeId;
                    if (!node)
                      {
                        throw std::runtime_error("makeCentroid(field): bad node found 2.0");
                      }
                    //double *  coord = m_eMesh.field_data(field, *node, null_u);
                    double *  coord = m_eMesh.field_data_inlined(field, *node);
                    if (coord)
                      {
                        for (int isp = 0; isp < spatialDim; isp++)
                          {
                            c_p[isp] += coord[isp]/dnpts;
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
        const stk::mesh::FieldVector & fields = m_eMesh.getFEM_meta_data()->get_fields();
        unsigned nfields = fields.size();
        //std::cout << "P[" << p_rank << "] info>    Number of fields = " << fields.size() << std::endl;
        for (unsigned ifld = 0; ifld < nfields; ifld++)
          {
            stk::mesh::FieldBase *field = fields[ifld];
            //std::cout << "P[" << m_eMesh.getRank() << "] field = " << field->name() << std::endl;
            makeCentroidField(element, needed_entity_rank, iSubDimOrd, field);
          }
      }

      /// do interpolation for all fields
      void interpolateFields()
      {
        const stk::mesh::FieldVector & fields = m_eMesh.getFEM_meta_data()->get_fields();
        unsigned nfields = fields.size();
        //std::cout << "P[" << p_rank << "] info>    Number of fields = " << fields.size() << std::endl;
        for (unsigned ifld = 0; ifld < nfields; ifld++)
          {
            stk::mesh::FieldBase *field = fields[ifld];
            //std::cout << "P[" << m_eMesh.getRank() << "] field = " << field->name() << std::endl;
            makeCentroid(field);
          }
      }


      /// check for adding new nodes to existing parts based on sub-entity part ownership

      void addToExistingParts(const stk::mesh::Entity& element,  stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd)
      {
        const std::vector< stk::mesh::Part * > & parts = m_eMesh.getFEM_meta_data()->get_parts();

        unsigned nparts = parts.size();

        //CHECK
        static SubDimCell_SDSEntityType subDimEntity;
        //subDimEntity.clear();
        getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
        static  SubDimCellData empty_SubDimCellData;
        SubDimCellData* nodeId_elementOwnderId_ptr = getFromMapPtr(subDimEntity);
        SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
        bool is_empty = nodeId_elementOwnderId_ptr == 0;

        if (is_empty)
          {
            throw std::runtime_error("addToExistingParts: no node found");
          }
        NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
        unsigned nidsz = nodeIds_onSE.size();

        for (unsigned i_nid = 0; i_nid < nidsz; i_nid++)
          {
            //stk::mesh::Entity * c_node = get_entity_node(*m_eMesh.getBulkData(), stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE[i_nid]);
            //stk::mesh::Entity * c_node = get_entity_node_Ia(*m_eMesh.getBulkData(), stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE, i_nid);
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

                //std::cout << "P[" << m_eMesh.getRank() << "] NodeRegistry::addToExistingParts Part[" << ipart << "]= " << part.name() << std::endl;
                //std::string part_name = part.name();

                // FIXME - is there a better way to determine if a part is one of the "standard" parts?
                if (stk::mesh::is_auto_declared_part(part))  //part_name[0] == '{')  //  is_auto_declared_part
                  continue;

                //std::cout << "P[" << p_rank << "] info>     Part[" << ipart << "]= " << part.name() 
                //              << " topology = " << (topology?CellTopology(topology).getName():"null")
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
                        //!!Entity * node = get_entity_node_II(*m_eMesh.getBulkData(),Node, nodeId);
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
                        m_eMesh.getBulkData()->change_entity_parts( *c_node, add_parts, remove_parts );
                        if (0)
                          {
                            std::cout << "P[" << m_eMesh.getRank() << "] adding node " << c_node->identifier() << " to   Part[" << ipart << "]= " << part.name() 
                                      << " topology = " << (topology ? CellTopology(topology).getName() : "null")
                                      << std::endl;
                          }
                      }

                  }
              }
          }
      }

      /// check for adding new nodes to existing parts based on sub-entity part ownership
      /// this version does it in bulk and thus avoids repeats on shared sub-dim entities
      void addToExistingPartsNew()
      {
        static std::vector<stk::mesh::Part*> add_parts(1, static_cast<stk::mesh::Part*>(0));
        static std::vector<stk::mesh::Part*> remove_parts;

        //std::cout << "tmp addToExistingPartsNew... " << std::endl;
        const std::vector< stk::mesh::Part * > & parts = m_eMesh.getFEM_meta_data()->get_parts();

        unsigned nparts = parts.size();
        for (unsigned ipart=0; ipart < nparts; ipart++)
          {
            stk::mesh::Part& part = *parts[ipart];

            //std::cout << "P[" << m_eMesh.getRank() << "] NodeRegistry::addToExistingParts Part[" << ipart << "]= " << part.name() << std::endl;
            //std::string part_name = part.name();

            // FIXME - is there a better way to determine if a part is one of the "standard" parts?
            if (stk::mesh::is_auto_declared_part(part)) //part_name[0] == '{')  // is_auto_declared_part
              continue;

            const CellTopologyData *const topology = stk::percept::PerceptMesh::get_cell_topology(part);
            const unsigned part_rank = part.primary_entity_rank();

            if (part_rank == stk::mesh::fem::FEMMetaData::NODE_RANK)
              {
                stk::mesh::Selector selector(part);

                //std::cout << "P[" << m_eMesh.getRank() << "] NodeRegistry::addToExistingPartsNew rank=Node = Part[" << ipart << "]= " << part.name() << std::endl;
                add_parts[0] = &part;

                SubDimCellToDataMap::iterator iter;
                //std::cout << "tmp m_cell_2_data_map.size() = " << m_cell_2_data_map.size() << std::endl;
                for (iter = m_cell_2_data_map.begin(); iter != m_cell_2_data_map.end(); ++iter)
                  {
                    const SubDimCell_SDSEntityType& subDimEntity = (*iter).first;
                    SubDimCellData& nodeId_elementOwnderId = (*iter).second;

                    NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
                    //std::cout << "P[" << p_rank << "] info>     Part[" << ipart << "]= " << part.name() 
                    //              << " topology = " << (topology?CellTopology(topology).getName():"null")
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
                              throw std::runtime_error("addToExistingParts: bad elem found 2");
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
                                std::cout << "addToExistingParts: " <<  nodeIds_onSE[i_nid] << " i_nid= " << i_nid << " nidsz= " << nidsz 
                                          << std::endl;
                                throw std::runtime_error("addToExistingParts: bad node found 0.3");
                              }

                            m_eMesh.getBulkData()->change_entity_parts( *c_node, add_parts, remove_parts );

                            if (0)
                              {
                                std::cout << "P[" << m_eMesh.getRank() << "] adding node " << c_node->identifier() << " to   Part[" << ipart << "]= " << part.name() 
                                          << " topology = " << (topology ? CellTopology(topology).getName() : "null")
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

      typedef bool (NodeRegistry::*ElementFunctionPrototype)( const stk::mesh::Entity& element, NeededEntityType& needed_entity_rank, unsigned iSubDimOrd);

      /// this is a helper method that loops over all sub-dimensional entities whose rank matches on of those in @param needed_entity_ranks
      ///    and registers that sub-dimensional entity as needing a new node.
      /// @param isGhost should be true if this element is a ghost, in which case this will call the appropriate method to set up for
      //     communications

      void //NodeRegistry::
      doForAllSubEntities(ElementFunctionPrototype function, const stk::mesh::Entity& element, vector<NeededEntityType>& needed_entity_ranks)
      {
        const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);
                
        CellTopology cell_topo(cell_topo_data);
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
                (this->*function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd);

              } // iSubDimOrd
          } // ineed_ent
      }


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

            //!
            unsigned erank = m_eMesh.element_rank();
            erank = stk::mesh::entity_rank(data.get<SDC_DATA_OWNING_ELEMENT_KEY>());
            stk::mesh::Entity * owning_element = get_entity_element(*m_eMesh.getBulkData(), erank, owning_elementId);
            //!

            if (!owning_element)
              {
                std::cout << "tmp owning_element = null, owning_elementId= " << owning_elementId 
                          << " nodeIds_onSE= " << nodeIds_onSE
                          << std::endl;
                throw std::logic_error("logic: hmmm #5.2");
              }
            if (!m_eMesh.isGhostElement(*owning_element))
              {
                //sz += 1;
                sz += nodeIds_onSE.size();
              }
          }
        return sz;
      }

      //========================================================================================================================
      // low-level interface
      // FIXME
      bool isParallelRun(unsigned size) { return true; }
      
      void checkDB()
      { 
        if (0)
          {
            unsigned sz=0;
            for (SubDimCellToDataMap::iterator cell_iter = m_cell_2_data_map.begin(); cell_iter != m_cell_2_data_map.end(); ++cell_iter)
              {
                SubDimCellData& data = (*cell_iter).second;
                stk::mesh::EntityId owning_elementId = stk::mesh::entity_id(data.get<SDC_DATA_OWNING_ELEMENT_KEY>());
                NodeIdsOnSubDimEntityType& nodeIds_onSE = data.get<SDC_DATA_GLOBAL_NODE_IDS>();

                stk::mesh::Entity * owning_element = m_eMesh.getBulkData()->get_entity(m_eMesh.element_rank(), owning_elementId);
                if (!owning_element)
                  throw std::logic_error("logic: hmmm #5.3");
                bool isGhost = m_eMesh.isGhostElement(*owning_element);
                if (!m_eMesh.isGhostElement(*owning_element))
                  {
                    ++sz;
                  }
                if (!isGhost)
                  std::cout << "P[" << m_eMesh.getRank() << "] owning_elementId = "  << owning_elementId << " isGhostElement = " << isGhost 
                            << " nodeId = " << nodeIds_onSE << std::endl;
              }
          }

      }

      /// allocate the send/recv buffers for all-to-all communication
      bool allocateBuffers()
      {
        stk::CommAll& comm_all = m_comm_all;
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
        stk::CommAll& comm_all = m_comm_all;
        //unsigned proc_size = comm_all.parallel_size();
        //unsigned proc_rank = comm_all.parallel_rank();

#if 0
        for (unsigned i_proc_rank = 0; i_proc_rank < proc_size; i_proc_rank++)
          {
            std::cout << "P[" << proc_rank << "] : i_proc_rank = " << i_proc_rank << " send buf size =  " 
                      <<   m_comm_all.send_buffer( i_proc_rank ).size() << " num in buf= " 
                      <<   m_comm_all.send_buffer( i_proc_rank ).size() / sizeof(CommDataType) <<  std::endl;
          }
#endif
        comm_all.communicate();

        stk::ParallelMachine pm = m_eMesh.getBulkData()->parallel();
        int failed = 0;
        stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );

        unpack();

      }

      void 
      unpack()
      {
        stk::CommAll& comm_all = m_comm_all;

        int failed = 0;
        std::string msg;

        stk::ParallelMachine pm = m_eMesh.getBulkData()->parallel();
        unsigned proc_size = m_eMesh.getBulkData()->parallel_size();
        unsigned proc_rank = comm_all.parallel_rank();

        vector<stk::mesh::EntityProc> nodes_to_ghost;

        if (proc_rank == 0)
          {
          }

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
            stk::mesh::Ghosting & ghosting = m_eMesh.getBulkData()->create_ghosting( std::string("new_nodes") );

            vector<stk::mesh::Entity*> receive;
            ghosting.receive_list( receive );
            m_eMesh.getBulkData()->change_ghosting( ghosting, nodes_to_ghost, receive);
          }

      }// unpack


      /// after registering all needed nodes, this method is used to request new nodes on this processor
      void createNewNodesInParallel()
      {
        unsigned num_nodes_needed = local_size();
        //std::cout << "P["<< m_eMesh.getRank() << "] num_nodes_needed= " << num_nodes_needed << std::endl;
        // FIXME
        // assert( bulk data is in modifiable mode)
        // create new entities on this proc
        vector<stk::mesh::Entity *> new_nodes;
        m_eMesh.createEntities( stk::mesh::fem::FEMMetaData::NODE_RANK, num_nodes_needed, new_nodes); 
      
        // set map values to new node id's
        unsigned inode=0;

        for (SubDimCellToDataMap::iterator cell_iter = m_cell_2_data_map.begin(); cell_iter != m_cell_2_data_map.end(); ++cell_iter)
          {
            SubDimCellData& data = (*cell_iter).second;
            stk::mesh::EntityId owning_elementId = stk::mesh::entity_id(data.get<SDC_DATA_OWNING_ELEMENT_KEY>());

            //!
            unsigned erank = m_eMesh.element_rank();
            erank = stk::mesh::entity_rank(data.get<SDC_DATA_OWNING_ELEMENT_KEY>());
            stk::mesh::Entity * owning_element = get_entity_element(*m_eMesh.getBulkData(), erank, owning_elementId);
            //!

            if (!owning_element)
              {
                throw std::logic_error("logic: hmmm #5.4");
              }
            if (!m_eMesh.isGhostElement(*owning_element))
              {
                VERIFY_OP(inode, < , num_nodes_needed, "UniformRefiner::doBreak() too many nodes");
                NodeIdsOnSubDimEntityType& nodeIds_onSE = data.get<SDC_DATA_GLOBAL_NODE_IDS>();
                if (nodeIds_onSE.m_entity_id_vector.size() != nodeIds_onSE.size())
                  {
                    throw std::logic_error("NodeRegistry:: createNewNodesInParallel logic err #0.0");
                  }

                for (unsigned ii = 0; ii < nodeIds_onSE.size(); ii++)
                  {
                    //nodeIds_onSE[ii] = new_nodes[inode]->identifier();
                    nodeIds_onSE[ii] = new_nodes[inode];
                    nodeIds_onSE.m_entity_id_vector[ii] = new_nodes[inode]->identifier();
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
        //stk::mesh::Entity * element = m_eMesh.getBulkData()->get_entity(m_eMesh.element_rank(), non_owning_elementId);

        //!
        unsigned erank = m_eMesh.element_rank();
        erank = non_owning_elementRank;
        stk::mesh::Entity * element = get_entity_element(*m_eMesh.getBulkData(), erank, non_owning_elementId);
        //!

        for (unsigned iid = 0; iid < nodeIds_onSE.size(); iid++)
          {
            //stk::mesh::Entity * node = get_entity_node_I(*m_eMesh.getBulkData(),stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE[iid]);
            //nodeIds_onSE.m_entity_vector[iid] = node;
            //stk::mesh::Entity * node = get_entity_node_Ia(*m_eMesh.getBulkData(),stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE, iid);
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

    private:
      percept::PerceptMesh& m_eMesh;
      stk::CommAll m_comm_all;
      SubDimCellToDataMap m_cell_2_data_map;
      ElementSideMap m_element_side_map;

      vector<stk::mesh::EntityProc> m_nodes_to_ghost;

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
