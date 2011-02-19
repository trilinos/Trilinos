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
#define NODE_REGISTRY_MAP_TYPE_PERCEPT_HASHTABLE 0
#define NODE_REGISTRY_MAP_TYPE_TEUCHOS_HASHTABLE 0
#define NODE_REGISTRY_MAP_TYPE_GOOGLE 0

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

#if NODE_REGISTRY_MAP_TYPE_TEUCHOS_HASHTABLE
#include <Teuchos_Hashtable.hpp>
#endif



#if NODE_REGISTRY_MAP_TYPE_GOOGLE
#include <google/sparse_hash_map>
#include <google/dense_hash_map>
#endif


namespace stk {
  namespace adapt {

    using namespace stk::mesh;
    using namespace stk::percept;
    using std::vector;
    using std::map;
    using std::set;

    typedef std::pair<EntityRank, unsigned> NeededEntityType;

    // using tuple here instead of pair to allow for future expansion

    // pair of node id and the owning element for a node on a sub-dimensional entity (like a face or edge)
    enum SubDimCellDataEnum {
      SDC_DATA_GLOBAL_NODE_IDS,  
      SDC_DATA_OWNING_ELEMENT_KEY
    };
    /// data on a sub-dim entity (global node ids on the entity, the owning element's id)
    //typedef EntityId NodeIdsOnSubDimEntityType;
    //typedef std::vector<EntityId> NodeIdsOnSubDimEntityType;

#define NODE_IDS_OLD 1
#if NODE_IDS_OLD
    struct NodeIdsOnSubDimEntityType : public std::vector<EntityId> 
    {
      typedef std::vector<EntityId> base_type;
      NodeIdsOnSubDimEntityType(unsigned sz=1, EntityId id=0 ) : base_type(sz,id) {}
      void pack(CommBuffer& buff) 
      { 
        buff.pack< unsigned > ( this->size() );
        for (unsigned ii = 0; ii < this->size(); ii++)
          {
            buff.pack<EntityId>( (*this)[ii] );
          }
      }
      void unpack(CommBuffer& buff) 
      { 
        unsigned sz;
        buff.unpack< unsigned > ( sz );
        this->resize( sz );
        for (unsigned ii = 0; ii < this->size(); ii++)
          {
            buff.unpack<EntityId>( (*this)[ii] );
          }
      }
    };
#else


#define MAX_NEW_NODE_ON_FACE 27u
    class NodeIdsOnSubDimEntityType : public stk::percept::NoMallocArray<EntityId, MAX_NEW_NODE_ON_FACE>
    {
    public:
      typedef  stk::percept::NoMallocArray<EntityId, MAX_NEW_NODE_ON_FACE> base_type;

      explicit NodeIdsOnSubDimEntityType()  : base_type() {}

      NodeIdsOnSubDimEntityType(unsigned sz ) : base_type(sz,0) {}

      NodeIdsOnSubDimEntityType(unsigned sz, EntityId id ) : base_type(sz,id) {}

      void pack(CommBuffer& buff) 
      { 
        buff.pack< unsigned > ( this->size() );
        for (unsigned ii = 0; ii < this->size(); ii++)
          {
            buff.pack<EntityId>( (*this)[ii] );
          }
      }
      void unpack(CommBuffer& buff) 
      { 
        unsigned sz;
        buff.unpack< unsigned > ( sz );
        this->resize( sz );
        for (unsigned ii = 0; ii < this->size(); ii++)
          {
            buff.unpack<EntityId>( (*this)[ii] );
          }
      }
    };
#endif

    inline std::ostream &operator<<(std::ostream& out, const boost::array<EntityId, 1>& arr)
    {
      out << arr[0];
      return out;
    }


    struct NodeIdsOnSubDimEntityType1 : public boost::array<EntityId, 1>
    {
      typedef boost::array<EntityId,1> base_type;
      //NodeIdsOnSubDimEntityType() : base_type() { (*this)[0] = 0u; }
      //NodeIdsOnSubDimEntityType(unsigned sz=1, EntityId id=0 ) : base_type(sz,id) {}
    };


    //typedef boost::tuple<NodeIdsOnSubDimEntityType, EntityId> SubDimCellData;
    typedef boost::tuple<NodeIdsOnSubDimEntityType, EntityKey> SubDimCellData;

    typedef SubDimCell<EntityId> SubDimCell_EntityId;

    inline std::ostream& operator<<(std::ostream& out,  SubDimCellData& val)
    {
      out << "SDC:: node ids= " << val.get<SDC_DATA_GLOBAL_NODE_IDS>() 
          << " owning element rank= " << stk::mesh::entity_rank(val.get<SDC_DATA_OWNING_ELEMENT_KEY>())
          << " owning element id= " << stk::mesh::entity_id(val.get<SDC_DATA_OWNING_ELEMENT_KEY>());
      return out;
    }

  }
}

#if NODE_REGISTRY_MAP_TYPE_TEUCHOS_HASHTABLE
namespace Teuchos
{
  //template <class T> int hashCode(const T& x);

  typedef stk::adapt::SubDimCell_EntityId SDCell_HashTable_Key;
  typedef stk::adapt::SubDimCellData      SDCell_HashTable_Value;
  typedef Teuchos::Hashtable<SDCell_HashTable_Key, SDCell_HashTable_Value> SDCell_HashTable;

  template <> 
  inline
  int hashCode(const SDCell_HashTable_Key& x)
  {
    return (int)x.getHash();
  }

}
#endif

#if NODE_REGISTRY_MAP_TYPE_PERCEPT_HASHTABLE
namespace stk {
  namespace adapt {
    typedef stk::adapt::SubDimCell_EntityId SDCell_HashTable_Key;
    typedef stk::adapt::SubDimCellData      SDCell_HashTable_Value;
    typedef percept::Hashtable<SDCell_HashTable_Key, SDCell_HashTable_Value> SDCell_HashTable;
  }
}

namespace Teuchos
{
  //template <class T> int hashCode(const T& x);

  template <> 
  inline
  int hashCode(const stk::adapt::SDCell_HashTable_Key& x)
  {
    return (int)x.getHash();
  }

}
#endif

namespace stk {
  namespace adapt {

#define NODE_REGISTRY_ENTITY_REPO_V0 0
#if NODE_REGISTRY_ENTITY_REPO_V0
    struct EntityPtr
    {
      Entity *ptr;
      EntityPtr() : ptr(0) {}
      EntityPtr(Entity* ptrin) : ptr(ptrin) {}
      EntityPtr( const EntityPtr & rhs ) : ptr( rhs.ptr ) {}
      EntityPtr & operator = ( const EntityPtr & rhs )
      { ptr = rhs.ptr ; return *this ; }
      //       EntityPtr & operator = (  EntityPtr * rhs )
      //       { ptr = rhs ; return *this ; }
      operator Entity *() { return ptr; }
    };
#else
    typedef Entity *EntityPtr;
#endif

    /// map of the node ids on a sub-dim entity to the data on the sub-dim entity

#if NODE_REGISTRY_MAP_TYPE_BOOST
#  ifdef STK_HAVE_TBB

    typedef tbb::scalable_allocator<std::pair<SubDimCell_EntityId const, SubDimCellData> > RegistryAllocator;
    typedef boost::unordered_map<SubDimCell_EntityId, SubDimCellData, my_fast_hash<EntityId,4>, my_fast_equal_to<EntityId,4>, RegistryAllocator > SubDimCellToDataMap;

#  else

    typedef boost::unordered_map<SubDimCell_EntityId, SubDimCellData, my_fast_hash<EntityId,4>, my_fast_equal_to<EntityId,4> > SubDimCellToDataMap;
    typedef boost::unordered_map<EntityId, EntityPtr > EntityRepo;

    typedef boost::tuple<const Entity *, EntityRank, unsigned> ElementSideTuple;

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

    typedef tbb::scalable_allocator<std::pair<SubDimCell_EntityId const, SubDimCellData> > RegistryAllocator;
    typedef google::sparse_hash_map<SubDimCell_EntityId, SubDimCellData, my_hash<EntityId,4>, my_equal_to<EntityId,4>, RegistryAllocator > SubDimCellToDataMap;
    typedef google::sparse_hash_map<EntityId, EntityPtr > EntityRepo;

#  else

    typedef google::sparse_hash_map<SubDimCell_EntityId, SubDimCellData, my_hash<EntityId,4>, my_equal_to<EntityId,4> > SubDimCellToDataMap;
    typedef boost::unordered_map<EntityId, EntityPtr > EntityRepo;

#  endif
#endif

#if NODE_REGISTRY_MAP_TYPE_TR1
    typedef tr1::unordered_map<SubDimCell_EntityId, SubDimCellData, my_hash<EntityId,4>, my_equal_to<EntityId,4> > SubDimCellToDataMap;
    typedef tr1::unordered_map<EntityId, EntityPtr > EntityRepo;
#endif

#if NODE_REGISTRY_MAP_TYPE_STD

#  if STK_ADAPT_SUBDIMCELL_USES_STL_SET
    typedef map<SubDimCell_EntityId, SubDimCellData> SubDimCellToDataMap;
#  endif

#  if STK_ADAPT_SUBDIMCELL_USES_STL_VECTOR
    typedef map<SubDimCell_EntityId, SubDimCellData, SubDimCell_compare<EntityId> > SubDimCellToDataMap;
#  endif

#  if STK_ADAPT_SUBDIMCELL_USES_NO_MALLOC_ARRAY
    typedef map<SubDimCell_EntityId, SubDimCellData, SubDimCell_compare<EntityId> > SubDimCellToDataMap;
#  endif

    typedef map<EntityId, EntityPtr > EntityRepo;
#endif

#if NODE_REGISTRY_MAP_TYPE_TEUCHOS_HASHTABLE
    class SubDimCellToDataMapType : public Teuchos::SDCell_HashTable
    {
    public:
      typedef Teuchos::SDCell_HashTable base_type;
      typedef Teuchos::SDCell_HashTable_Key Key;
      typedef Teuchos::SDCell_HashTable_Value Value;

      SubDimCellToDataMapType(int capacity=101, double rehashDensity = 0.8) : base_type(capacity, rehashDensity) {}

      //Value& operator[](Key& key) { return *const_cast<Value *>(&(get(key))); }
      //const Value& operator[](Key& key) const { return get(key); }
      
    };

    typedef SubDimCellToDataMapType SubDimCellToDataMap;
    //typedef SubDimCellToDataMapType TSubDimCellToDataMap;

    typedef boost::unordered_map<EntityId, EntityPtr > EntityRepo;

#endif

#if NODE_REGISTRY_MAP_TYPE_PERCEPT_HASHTABLE
    typedef SDCell_HashTable SubDimCellToDataMap;
    typedef boost::unordered_map<EntityId, EntityPtr > EntityRepo;
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
    //typedef boost::tuple<EntityRank, unsigned, EntityId, NodeIdsOnSubDimEntityType> CommDataType;
    //!typedef boost::tuple<EntityRank, unsigned, EntityId> CommDataType;
    typedef boost::tuple<EntityRank, unsigned, EntityKey> CommDataType;

    enum {
      MODE_SERIAL,
      MODE_BUFFER_SIZING,
      MODE_BUFFERS_ALLOCD,
      MODE_SEND_DONE
      // after unpack, reset to MODE_SERIAL, i.e. it's cyclical
    };



    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    class NodeRegistry 
    {
    private:

    public:
      //========================================================================================================================
      // high-level interface
      //NodeRegistry(percept::PerceptMesh& eMesh) : m_eMesh(eMesh), m_comm_all(eMesh.getBulkData()->parallel()), m_gee_cnt(0), m_gen_cnt(0),
      //m_entity_repo(stk::mesh::EntityRankEnd)

      NodeRegistry(percept::PerceptMesh& eMesh) : m_eMesh(eMesh), m_comm_all(eMesh.getBulkData()->parallel()),
                                                  // why does this cause failures? 
                                                  //m_cell_2_data_map(eMesh.getNumberElements()*8u),
                                                  m_gee_cnt(0), m_gen_cnt(0),
                                                  m_entity_repo(stk::mesh::EntityRankEnd)
      {
#if NODE_REGISTRY_MAP_TYPE_GOOGLE
        //SubDimCell_EntityId empty_key;
        //empty_key.insert( std::numeric_limits<EntityId>::max() );
        //m_cell_2_data_map.set_empty_key(empty_key);

        SubDimCell_EntityId deleted_key;
        deleted_key.insert( std::numeric_limits<EntityId>::max() - 1u ); 
        m_cell_2_data_map.set_deleted_key(deleted_key);
#endif
      }


      void initialize() //stk::CommAll& comm_all)
      {
        //m_comm_all = &comm_all;
#if NODE_REGISTRY_MAP_TYPE_TEUCHOS_HASHTABLE || NODE_REGISTRY_MAP_TYPE_PERCEPT_HASHTABLE
        //m_cell_2_data_map.clear();
#else
        m_cell_2_data_map.clear();
#endif        
        for (unsigned i = 0; i < mesh::EntityRankEnd; i++) m_entity_repo[i].clear();
      }

      void //NodeRegistry::
      beginRegistration()
      {
      }

      void //NodeRegistry::
      endRegistration()
      {
        //putInESMap();

        m_eMesh.getBulkData()->modification_begin();  
        this->createNewNodesInParallel(); 
        m_nodes_to_ghost.resize(0);

#if STK_ADAPT_NODEREGISTRY_DO_REHASH
        m_cell_2_data_map.rehash(m_cell_2_data_map.size());
#endif
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
      }

      void //NodeRegistry::
      endCheckForRemote()
      {
        stk::ParallelMachine pm = m_eMesh.getBulkData()->parallel();
        int failed = 0;
        stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );

        this->allocateBuffers();  

#if STK_ADAPT_NODEREGISTRY_DO_REHASH
        m_cell_2_data_map.rehash(m_cell_2_data_map.size());
#endif


      }

      void //NodeRegistry::
      beginGetFromRemote()
      {

      }
      void //NodeRegistry::
      endGetFromRemote()
      {
        stk::ParallelMachine pm = m_eMesh.getBulkData()->parallel();
        int failed = 0;
        stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );

        this->communicate();  

        failed = 0;
        stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );

        if (1)
          {
            Ghosting & ghosting = m_eMesh.getBulkData()->create_ghosting( std::string("new_nodes") );

            vector<Entity*> receive;

            ghosting.receive_list( receive );

            m_eMesh.getBulkData()->change_ghosting( ghosting, m_nodes_to_ghost, receive);

          }

        failed = 0;
        stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );

        m_eMesh.getBulkData()->modification_end();  
      }


      void putInESMap(const Entity& element, EntityRank& needed_entity_rank, unsigned iSubDimOrd, SubDimCellData& data) 
      {
        ElementSideTuple est(&element, needed_entity_rank, iSubDimOrd);
        m_element_side_map[est] = data;
      }

      SubDimCellData& getFromESMap(const Entity& element, EntityRank& needed_entity_rank, unsigned iSubDimOrd)
      {
        ElementSideTuple est(&element, needed_entity_rank, iSubDimOrd);
        return m_element_side_map[est];
      }

      SubDimCellData& getFromESMap1(SubDimCell_EntityId& subDimEntity, const Entity& element, EntityRank& needed_entity_rank, unsigned iSubDimOrd)
      {
#if 0
        return getFromMap(subDimEntity);
#else
        static SubDimCellData empty_SubDimCellData;
        //ElementSideTuple est(&element, needed_entity_rank, iSubDimOrd);
        SubDimCellData& get1 = getFromESMap(element, needed_entity_rank, iSubDimOrd);
        if (get1 == empty_SubDimCellData)
          {
            SubDimCellData& get2 = getFromMap(subDimEntity);
            if (get2 == empty_SubDimCellData)
               return get2;
             else
               putInESMap(element, needed_entity_rank, iSubDimOrd, get2);
             return get2;
          }
        else
          {
            return get1;
          }
#endif
      }

      /// Register the need for a new node on the sub-dimensional entity @param subDimEntity on element @param element.
      /// If the element is a ghost element, the entity is still registered: the locality/ownership of the new entity
      /// can be determined by the locality of the element (ghost or not).
      bool registerNeedNewNode(const Entity& element, NeededEntityType& needed_entity_rank, unsigned iSubDimOrd) // SubDimCell_EntityId* subDimEntity=0, )
      {
        SubDimCell_EntityId subDimEntity;
        getSubDimEntity(subDimEntity, element, needed_entity_rank.first, iSubDimOrd);
        
#if NODE_IDS_OLD
        static SubDimCellData empty_SubDimCellData;
#else
        static SubDimCellData empty_SubDimCellData(MAX_NEW_NODE_ON_FACE,0u);
#endif

#if NODE_REGISTRY_MAP_TYPE_PERCEPT_HASHTABLE
        bool is_empty = !m_cell_2_data_map.containsKey(subDimEntity);
        SubDimCellData *nid1 = &empty_SubDimCellData;
        if (!is_empty) nid1 = &m_cell_2_data_map[subDimEntity];
        SubDimCellData& nodeId_elementOwnderId = *nid1;
#else
        SubDimCellData* nodeId_elementOwnderId_ptr = getFromMapPtr(subDimEntity);
        SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
        bool is_empty = nodeId_elementOwnderId_ptr == 0;
#endif

        // if empty or if my id is the smallest, make this element the owner
        //if (is_empty || element.identifier() < nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>() )
//         bool should_put_in = 
//           (element.identifier()  < stk::mesh::entity_id(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>()))
//           && (element.entity_rank() > stk::mesh::entity_rank(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>()));
        bool should_put_in = 
          (element.identifier()  < stk::mesh::entity_id(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>()))
          || (element.entity_rank() > stk::mesh::entity_rank(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>()));

        // once it's in, the assertion should be:
        // owning_elementId < non_owning_elementId && owning_elementRank >= non_owning_elementRank
        if (is_empty || should_put_in)
          {
#if NODE_REGISTRY_MAP_TYPE_TEUCHOS_HASHTABLE || NODE_REGISTRY_MAP_TYPE_PERCEPT_HASHTABLE
            NodeIdsOnSubDimEntityType nids(needed_entity_rank.second, 0u);
            // new SubDimCellData SDC_DATA_OWNING_ELEMENT_KEY
            m_cell_2_data_map.put(subDimEntity, SubDimCellData(nids, EntityKey(element.entity_rank(), element.identifier()) ) );
#else
            // new SubDimCellData SDC_DATA_OWNING_ELEMENT_KEY
            SubDimCellData data(NodeIdsOnSubDimEntityType(needed_entity_rank.second, 0u), EntityKey(element.entity_rank(), element.identifier()) );
            putInMap(subDimEntity,  data);
#endif

            return true;
          }
        return false;
      }

      /// check the newly registered node from the registry, which does one of three things, depending on what mode we are in:
      ///   1. counts buffer in prep for sending (just does a pack)
      ///   2. packs the buffer (after buffers are alloc'd)
      ///   3. returns the new node after all communications are done
      bool checkForRemote(const Entity& element, NeededEntityType& needed_entity_rank, unsigned iSubDimOrd)
      {
        EXCEPTWATCH;
        static SubDimCellData empty_SubDimCellData;
        static CommDataType buffer_entry;

        bool isGhost = m_eMesh.isGhostElement(element);

        if (!isGhost) return true;

        SubDimCell_EntityId subDimEntity;
        getSubDimEntity(subDimEntity, element, needed_entity_rank.first, iSubDimOrd);

        stk::CommAll& comm_all = m_comm_all;
        unsigned proc_size = comm_all.parallel_size();
        unsigned proc_rank = comm_all.parallel_rank();
        unsigned owner_proc_rank = element.owner_rank();


#if NODE_REGISTRY_MAP_TYPE_PERCEPT_HASHTABLE
        bool is_empty = !m_cell_2_data_map.containsKey(subDimEntity);
        SubDimCellData *nid1 = &empty_SubDimCellData;
        if (!is_empty) nid1 = &m_cell_2_data_map[subDimEntity];
        SubDimCellData& nodeId_elementOwnderId = *nid1;
#else
        SubDimCellData* nodeId_elementOwnderId_ptr = getFromMapPtr(subDimEntity);
        SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
        bool is_empty = nodeId_elementOwnderId_ptr == 0;

#endif

        if (is_empty)
          {
            std::cout << "element= " << element 
                      << " needed_entity_rank= " << needed_entity_rank.first<< " " << needed_entity_rank.second << std::endl;
            std::cout << "subDimEntity= " << subDimEntity << std::endl;
            std::cout << "nodeId_elementOwnderId= " << nodeId_elementOwnderId << std::endl;
            std::cout << "empty_SubDimCellData= " << empty_SubDimCellData << std::endl;
            throw std::logic_error("hmm....");
            return false;
          }
        else
          {
                                        
            unsigned owning_elementId = stk::mesh::entity_id(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>());
            unsigned owning_elementRank = stk::mesh::entity_rank(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>());

#if NODE_REGISTRY_MAP_TYPE_PERCEPT_HASHTABLE
            if (!owning_elementId)
              {
                is_empty = !m_cell_2_data_map.containsKey(subDimEntity); // true
                std::cout << "error: owning_elementId = 0" << std::endl;
              }
#endif

            NodeIdsOnSubDimEntityType nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();

            // error check
            bool isNotOK = (element.identifier() < owning_elementId) && (element.entity_rank() > owning_elementRank) ;

            //if ( element.identifier() < owning_elementId )
            if ( isNotOK )
              {
                std::cout << "P[" << proc_rank << "] elem id = " << element.identifier() 
                          << " nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>() = " 
                          << owning_elementId
                          << std::endl;
                throw std::logic_error("logic: in getNewNode, owning element info is wrong"); 
              }

            //!
            unsigned erank = mesh::Element;
            //erank = element.entity_rank();
            erank = owning_elementRank;
            //VERIFY_OP(erank, <=, owning_elementRank , "erank...");
            Entity * owning_element = get_entity_element(*m_eMesh.getBulkData(), erank, owning_elementId);
            //!

            if (!owning_element)
              throw std::logic_error("logic: hmmm #5");

            bool owning_element_is_ghost = m_eMesh.isGhostElement(*owning_element);

            // if this element is a ghost, and the owning element of the node is not a ghost, send info
            //   to ghost element's owner proc
            if (!owning_element_is_ghost && isGhost)
              {
                buffer_entry = CommDataType(
                                            needed_entity_rank.first,
                                            iSubDimOrd, 
                                            //!!element.identifier()
                                            element.key()
                                            //,nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>() 
                                            );

                unsigned nidsz = nodeIds_onSE.size();
                for (unsigned iid = 0; iid < nidsz; iid++)
                  {
                    if (nodeIds_onSE[iid] == 0)
                      {
                        throw std::logic_error("logic: hmmm #5.0");
                      }
                    //Entity * new_node = m_eMesh.getBulkData()->get_entity(Node, nodeIds_onSE[iid]);
                    Entity * new_node = get_entity_node(*m_eMesh.getBulkData(), Node, nodeIds_onSE[iid]);
                    if (!new_node)
                      throw std::logic_error("logic: hmmm #5.1");

                    m_nodes_to_ghost.push_back( EntityProc(new_node, owner_proc_rank) );
                  }

                if (isParallelRun(proc_size))
                  { 
                    //std::cout << "P[" << proc_rank << "] : pack " << buffer_entry << " owner_proc_rank= " << owner_proc_rank << std::endl;
                    m_comm_all.send_buffer( owner_proc_rank ).pack< CommDataType > (buffer_entry);
                    NodeIdsOnSubDimEntityType& nids = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
                    //m_comm_all.send_buffer( owner_proc_rank ).pack< NodeIdsOnSubDimEntityType > (nids);
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

      bool getFromRemote(const Entity& element, NeededEntityType& needed_entity_rank, unsigned iSubDimOrd)
      {
        return checkForRemote(element, needed_entity_rank, iSubDimOrd);
      }



      
#if !NODE_REGISTRY_ENTITY_REPO_V0
      inline Entity* get_entity_using_find(EntityRank& rank, const EntityId& id) const
      {
        const EntityRepo::const_iterator i = m_entity_repo[rank].find( id );
        return i != m_entity_repo[rank].end() ? i->second : NULL ;
      }
#endif

      inline Entity *get_entity(BulkData& bulk, EntityRank rank, EntityId id)
      {

#if STK_ADAPT_NODEREGISTRY_USE_ENTITY_REPO
#if NODE_REGISTRY_ENTITY_REPO_V0
        Entity* entity = m_entity_repo[rank][id];
#else
        Entity* entity = get_entity_using_find(rank, id);
#endif
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

      inline Entity *get_entity_element(BulkData& bulk, EntityRank rank, EntityId id)
      {
        //m_gee_cnt++;
        return get_entity(bulk, rank, id);
      }
      inline Entity *get_entity_node(BulkData& bulk, EntityRank rank, EntityId id)
      {
        //m_gen_cnt++;
        return get_entity(bulk, rank, id);
      }


      NodeIdsOnSubDimEntityType& getNewNodesOnSubDimEntity(const Entity& element,  EntityRank& needed_entity_rank, unsigned iSubDimOrd)
      {
        EXCEPTWATCH;
        static SubDimCell_EntityId subDimEntity;
        getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
        static  SubDimCellData empty_SubDimCellData;

#if NODE_REGISTRY_MAP_TYPE_PERCEPT_HASHTABLE
        bool is_empty = !m_cell_2_data_map.containsKey(subDimEntity);
        SubDimCellData *nid1 = &empty_SubDimCellData;
        if (!is_empty) nid1 = const_cast<SubDimCellData *>(&m_cell_2_data_map.get(subDimEntity));
        SubDimCellData& nodeId_elementOwnderId = *nid1;
#else
        SubDimCellData* nodeId_elementOwnderId_ptr = getFromMapPtr(subDimEntity);
        SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
        bool is_empty = nodeId_elementOwnderId_ptr == 0;
#endif

        if (is_empty)
          {
            const CellTopologyData * const cell_topo_data = get_cell_topology(element);
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
      void makeCentroid(const Entity& element,  EntityRank needed_entity_rank, unsigned iSubDimOrd)
      {
        makeCentroid(element, needed_entity_rank, iSubDimOrd, m_eMesh.getCoordinatesField());
      }

      void makeCentroid(const Entity& element,  EntityRank needed_entity_rank, unsigned iSubDimOrd, FieldBase *field)
      {
        EXCEPTWATCH;

        int spatialDim = m_eMesh.getSpatialDim();
        EntityRank field_rank = mesh::Node;
        {
          unsigned nfr = field->restrictions().size();
          //if (printInfo) std::cout << "P[" << p_rank << "] info>    number of field restrictions= " << nfr << std::endl;
          for (unsigned ifr = 0; ifr < nfr; ifr++)
            {
              const FieldRestriction& fr = field->restrictions()[ifr];
              //mesh::Part& frpart = metaData.get_part(fr.ordinal());
              field_rank = fr.type();
              spatialDim = fr.stride[0] ;
            }
        }

        if (field_rank != mesh::Node)
          {
            return;
          }

        unsigned *null_u = 0;
        static SubDimCell_EntityId subDimEntity;
        //subDimEntity.clear();
        getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
        static SubDimCellData empty_SubDimCellData;

#if NODE_REGISTRY_MAP_TYPE_PERCEPT_HASHTABLE
        bool is_empty = !m_cell_2_data_map.containsKey(subDimEntity);
        SubDimCellData *nid1 = &empty_SubDimCellData;
        if (!is_empty) nid1 = &m_cell_2_data_map[subDimEntity];
        SubDimCellData& nodeId_elementOwnderId = *nid1;
#else
        SubDimCellData* nodeId_elementOwnderId_ptr = getFromMapPtr(subDimEntity);
        SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
        bool is_empty = nodeId_elementOwnderId_ptr == 0;

#endif

        if (is_empty)
          {
#if 0
            std::cout << "NodeRegistry::makeCentroid: no node found, subDimEntity= " << subDimEntity 
                      << " element= " << element 
                      << " element.entity_rank() = " << element.entity_rank()
                      << " needed_entity_rank= " << needed_entity_rank
                      << " iSubDimOrd= " << iSubDimOrd << std::endl;
#endif
            throw std::runtime_error("makeCentroid: no node found");
          }
        NodeIdsOnSubDimEntityType nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
        if (nodeIds_onSE.size() != 1)
          throw std::runtime_error("makeCentroid not ready for multiple nodes");
        //Entity * c_node = m_eMesh.getBulkData()->get_entity(Node, nodeIds_onSE[0]);
        Entity * c_node = get_entity_node(*m_eMesh.getBulkData(),Node, nodeIds_onSE[0]);

        if (!c_node)
          {
            throw std::runtime_error("makeCentroid: bad node found 0");
          }

        //std::vector<double> c_p(spatialDim, 0.0);
        double c_p[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        if (needed_entity_rank == mesh::Element)
          {
            const mesh::PairIterRelation elem_nodes = element.relations(Node);
            unsigned npts = elem_nodes.size();
            double dnpts = elem_nodes.size();
            for (unsigned ipts = 0; ipts < npts; ipts++)
              {
                //EntityId nodeId = elem_nodes[ipts];
                Entity * node = elem_nodes[ipts].entity(); //m_eMesh.getBulkData()->get_entity(Node, nodeId);
                if (!node)
                  {
                    throw std::runtime_error("makeCentroid: bad node found 1");
                  }
                //double * const coord = stk::mesh::field_data( *field , *node );
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
        else
          {
            double dnpts = subDimEntity.size();
            for (SubDimCell_EntityId::iterator ids = subDimEntity.begin(); ids != subDimEntity.end(); ids++)
              {
                EntityId nodeId = *ids;
                //Entity * node = m_eMesh.getBulkData()->get_entity(Node, nodeId);
                Entity * node = get_entity_node(*m_eMesh.getBulkData(),Node, nodeId);
                if (!node)
                  {
                    throw std::runtime_error("makeCentroid: bad node found 2");
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
                //coord_str += toString(c_coord[isp])+ " ";
              }
          }
        //std::cout << "P[" << m_eMesh.getRank() << "] needed_entity_rank= " << needed_entity_rank << " coord= " << coord_str << std::endl;
        
      }

      /// makes coordinates of this new node be the centroid of its sub entity - this version does it for all new nodes
      void makeCentroid(FieldBase *field)
      {
        EXCEPTWATCH;
        unsigned *null_u = 0;

        int spatialDim = m_eMesh.getSpatialDim();
        EntityRank field_rank = mesh::Node;
        {
          EXCEPTWATCH;
          unsigned nfr = field->restrictions().size();
          //if (printInfo) std::cout << "P[" << p_rank << "] info>    number of field restrictions= " << nfr << std::endl;
          for (unsigned ifr = 0; ifr < nfr; ifr++)
            {
              const FieldRestriction& fr = field->restrictions()[ifr];
              //mesh::Part& frpart = metaData.get_part(fr.ordinal());
              field_rank = fr.type();
              spatialDim = fr.stride[0] ;
            }
        }
        if (field_rank != mesh::Node)
          {
            return;
          }

        if (!field)
          throw std::runtime_error("NodeRegistry::makeCentroid(field): field is null");


        SubDimCellToDataMap::iterator iter;

#if NODE_REGISTRY_MAP_TYPE_PERCEPT_HASHTABLE
        Teuchos::Array<SubDimCell_EntityId> keys;
        Teuchos::Array<SubDimCellData> values;
        m_cell_2_data_map.arrayify(keys, values);
        for (int itkv = 0; itkv < keys.size(); itkv++)
          {
            EXCEPTWATCH;
            SubDimCell_EntityId& subDimEntity = keys[itkv];
            SubDimCellData& nodeId_elementOwnderId = values[itkv];
#else
        for (iter = m_cell_2_data_map.begin(); iter != m_cell_2_data_map.end(); ++iter)
          {
            EXCEPTWATCH;
            const SubDimCell_EntityId& subDimEntity = (*iter).first;
            SubDimCellData& nodeId_elementOwnderId = (*iter).second;
#endif
          
            NodeIdsOnSubDimEntityType nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
            static const SubDimCellData empty_SubDimCellData;

#if NODE_REGISTRY_MAP_TYPE_PERCEPT_HASHTABLE
            bool is_empty = !m_cell_2_data_map.containsKey(subDimEntity);
#else
            bool is_empty = (nodeId_elementOwnderId == empty_SubDimCellData);
#endif

            if (is_empty)
              {
                throw new std::runtime_error("makeCentroid empty cell found");
              }

            if (nodeIds_onSE.size() != 1)
              {
                continue;
              }

            EntityRank needed_entity_rank = mesh::Node;
            // SPECIAL CASE
            // SPECIAL CASE
            // SPECIAL CASE
            if (subDimEntity.size() == 1)
              {
                needed_entity_rank = mesh::Element;
              }

            if (nodeIds_onSE[0] == 0)  
              {
                continue; 
              }

            Entity * c_node = get_entity_node(*m_eMesh.getBulkData(), mesh::Node, nodeIds_onSE[0]);

            if (!c_node)
              {
                throw std::runtime_error("makeCentroid: bad node found 0.0");
              }

            double c_p[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

            if (needed_entity_rank == mesh::Element)
              {
                EXCEPTWATCH;
                Entity *element_p = 0; 
                {
                  EntityId elementId = *subDimEntity.begin();
                  element_p = get_entity_element(*m_eMesh.getBulkData(), mesh::Element, elementId);
                  if (!element_p)
                    {
                      throw std::runtime_error("makeCentroid: bad elem found 2");
                    }
                }

                Entity& element = *element_p;
                bool element_is_ghost = m_eMesh.isGhostElement(element);
                if (element_is_ghost)
                  {
                    //std::cout << "tmp found ghost" << std::endl;
                  }
                else
                  {
                    const mesh::PairIterRelation elem_nodes = element.relations(Node);
                    unsigned npts = elem_nodes.size();
                    double dnpts = elem_nodes.size();
                    for (unsigned ipts = 0; ipts < npts; ipts++)
                      {
                        Entity * node = elem_nodes[ipts].entity(); 
                        if (!node)
                          {
                            throw std::runtime_error("makeCentroid: bad node found 1.0");
                          }

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
              }
            else
              {
                EXCEPTWATCH;
                double dnpts = subDimEntity.size();

                for (SubDimCell_EntityId::const_iterator ids = subDimEntity.begin(); ids != subDimEntity.end(); ++ids)
                  {
                    EntityId nodeId = *ids;

                    Entity * node = get_entity_node(*m_eMesh.getBulkData(), mesh::Node, nodeId);
                    if (!node)
                      {
                        throw std::runtime_error("makeCentroid: bad node found 2.0");
                      }
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

            // set coords
            {
              EXCEPTWATCH;

              double *  c_coord = m_eMesh.field_data(field, *c_node, null_u);

              if (c_coord)
                {
                  for (int isp = 0; isp < spatialDim; isp++)
                    {
                      c_coord[isp] = c_p[isp];
                    }
                }
            }
          }
      }
  
      /// do interpolation for all fields
      void interpolateFields(const Entity& element,  EntityRank needed_entity_rank, unsigned iSubDimOrd)
      {
        const FieldVector & fields = m_eMesh.getMetaData()->get_fields();
        unsigned nfields = fields.size();
        //std::cout << "P[" << p_rank << "] info>    Number of fields = " << fields.size() << std::endl;
        for (unsigned ifld = 0; ifld < nfields; ifld++)
          {
            FieldBase *field = fields[ifld];
            //std::cout << "P[" << m_eMesh.getRank() << "] field = " << field->name() << std::endl;
            makeCentroid(element, needed_entity_rank, iSubDimOrd, field);
          }
      }

      /// do interpolation for all fields
      void interpolateFields()
      {
        const FieldVector & fields = m_eMesh.getMetaData()->get_fields();
        unsigned nfields = fields.size();
        //std::cout << "P[" << p_rank << "] info>    Number of fields = " << fields.size() << std::endl;
        for (unsigned ifld = 0; ifld < nfields; ifld++)
          {
            FieldBase *field = fields[ifld];
            //std::cout << "P[" << m_eMesh.getRank() << "] field = " << field->name() << std::endl;
            makeCentroid(field);
          }
      }


      /// check for adding new nodes to existing parts based on sub-entity part ownership

      void addToExistingParts(const Entity& element,  EntityRank needed_entity_rank, unsigned iSubDimOrd)
      {
        const std::vector< stk::mesh::Part * > & parts = m_eMesh.getMetaData()->get_parts();

        unsigned nparts = parts.size();

        //CHECK
        static SubDimCell_EntityId subDimEntity;
        //subDimEntity.clear();
        getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
        static  SubDimCellData empty_SubDimCellData;
#if NODE_REGISTRY_MAP_TYPE_PERCEPT_HASHTABLE
        bool is_empty = !m_cell_2_data_map.containsKey(subDimEntity);
        SubDimCellData *nid1 = &empty_SubDimCellData;
        if (!is_empty) nid1 = &m_cell_2_data_map[subDimEntity];
        SubDimCellData& nodeId_elementOwnderId = *nid1;
#else
        SubDimCellData* nodeId_elementOwnderId_ptr = getFromMapPtr(subDimEntity);
        SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
        bool is_empty = nodeId_elementOwnderId_ptr == 0;

#endif

        if (is_empty)
          {
            throw std::runtime_error("addToExistingParts: no node found");
          }
        NodeIdsOnSubDimEntityType nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
        unsigned nidsz=nodeIds_onSE.size();

        //if (nodeIds_onSE.size() != 1)
        //  throw std::runtime_error("addToExistingParts: not ready for multiple nodes");
        for (unsigned i_nid = 0; i_nid < nidsz; i_nid++)
          {
            Entity * c_node = get_entity_node(*m_eMesh.getBulkData(), Node, nodeIds_onSE[i_nid]);

            if (!c_node)
              {
                std::cout << "addToExistingParts: " <<  nodeIds_onSE[i_nid] << " i_nid= " << i_nid << " nidsz= " << nidsz 
                          << " needed_entity_rank= " << needed_entity_rank << " iSubDimOrd= " << iSubDimOrd << std::endl;
                throw std::runtime_error("addToExistingParts: bad node found 0.1");
              }

            for (unsigned ipart=0; ipart < nparts; ipart++)
              {
                Part& part = *parts[ipart];
                mesh::Selector selector(part);

                //std::cout << "P[" << m_eMesh.getRank() << "] NodeRegistry::addToExistingParts Part[" << ipart << "]= " << part.name() << std::endl;
                std::string part_name = part.name();

                // FIXME - is there a better way to determine if a part is one of the "standard" parts?
                if (part_name[0] == '{')
                  continue;

                //std::cout << "P[" << p_rank << "] info>     Part[" << ipart << "]= " << part.name() 
                //              << " topology = " << (topology?CellTopology(topology).getName():"null")
                //              << std::endl;


                bool found = true;
                if (needed_entity_rank == mesh::Element)
                  {
                    const mesh::PairIterRelation elem_nodes = element.relations(Node);
                    unsigned npts = elem_nodes.size();
                    for (unsigned ipts = 0; ipts < npts; ipts++)
                      {
                        Entity * node = elem_nodes[ipts].entity(); 
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
                    for (SubDimCell_EntityId::iterator ids = subDimEntity.begin(); ids != subDimEntity.end(); ++ids)
                      {
                        EntityId nodeId = *ids;
                        Entity * node = get_entity_node(*m_eMesh.getBulkData(),Node, nodeId);
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
                    const CellTopologyData *const topology = stk::mesh::get_cell_topology(part);
                    const unsigned part_rank = part.primary_entity_rank();

                    //if (!topology)
                    if (part_rank == stk::mesh::Node)
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
      void addToExistingParts()
      {
        const std::vector< stk::mesh::Part * > & parts = m_eMesh.getMetaData()->get_parts();

        unsigned nparts = parts.size();
        for (unsigned ipart=0; ipart < nparts; ipart++)
          {
            Part& part = *parts[ipart];
            mesh::Selector selector(part);

            //std::cout << "P[" << m_eMesh.getRank() << "] NodeRegistry::addToExistingParts Part[" << ipart << "]= " << part.name() << std::endl;
            std::string part_name = part.name();

            // FIXME - is there a better way to determine if a part is one of the "standard" parts?
            if (part_name[0] == '{')
              continue;

            const CellTopologyData *const topology = stk::mesh::get_cell_topology(part);
            const unsigned part_rank = part.primary_entity_rank();

            if (part_rank == stk::mesh::Node)
              {
                std::vector<stk::mesh::Part*> add_parts(1, &part);
                std::vector<stk::mesh::Part*> remove_parts;

#if NODE_REGISTRY_MAP_TYPE_PERCEPT_HASHTABLE
                Teuchos::Array<SubDimCell_EntityId> keys;
                Teuchos::Array<SubDimCellData> values;
                m_cell_2_data_map.arrayify(keys, values);
                for (int itkv = 0; itkv < keys.size(); itkv++)
                  {
                    EXCEPTWATCH;
                    //const SubDimCell_EntityId& subDimEntity = (*iter).first();
                    //SubDimCellData& nodeId_elementOwnderId = (*iter).second();
                    SubDimCell_EntityId& subDimEntity = keys[itkv];
                    SubDimCellData& nodeId_elementOwnderId = values[itkv];
#else
                SubDimCellToDataMap::iterator iter;

                for (iter = m_cell_2_data_map.begin(); iter != m_cell_2_data_map.end(); ++iter)
                  {
                    const SubDimCell_EntityId& subDimEntity = (*iter).first;
                    SubDimCellData& nodeId_elementOwnderId = (*iter).second;
#endif

                    NodeIdsOnSubDimEntityType nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
                    unsigned nidsz=nodeIds_onSE.size();

                    for (unsigned i_nid = 0; i_nid < nidsz; i_nid++)
                      {
                        Entity * c_node = 0;

                        c_node = get_entity_node(*m_eMesh.getBulkData(), Node, nodeIds_onSE[i_nid]);

                        if (!c_node)
                          {
                            std::cout << "addToExistingParts: " <<  nodeIds_onSE[i_nid] << " i_nid= " << i_nid << " nidsz= " << nidsz 
                                      << std::endl;
                            throw std::runtime_error("addToExistingParts: bad node found 0.3");
                          }

                        //std::cout << "P[" << p_rank << "] info>     Part[" << ipart << "]= " << part.name() 
                        //              << " topology = " << (topology?CellTopology(topology).getName():"null")
                        //              << std::endl;

                        bool found = true;
                        EntityRank needed_entity_rank = mesh::Node;
                        //
                        // SPECIAL CASE
                        // SPECIAL CASE
                        // SPECIAL CASE
                        //
                        if( subDimEntity.size() == 1)
                          {
                            needed_entity_rank = mesh::Element;
                          }

                        if (needed_entity_rank == mesh::Element)
                          {
                            Entity *element_p = 0; 
                            {
                              EntityId elementId = *subDimEntity.begin();
                              element_p = get_entity_element(*m_eMesh.getBulkData(), mesh::Element, elementId);
                              if (!element_p)
                                {
                                  throw std::runtime_error("addToExistingParts: bad elem found 2");
                                }
                            }

                            Entity& element = *element_p;

                            const mesh::PairIterRelation elem_nodes = element.relations(Node);
                            unsigned npts = elem_nodes.size();
                            for (unsigned ipts = 0; ipts < npts; ipts++)
                              {
                                Entity * node = elem_nodes[ipts].entity(); 
                                if (!node)
                                  {
                                    throw std::runtime_error("addToExistingParts: bad node found 1.3");
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
                            for (SubDimCell_EntityId::const_iterator ids = subDimEntity.begin(); ids != subDimEntity.end(); ++ids)
                              {
                                EntityId nodeId = *ids;
                                Entity * node = get_entity_node(*m_eMesh.getBulkData(), Node, nodeId);
                                if (!node)
                                  {
                                    throw std::runtime_error("addToExistingParts: bad node found 2.3");
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
      }

      SubDimCellData& getNewNodeAndOwningElement(SubDimCell_EntityId& subDimEntity)
      {
        return m_cell_2_data_map[subDimEntity];
      }


#define NODE_REGISTRY_MAP_ACCESSORS_INLINED 1

      SubDimCellData * getFromMapPtr(const SubDimCell_EntityId& subDimEntity) const
#if NODE_REGISTRY_MAP_ACCESSORS_INLINED
        {
          const SubDimCellToDataMap::const_iterator i = m_cell_2_data_map.find( subDimEntity );
          //return *(const_cast<SubDimCellData *>(&(i->second)));

          return i != m_cell_2_data_map.end() ? (const_cast<SubDimCellData *>(&(i->second))) : 0 ;

          //return m_cell_2_data_map[subDimEntity];
        }
#else
      ;
#endif

      SubDimCellData& getFromMap(const SubDimCell_EntityId& subDimEntity) const
#if NODE_REGISTRY_MAP_ACCESSORS_INLINED
        {
          const SubDimCellToDataMap::const_iterator i = m_cell_2_data_map.find( subDimEntity );
          return * (const_cast<SubDimCellData *>(&(i->second)));

          //return i != m_entities.end() ? i->second : NULL ;

          //return m_cell_2_data_map[subDimEntity];
        }
#else
      ;
#endif
      void putInMap(SubDimCell_EntityId& subDimEntity, SubDimCellData& data)
#if NODE_REGISTRY_MAP_ACCESSORS_INLINED
      {
        m_cell_2_data_map[subDimEntity] = data;
      }
#else
      ;
#endif

      typedef bool (NodeRegistry::*ElementFunctionPrototype)( const Entity& element, NeededEntityType& needed_entity_rank, unsigned iSubDimOrd);

      /// this is a helper method that loops over all sub-dimensional entities whose rank matches on of those in @param needed_entity_ranks
      ///    and registers that sub-dimensional entity as needing a new node.
      /// @param isGhost should be true if this element is a ghost, in which case this will call the appropriate method to set up for
      //     communications

      void //NodeRegistry::
      doForAllSubEntities(ElementFunctionPrototype function, const Entity& element, vector<NeededEntityType>& needed_entity_ranks)
      {
        const CellTopologyData * const cell_topo_data = get_cell_topology(element);
                
        CellTopology cell_topo(cell_topo_data);
        const mesh::PairIterRelation elem_nodes = element.relations(Node);

        for (unsigned ineed_ent=0; ineed_ent < needed_entity_ranks.size(); ineed_ent++)
          {
            unsigned numSubDimNeededEntities = 0;
            EntityRank needed_entity_rank = needed_entity_ranks[ineed_ent].first;

            if (needed_entity_rank == Edge)
              {
                numSubDimNeededEntities = cell_topo_data->edge_count;
              }
            else if (needed_entity_rank == Face)
              {
                numSubDimNeededEntities = cell_topo_data->side_count;
              }
            else if (needed_entity_rank == mesh::Element)
              {
                numSubDimNeededEntities = 1;
              }

            for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
              {
                /// note: at this level of granularity we can do single edge refinement, hanging nodes, etc.
                //SubDimCell_EntityId subDimEntity;
                //getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
                (this->*function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd);

              } // iSubDimOrd
          } // ineed_ent
      }


      /// fill 
      ///    @param subDimEntity with the EntityId's of 
      ///    the ordinal @param iSubDimOrd sub-dimensional entity of
      ///    @param element of rank
      ///    @param needed_entity_rank
      ///
      void //NodeRegistry::
      getSubDimEntity(SubDimCell_EntityId& subDimEntity, const Entity& element, EntityRank needed_entity_rank, unsigned iSubDimOrd)
      {
        subDimEntity.clear();
        // in the case of elements, we don't share any nodes so we just make a map of element id to node
        if (needed_entity_rank == mesh::Element)
          {
            subDimEntity.insert(element.identifier());
            return;
          }

        const CellTopologyData * const cell_topo_data = get_cell_topology(element);
                
        //CellTopology cell_topo(cell_topo_data);
        const mesh::PairIterRelation elem_nodes = element.relations(Node);

        const unsigned *  inodes = 0;
        unsigned nSubDimNodes = 0;
        static const unsigned edge_nodes_2[2] = {0,1};
        static const unsigned face_nodes_3[3] = {0,1,2};
        static const unsigned face_nodes_4[4] = {0,1,2,3};

        // special case for faces in 3D
        if (needed_entity_rank == Face && needed_entity_rank == element.entity_rank())
          {
            nSubDimNodes = cell_topo_data->vertex_count;

            // note, some cells have sides with both 3 and 4 nodes (pyramid, prism)
            if (nSubDimNodes ==3 )
              inodes = face_nodes_3;
            else
              inodes = face_nodes_4;

          }
        // special case for edges in 2D
        else if (needed_entity_rank == Edge && needed_entity_rank == element.entity_rank())
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
        else if (needed_entity_rank == Edge)
          {
            inodes = cell_topo_data->edge[iSubDimOrd].node;
            nSubDimNodes = 2;
          }
        else if (needed_entity_rank == Face)
          {
            nSubDimNodes = cell_topo_data->side[iSubDimOrd].topology->vertex_count;
            // note, some cells have sides with both 3 and 4 nodes (pyramid, prism)
            inodes = cell_topo_data->side[iSubDimOrd].node;
          }

        //subDimEntity.reserve(nSubDimNodes);
        for (unsigned jnode = 0; jnode < nSubDimNodes; jnode++)
          {
            subDimEntity.insert(elem_nodes[inodes[jnode]].entity()->identifier());
          }

      }


      // FIXME
      unsigned total_size() { 
        //throw std::runtime_error("not ready");
        //return m_cell_2_data_map.size(); 
        unsigned sz=0;
#if NODE_REGISTRY_MAP_TYPE_PERCEPT_HASHTABLE
        Teuchos::Array<SubDimCell_EntityId> keys;
        Teuchos::Array<SubDimCellData> values;
        m_cell_2_data_map.arrayify(keys, values);
        for (int itkv = 0; itkv < keys.size(); itkv++)
          {
            SubDimCellData& data = values[itkv];
#else

        for (SubDimCellToDataMap::iterator cell_iter = m_cell_2_data_map.begin(); cell_iter != m_cell_2_data_map.end(); ++cell_iter)
          {
            SubDimCellData& data = (*cell_iter).second;
#endif
            NodeIdsOnSubDimEntityType& nodeIds_onSE = data.get<SDC_DATA_GLOBAL_NODE_IDS>();

            sz += nodeIds_onSE.size();
          }
        return sz;
      }

      unsigned local_size() 
      { 
        unsigned sz=0;
#if NODE_REGISTRY_MAP_TYPE_PERCEPT_HASHTABLE
        Teuchos::Array<SubDimCell_EntityId> keys;
        Teuchos::Array<SubDimCellData> values;
        m_cell_2_data_map.arrayify(keys, values);
        for (int itkv = 0; itkv < keys.size(); itkv++)
          {
            SubDimCell_EntityId& subDimEntity = keys[itkv];
            SubDimCellData& data = values[itkv];
#else
        for (SubDimCellToDataMap::iterator cell_iter = m_cell_2_data_map.begin(); cell_iter != m_cell_2_data_map.end(); ++cell_iter)
          {
            SubDimCellData& data = (*cell_iter).second;
#endif
            EntityId owning_elementId = stk::mesh::entity_id(data.get<SDC_DATA_OWNING_ELEMENT_KEY>());

#if NODE_REGISTRY_MAP_TYPE_PERCEPT_HASHTABLE
            if (!owning_elementId)
              {
                bool inTable = m_cell_2_data_map.containsKey(subDimEntity);
                
                std::cout << "tmp ERROR: owning_elementId= " << owning_elementId << " key= " << subDimEntity << " inTable= " << inTable << std::endl;
              }
#endif

            NodeIdsOnSubDimEntityType& nodeIds_onSE = data.get<SDC_DATA_GLOBAL_NODE_IDS>();

            //!
            unsigned erank = mesh::Element;
            erank = stk::mesh::entity_rank(data.get<SDC_DATA_OWNING_ELEMENT_KEY>());
            Entity * owning_element = get_entity_element(*m_eMesh.getBulkData(), erank, owning_elementId);
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
#if NODE_REGISTRY_MAP_TYPE_PERCEPT_HASHTABLE
            Teuchos::Array<SubDimCell_EntityId> keys;
            Teuchos::Array<SubDimCellData> values;
            m_cell_2_data_map.arrayify(keys, values);
            for (int itkv = 0; itkv < keys.size(); itkv++)
              {
                SubDimCellData& data = values[itkv];
#else
            for (SubDimCellToDataMap::iterator cell_iter = m_cell_2_data_map.begin(); cell_iter != m_cell_2_data_map.end(); ++cell_iter)
              {
                SubDimCellData& data = (*cell_iter).second;
#endif
                EntityId owning_elementId = stk::mesh::entity_id(data.get<SDC_DATA_OWNING_ELEMENT_KEY>());
                NodeIdsOnSubDimEntityType& nodeIds_onSE = data.get<SDC_DATA_GLOBAL_NODE_IDS>();

                Entity * owning_element = m_eMesh.getBulkData()->get_entity(mesh::Element, owning_elementId);
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

        vector<EntityProc> nodes_to_ghost;

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
                    // typedef boost::tuple::tuple<EntityRank, unsigned, EntityId, EntityId> CommDataType;


                    recv_buffer.unpack< CommDataType >( buffer_entry );
                    nodeIds_onSE.unpack(recv_buffer);
                    //recv_buffer.unpack< NodeIdsOnSubDimEntityType > (nodeIds_onSE);


                    //std::cout << "P[" << proc_rank << "] unpack for buffer from proc= " << from_proc << " " << buffer_entry << std::endl;
                    createNodeAndConnect(buffer_entry, nodeIds_onSE, from_proc, nodes_to_ghost);
                  }
              }
          }
        catch ( std::exception &x )
          {
            failed = 1;
            msg = x.what();
          }

        stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );
        if ( failed )
          {
            throw std::runtime_error( msg );
          }

        if (nodes_to_ghost.size())
          {
            Ghosting & ghosting = m_eMesh.getBulkData()->create_ghosting( std::string("new_nodes") );

            vector<Entity*> receive;
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
        vector<Entity *> new_nodes;
        m_eMesh.createEntities( Node, num_nodes_needed, new_nodes); 
      
        // set map values to new node id's
        unsigned inode=0;

#if NODE_REGISTRY_MAP_TYPE_PERCEPT_HASHTABLE
        Teuchos::Array<SubDimCell_EntityId> keys;
        Teuchos::Array<SubDimCellData> values;
        m_cell_2_data_map.arrayify(keys, values);
        for (int itkv = 0; itkv < keys.size(); itkv++)
          {
            SubDimCellData& data = values[itkv];
#else
        for (SubDimCellToDataMap::iterator cell_iter = m_cell_2_data_map.begin(); cell_iter != m_cell_2_data_map.end(); ++cell_iter)
          {
            SubDimCellData& data = (*cell_iter).second;
#endif
            EntityId owning_elementId = stk::mesh::entity_id(data.get<SDC_DATA_OWNING_ELEMENT_KEY>());

            //!
            unsigned erank = mesh::Element;
            erank = stk::mesh::entity_rank(data.get<SDC_DATA_OWNING_ELEMENT_KEY>());
            Entity * owning_element = get_entity_element(*m_eMesh.getBulkData(), erank, owning_elementId);
            //!

            if (!owning_element)
              {
                throw std::logic_error("logic: hmmm #5.4");
              }
            if (!m_eMesh.isGhostElement(*owning_element))
              {
                VERIFY_OP(inode, < , num_nodes_needed, "UniformRefiner::doBreak() too many nodes");
                NodeIdsOnSubDimEntityType& nodeIds_onSE = data.get<SDC_DATA_GLOBAL_NODE_IDS>();
                for (unsigned ii = 0; ii < nodeIds_onSE.size(); ii++)
                  {
                    nodeIds_onSE[ii] = new_nodes[inode]->identifier();
                    inode++;
                  }
                //data.get<SDC_DATA_GLOBAL_NODE_IDS>()[0] = new_nodes[inode]->identifier();
              }
          }

      }


      /// unpacks the incoming information in @param buffer_entry and adds that information to my local node registry
      /// (i.e. the map of sub-dimensional entity to global node id is updated)
      void 
      createNodeAndConnect(CommDataType& buffer_entry, NodeIdsOnSubDimEntityType nodeIds_onSE, unsigned from_proc, vector<EntityProc>& nodes_to_ghost)
      {
        //EntityId&                  non_owning_elementId = buffer_entry.get<CDT_NON_OWNING_ELEMENT_KEY>();

        EntityRank& needed_entity_rank                    = buffer_entry.get<CDT_NEEDED_ENTITY_RANK>();
        unsigned    iSubDimOrd                            = buffer_entry.get<CDT_SUB_DIM_ENTITY_ORDINAL>();
        EntityKey&  non_owning_elementKey                 = buffer_entry.get<CDT_NON_OWNING_ELEMENT_KEY>();
        EntityId    non_owning_elementId                  = stk::mesh::entity_id(non_owning_elementKey);
        EntityRank  non_owning_elementRank                = stk::mesh::entity_rank(non_owning_elementKey);

        // create a new relation here?  no, we are going to delete this element, so we just register that the new node is attached to 
        //Entity * element = m_eMesh.getBulkData()->get_entity(mesh::Element, non_owning_elementId);

        //!
        unsigned erank = mesh::Element;
        erank = non_owning_elementRank;
        Entity * element = get_entity_element(*m_eMesh.getBulkData(), erank, non_owning_elementId);
        //!

        for (unsigned iid = 0; iid < nodeIds_onSE.size(); iid++)
          {
            Entity * node = get_entity_node(*m_eMesh.getBulkData(),Node, nodeIds_onSE[iid]);
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
                  
        SubDimCell_EntityId subDimEntity;
        getSubDimEntity(subDimEntity, *element, needed_entity_rank, iSubDimOrd);
        SubDimCellData& subDimCellData = getNewNodeAndOwningElement(subDimEntity);
        // assert it is empty?

        subDimCellData.get<SDC_DATA_GLOBAL_NODE_IDS>() = nodeIds_onSE;

#ifndef NDEBUG
        EntityId owning_element_id = stk::mesh::entity_id(subDimCellData.get<SDC_DATA_OWNING_ELEMENT_KEY>());
        EntityRank owning_element_rank = stk::mesh::entity_rank(subDimCellData.get<SDC_DATA_OWNING_ELEMENT_KEY>());
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

#if NODE_REGISTRY_MAP_TYPE_TEUCHOS_HASHTABLE
      //TSubDimCellToDataMap m_cell_2_data_map_t;
#endif

      vector<EntityProc> m_nodes_to_ghost;

    public:
      int m_gee_cnt;
      int m_gen_cnt;
      //const CellTopologyData * const m_cell_topo_data;
      std::vector<EntityRepo> m_entity_repo;
    };

  }
}
#endif
