#ifndef stk_adapt_NodeRegistry_hpp
#define stk_adapt_NodeRegistry_hpp

#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <cmath>
#include <utility>
#include <math.h>
#include <map>
#include <set>
#include <vector>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

#include <stk_percept/stk_mesh.hpp>
#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/Util.hpp>

#include <boost/array.hpp>

#include <boost/tuple/tuple_io.hpp>
#include <boost/tuple/tuple_comparison.hpp>

/// define only one of these to be 1
#define NODE_REGISTRY_MAP_TYPE_BOOST 1
#define NODE_REGISTRY_MAP_TYPE_TR1 0
#define NODE_REGISTRY_MAP_TYPE_STD 0

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
    // FIXME - prefix names with SDCD_
    enum SubDimCellDataEnum {
      GLOBAL_NODE_IDS,  
      OWNING_ELEMENT_ID
    };
    /// data on a sub-dim entity (global node ids on the entity, the owning element's id)
    //typedef EntityId NodeIdsOnSubDimEntityType;
    //typedef std::vector<EntityId> NodeIdsOnSubDimEntityType;

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


    typedef boost::tuple<NodeIdsOnSubDimEntityType, EntityId> SubDimCellData;

    typedef SubDimCell<EntityId> SubDimCell_EntityId;

    /// map of the node ids on a sub-dim entity to the data on the sub-dim entity
#if NODE_REGISTRY_MAP_TYPE_BOOST

    typedef boost::unordered_map<SubDimCell_EntityId, SubDimCellData, my_hash<EntityId,4>, my_equal_to<EntityId,4> > SubDimCellToDataMap;
#endif

#if NODE_REGISTRY_MAP_TYPE_TR1
    typedef tr1::unordered_map<SubDimCell_EntityId, SubDimCellData, my_hash<EntityId,4>, my_equal_to<EntityId,4> > SubDimCellToDataMap;
#endif

#if NODE_REGISTRY_MAP_TYPE_STD
#if SUBDIMCELL_USES_STL_SET
    typedef map<SubDimCell_EntityId, SubDimCellData> SubDimCellToDataMap;
#else
    typedef map<SubDimCell_EntityId, SubDimCellData, SubDimCell_compare<EntityId> > SubDimCellToDataMap;
#endif
#endif

    // Rank of sub-dim cells needing new nodes, which sub-dim entity, one non-owning element identifier, nodeId_elementOwnderId.first 
    // FIXME - consider using bitfields to pack first two entries into a short - does this save anything on buffer size?
    // FIXME - prefix names with CDT_
    enum CommDataTypeEnum {
      NEEDED_ENTITY_RANK,
      SUB_DIM_ENTITY_ORDINAL,
      NON_OWNING_ELEMENT_ID
      //,NEW_NODE_IDS
    };
    enum 
      {
        NEW_NODE_IDS
      };
    //typedef boost::tuple<EntityRank, unsigned, EntityId, NodeIdsOnSubDimEntityType> CommDataType;
    typedef boost::tuple<EntityRank, unsigned, EntityId> CommDataType;

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
      NodeRegistry(percept::PerceptMesh& eMesh) : m_eMesh(eMesh), m_comm_all(eMesh.getBulkData()->parallel())
      {

      }

      void initialize() //stk::CommAll& comm_all)
      {
        //m_comm_all = &comm_all;
        m_cell_2_data_map.clear();
      }

      void //NodeRegistry::
      beginRegistration()
      {
      }

      void //NodeRegistry::
      endRegistration()
      {
        m_eMesh.getBulkData()->modification_begin();  
        this->createNewNodesInParallel(); 
        m_nodes_to_ghost.resize(0);
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

      /// Register the need for a new node on the sub-dimensional entity @param subDimEntity on element @param element.
      /// If the element is a ghost element, the entity is still registered: the locality/ownership of the new entity
      /// can be determined by the locality of the element (ghost or not).
      bool registerNeedNewNode(const Entity& element, NeededEntityType& needed_entity_rank, unsigned iSubDimOrd) // SubDimCell_EntityId* subDimEntity=0, )
      {
        SubDimCell_EntityId subDimEntity;
        getSubDimEntity(subDimEntity, element, needed_entity_rank.first, iSubDimOrd);
        
        static SubDimCellData empty_SubDimCellData;

        SubDimCellData& nodeId_elementOwnderId = m_cell_2_data_map[subDimEntity];

        // if empty or if my id is the smallest, make this element the owner
        if (nodeId_elementOwnderId == empty_SubDimCellData ||
            element.identifier() < nodeId_elementOwnderId.get<OWNING_ELEMENT_ID>() )
          {
            //!NodeIdsOnSubDimEntityType nids(1u, 0u);
            NodeIdsOnSubDimEntityType nids(needed_entity_rank.second, 0u);
            m_cell_2_data_map[subDimEntity] = SubDimCellData(nids, element.identifier());

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
        static const SubDimCellData empty_SubDimCellData;
        static CommDataType buffer_entry;

        SubDimCell_EntityId subDimEntity;
        getSubDimEntity(subDimEntity, element, needed_entity_rank.first, iSubDimOrd);

        SubDimCellData& nodeId_elementOwnderId = m_cell_2_data_map[subDimEntity];
        stk::CommAll& comm_all = m_comm_all;
        unsigned proc_size = comm_all.parallel_size();
        unsigned proc_rank = comm_all.parallel_rank();
        unsigned owner_proc_rank = element.owner_rank();

        bool isGhost = m_eMesh.isGhostElement(element);

        if (nodeId_elementOwnderId == empty_SubDimCellData)
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
            unsigned owning_elementId = nodeId_elementOwnderId.get<OWNING_ELEMENT_ID>();
            NodeIdsOnSubDimEntityType nodeIds_onSE = nodeId_elementOwnderId.get<GLOBAL_NODE_IDS>();

            // error check
            if ( element.identifier() < owning_elementId )
              {
                std::cout << "P[" << proc_rank << "] elem id = " << element.identifier() 
                          << " nodeId_elementOwnderId.get<OWNING_ELEMENT_ID>() = " 
                          << owning_elementId
                          << std::endl;
                throw std::logic_error("logic: in getNewNode, owning element info is wrong"); 
              }

            Entity * owning_element = m_eMesh.getBulkData()->get_entity(mesh::Element, owning_elementId);
            if (!owning_element)
              throw std::logic_error("logic: hmmm #5");

            //stk::CommBuffer & send_buffer = comm_all.send_buffer( owner_proc_rank );

            bool owning_element_is_ghost = m_eMesh.isGhostElement(*owning_element);

            // if this element is a ghost, and the owning element of the node is not a ghost, send info
            //   to ghost element's owner proc
            if (!owning_element_is_ghost && isGhost)
              {
                buffer_entry = CommDataType(
                                            needed_entity_rank.first,
                                            iSubDimOrd, 
                                            element.identifier()
                                            //,nodeId_elementOwnderId.get<GLOBAL_NODE_IDS>() 
                                            );

                unsigned nidsz = nodeIds_onSE.size();
                for (unsigned iid = 0; iid < nidsz; iid++)
                  {
                    if (nodeIds_onSE[iid] == 0)
                      {
                        throw std::logic_error("logic: hmmm #5.0");
                      }
                    Entity * new_node = m_eMesh.getBulkData()->get_entity(Node, nodeIds_onSE[iid]);
                    if (!new_node)
                      throw std::logic_error("logic: hmmm #5.1");

                    m_nodes_to_ghost.push_back( EntityProc(new_node, owner_proc_rank) );
                  }

                if (isParallelRun(proc_size))
                  { 
                    //std::cout << "P[" << proc_rank << "] : pack " << buffer_entry << " owner_proc_rank= " << owner_proc_rank << std::endl;
                    m_comm_all.send_buffer( owner_proc_rank ).pack< CommDataType > (buffer_entry);
                    NodeIdsOnSubDimEntityType& nids = nodeId_elementOwnderId.get<GLOBAL_NODE_IDS>();
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

      NodeIdsOnSubDimEntityType getNewNodesOnSubDimEntity(const Entity& element,  EntityRank& needed_entity_rank, unsigned iSubDimOrd)
      {
        EXCEPTWATCH;
        SubDimCell_EntityId subDimEntity;
        getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
        static const SubDimCellData empty_SubDimCellData;
        SubDimCellData& nodeId_elementOwnderId = m_cell_2_data_map[subDimEntity];
        if (nodeId_elementOwnderId == empty_SubDimCellData)
          {
            std::cout << "NodeRegistry::getNewNodesOnSubDimEntity: no node found, subDimEntity= " << subDimEntity 
                      << " element= " << element 
                      << " element.entity_rank() = " << element.entity_rank()
                      << " needed_entity_rank= " << needed_entity_rank
                      << " iSubDimOrd= " << iSubDimOrd << std::endl;
            throw std::runtime_error("NodeRegistry::getNewNodesOnSubDimEntity: no node found");

            //return 0;
          }
        NodeIdsOnSubDimEntityType nodeId = nodeId_elementOwnderId.get<GLOBAL_NODE_IDS>();
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
        unsigned *null_u = 0;
        SubDimCell_EntityId subDimEntity;
        getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
        static const SubDimCellData empty_SubDimCellData;
        SubDimCellData& nodeId_elementOwnderId = m_cell_2_data_map[subDimEntity];
        if (nodeId_elementOwnderId == empty_SubDimCellData)
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
        NodeIdsOnSubDimEntityType nodeIds_onSE = nodeId_elementOwnderId.get<GLOBAL_NODE_IDS>();
        if (nodeIds_onSE.size() != 1)
          throw std::runtime_error("makeCentroid not ready for multiple nodes");
        Entity * c_node = m_eMesh.getBulkData()->get_entity(Node, nodeIds_onSE[0]);

        if (!c_node)
          {
            throw std::runtime_error("makeCentroid: bad node found 0");
          }

        int spatialDim = m_eMesh.getSpatialDim();
        {
          unsigned nfr = field->restrictions().size();
          //if (printInfo) std::cout << "P[" << p_rank << "] info>    number of field restrictions= " << nfr << std::endl;
          for (unsigned ifr = 0; ifr < nfr; ifr++)
            {
              const FieldRestriction& fr = field->restrictions()[ifr];
              //mesh::Part& frpart = m_eMesh.getMetaData()->get_part(fr.ordinal());
              spatialDim = fr.stride[0] ;
            }
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
                Entity * node = m_eMesh.getBulkData()->get_entity(Node, nodeId);
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

      /// check for adding new nodes to existing parts based on sub-entity part ownership
      void addToExistingParts(const Entity& element,  EntityRank needed_entity_rank, unsigned iSubDimOrd)
      {
        const std::vector< stk::mesh::Part * > & parts = m_eMesh.getMetaData()->get_parts();

        unsigned nparts = parts.size();

        SubDimCell_EntityId subDimEntity;
        getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
        static const SubDimCellData empty_SubDimCellData;
        SubDimCellData& nodeId_elementOwnderId = m_cell_2_data_map[subDimEntity];
        if (nodeId_elementOwnderId == empty_SubDimCellData)
          {
            throw std::runtime_error("addToExistingParts: no node found");
          }
        NodeIdsOnSubDimEntityType nodeIds_onSE = nodeId_elementOwnderId.get<GLOBAL_NODE_IDS>();
        unsigned nidsz=nodeIds_onSE.size();

        //if (nodeIds_onSE.size() != 1)
        //  throw std::runtime_error("addToExistingParts: not ready for multiple nodes");
        for (unsigned i_nid = 0; i_nid < nidsz; i_nid++)
          {
            Entity * c_node = m_eMesh.getBulkData()->get_entity(Node, nodeIds_onSE[i_nid]);

            if (!c_node)
              {
                std::cout << "addToExistingParts: " <<  nodeIds_onSE[i_nid] << " i_nid= " << i_nid << " nidsz= " << nidsz 
                          << " needed_entity_rank= " << needed_entity_rank << " iSubDimOrd= " << iSubDimOrd << std::endl;
                throw std::runtime_error("addToExistingParts: bad node found 0");
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
                            throw std::runtime_error("addToExistingParts: bad node found 1");
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
                    for (SubDimCell_EntityId::iterator ids = subDimEntity.begin(); ids != subDimEntity.end(); ids++)
                      {
                        EntityId nodeId = *ids;
                        Entity * node = m_eMesh.getBulkData()->get_entity(Node, nodeId);
                        if (!node)
                          {
                            throw std::runtime_error("addToExistingParts: bad node found 2");
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

      SubDimCellData& getNewNodeAndOwningElement(SubDimCell_EntityId& subDimEntity)
      {
        return m_cell_2_data_map[subDimEntity];
      }

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
                
        CellTopology cell_topo(cell_topo_data);
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
        for (SubDimCellToDataMap::iterator cell_iter = m_cell_2_data_map.begin(); cell_iter != m_cell_2_data_map.end(); cell_iter++)
          {
            SubDimCellData& data = (*cell_iter).second;
            //EntityId& owning_elementId = data.get<OWNING_ELEMENT_ID>();
            NodeIdsOnSubDimEntityType& nodeIds_onSE = data.get<GLOBAL_NODE_IDS>();

            sz += nodeIds_onSE.size();
          }
        return sz;
      }
      unsigned local_size() 
      { 
        unsigned sz=0;
        for (SubDimCellToDataMap::iterator cell_iter = m_cell_2_data_map.begin(); cell_iter != m_cell_2_data_map.end(); cell_iter++)
          {
            SubDimCellData& data = (*cell_iter).second;
            EntityId& owning_elementId = data.get<OWNING_ELEMENT_ID>();
            NodeIdsOnSubDimEntityType& nodeIds_onSE = data.get<GLOBAL_NODE_IDS>();

            Entity * owning_element = m_eMesh.getBulkData()->get_entity(mesh::Element, owning_elementId);
            if (!owning_element)
              throw std::logic_error("logic: hmmm #5.1");
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
            for (SubDimCellToDataMap::iterator cell_iter = m_cell_2_data_map.begin(); cell_iter != m_cell_2_data_map.end(); cell_iter++)
              {
                SubDimCellData& data = (*cell_iter).second;
                EntityId& owning_elementId = data.get<OWNING_ELEMENT_ID>();
                NodeIdsOnSubDimEntityType& nodeIds_onSE = data.get<GLOBAL_NODE_IDS>();

                Entity * owning_element = m_eMesh.getBulkData()->get_entity(mesh::Element, owning_elementId);
                if (!owning_element)
                  throw std::logic_error("logic: hmmm #5.1");
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
        for (SubDimCellToDataMap::iterator cell_iter = m_cell_2_data_map.begin(); cell_iter != m_cell_2_data_map.end(); cell_iter++)
          {
            SubDimCellData& data = (*cell_iter).second;
            EntityId& owning_elementId = data.get<OWNING_ELEMENT_ID>();

            Entity * owning_element = m_eMesh.getBulkData()->get_entity(mesh::Element, owning_elementId);
            if (!owning_element)
              {
                throw std::logic_error("logic: hmmm #5.2");
              }
            if (!m_eMesh.isGhostElement(*owning_element))
              {
                VERIFY_OP(inode, < , num_nodes_needed, "UniformRefiner::doBreak() too many nodes");
                NodeIdsOnSubDimEntityType& nodeIds_onSE = data.get<GLOBAL_NODE_IDS>();
                for (unsigned ii = 0; ii < nodeIds_onSE.size(); ii++)
                  {
                    nodeIds_onSE[ii] = new_nodes[inode]->identifier();
                    inode++;
                  }
                //data.get<GLOBAL_NODE_IDS>()[0] = new_nodes[inode]->identifier();
              }
          }

      }


      /// unpacks the incoming information in @param buffer_entry and adds that information to my local node registry
      /// (i.e. the map of sub-dimensional entity to global node id is updated)
      void 
      createNodeAndConnect(CommDataType& buffer_entry, NodeIdsOnSubDimEntityType nodeIds_onSE, unsigned from_proc, vector<EntityProc>& nodes_to_ghost)
      {
        EntityRank&                needed_entity_rank   = buffer_entry.get<NEEDED_ENTITY_RANK>();
        unsigned                   iSubDimOrd           = buffer_entry.get<SUB_DIM_ENTITY_ORDINAL>();
        EntityId&                  non_owning_elementId = buffer_entry.get<NON_OWNING_ELEMENT_ID>();

        // create a new relation here?  no, we are going to delete this element, so we just register that the new node is attached to 
        Entity * element = m_eMesh.getBulkData()->get_entity(mesh::Element, non_owning_elementId);
        
        for (unsigned iid = 0; iid < nodeIds_onSE.size(); iid++)
          {
            Entity * node = m_eMesh.getBulkData()->get_entity(Node, nodeIds_onSE[iid]);
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

        subDimCellData.get<GLOBAL_NODE_IDS>() = nodeIds_onSE;
        EntityId& owning_element_id = subDimCellData.get<OWNING_ELEMENT_ID>();
        VERIFY_OP(owning_element_id, !=, non_owning_elementId, "createNodeAndConnect:: bad elem ids");
        VERIFY_OP(owning_element_id, < , non_owning_elementId, "createNodeAndConnect:: bad elem ids 2");

      }

    private:
      percept::PerceptMesh& m_eMesh;
      stk::CommAll m_comm_all;
      SubDimCellToDataMap m_cell_2_data_map;

      vector<EntityProc> m_nodes_to_ghost;
    };

  }
}
#endif
