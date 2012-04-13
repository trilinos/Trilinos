#include <stk_adapt/NodeRegistry.hpp>

#if 0 && NODE_REGISTRY_MAP_TYPE_TEUCHOS_HASHTABLE
namespace Teuchos
{
  //template <class T> int hashCode(const T& x);

  template <> int hashCode(const SDCell_HashTable_Key& x)
  {
    return (int)x.getHash();
  }

}
#endif

namespace stk {
  namespace adapt {

#if !NODE_REGISTRY_MAP_ACCESSORS_INLINED
    SubDimCellData& NodeRegistry::getFromMap(SubDimCell_EntityId& subDimEntity)
    {
      //ftest(subDimEntity);
      return m_cell_2_data_map[subDimEntity];
    }
    void NodeRegistry::putInMap(SubDimCell_EntityId& subDimEntity, SubDimCellData& data)
    {
      //m_cell_2_data_map.insert(std::subDimEntity, data);
      //SubDimCellData& dataInMap = m_cell_2_data_map[subDimEntity];
      //dataInMap = data;
      //ftest(subDimEntity);
      m_cell_2_data_map[subDimEntity] = data;
    }
#endif

      /// fill
      ///    @param subDimEntity with the stk::mesh::EntityId's of
      ///    the ordinal @param iSubDimOrd sub-dimensional entity of
      ///    @param element of rank
      ///    @param needed_entity_rank
      ///
    void NodeRegistry::
      noInline_getSubDimEntity(SubDimCell_SDSEntityType& subDimEntity, const stk::mesh::Entity& element, stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd)
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

#if STK_ADAPT_HAVE_YAML_CPP
#define DEBUG_YAML 0

    void NodeRegistry::serialize_write(YAML::Emitter& emitter, std::string msg)
    {
      SubDimCellToDataMap::iterator iter;
      SubDimCellToDataMap& map = m_cell_2_data_map;
      //std::cout << msg << " tmp serialize_write map size: " << map.size() << std::endl;

      if (0) emitter << YAML::Anchor("NodeRegistry::map");   YAML_ERRCHECK;
      //emitter << YAML::Flow;      YAML_ERRCHECK;
      emitter << YAML::BeginMap;      YAML_ERRCHECK;

      // key.serialized = { nodeid_0,... : set<EntityId> }
      // value.serialized = { {new_nid0, new_nid1,...}:vector<EntityId>, {elem_own[rank, ele_id]:EntityKey} }

      int jj=0;
      for (iter = map.begin(); iter != map.end(); ++iter)
        {
          const SubDimCell_SDSEntityType& subDimEntity = (*iter).first;
          SubDimCellData& nodeId_elementOwnderId = (*iter).second;
            
          //emitter << YAML::Key << subDimEntity;
          emitter << YAML::Key;       YAML_ERRCHECK;
          emitter << YAML::Flow;      YAML_ERRCHECK;
          emitter << YAML::BeginSeq;  YAML_ERRCHECK;
          if (0) emitter << YAML::Anchor(std::string("seq")+boost::lexical_cast<std::string>(jj++));   YAML_ERRCHECK;
          for (unsigned k=0; k < subDimEntity.size(); k++)
            {
              //std::cout << " " << subDimEntity[k]->identifier() << " ";
              emitter << subDimEntity[k]->identifier();          YAML_ERRCHECK;
            }
          emitter << YAML::EndSeq;      YAML_ERRCHECK;

          NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
          stk::mesh::EntityKey& value_entity_key = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>();

          //emitter << YAML::Scalar << nodeIds_onSE.size()
          emitter << YAML::Value;          YAML_ERRCHECK;
          emitter << YAML::Flow;           YAML_ERRCHECK;
          emitter << YAML::BeginSeq;       YAML_ERRCHECK;
          emitter << stk::mesh::entity_rank(value_entity_key);
          emitter << stk::mesh::entity_id(value_entity_key);
          for (unsigned ii = 0; ii < nodeIds_onSE.size(); ii++)
            {
              //emitter << (int)nodeIds_onSE[ii]->identifier();      YAML_ERRCHECK;
              emitter << (int)nodeIds_onSE.m_entity_id_vector[ii];      YAML_ERRCHECK;
            }
          emitter << YAML::EndSeq;     YAML_ERRCHECK;
        }

      emitter << YAML::EndMap;    YAML_ERRCHECK;

    }

    void NodeRegistry::serialize_read(std::ifstream& file_in,  std::string msg, bool force_have_node)
    {
      m_eMesh.getBulkData()->modification_begin();

      YAML::Parser parser(file_in);
      YAML::Node doc;
      if (DEBUG_YAML)
        std::cout 
          << "\n serialize_read..." 
          << " YAML::NodeType::Null= " <<  YAML::NodeType::Null
          << " YAML::NodeType::Scalar= " <<  YAML::NodeType::Scalar
          << " YAML::NodeType::Sequence= " <<  YAML::NodeType::Sequence
          << " YAML::NodeType::Map= " <<  YAML::NodeType::Map 
          << std::endl;

      SubDimCellToDataMap& map = m_cell_2_data_map;
      //std::cout << msg << " tmp serialize_read map size: " << map.size() << std::endl;

      try {
        while(parser.GetNextDocument(doc)) {
          if (DEBUG_YAML) std::cout << "s_r doc.Type() = " << doc.Type() << " doc.Tag()= " << doc.Tag() << " doc.size= " << doc.size() << std::endl;
          if (doc.Type() == YAML::NodeType::Map)
            {
              for(YAML::Iterator it=doc.begin();it!=doc.end();++it) {
                typedef stk::mesh::EntityId SDSEntityType_ID;
                //typedef stk::mesh::Entity * SDSEntityType;
                SDSEntityType_ID key_quantum;
                //typedef SubDimCell<SDSEntityType> SubDimCell_SDSEntityType;
                SubDimCell_SDSEntityType key; // subDimEntity = (*iter).first;

                //struct NodeIdsOnSubDimEntityType : public std::vector<NodeIdsOnSubDimEntityTypeQuantum>
                // { 
                //     typedef std::vector<stk::mesh::EntityId> entity_id_vector_type;

                stk::mesh::EntityId value_tuple_0_quantum;
                NodeIdsOnSubDimEntityType value_tuple_0;
                //stk::mesh::EntityKey::raw_key_type value_tuple_1;

                //typedef boost::tuple<NodeIdsOnSubDimEntityType, stk::mesh::EntityKey> SubDimCellData;
                SubDimCellData value; // nodeId_elementOwnderId = (*iter).second;
                NodeIdsOnSubDimEntityType& nodeIds_onSE = value.get<SDC_DATA_GLOBAL_NODE_IDS>();
                nodeIds_onSE.resize(0);
                stk::mesh::EntityKey& value_entity_key = value.get<SDC_DATA_OWNING_ELEMENT_KEY>();
                // value = { {new_node0, new_node1,...}:[vector<Entity*>,vector<EntityId>], {elem_own[rank, ele_id]:EntityKey} }
                // key = { nodePtr_0,... : set<Entity *> }
                // value.serialized = { {new_nid0, new_nid1,...}:vector<EntityId>, {elem_own[rank, ele_id]:EntityKey} }
                // key.serialized = { nodeid_0,... : set<EntityId> }
                

                //if (DEBUG_YAML) std::cout << "it.first().Type() = " << it.first().Type() << " it.first().Tag()= " << it.first().Tag() << std::endl;
                //if (DEBUG_YAML) std::cout << "it.second().Type() = " << it.second().Type() << " it.second().Tag()= " << it.second().Tag() << std::endl;
                const YAML::Node& keySeq = it.first();
                for(YAML::Iterator itk=keySeq.begin();itk!=keySeq.end();++itk) {
                  *itk >> key_quantum;
                  if (DEBUG_YAML) std::cout << "s_r key_quantum= " << key_quantum << std::endl;
                  SDSEntityType node = m_eMesh.getBulkData()->get_entity(0, key_quantum);
                  //key.insert(const_cast<stk::mesh::Entity*>(&element) );
                  if (!node)
                    {
                      if (force_have_node)
                        throw std::runtime_error("NodeRegistry::serialize_read: null node returned from get_entity");
                      else
                        {
                          stk::mesh::PartVector parts(1, &m_eMesh.getFEM_meta_data()->universal_part());
                          node = &m_eMesh.getBulkData()->declare_entity(0, static_cast<stk::mesh::EntityId>(key_quantum), parts);
                        }
                    }
                  
                  key.insert( node );

                }
              
                int iseq=0;
                const YAML::Node& valSeq = it.second();
                stk::mesh::EntityRank rank;
                stk::mesh::EntityKey::raw_key_type id;
                for(YAML::Iterator itv=valSeq.begin();itv!=valSeq.end();++itv,++iseq) {
                  if (iseq == 0)
                    {
                      *itv >> rank;
                    }
                  else if (iseq == 1)
                    {
                      *itv >> id;
                      stk::mesh::EntityKey entityKey(rank,id);
                      if (DEBUG_YAML) std::cout << "s_r value_tuple_1= " << rank << " " << id << std::endl;
                      value_entity_key = stk::mesh::EntityKey(rank,id);
                      if (DEBUG_YAML) std::cout << "s_r owning element rank= " << stk::mesh::entity_rank(value.get<SDC_DATA_OWNING_ELEMENT_KEY>())
                                                << " owning element id= " << stk::mesh::entity_id(value.get<SDC_DATA_OWNING_ELEMENT_KEY>()) 
                                                << std::endl;
                    }
                  else
                    {
                      *itv >> value_tuple_0_quantum;

                      //stk::mesh::EntityId owning_elementId = stk::mesh::entity_id(data.get<SDC_DATA_OWNING_ELEMENT_KEY>());
                      nodeIds_onSE.m_entity_id_vector.push_back(value_tuple_0_quantum);
                      stk::mesh::Entity *entity = m_eMesh.getBulkData()->get_entity(0, value_tuple_0_quantum);
                      if (!entity)
                        {
                          if (force_have_node)
                            throw std::runtime_error("NodeRegistry::serialize_read: null node returned from get_entity 2");
                          else
                            {
                              stk::mesh::PartVector parts(1, &m_eMesh.getFEM_meta_data()->universal_part());
                              entity = &m_eMesh.getBulkData()->declare_entity(0, static_cast<stk::mesh::EntityId>(value_tuple_0_quantum), 
                                                                              parts);
                            }
                        }

                      nodeIds_onSE.push_back(entity);

                      if (DEBUG_YAML) std::cout << "s_r value_tuple_0_quantum= " << value_tuple_0_quantum << " entity= " << entity 
                                                <<  " len0= " << nodeIds_onSE.size() << " len1= " << nodeIds_onSE.m_entity_id_vector.size() << std::endl;
                    }
                }
              
                map[key] = value;
              }
#if 0
              const YAML::Node& node = *it;
              std::cout << "it= " << &node << std::endl;
              switch (it->Type())
                {
                case YAML::NodeType::Null:
                  break;
                case YAML::NodeType::Scalar:
                  break;
                case YAML::NodeType::Sequence:
                  break;
                case YAML::NodeType::Map:
                  break;
              
                }
#endif
            }
        }
      }
      catch(YAML::ParserException& e) {
        std::cout << e.what() << "\n";
        throw std::runtime_error(std::string("yaml parsing error: ")+e.what());
      }

      m_eMesh.getBulkData()->modification_end();
    }

#endif


  }
}
