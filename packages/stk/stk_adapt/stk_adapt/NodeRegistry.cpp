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

#if STK_ADAPT_HAVE_YAML_CPP
#define DEBUG_YAML 1
#define YAML_ERRCHECK do { if (DEBUG_YAML && !emitter.good()) { std::cout << "Emitter error: " << __FILE__ << ":" << __LINE__ << " emitter.good()= " \
                                                                          << emitter.good() << " Error Message: " << emitter.GetLastError() << std::endl; return;} } while(0)

    void NodeRegistry::serialize_write(YAML::Emitter& emitter, std::string msg)
    {
      SubDimCellToDataMap::iterator iter;
      SubDimCellToDataMap& map = m_cell_2_data_map;
      std::cout << msg << " tmp serialize_write map size: " << map.size() << std::endl;

      if (0) emitter << YAML::Anchor("NodeRegistry::map");   YAML_ERRCHECK;
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
          emitter << YAML::BeginSeq;  YAML_ERRCHECK;
          if (0) emitter << YAML::Anchor(std::string("seq")+boost::lexical_cast<std::string>(jj++));   YAML_ERRCHECK;
          for (unsigned k=0; k < subDimEntity.size(); k++)
            {
              //std::cout << " " << subDimEntity[k]->identifier() << " ";
              emitter << subDimEntity[k]->identifier();          YAML_ERRCHECK;
            }
          emitter << YAML::EndSeq;      YAML_ERRCHECK;

          emitter << YAML::Value;          YAML_ERRCHECK;
          NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
          stk::mesh::EntityKey& value_entity_key = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>();

          //emitter << YAML::Scalar << nodeIds_onSE.size()
          emitter << YAML::BeginSeq;       YAML_ERRCHECK;
          emitter << stk::mesh::entity_rank(value_entity_key);
          emitter << stk::mesh::entity_id(value_entity_key);
          for (unsigned ii = 0; ii < nodeIds_onSE.size(); ii++)
            {
              emitter << (int)nodeIds_onSE[ii]->identifier();      YAML_ERRCHECK;
            }
          emitter << YAML::EndSeq;     YAML_ERRCHECK;
        }

      emitter << YAML::EndMap;    YAML_ERRCHECK;

    }

    void NodeRegistry::serialize_read(std::ifstream& file_in, std::string msg)
    {
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
      std::cout << msg << " tmp serialize_read map size: " << map.size() << std::endl;

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
                  if (!node) throw std::runtime_error("null node");
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

    }

#endif


  }
}
