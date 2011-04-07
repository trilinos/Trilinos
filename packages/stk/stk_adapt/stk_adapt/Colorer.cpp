#include <stk_percept/function/ElementOp.hpp>

#include <stk_adapt/Colorer.hpp>



namespace stk {
  namespace adapt {

    using namespace std;

    template<typename STD_Set, typename Key > bool contains(STD_Set& set, Key key) { return set.find(key) != set.end(); }

    Colorer::Colorer(std::vector< ColorerSetType >& element_colors, std::vector<stk::mesh::EntityRank> ranks ) : m_element_colors(element_colors), m_entityRanks()
      {
        //stk::mesh::EntityRank ranks[2] = {Face, Element};
        if (ranks.size())
          {
            m_entityRanks = ranks;
          }
        else
          {
            //             m_entityRanks.push_back(stk_mesh_Face);
            //             m_entityRanks.push_back(stk_mesh_Element);
            throw std::runtime_error("Colorer:: you must pass in non-zero length ranks");
          }
      }

    Colorer::Colorer(std::vector<stk::mesh::EntityRank> ranks ) : m_element_colors(m_element_colors_internal), m_entityRanks()
    {
      //stk::mesh::EntityRank ranks[2] = {Face, Element};
      if (ranks.size())
        {
          m_entityRanks = ranks;
        }
      else
        {
//           m_entityRanks.push_back(stk_mesh_Face);
//           m_entityRanks.push_back(stk_mesh_Element);
            throw std::runtime_error("Colorer:: you must pass in non-zero length ranks");
        }
    }

    std::vector< ColorerSetType >& Colorer::
    getElementColors() { return m_element_colors; }

    void Colorer::
    color(percept::PerceptMesh& eMesh, unsigned * elementType,  stk::mesh::PartVector* fromParts, stk::mesh::FieldBase *element_color_field)
    {
      const unsigned MAX_COLORS=1000;
      vector< ColorerNodeSetType > node_colors(MAX_COLORS+1); 
      m_element_colors = vector< ColorerSetType > (MAX_COLORS+1);
      ColorerElementSetType all_elements; 

      mesh::Selector selector(eMesh.getFEM_meta_data()->universal_part());
      if (fromParts) 
        {
          if (0)
            {
              //std::cout << "tmp Colorer::color fromParts= " << *fromParts << std::endl;
              //std::cout << "tmp Colorer::color elementType= " << *elementType << std::endl;
              for (unsigned i_part = 0; i_part < fromParts->size(); i_part++)
                {
                  std::cout << "tmp Colorer::color i_part = " << i_part << " fromParts= " << (*fromParts)[i_part]->name() << std::endl;
                }
            }

          selector = mesh::selectUnion(*fromParts);
        }

      stk::mesh::BulkData& bulkData = *eMesh.getBulkData();
      int ncolor = 0;
      int nelem = 0;
      for (unsigned icolor = 0; icolor < MAX_COLORS; icolor++)
        {
          int num_colored_this_pass = 0;
          for (unsigned irank = 0; irank < m_entityRanks.size(); irank++)
            {
              const vector<stk::mesh::Bucket*> & buckets = bulkData.buckets( m_entityRanks[irank] );
              for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
                {
                  if (selector(**k))  
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    bool doThisBucket = true;
                    const CellTopologyData * const bucket_cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(bucket);
                    shards::CellTopology topo(bucket_cell_topo_data);
                    if (elementType && (topo.getKey() != *elementType))
                      {
                        doThisBucket = false;
                      }

                    if (0 && doThisBucket)
                      {
                        std::cout << "tmp color = " << icolor << " bucket topo name= " << topo.getName() << " key= " << topo.getKey() 
                                  << " elementType= " << (elementType?  *elementType : 0) << " doThisBucket= " << doThisBucket << std::endl;
                      }

                    if (doThisBucket)
                      {
                        const unsigned num_elements_in_bucket = bucket.size();
                        nelem += num_elements_in_bucket;
                
                        for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                          {
                            stk::mesh::Entity& element = bucket[iElement];

                            if (0)
                              std::cout << "tmp color = " << icolor << " bucket topo name= " << topo.getName() << " key= " << topo.getKey() 
                                        << " elementId = " << element.identifier() << " element = " << element << std::endl;

                            stk::mesh::EntityId elem_id = element.identifier();
                            if (contains(all_elements, elem_id))
                              continue;

                            const stk::mesh::PairIterRelation elem_nodes = element.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );  //! check for reference
                            unsigned num_node = elem_nodes.size(); 
                            bool none_in_this_color = true;
                            static std::vector<stk::mesh::EntityId> node_ids(100);
                            node_ids.reserve(num_node);
                            for (unsigned inode=0; inode < num_node; inode++)
                              {
                                stk::mesh::Entity & node = *elem_nodes[ inode ].entity();
                                stk::mesh::EntityId nid = node.identifier();
                                node_ids[inode] = nid;
                                if (contains(node_colors[icolor], nid))
                                  {
                                    none_in_this_color = false;
                                    break;
                                  }
                              }
                            if (none_in_this_color)
                              {
                                ++num_colored_this_pass;
                                if (element_color_field)
                                  {
                                    double *fdata = stk::mesh::field_data( *static_cast<const percept::ScalarFieldType *>(element_color_field) , element );
                                    fdata[0] = double(icolor);
                                  }
#if STK_ADAPT_COLORER_SET_TYPE_USE_VECTOR
                                m_element_colors[icolor].push_back(&element);
#else
                                m_element_colors[icolor].insert(&element);
#endif
                                all_elements.insert(elem_id);
                                for (unsigned inode=0; inode < num_node; inode++)
                                  {
                                    node_colors[icolor].insert(node_ids[inode]);
                                  }
                              }
                          }  // elements in bucket
                      } // doThisBucket
                  } // selection
                } // buckets
            } // irank
          if (0 == num_colored_this_pass)
            {
              break;
            }
          ++ncolor;
          if (ncolor == MAX_COLORS-1)
            {
              throw std::runtime_error("broken algorithm in mesh colorer");
            }
        } // icolor

      //std::cout << "tmp ncolor = " << ncolor << " nelem= " << nelem << std::endl;

      m_element_colors.resize(ncolor);
    }

  } // namespace adapt
} // namespace stk
