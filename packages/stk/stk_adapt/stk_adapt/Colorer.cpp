#include <stk_percept/function/ElementOp.hpp>

#include <stk_adapt/Colorer.hpp>



namespace stk {
  namespace adapt {

    using namespace std;
    using namespace mesh;

#if 0
    class DoElementColor : public ElementOp
    {
    public:
      virtual bool operator()(const Entity& element, FieldBase *field,  const mesh::BulkData& bulkData)
      {
        
      }
      virtual void init_elementOp(){}
      virtual void fini_elementOp(){}
    };
#endif



    template<typename STD_Set, typename Key > bool contains(STD_Set& set, Key key) { return set.find(key) != set.end(); }

    std::vector< ColorerSetType >& Colorer::
    getElementColors() { return m_element_colors; }

    void Colorer::
    color(percept::PerceptMesh& eMesh, unsigned * elementType, FieldBase *element_color_field)
    {
      const unsigned MAX_COLORS=100;
      vector< ColorerSetType > node_colors(MAX_COLORS+1); 
      m_element_colors = vector< ColorerSetType > (MAX_COLORS+1);
      ColorerSetType all_elements; 

      BulkData& bulkData = *eMesh.getBulkData();
      int ncolor=0;
      int nelem = 0;
      //for (unsigned icolor = 1; icolor <= MAX_COLORS; icolor++)
      for (unsigned icolor = 0; icolor < MAX_COLORS; icolor++)
        {
          for (unsigned irank = 0; irank < m_entityRanks.size(); irank++)
            {
              const vector<Bucket*> & buckets = bulkData.buckets( m_entityRanks[irank] );
              int num_colored_this_pass = 0;
              for ( vector<Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
                {
                  //if (select_owned(**k))  // this is where we do part selection
                  {
                    Bucket & bucket = **k ;
                    const unsigned num_elements_in_bucket = bucket.size();
                    nelem += num_elements_in_bucket;

                    if (0)
                      {
                        const CellTopologyData * const bucket_cell_topo_data = stk::mesh::get_cell_topology(bucket);
                        shards::CellTopology topo(bucket_cell_topo_data);
                        std::cout << "topo name= " << topo.getName() << std::endl;
                      }
                
                    for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                      {
                        Entity& element = bucket[iElement];
                        if (contains(all_elements, element.identifier()))
                          continue;

                        const PairIterRelation elem_nodes = element.relations( mesh::Node );  //! check for reference
                        unsigned num_node = elem_nodes.size(); 
                        bool none_in_this_color = true;
                        for (unsigned inode=0; inode < num_node; inode++)
                          {
                            Entity & node = *elem_nodes[ inode ].entity();

                            if (contains(node_colors[icolor], node.identifier()))
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
                            m_element_colors[icolor].insert(element.identifier());
                            all_elements.insert(element.identifier());
                            for (unsigned inode=0; inode < num_node; inode++)
                              {
                                Entity & node = *elem_nodes[ inode ].entity();
                                node_colors[icolor].insert(node.identifier());
                              }
                          }
                      }  // elements in bucket
                  } // selection
                } // buckets
              if (0 == num_colored_this_pass)
                {
                  break;
                }
              ++ncolor;
              if (ncolor == MAX_COLORS-1)
                {
                  throw std::runtime_error("broken algorithm in mesh colorer");
                }
            } // irank
        } // icolor

      //std::cout << "tmp ncolor = " << ncolor << " nelem= " << nelem << std::endl;

      m_element_colors.resize(ncolor);
    }

  } // namespace adapt
} // namespace stk
