// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/function/ElementOp.hpp>

#include <adapt/Colorer.hpp>
#include <percept/PerceptMesh.hpp>


  namespace percept {

    template<typename STD_Set, typename Key > bool contains(STD_Set& set, Key key) { return set.find(key) != set.end(); }

    Colorer::Colorer(std::vector< ColorerSetType >& element_colors, std::vector<stk::mesh::EntityRank> ranks ) : m_element_colors(element_colors), m_entityRanks(),
                                                                                                                 m_noColoring(true)
    {
      if (ranks.size())
        {
          m_entityRanks = ranks;
        }
      else
        {
          throw std::runtime_error("Colorer:: you must pass in non-zero length ranks");
        }
    }

    Colorer::Colorer(std::vector<stk::mesh::EntityRank> ranks ) : m_element_colors(m_element_colors_internal), m_entityRanks(), m_noColoring(true)
    {
      if (ranks.size())
        {
          m_entityRanks = ranks;
        }
      else
        {
          throw std::runtime_error("Colorer:: you must pass in non-zero length ranks");
        }
    }

    void Colorer::setNoColoring(bool no_coloring) { m_noColoring = no_coloring; }
    bool Colorer::getNoColoring() { return m_noColoring; }

    std::vector< ColorerSetType >& Colorer::
    getElementColors() { return m_element_colors; }

    void Colorer::
    color(percept::PerceptMesh& eMesh, unsigned * elementType,  stk::mesh::PartVector* fromParts, stk::mesh::FieldBase *element_color_field)
    {
      const unsigned MAX_COLORS=1000;
      std::vector< ColorerNodeSetType > node_colors(MAX_COLORS+1);
      ColorerElementSetType all_elements;

      stk::mesh::Selector selector(eMesh.get_fem_meta_data()->universal_part());
      if (fromParts)
        {
          if (0)
            {
              std::cout << "tmp Colorer::color fromParts= " << *fromParts << std::endl;
              std::cout << "tmp Colorer::color elementType= " << *elementType << std::endl;
              for (unsigned i_part = 0; i_part < fromParts->size(); i_part++)
                {
                  std::cout << "tmp Colorer::color i_part = " << i_part << " fromParts= " << (*fromParts)[i_part]->name() << std::endl;
                }
            }

          selector = stk::mesh::selectUnion(*fromParts);
        }

      stk::mesh::BulkData& bulkData = *eMesh.get_bulk_data();
      unsigned ncolor = 0;
      // int nelem = 0;
      unsigned num_max_colors = MAX_COLORS;
      if (m_noColoring)
        num_max_colors = 1;

      m_element_colors = std::vector< ColorerSetType > (num_max_colors+1);

      for (unsigned icolor = 0; icolor < num_max_colors; icolor++)
        {
          int num_colored_this_pass = 0;
          for (unsigned irank = 0; irank < m_entityRanks.size(); irank++)
            {
              const stk::mesh::BucketVector & buckets = bulkData.buckets( m_entityRanks[irank] );
              for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
                {
                  if (selector(**k))
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    bool doThisBucket = true;
                    const CellTopologyData * const bucket_cell_topo_data = eMesh.get_cell_topology(bucket);
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
                        // nelem += num_elements_in_bucket;

                        for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                          {
                            stk::mesh::Entity element = bucket[iElement];

                            if (0)
                              std::cout << "tmp color = " << icolor << " bucket topo name= " << topo.getName() << " key= " << topo.getKey()
                                        << " elementId = " << eMesh.identifier(element) << " element = " << element << std::endl;

                            stk::mesh::EntityId elem_id = eMesh.identifier(element);

                            if (!m_noColoring && contains(all_elements, elem_id))
                              continue;

                            bool none_in_this_color = true;
                            static std::vector<stk::mesh::EntityId> node_ids(100);
                            unsigned num_node = 0;

                            if (!m_noColoring)
                              {
                                const percept::MyPairIterRelation elem_nodes (eMesh, element, stk::topology::NODE_RANK );
                                num_node = elem_nodes.size();
                                node_ids.reserve(num_node);
                                for (unsigned inode=0; inode < num_node; inode++)
                                  {
                                    stk::mesh::Entity node = elem_nodes[ inode ].entity();
                                    stk::mesh::EntityId nid = eMesh.identifier(node);
                                    node_ids[inode] = nid;
                                    if (contains(node_colors[icolor], nid))
                                      {
                                        none_in_this_color = false;
                                        break;
                                      }
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
                                m_element_colors[icolor].push_back(element);
                                if (!m_noColoring)
                                  {
                                    all_elements.insert(elem_id);
                                    for (unsigned inode=0; inode < num_node; inode++)
                                      {
                                        node_colors[icolor].insert(node_ids[inode]);
                                      }
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
          if (ncolor == num_max_colors-1)
            {
              throw std::runtime_error("broken algorithm in mesh colorer");
            }
        } // icolor

      //std::cout << "tmp ncolor = " << ncolor << " nelem= " << nelem << std::endl;

      m_element_colors.resize(ncolor);
    }

  } // namespace percept
