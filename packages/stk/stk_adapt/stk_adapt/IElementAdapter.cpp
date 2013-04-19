#include <stk_adapt/IElementAdapter.hpp>

namespace stk {
  namespace adapt {

    void IElementAdapter::
    refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity element,
                                            vector<NeededEntityType>& needed_entity_ranks)
    {
      const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);

      CellTopology cell_topo(cell_topo_data);
      const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::MetaData::NODE_RANK);

      //VectorFieldType* coordField = m_eMesh.get_coordinates_field();

      int markInfo = markElement(element);

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
          else if (needed_entity_rank == stk::mesh::MetaData::ELEMENT_RANK)
            {
              numSubDimNeededEntities = 1;
            }

          bool needNodes = ( markInfo & DO_REFINE);
            {
              for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                {
                  (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, needNodes);
                } // iSubDimOrd
            }
        } // ineed_ent
    }


    ElementUnrefineCollection IElementAdapter::buildUnrefineList()
    {
      ElementUnrefineCollection elements_to_unref;

      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( stk::mesh::MetaData::ELEMENT_RANK );

      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          {
            stk::mesh::Bucket & bucket = **k ;

            const unsigned num_entity_in_bucket = bucket.size();
            for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
              {
                stk::mesh::Entity element = bucket[ientity];

                // FIXME
                // skip elements that are already a parent (if there's no family tree yet, it's not a parent, so avoid throwing an error is isParentElement)
                const bool check_for_family_tree = false;
                bool isParent = m_eMesh.isParentElement(element, check_for_family_tree);

                if (isParent)
                  continue;

                const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::MetaData::NODE_RANK);

                if (elem_nodes.size())
                  {
                    int markInfo = markElement(element);
                    if (markInfo & DO_UNREFINE)
                      {
                        elements_to_unref.insert(element);
                      }
                  }
              }
          }
        }

      return elements_to_unref;
    }

  }
}

