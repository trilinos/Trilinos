#ifndef stk_adapt_IElementAdapter_hpp
#define stk_adapt_IElementAdapter_hpp

#include <stk_adapt/IAdapter.hpp>

namespace stk_classic {
  namespace adapt {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     * An IElementAdapter is an abstract base class for derived classes that are required to overload the mark method,
     *   which supplies the derived class with the element to be marked for refine, unrefine, or both (@see IAdapter::AdaptInstruction)
     */
    class IElementAdapter : public IAdapter
    {
    public:
      IElementAdapter(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk_classic::mesh::FieldBase *proc_rank_field=0);

      virtual ElementUnrefineCollection  buildUnrefineList() ;

    protected:

      /// Client supplies this method - given an element return instruction on what to do to the element:
      ///    DO_NOTHING (nothing), DO_REFINE (refine), DO_UNREFINE
      virtual int mark(const stk_classic::mesh::Entity& element) = 0;

      virtual void
      refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk_classic::mesh::Entity& element, 
                                              vector<NeededEntityType>& needed_entity_ranks);

    };

    // This is a very specialized test that is used in unit testing only (see unit_localRefiner/break_tri_to_tri_N_3 in UnitTestLocalRefiner.cpp)

    IElementAdapter::IElementAdapter(percept::PerceptMesh& eMesh, UniformRefinerPatternBase &  bp, stk_classic::mesh::FieldBase *proc_rank_field) : 
      IAdapter(eMesh, bp, proc_rank_field)
    {
    }

    void IElementAdapter::
    refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk_classic::mesh::Entity& element, 
                                            vector<NeededEntityType>& needed_entity_ranks)
    {
      const CellTopologyData * const cell_topo_data = stk_classic::percept::PerceptMesh::get_cell_topology(element);
                
      CellTopology cell_topo(cell_topo_data);
      const mesh::PairIterRelation elem_nodes = element.relations(stk_classic::mesh::fem::FEMMetaData::NODE_RANK);

      //VectorFieldType* coordField = m_eMesh.get_coordinates_field();

      int markInfo = mark(element);

      for (unsigned ineed_ent=0; ineed_ent < needed_entity_ranks.size(); ineed_ent++)
        {
          unsigned numSubDimNeededEntities = 0;
          stk_classic::mesh::EntityRank needed_entity_rank = needed_entity_ranks[ineed_ent].first;

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

      const vector<stk_classic::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.element_rank() );

      for ( vector<stk_classic::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
        {
          {
            stk_classic::mesh::Bucket & bucket = **k ;

            const unsigned num_entity_in_bucket = bucket.size();
            for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
              {
                stk_classic::mesh::Entity& element = bucket[ientity];

                // FIXME
                // skip elements that are already a parent (if there's no family tree yet, it's not a parent, so avoid throwing an error is isParentElement)
                const bool check_for_family_tree = false;  
                bool isParent = m_eMesh.isParentElement(element, check_for_family_tree);
              
                if (isParent)
                  continue;
                
                const mesh::PairIterRelation elem_nodes = element.relations(stk_classic::mesh::fem::FEMMetaData::NODE_RANK);

                if (elem_nodes.size() && m_eMesh.isChildWithoutNieces(element, false))
                  {
                    int markInfo = mark(element);
                    if (markInfo & DO_UNREFINE)
                      {
                        elements_to_unref.insert(&element);
                      }
                  }
              }
          }
        }

      return elements_to_unref;
    }

  }
}
#endif
