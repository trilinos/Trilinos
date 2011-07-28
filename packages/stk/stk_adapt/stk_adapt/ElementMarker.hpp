#ifndef stk_adapt_ElementMarker_hpp
#define stk_adapt_ElementMarker_hpp

#include <stk_adapt/Marker.hpp>

namespace stk {
  namespace adapt {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     * An ElementMarker is an abstract base class for derived classes that are required to overload the mark method,
     *   which supplies the derived class with the element to be marked, and markUnrefine method specifying elements to
     *   unrefine.
     */
    class ElementMarker : public Marker
    {
    public:
      ElementMarker(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0);

      virtual ElementUnrefineCollection  buildUnrefineList() ;

    protected:

      /// Client supplies this method - given an element return instruction on what to do to the element:
      ///    0 (nothing), 1 (refine)

      virtual int mark(const stk::mesh::Entity& element) = 0;

      /// Client supplies this method - given an element return instruction on what to do to the element:
      ///    -1 (unrefine), 0 (nothing)

      virtual int markUnrefine(const stk::mesh::Entity& element) = 0;

      virtual void
      apply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity& element, 
                                              vector<NeededEntityType>& needed_entity_ranks);

    };

    // This is a very specialized test that is used in unit testing only (see unit_localRefiner/break_tri_to_tri_N_3 in UnitTestLocalRefiner.cpp)

    ElementMarker::ElementMarker(percept::PerceptMesh& eMesh, UniformRefinerPatternBase &  bp, stk::mesh::FieldBase *proc_rank_field) : 
      Marker(eMesh, bp, proc_rank_field)
    {
    }

    void ElementMarker::
    apply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity& element, 
                                            vector<NeededEntityType>& needed_entity_ranks)
    {
      const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);
                
      CellTopology cell_topo(cell_topo_data);
      const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);

      //VectorFieldType* coordField = m_eMesh.getCoordinatesField();

      int markInfo = mark(element);
#if 0
      if (markInfo <= 0)
        return;
#endif


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

          bool needNodes = (1 == markInfo);
            {
              for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                {
                  (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, needNodes);
                } // iSubDimOrd
            }
        } // ineed_ent
    }


    ElementUnrefineCollection ElementMarker::buildUnrefineList()
    {
      ElementUnrefineCollection elements_to_unref;

      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.getBulkData()->buckets( m_eMesh.element_rank() );

      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
        {
          {
            stk::mesh::Bucket & bucket = **k ;

            const unsigned num_entity_in_bucket = bucket.size();
            for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
              {
                stk::mesh::Entity& element = bucket[ientity];

                // FIXME
                // skip elements that are already a parent (if there's no family tree yet, it's not a parent, so avoid throwing an error is isParentElement)
                const bool check_for_family_tree = false;  
                bool isParent = m_eMesh.isParentElement(element, check_for_family_tree);
              
                if (isParent)
                  continue;
                
                const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);

                if (elem_nodes.size() && m_eMesh.isChildWithoutNieces(element, false))
                  {
                    int markInfo = markUnrefine(element);
                    if (markInfo < 0)
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
