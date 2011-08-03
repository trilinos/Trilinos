#ifndef stk_adapt_PredicateBasedMarker_hpp
#define stk_adapt_PredicateBasedMarker_hpp

#include <functional>

#include <stk_adapt/Marker.hpp>

namespace stk {
  namespace adapt {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     *  Predicate-based marker
     */
    typedef std::unary_function<stk::mesh::Entity& , bool> MarkerPredicateFunctor;

    template<class RefinePredicate, class UnrefinePredicate>
    class PredicateBasedMarker : public Marker
    {
      RefinePredicate m_predicate_refine;
      UnrefinePredicate m_predicate_unrefine;

    public:

      PredicateBasedMarker(RefinePredicate predicate_refine, UnrefinePredicate predicate_unrefine,
                           percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0) :
        Marker(eMesh, bp, proc_rank_field), m_predicate_refine(predicate_refine), m_predicate_unrefine(predicate_unrefine)
      {
      }

      virtual ElementUnrefineCollection  buildUnrefineList()
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

                      if (m_predicate_unrefine(element))
                        elements_to_unref.insert(&element);
                    }
                }
            }
          }

        return elements_to_unref;
      }

    protected:

      virtual void 
      apply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity& element, 
            vector<NeededEntityType>& needed_entity_ranks)
      {
        const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);
                
        CellTopology cell_topo(cell_topo_data);
        const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);

        //VectorFieldType* coordField = m_eMesh.getCoordinatesField();

        bool markInfo = m_predicate_refine(element);

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

            for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
              {
                (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, markInfo);
              }
          }
      }

    };



  }
}
#endif
