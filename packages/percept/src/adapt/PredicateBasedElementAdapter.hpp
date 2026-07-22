// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_PredicateBasedElementAdapter_hpp
#define adapt_PredicateBasedElementAdapter_hpp

#include <functional>
#include <cmath>

#include <adapt/IAdapter.hpp>

  namespace percept {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     *  Predicate-based marker
     *
     *  The functor @class RefinePredicate should supply an operator() that returns an entry from AdaptInstruction,
     *    either to do nothing, refine, unrefine, or both refine & unrefine (useful for unit testing, etc.)
     */
    typedef std::function<bool(stk::mesh::Entity)> AdapterPredicateFunctor;

    /// This class (and derived classes) supports basic element-based marking, 
    ///   a flavor of quality-improved element-based marking,
    ///   and a high-quality transition-element approach
    /// 1. Basic element-based marking - elements are marked and isotropically refined; any remaining
    ///    marks on neighboring elements are used to refine elements, which can lead to poor quality elements
    ///    This is the default for this class
    /// 2. Quality-improved element-based marking - elements are marked isotropically, and iterations are
    ///    performed to mark neighboring elements such that longest edges get broken to improve quality.
    ///    To use: derived class QualityElementAdapter
    /// 3. High-quality transition-element approach: elements are marked isotropically, and transition
    ///    elements are then created.  Each iteration enforces the two-to-one rule so that no two
    ///    elements are more than one refine level different (@see ElementRefinePredicate).
    ///    At the start of each refine or unrefine pass, all transition elements are removed, then
    ///    refine/unrefine done, then transition elements are recreated from the remaining marks.
    ///    This algorithm is based on the hanging-node approach.
    ///    To use: HanginNodeAdapter and use one of the RefinerPattern* classes
    ///    with TransitionElement in the name.

    template<class RefinePredicate>
    class PredicateBasedElementAdapter : public IAdapter
    {
    protected:
      RefinePredicate& m_predicate_refine;
      stk::mesh::FieldBase *m_coordField;

    public:

      PredicateBasedElementAdapter(RefinePredicate& predicate_refine,
                           percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0) :
        IAdapter(eMesh, bp, proc_rank_field), m_predicate_refine(predicate_refine)
      {
      }

      RefinePredicate& getRefinePredicate() { return m_predicate_refine; }

      virtual void buildUnrefineList(ElementUnrefineCollection& elements_to_unref)
      {
        //ElementUnrefineCollection elements_to_unref(*m_eMesh.get_bulk_data());
        elements_to_unref.clear();

        const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( stk::topology::ELEMENT_RANK );

        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
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

                  const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);

                  if (elem_nodes.size())
                    {
                      if (m_predicate_refine(element) & DO_UNREFINE)
                        elements_to_unref.insert(element);
                    }
                }
            }
          }

        //return elements_to_unref;
      }

    protected:

      virtual void
      refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity element,
            vector<NeededEntityType>& needed_entity_ranks, const CellTopologyData * const bucket_topo_data)
      {
        const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);

        CellTopology cell_topo(cell_topo_data);

        bool markInfo = (m_predicate_refine(element) & DO_REFINE);

        for (unsigned ineed_ent=0; ineed_ent < needed_entity_ranks.size(); ineed_ent++)
          {
            unsigned numSubDimNeededEntities = 0;
            stk::mesh::EntityRank needed_entity_rank = needed_entity_ranks[ineed_ent].first;

            bool markInfo1 = markInfo;
            if (needed_entity_rank == m_eMesh.edge_rank())
              {
                numSubDimNeededEntities = cell_topo_data->edge_count;
              }
            else if (needed_entity_rank == m_eMesh.face_rank())
              {
                numSubDimNeededEntities = cell_topo_data->side_count;
                if (m_predicate_refine.m_mark_centroid_always && m_eMesh.get_spatial_dim() == 3)
                  markInfo1 = true;
              }
            else if (needed_entity_rank == stk::topology::ELEMENT_RANK)
              {
                numSubDimNeededEntities = 1;
                if (m_predicate_refine.m_mark_centroid_always)
                  markInfo1 = true;
              }

            // setup to mark all edges
            std::vector<bool> markInfoVec(numSubDimNeededEntities, markInfo1);

            // allow for overload here for special cases
            modify_marks(element, needed_entity_rank, markInfoVec, markInfo1);

            if (needed_entity_ranks[ineed_ent].third.size())
              {
                VERIFY_OP_ON(needed_entity_ranks[ineed_ent].third.size(), ==, numSubDimNeededEntities, "bad size");
                for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                  {
                    if (!needed_entity_ranks[ineed_ent].third[iSubDimOrd])
                      markInfoVec[iSubDimOrd] = false;
                  }
              }
            for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
              {
                (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, markInfoVec[iSubDimOrd], bucket_topo_data);
              }
          }
      }


      void getEdgeMarkInfo(stk::mesh::Entity element, stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd, bool& is_marked, bool& is_not_marked)
      {
        is_marked = false;
        is_not_marked = false;

        SubDimCell_SDCEntityType subDimEntity(&m_eMesh);
        getNodeRegistry().getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);

        static SubDimCellData new_SubDimCellData;
        static SubDimCellData empty_SubDimCellData;

        SubDimCellData* nodeId_elementOwnderId_ptr = getNodeRegistry().getFromMapPtr(subDimEntity);
        SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
        bool is_empty = nodeId_elementOwnderId_ptr == 0;
        if (!is_empty)
          {
            NodeIdsOnSubDimEntityType& nodeIds_onSE = std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnderId);
            if (nodeIds_onSE.size())
              {
                unsigned mark = nodeIds_onSE.m_mark;
                is_marked = mark & NodeRegistry::NR_MARK;
                is_not_marked = mark & NodeRegistry::NR_MARK_NONE;
              }
          }

      }

    };


  }

#endif
