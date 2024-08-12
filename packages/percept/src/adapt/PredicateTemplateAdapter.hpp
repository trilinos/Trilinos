// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_PredicateTemplateAdapter_hpp
#define adapt_PredicateTemplateAdapter_hpp

#include <functional>
#include <cmath>

#include <adapt/IAdapter.hpp>
#include <stk_mesh/base/MeshUtils.hpp>

  namespace percept {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     *  Predicate-based marker, using templates
     *
     *  The functor @class RefinePredicate should supply an operator() that returns an entry from AdaptInstruction,
     *    either to do nothing, refine, unrefine, or both refine & unrefine (useful for unit testing, etc.)
     *
     *  Template algorithm (see Mavriplis [2000])
     *  1. mark only non-transition elements (NTE) - one iter only
     *  2. mark all edges of parents of transition elements (TE), for TE's
     *       that are requested to be refined, or that have a neighboring NTE marked
     *       for refine - multiple iters
     *  3. remove all TE's
     *  4. upgrade all marks iteratively until the template
     *       requirements are satisfied (a mark pattern matches a template) - multiple iters
     *  5. refine recursively using similar algorithm to unrefinement remesh
     *
     *  Note: the steps above are now embedded in Stage_2_Mark_TE_Parents, so we only have 2 stages now.
     */
    typedef std::function<bool(stk::mesh::Entity)> AdapterPredicateFunctor;
    enum PTA_Stage {
      Stage_None,
      Stage_1_Mark_NTE,
      Stage_2_Mark_TE_Parents
    };

    template<class RefinePredicate>
    class PredicateTemplateAdapter : public IAdapter
    {
      RefinePredicate& m_predicate_refine;
      bool m_useQualityMarking;
      stk::mesh::FieldBase *m_coordField;

      int m_stage;  // stage of the algorithm corresponding to comments above
      int m_iter;
      int m_num_registration_loops;
      bool m_changed;
      bool m_replace_pass;
      bool m_debug;
      TransitionElementType *m_transition_element_field;

    public:

      PredicateTemplateAdapter(RefinePredicate& predicate_refine,
                           percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0) :
        IAdapter(eMesh, bp, proc_rank_field), m_predicate_refine(predicate_refine), m_useQualityMarking(false),
        m_stage(Stage_None),
        m_iter(0), m_num_registration_loops(0),
        m_changed(false),
        m_replace_pass(false),
        m_debug(false)
      {
        m_transition_element_field = m_eMesh.get_transition_element_field();
        if (!m_transition_element_field)
          throw std::runtime_error("can't use PredicateTemplateAdapter without a transition_element field");
      }

      RefinePredicate& getRefinePredicate() { return m_predicate_refine; }
      void setUseQualityMarking(bool val) { m_useQualityMarking = val; }
      bool getUseQualityMarking() { return m_useQualityMarking; }

      virtual void
      doBreak(int num_registration_loops)
      {
        EXCEPTWATCH;
        // stk::diag::setEnabledTimerMetricsMask(stk::diag::METRICS_CPU_TIME | stk::diag::METRICS_WALL_TIME);
        // static stk::diag::Timer timerAdapt_("Adapt", Refiner::rootTimer());
        // stk::diag::TimeBlock tbTimerAdapt_(timerAdapt_);

        initializeRefine();

        int nq_iter_single = 1, nq_iter_multi = 10;

        m_stage = Stage_1_Mark_NTE;
        doMark(nq_iter_single);

        m_stage = Stage_2_Mark_TE_Parents;
        doMark(nq_iter_multi);

        // done as post-step of previous stage
        //m_stage = Stage_3_Remove_TE;
        //doMark(nq_iter_single);

        //m_stage = Stage_4_UpgradeMarks;
        //doMark(nq_iter_multi);

        //m_stage = Stage_5_Refine;
        //doRefine();
        //doRefine();
      }

      int isTransitionElement(stk::mesh::Entity element)
      {
        int *transition_element = stk::mesh::field_data( *m_transition_element_field , element );
        return transition_element[0];
      }
      void setTransitionElement(stk::mesh::Entity element, int val)
      {
        int *transition_element = stk::mesh::field_data( *m_transition_element_field , element );
        transition_element[0] = val;
      }

      virtual void
      preMark(int iter, int num_registration_loops)
      {
        m_iter = iter;
        m_num_registration_loops = num_registration_loops;

        if (m_iter == 0 && (m_stage == Stage_2_Mark_TE_Parents))
          {
            m_changed = true;
          }
      }

      virtual bool
      postMark(int iter, int num_registration_loops)
      {
        m_eMesh.initializeIdServer();

        bool return_val_break_registration_loop = false;
        if (iter == num_registration_loops - 1 && m_stage == Stage_2_Mark_TE_Parents)  // FIXME
          {
            // remove all transition elements in preparation for upgrading marks
            const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(stk::topology::ELEMENT_RANK + 1u);

            SetOfEntities family_trees_to_be_removed(*m_eMesh.get_bulk_data());
            SetOfEntities te_to_be_removed(*m_eMesh.get_bulk_data());
            SetOfEntities te_to_be_removed_with_ghosts(*m_eMesh.get_bulk_data());

            m_eMesh.get_bulk_data()->modification_begin();

            //mesh::Selector on_locally_owned_part =  ( m_eMesh.get_fem_meta_data()->locally_owned_part() );
            const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.element_rank() );

            for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                stk::mesh::Bucket & bucket = **k ;
                //if (on_locally_owned_part(bucket))
                  {
                    const unsigned num_entity_in_bucket = bucket.size();
                    for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                      {
                        stk::mesh::Entity element = bucket[ientity];
                        bool isGhostElement = m_eMesh.isGhostElement(element);
                        bool isTE = isTransitionElement(element);
                        if (isTE)
                          {
                            percept::MyPairIterRelation te_to_family_tree_relations (m_eMesh, element, FAMILY_TREE_RANK);

                            // look for level 0 only - these are children with no children
                            unsigned child_ft_level_0 = m_eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, element);

                            stk::mesh::Entity family_tree = te_to_family_tree_relations[child_ft_level_0].entity();
                            percept::MyPairIterRelation family_tree_relations (m_eMesh, family_tree,
                                                                               stk::topology::ELEMENT_RANK);

                            stk::mesh::Entity parent = family_tree_relations[FAMILY_TREE_PARENT].entity();
                            if (!m_eMesh.is_valid(parent))
                              {
                                throw std::logic_error("PredicateTemplateAdapter parent == null");
                              }
                            family_trees_to_be_removed.insert(family_tree);
                            if (!isGhostElement)
                              te_to_be_removed.insert(element);
                            te_to_be_removed_with_ghosts.insert(element);
                          }
                      }
                  }
              }

            // FIXME - side sets....
            if (m_debug)
              std::cout << "tmp srk postMark(Stage_2_Mark_TE_Parents) family_trees_to_be_removed.size() = " << family_trees_to_be_removed.size()
                        << " te_to_be_removed.size= " << te_to_be_removed.size() << std::endl;

            //m_nodeRegistry->clear_element_owner_data(children_to_be_removed_with_ghosts);
            //m_nodeRegistry->clear_element_owner_data(te_to_be_removed);

            m_replace_pass = true;
            replaceNodeRegistryOwnership(te_to_be_removed_with_ghosts, m_eMesh.element_rank());
            m_replace_pass = false;

            removeFamilyTrees(family_trees_to_be_removed);
            removeChildElements(te_to_be_removed);

            {
              SetOfEntities elements_to_be_remeshed(*m_eMesh.get_bulk_data());
              const stk::mesh::BucketVector & buckets1 = m_eMesh.get_bulk_data()->buckets( m_eMesh.element_rank() );
              for ( stk::mesh::BucketVector::const_iterator k = buckets1.begin() ; k != buckets1.end() ; ++k )
                {
                  stk::mesh::Bucket & bucket = **k ;
                  //if (on_locally_owned_part(bucket))
                  {
                    const unsigned num_entity_in_bucket = bucket.size();
                    for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                      {
                        stk::mesh::Entity element = bucket[ientity];
                        bool isGhostElement = m_eMesh.isGhostElement(element);
                        if (!isGhostElement && m_eMesh.numChildren(element) == 0)
                          {
                            elements_to_be_remeshed.insert(element);
                          }
                      }
                  }
                }
              int hint = elements_to_be_remeshed.size()*10;
              m_replace_pass = true;
              remeshElements(elements_to_be_remeshed, m_eMesh.element_rank(), hint);
              m_replace_pass = false;


              removeEmptyElements();
              set_active_part();
              fix_side_sets_2();
              //stk::mesh::fixup_ghosted_to_shared_nodes(*m_eMesh.get_bulk_data());
              m_eMesh.get_bulk_data()->modification_end();
              m_nodeRegistry->clear_element_owner_data_phase_2();
              m_eMesh.get_bulk_data()->modification_begin();
              set_active_part();

              m_nodeRegistry->addToExistingPartsNew();
              m_nodeRegistry->prolongate(m_eMesh.get_coordinates_field());
              m_nodeRegistry->prolongateFields();

              //stk::mesh::fixup_ghosted_to_shared_nodes(*m_eMesh.get_bulk_data());
              m_eMesh.get_bulk_data()->modification_end();
            }
            return_val_break_registration_loop = false;
          }

        return return_val_break_registration_loop;
      }

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
        //const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);

        //CoordinatesFieldType* coordField = m_eMesh.get_coordinates_field();

        bool markInfo = (m_predicate_refine(element) & DO_REFINE);
        bool isTE = isTransitionElement(element);

        if (isTE || markInfo)
          {
            if (0 && m_debug)
              std::cout << "tmp srk refineMethodApply(-1 ... " << m_stage << ", " << m_iter << " ) isTE = " << isTE << " markInfo= " << markInfo << " for elem= "  << m_eMesh.identifier(element) << std::endl;
          }

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
            else if (needed_entity_rank == stk::topology::ELEMENT_RANK)
              {
                numSubDimNeededEntities = 1;
              }

            if (m_replace_pass)
              {
                for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                  {
                    (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, false,bucket_topo_data);
                  }
              }
            else if (m_stage == Stage_1_Mark_NTE)
              {
                // only mark NTE
                bool markInfo1 = markInfo;
                if (isTE)
                  {
                    markInfo1 = false;
                  }
                std::vector<bool> markInfoVec(numSubDimNeededEntities, markInfo1);

                if (isTE || markInfo)
                  {
                    if (m_debug)
                      std::cout << "tmp srk refineMethodApply(Stage_1_Mark_NTE ... " << m_stage << ", " << m_iter << " [" << m_num_registration_loops << "] ) isTE = " << isTE << " markInfo= " << markInfo << " for elem= "  << m_eMesh.identifier(element) << std::endl;
                  }

                for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                  {
                    (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, markInfoVec[iSubDimOrd],bucket_topo_data);
                  }
              }
            else if (m_stage == Stage_2_Mark_TE_Parents)
              {
                // mark parents of TE's
                bool markInfo1 = false;
                if (isTE)
                  {
                    bool at_least_one_edge_marked = false;
                    for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                      {
                        bool is_marked = false;
                        bool is_not_marked = false;
                        getEdgeMarkInfo(element, needed_entity_rank, iSubDimOrd, is_marked, is_not_marked);
                        if (is_marked)
                          {
                            at_least_one_edge_marked = true;
                            break;
                          }
                      }
                    stk::mesh::Entity parent = m_eMesh.getParent(element, true);
                    VERIFY_OP_ON(parent, !=, stk::mesh::Entity(), "bad entity");
                    if (m_debug)
                      std::cout << "tmp srk refineMethodApply(Stage_2_Mark_TE_Parents ... " << m_stage << ", " << m_iter
                                << " [" << m_num_registration_loops << "]   ) found transition element = " << m_eMesh.identifier(element)
                                << " at_least_one_edge_marked= " << at_least_one_edge_marked
                                << " markInfo = " << markInfo
                                << " parent = " << m_eMesh.identifier(parent)
                                << std::endl;

                    if (markInfo || at_least_one_edge_marked)
                      {
                        markInfo1 = true;
                        std::vector<bool> markInfoVec(numSubDimNeededEntities, markInfo1);
                        for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                          {
                            (m_nodeRegistry ->* function)(parent, needed_entity_ranks[ineed_ent], iSubDimOrd, markInfoVec[iSubDimOrd],bucket_topo_data);
                          }
                      }
                  }
                else
                  {
                    // upgrade marks
                    unsigned nmarks=0;
                    std::vector<bool> markInfoVec(numSubDimNeededEntities, false);
                    std::vector<bool> markInfoVecAll(numSubDimNeededEntities, true);
                    for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                      {
                        bool is_marked = false;
                        bool is_not_marked = false;
                        getEdgeMarkInfo(element, needed_entity_rank, iSubDimOrd, is_marked, is_not_marked);
                        if (is_marked)
                          {
                            ++nmarks;
                            markInfoVec[iSubDimOrd] = true;
                          }
                      }

                    if (nmarks && m_debug)
                      std::cout << "tmp srk refineMethodApply(Stage_2_Mark_TE_Parents ... " << m_stage << ", " << m_iter 
                                << " [" << m_num_registration_loops << "] ) nmarks= " << nmarks 
                                << " for elem= " << m_eMesh.identifier(element)
                                << std::endl;

                    // FIXME for 3D
                    if (nmarks == 1 || nmarks == 3)
                      {
                        // ok
                      }
                    else if (nmarks == 2)
                      {
                        markInfoVec = markInfoVecAll;
                        if (m_debug)
                          std::cout << "tmp srk found nmarks=2  for elem= " << m_eMesh.identifier(element) << std::endl;
                      }
                    for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                      {
                        (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, markInfoVec[iSubDimOrd],bucket_topo_data);
                      }
                  }
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
