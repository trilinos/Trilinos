#include <exception>
#include <fstream>
#include <set>
#include <typeinfo>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <stk_adapt/Refiner.hpp>

#include <stk_percept/MeshUtil.hpp>

#define DEBUG_UNREF 0
#define DEBUG_UNREF_1 0
#define DEBUG_UNREF_2 0
#define DEBUG_UNREF_3 0

namespace stk {
  namespace adapt {
    using namespace std;
    using namespace percept;

    // ====================================================================================================
    // ====================================================================================================
    // ====================================================================================================
    bool use_idServer = true;

    void Refiner::
    replaceNodeRegistryOwnership(ElementUnrefineCollection& elements_to_delete, stk::mesh::EntityRank rank)
    {
      // assumes clear_element_owner_data has been invoked
      for (unsigned irank=0; irank < m_ranks.size(); irank++)
        {
          if (rank == m_ranks[irank])
            {
              //unsigned num_new_elem_during_remesh = 0;
              vector<NeededEntityType> needed_entity_ranks;
              m_breakPattern[irank]->fillNeededEntities(needed_entity_ranks);

              const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( rank );

              for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
                {
                  {
                    stk::mesh::Bucket & bucket = **k ;
                    const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(bucket);

                    // this can happen for empty elements (created as a pool to be used for creating new refined elems)
                    if (!cell_topo_data) continue;

                    CellTopology cell_topo(cell_topo_data);
                    unsigned elementType = cell_topo.getKey();
                    unsigned bpElementType = m_breakPattern[irank]->getFromTypeKey();
                    if (elementType == bpElementType)
                      {
                        const unsigned num_entity_in_bucket = bucket.size();
                        for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                          {
                            stk::mesh::Entity element = bucket[ientity];
                            if (elements_to_delete.find(element) == elements_to_delete.end())
                              {
                                refineMethodApply(&NodeRegistry::replaceElementOwnership, element, needed_entity_ranks);
                              }
                          }
                      }
                  }
                }
            }
        }
    }


    // get set of children with no nieces that are to be removed; put their parents in a list;
    //   put the children's family trees in a list to be removed
    void Refiner::getChildrenToBeRemoved(ElementUnrefineCollection& elements_to_unref,
                                         SetOfEntities& children_to_be_removed, SetOfEntities& children_to_be_removed_with_ghosts,
                                         SetOfEntities& copied_children_to_be_removed_NOT_USED,
                                         SetOfEntities& family_trees_to_be_removed,
                                         SetOfEntities& parent_elements)
    {
      const unsigned FAMILY_TREE_RANK = stk::mesh::MetaData::ELEMENT_RANK + 1u;

      unsigned nchild_removed = 0;

      for (ElementUnrefineCollection::iterator u_iter = elements_to_unref.begin();
           u_iter != elements_to_unref.end(); ++u_iter)
        {
          stk::mesh::Entity element_p = *u_iter;
          bool isGhostElement = m_eMesh.isGhostElement(element_p);
          bool isParent = m_eMesh.hasFamilyTree(element_p) && m_eMesh.isParentElement(element_p);

          if (isParent)
            {
              throw std::logic_error(" found parent where child expected");
              continue;
            }

          percept::MyPairIterRelation child_to_family_tree_relations (m_eMesh, element_p,FAMILY_TREE_RANK);

          // look for level 0 only - these are children with no children
          unsigned child_ft_level_0 = m_eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, element_p);

          stk::mesh::Entity family_tree = child_to_family_tree_relations[child_ft_level_0].entity();
          percept::MyPairIterRelation family_tree_relations (m_eMesh, family_tree,stk::mesh::MetaData::ELEMENT_RANK);
          if (family_tree_relations.size() == 0)
            {
              throw std::logic_error("Refiner::unrefineTheseElements family_tree_relations.size() == 0");
            }

          stk::mesh::Entity parent = family_tree_relations[FAMILY_TREE_PARENT].entity();
          if (!m_eMesh.is_valid(parent))
            {
              throw std::logic_error("Refiner::unrefineTheseElements parent == null");
            }
          parent_elements.insert(parent);

          family_trees_to_be_removed.insert(family_tree);

          ++nchild_removed;

          if (!isGhostElement)
            {
              children_to_be_removed.insert( element_p );
            }
          children_to_be_removed_with_ghosts.insert( element_p );

        }
#if DEBUG_UNREF
      std::cout << "tmp nchild_removed=: " << nchild_removed << std::endl;
#endif

#if DEBUG_UNREF_1
      std::cout << "tmp children_to_be_removed size = " << children_to_be_removed.size() << std::endl;
      std::cout << "tmp nchild_removed=: " << nchild_removed << std::endl;
#endif

    }

    void Refiner::removeFamilyTrees(SetOfEntities& family_trees_to_be_removed)
    {
      for(SetOfEntities::iterator family_tree_it = family_trees_to_be_removed.begin();
          family_tree_it != family_trees_to_be_removed.end(); ++family_tree_it)
        {
          stk::mesh::Entity family_tree = *family_tree_it;
          if ( ! m_eMesh.get_bulk_data()->destroy_entity( family_tree ) )
            {
              throw std::logic_error("Refiner::removeFamilyTrees couldn't remove element, destroy_entity returned false for family_tree.");
            }
        }
    }

    void Refiner::removeDeletedNodes(NodeSetType& deleted_nodes)
    {
      //std::cout << "P["<< m_eMesh.get_rank() << "] removeDeletedNodes deleted_nodes.size()= " << deleted_nodes.size() << std::endl;
      for(SetOfEntities::iterator node_it = deleted_nodes.begin();
          node_it != deleted_nodes.end(); ++node_it)
        {
          stk::mesh::Entity node = *node_it;
          if ( ! m_eMesh.get_bulk_data()->destroy_entity( node ) )
            {
              //throw std::logic_error("Refiner::removeDeletedNodes couldn't remove node, destroy_entity returned false for node.");
            }
        }
    }


    void Refiner::getSideElemsToBeRemoved(NodeSetType& deleted_nodes,
                                          SetOfEntities& children_to_be_removed, SetOfEntities& side_elem_set_to_be_removed, SetOfEntities& family_trees_to_be_removed, SetOfEntities& parent_side_elements)
    {
      if (getIgnoreSideSets()) return;

      const unsigned FAMILY_TREE_RANK = stk::mesh::MetaData::ELEMENT_RANK + 1u;

      for(SetOfEntities::iterator child_it = children_to_be_removed.begin();
          child_it != children_to_be_removed.end(); ++child_it)
        {
          stk::mesh::Entity child = *child_it;

          if (m_eMesh.hasFamilyTree(child) && !m_eMesh.isChildWithoutNieces(child))
            throw std::logic_error("error 34");

          // add sideset elements to list to be removed (and their family tree info)
          percept::MyPairIterRelation side_relations (m_eMesh, child,m_eMesh.side_rank());
          for (unsigned jside = 0; jside < side_relations.size(); jside++)
            {
              stk::mesh::Entity side_element = side_relations[jside].entity();
              // FIXME err check
              if (m_eMesh.hasFamilyTree(side_element) && !m_eMesh.isChildWithoutNieces(side_element))
                {
                  throw std::logic_error("error 35");
                }
              //if (m_eMesh.hasFamilyTree(*side_element) && m_eMesh.isChildWithoutNieces(*side_element, false))
              if (m_eMesh.hasFamilyTree(side_element))
              {
                unsigned side_elem_child_ft_level_0 = m_eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, side_element);

                percept::MyPairIterRelation side_element_to_family_tree_relations (m_eMesh, side_element,FAMILY_TREE_RANK);
                stk::mesh::Entity family_tree = side_element_to_family_tree_relations[side_elem_child_ft_level_0].entity();
                percept::MyPairIterRelation family_tree_relations (m_eMesh, family_tree,m_eMesh.entity_rank(side_element));
                if (family_tree_relations.size() == 0)
                  {
                    throw std::logic_error("Refiner::getSideElemsToBeRemoved family_tree_relations.size() == 0 [1]");
                  }

                for (unsigned ichild=1; ichild < family_tree_relations.size(); ichild++)
                  {
                    stk::mesh::Entity side_elem_sibling = family_tree_relations[ichild].entity();

                    //std::cout << "tmp unref side element id= " << side_element->identifier() << std::endl;
                    side_elem_set_to_be_removed.insert(side_elem_sibling);
                  }

                stk::mesh::Entity parent = family_tree_relations[FAMILY_TREE_PARENT].entity();
                parent_side_elements.insert(parent);
                family_trees_to_be_removed.insert(family_tree);
              }
            }
        }
    }

    void Refiner::getSideParentsToBeRemeshed(SetOfEntities& parents_to_be_remeshed, SetOfEntities& parent_side_elements)
    {
      if (getIgnoreSideSets()) return;

      //const unsigned FAMILY_TREE_RANK = stk::mesh::MetaData::ELEMENT_RANK + 1u;

      for(SetOfEntities::iterator parent_it = parents_to_be_remeshed.begin();
          parent_it != parents_to_be_remeshed.end(); ++parent_it)
        {
          stk::mesh::Entity parent = *parent_it;

          //           if (m_eMesh.hasFamilyTree(parent) && !m_eMesh.isChildWithoutNieces(child))
          //             throw std::logic_error("error 34");

          // add sideset elements to list to be removed (and their family tree info)
          percept::MyPairIterRelation side_relations (m_eMesh, parent, m_eMesh.side_rank());
          for (unsigned jside = 0; jside < side_relations.size(); jside++)
            {
              stk::mesh::Entity side_element = side_relations[jside].entity();

              if (0 == m_eMesh.numChildren(side_element))
              {
                parent_side_elements.insert(side_element);
                //family_trees_to_be_removed.insert(family_tree);
              }
            }
        }
    }



    void Refiner::removeChildElements(SetOfEntities& children_to_be_removed, ElementUnrefineCollection *elements_to_unref_0)
    {
      const unsigned FAMILY_TREE_RANK = stk::mesh::MetaData::ELEMENT_RANK + 1u;
      for(SetOfEntities::iterator child_it = children_to_be_removed.begin();
          child_it != children_to_be_removed.end(); ++child_it)
        {
          stk::mesh::Entity child = *child_it;

          if (m_eMesh.isGhostElement(child))
          {
              throw std::logic_error("Refiner::removeChildElements couldn't remove element, Ghost is true.");
          }
          if ( ! m_eMesh.get_bulk_data()->destroy_entity( child ) )
            {
              CellTopology cell_topo(m_eMesh.get_cell_topology(child));

              //const percept::MyPairIterRelation elem_relations ( child->relations(child->entity_rank(m_eMesh,)+1);
              const percept::MyPairIterRelation child_to_ft_relations (m_eMesh, child,FAMILY_TREE_RANK);
#if DEBUG_UNREF
              std::cout << "tmp Refiner::removeChildElements couldn't remove element  cell= " << cell_topo.getName() << std::endl;
              std::cout << "tmp child_to_ft_relations.size() = " << child_to_ft_relations.size() << std::endl;
              //std::cout << "tmp ft_id loc, outerloop= " << child_to_family_tree_relations[0].entity()->identifier() << " " << family_tree_id << std::endl;

              m_eMesh.print_entity(std::cout, child);
#endif

              throw std::logic_error("Refiner::removeChildElements couldn't remove element, destroy_entity returned false.");
            }
          if (elements_to_unref_0)
            {
              elements_to_unref_0->erase(child);
            }
        }
    }

    void Refiner::removeSideElements(SetOfEntities& side_elem_set_to_be_removed, SetOfEntities& elements_to_be_deleted)
    {
      for(SetOfEntities::iterator side_elem_it = side_elem_set_to_be_removed.begin();
          side_elem_it != side_elem_set_to_be_removed.end(); ++side_elem_it)
        {
          stk::mesh::Entity side_elem = *side_elem_it;
          bool del =  m_eMesh.get_bulk_data()->destroy_entity( side_elem );
          if ( ! del )
            {
              std::cout << "Refiner::removeSideElements couldn't remove side element, destroy_entity returned false." << side_elem << std::endl;
              del =  m_eMesh.get_bulk_data()->destroy_entity( side_elem );
              if (1)
                {
                  percept::MyPairIterRelation rels (m_eMesh, side_elem,stk::mesh::MetaData::ELEMENT_RANK);
                  for (unsigned irels=0; irels < rels.size(); irels++)
                    {
                      bool in_del = (elements_to_be_deleted.find(rels[irels].entity()) != elements_to_be_deleted.end());
                      std::cout << "Refiner::removeSideElements found element: in del list= " << in_del << std::endl;
                    }
                }
              throw std::logic_error("Refiner::removeSideElements couldn't remove side element, destroy_entity returned false.");
            }
        }
    }


    // for the given set of parent elements, remesh the interior after the children have been removed and any deleted nodes have been removed from the NodeRegistry
    void Refiner::remesh(SetOfEntities& parent_elements)
    {
      //
      // FIXME refactor to a generic function operating on a collection of elements; incorporate with the doBreak() calls above
      //
      // remesh

#if DEBUG_UNREF
      //std::cout << "tmp remesh:: parent_elements.size() [elements to be remeshed] = " << parent_elements.size() << std::endl;
#endif

      // FIXME for performance
      //static NewSubEntityNodesType s_new_sub_entity_nodes(stk::percept::EntityRankEnd);
      NewSubEntityNodesType s_new_sub_entity_nodes(stk::percept::EntityRankEnd);

      NewSubEntityNodesType& new_sub_entity_nodes = s_new_sub_entity_nodes;

      for (unsigned irank=0; irank < m_ranks.size(); irank++)
        {
          unsigned num_new_elem_during_remesh = 0;
          vector<NeededEntityType> needed_entity_ranks;
          m_breakPattern[irank]->fillNeededEntities(needed_entity_ranks);

          vector<stk::mesh::Entity> new_elements;

          // count num new elements needed on this proc (served by UniformRefinerPattern)
          unsigned num_elem_not_ghost = 0u;

          for (ElementUnrefineCollection::iterator p_iter = parent_elements.begin();
               p_iter != parent_elements.end(); ++p_iter)
            {
              stk::mesh::Entity parent = *p_iter;

              stk::mesh::Entity element = parent;

              const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
              CellTopology cell_topo(cell_topo_data);
              unsigned elementType = cell_topo.getKey();
              unsigned bpElementType = m_breakPattern[irank]->getFromTypeKey();
              if (elementType == bpElementType)
                {
                  if (!m_eMesh.isGhostElement(element))
                    {
                      ++num_elem_not_ghost;
                    }
                }
            }

          unsigned num_elem_needed = num_elem_not_ghost * m_breakPattern[irank]->getNumNewElemPerElem();

          //std::cout << "tmp remesh::rank[" << irank << "] = " << m_ranks[irank] << " , num_elem_needed= " << num_elem_needed << std::endl;

#if DEBUG_UNREF
          //std::cout << "tmp remesh::rank[irank], num_elem_needed= " << m_ranks[irank] << " " << num_elem_needed << std::endl;
#endif


          // create new entities on this proc
          new_elements.resize(0);
          if (m_eMesh.getEntityPool().size())
            {
              if (!m_eMesh.getEntitiesFromPool(m_ranks[irank], num_elem_needed, new_elements))
                {
                  throw std::logic_error("entity pool deplenished");
                }
            }
          else if (use_idServer)
            {
              new_elements.resize(num_elem_needed);
              stk::mesh::PartVector empty ;

              for (unsigned ii=0; ii < num_elem_needed; ii++)
                {
                  stk::mesh::EntityId new_id = m_eMesh.getNextId(m_ranks[irank]);
                  new_elements[ii] = m_eMesh.get_bulk_data()->declare_entity( m_ranks[irank], new_id, empty );
                }
            }
          else
            {
              m_eMesh.createEntities( m_ranks[irank], num_elem_needed, new_elements);
            }
          vector<stk::mesh::Entity>::iterator element_pool_it = new_elements.begin();

          // FIXME - we could directly call this with a refactor to change elementColors passed in here as a generic collection + checking for element Type
          //
          //createElementsAndNodesAndConnectLocal(m_ranks[irank], m_breakPattern[irank], elementColors, needed_entity_ranks, new_elements);

          for (ElementUnrefineCollection::iterator p_iter = parent_elements.begin();
               p_iter != parent_elements.end(); ++p_iter)
            {
              stk::mesh::Entity parent_p = *p_iter;
              const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(parent_p);
              CellTopology cell_topo(cell_topo_data);
              unsigned elementType = cell_topo.getKey();
              unsigned bpElementType = m_breakPattern[irank]->getFromTypeKey();
              if (elementType == bpElementType)
                {
                  stk::mesh::Entity parent = parent_p;

                  if (!m_eMesh.isGhostElement(parent))
                    {
#if DEBUG_UNREF_2
                      //std::cout << "P["<< m_eMesh.get_rank() << "] eMesh.owner_rank(parent) = " << eMesh.owner_rank(parent) << std::endl;
                      std::cout << "tmp Parent to be remeshed = ";
                      m_eMesh.print_entity(std::cout, parent);
#endif
                      if (createNewNeededNodeIds(cell_topo_data, parent, needed_entity_ranks, new_sub_entity_nodes))
                        {
                          //std::cout << "typeid= " << typeid(*breakPattern).name() << std::endl;
                          throw std::logic_error("remesh:: createNewNeededNodeIds failed");
                        }
                      vector<stk::mesh::Entity>::iterator element_pool_it_b4 = element_pool_it;
                      m_breakPattern[irank]->createNewElements(m_eMesh, *m_nodeRegistry, parent, new_sub_entity_nodes, element_pool_it, m_proc_rank_field);
                      vector<stk::mesh::Entity>::iterator element_pool_it_af = element_pool_it;
                      num_new_elem_during_remesh += (element_pool_it_af - element_pool_it_b4);
                    }
                }
            }
#if DEBUG_UNREF
          //std::cout << "tmp remesh:: nchild elements during remesh for rank[irank] = " << m_ranks[irank] << " " << num_new_elem_during_remesh << std::endl;
#endif

        } // irank

    }


    void
    Refiner::
    unrefineAll()
    {
      ElementUnrefineCollection elements_to_unref(*m_eMesh.get_bulk_data());

      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( stk::mesh::MetaData::ELEMENT_RANK );

      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          //if (removePartSelector(**k))
          {
            stk::mesh::Bucket & bucket = **k ;

            const unsigned num_entity_in_bucket = bucket.size();
            for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
              {
                stk::mesh::Entity element = bucket[ientity];
                // FIXME
                // skip elements that are already a parent (if there's no family tree yet, it's not a parent, so avoid throwing an error in isParentElement)
                const bool check_for_family_tree = false;
                bool isParent = m_eMesh.isParentElement(element, check_for_family_tree);

                if (isParent)
                  continue;

                const percept::MyPairIterRelation elem_nodes (m_eMesh, element,stk::mesh::MetaData::NODE_RANK);

                if (elem_nodes.size() && (m_eMesh.hasFamilyTree(element) && m_eMesh.isChildWithoutNieces(element, false) ) )
                  {
                    bool elementIsGhost = m_eMesh.isGhostElement(element);

                    if (!elementIsGhost)
                      {
                        elements_to_unref.insert(element);
                      }
                  }
              }
          }
        }
      unrefineTheseElements(elements_to_unref);
    }

    void Refiner::
    filterUnrefSet(ElementUnrefineCollection& elements_to_unref)
    {
      int print_filter_info = DEBUG_UNREF_1;
      if (print_filter_info)  std::cout << "P["<< m_eMesh.get_rank() << "] filterUnrefSet: initial set size = " << elements_to_unref.size() << std::endl;

      const unsigned FAMILY_TREE_RANK = stk::mesh::MetaData::ELEMENT_RANK + 1u;
      ElementUnrefineCollection elements_to_unref_copy(*m_eMesh.get_bulk_data());

      typedef std::set<stk::mesh::Entity, stk::mesh::EntityLess> SetOfEntities;

      int num_is_parent = 0;
      int num_elem_nodes_0 = 0;
      int num_has_nieces = 0;
      int num_all_siblings_in_unref_set = 0;

      for (ElementUnrefineCollection::iterator u_iter = elements_to_unref.begin();
           u_iter != elements_to_unref.end(); ++u_iter)
        {
          stk::mesh::Entity element = *u_iter;

          const bool check_for_family_tree = false;
          bool isParent = m_eMesh.isParentElement(element, check_for_family_tree);

          if (isParent)
            {
              ++num_is_parent;
              continue;
            }

          const percept::MyPairIterRelation elem_nodes (m_eMesh, element,stk::mesh::MetaData::NODE_RANK);
          int elem_nodes_size = elem_nodes.size();
          bool has_no_nieces = m_eMesh.hasFamilyTree(element) && m_eMesh.isChildWithoutNieces(element, false);

          if (!elem_nodes_size) num_elem_nodes_0++;
          if (!has_no_nieces) {
            num_has_nieces++;
          }
          if (elem_nodes_size && has_no_nieces)
            {
              percept::MyPairIterRelation child_to_family_tree_relations (m_eMesh, element,  FAMILY_TREE_RANK);

              // look for level 0 only - these are children with no children
              unsigned child_ft_level_0 = m_eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, element);

              stk::mesh::Entity family_tree = child_to_family_tree_relations[child_ft_level_0].entity();
              percept::MyPairIterRelation family_tree_relations (m_eMesh, family_tree,  stk::mesh::MetaData::ELEMENT_RANK);
              if (family_tree_relations.size() == 0)
                {
                  throw std::logic_error("Refiner::filterUnrefSet family_tree_relations.size() == 0");
                }

              SetOfEntities side_elem_set(*m_eMesh.get_bulk_data());

              bool all_siblings_in_unref_set = true;
              bool all_side_sets_ok = true;

              // siblings
              for (unsigned ichild=1; ichild < family_tree_relations.size(); ichild++)
                {
                  stk::mesh::Entity child = family_tree_relations[ichild].entity();
                  if (m_eMesh.isParentElement(child))
                    {
                      throw std::logic_error("Refiner::filterUnrefSet isParentElement not expected");
                    }
                  bool in_unref_set = elements_to_unref.find(child) != elements_to_unref.end();
                  if (!in_unref_set)
                    {
                      all_siblings_in_unref_set = false;
                      break;
                    }

                  {
                    percept::MyPairIterRelation side_relations (m_eMesh, child,  m_eMesh.side_rank());
                    for (unsigned jside = 0; jside < side_relations.size(); jside++)
                      {
                        stk::mesh::Entity side_element = side_relations[jside].entity();
                        side_elem_set.insert(side_element);

                        if (0)
                        {
                          percept::MyPairIterRelation side_elem_to_family_tree_relations (m_eMesh, side_element,  FAMILY_TREE_RANK);
                            stk::mesh::Entity side_elem_family_tree_0 = side_elem_to_family_tree_relations[0].entity();
                            std::cout << "side_elem_family_tree_0= " << m_eMesh.identifier(side_elem_family_tree_0) << std::endl;
                        }

                        // FIXME if (!m_eMesh.hasFamilyTree(*side_element) || !m_eMesh.isChildWithoutNieces(*side_element))
                        if (m_eMesh.hasFamilyTree(side_element) && !m_eMesh.isChildWithoutNieces(side_element))
                          {
                            //std::cout << "error 35" << std::endl;
                            all_side_sets_ok=false;
                            break;
                          }
                      }
                  }
                }

              if (all_side_sets_ok)
                {
                  for (unsigned ichild=1; ichild < family_tree_relations.size(); ichild++)
                    {
                      stk::mesh::Entity child = family_tree_relations[ichild].entity();

                      {
                        percept::MyPairIterRelation side_relations (m_eMesh, child,  m_eMesh.side_rank());
                        for (unsigned jside = 0; jside < side_relations.size(); jside++)
                          {
                            stk::mesh::Entity side_element = side_relations[jside].entity();
                            //side_elem_set.insert(side_element);

                            if (!m_eMesh.hasFamilyTree(side_element)) continue;

                            percept::MyPairIterRelation side_elem_to_family_tree_relations (m_eMesh, side_element,  FAMILY_TREE_RANK);

                            // look for level 0 only - these are children with no children
                            unsigned side_elem_ft_level_0 = m_eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, side_element);

                            stk::mesh::Entity side_elem_family_tree = side_elem_to_family_tree_relations[side_elem_ft_level_0].entity();
                            percept::MyPairIterRelation side_elem_family_tree_relations (m_eMesh, side_elem_family_tree,  m_eMesh.side_rank());
                            for (unsigned ise_child=1; ise_child < side_elem_family_tree_relations.size(); ise_child++)
                              {
                                stk::mesh::Entity se_sibling = side_elem_family_tree_relations[ise_child].entity();
                                bool in_set = (side_elem_set.find(se_sibling) != side_elem_set.end());
                                if (!in_set)
                                  {
                                    all_side_sets_ok = false;
                                    break;
                                  }
                              }
                          }
                      }
                    }
                }

              all_side_sets_ok=true;
              //all_siblings_in_unref_set=true;
              if (!all_siblings_in_unref_set) {
                ++num_all_siblings_in_unref_set;
              }
              if (all_siblings_in_unref_set && all_side_sets_ok)
                {
                  for (unsigned ichild=1; ichild < family_tree_relations.size(); ichild++)
                    {
                      stk::mesh::Entity child = family_tree_relations[ichild].entity();
                      elements_to_unref_copy.insert(child);
                    }
                }
            }
        }
      if (print_filter_info) std::cout << "tmp filterUnrefSet::elements_to_unref.size = " << elements_to_unref.size()
                                       << " filtered size= " << elements_to_unref_copy.size()
                                       << " num_all_siblings_in_unref_set= " << num_all_siblings_in_unref_set
                                       << " num_has_nieces= " << num_has_nieces
                                       << " num_elem_nodes_0= " << num_elem_nodes_0
                                       << " num_is_parent= " << num_is_parent
                                       << std::endl;
      elements_to_unref = elements_to_unref_copy;

    }

    void Refiner::
    getKeptNodes(NodeSetType& kept_nodes, ElementUnrefineCollection& elements_to_unref)
    {
      // mark kept nodes
      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( stk::mesh::MetaData::ELEMENT_RANK );

      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          {
            stk::mesh::Bucket & bucket = **k ;

            const unsigned num_entity_in_bucket = bucket.size();
            for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
              {
                stk::mesh::Entity element = bucket[ientity];

                const percept::MyPairIterRelation elem_nodes (m_eMesh, element,  stk::mesh::MetaData::NODE_RANK);

                bool tt = (elem_nodes.size() && m_eMesh.isLeafElement(element));
                tt =  !m_eMesh.hasFamilyTree(element) || m_eMesh.numChildren(element) == 0 || !m_eMesh.hasGrandChildren(element, true);

                if (tt)
                  {
                    bool in_unref_set = elements_to_unref.find( element ) != elements_to_unref.end();
                    //bool isGhostElement = m_eMesh.isGhostElement(element);
                    //if (!in_unref_set && !isGhostElement)
                    if (!in_unref_set)
                      {
                        for (unsigned inode=0; inode < elem_nodes.size(); inode++)
                          {
                            stk::mesh::Entity node = elem_nodes[inode].entity();
                            kept_nodes.insert(node);
#if DEBUG_UNREF_2
                            std::cout << "tmp kept node: " << *node << " ";
                            m_eMesh.print_entity(std::cout, *node);
#endif

                          }
                      }
                  }
              }
          }
        }

#if DEBUG_UNREF
      std::cout << "getKeptNodes: elements_to_unref.size= " << elements_to_unref.size() << " kept_nodes.size= " << kept_nodes.size() << std::endl;
#endif


    }

    void Refiner::
    getDeletedNodes(NodeSetType& deleted_nodes, const NodeSetType& kept_nodes, ElementUnrefineCollection& elements_to_unref)
    {
      // mark deleted nodes (nodes in list of elements_to_unref
      unsigned count_candidate_elements = 0;
      unsigned count_in_kept_nodes_set = 0;
      unsigned count_non_candidate_elements = 0;
      unsigned count_not_in_kept_nodes_set = 0;
      for (ElementUnrefineCollection::iterator u_iter = elements_to_unref.begin();
           u_iter != elements_to_unref.end(); ++u_iter)
        {
          stk::mesh::Entity element = *u_iter;

          const percept::MyPairIterRelation elem_nodes (m_eMesh, element,  stk::mesh::MetaData::NODE_RANK);

          if (!m_eMesh.isGhostElement(element) && ((elem_nodes.size() && (!m_eMesh.hasFamilyTree(element) || !m_eMesh.isParentElement(element)))))
            //if (m_eMesh.isChildElement(*element) && !m_eMesh.isGhostElement(*element))
            {
              ++count_candidate_elements;

              for (unsigned inode=0; inode < elem_nodes.size(); inode++)
                {
                  stk::mesh::Entity node = elem_nodes[inode].entity();
                  bool in_kept_nodes_set = kept_nodes.find( node ) != kept_nodes.end();
                  if (in_kept_nodes_set) count_in_kept_nodes_set++;
                  if (!in_kept_nodes_set)
                    {
                      ++count_not_in_kept_nodes_set;
                      deleted_nodes.insert(node);
#if DEBUG_UNREF_2
                      std::cout << "tmp deleted node: " << *node << " ";
                      m_eMesh.print_entity(std::cout, *node);
#endif
                    }
                }
            }
          else
            {
              ++count_non_candidate_elements;
            }
        }
      if (DEBUG_UNREF_1) std::cout << "tmp kept_nodes size= " << kept_nodes.size() << " deleted_nodes size= " << deleted_nodes.size()
                                   << " count_candidate_elements= " << count_candidate_elements
                                   << " count_in_kept_nodes_set= " << count_in_kept_nodes_set
                                   << " count_non_candidate_elements= " << count_non_candidate_elements
                                   << " count_not_in_kept_nodes_set= " << count_not_in_kept_nodes_set
                                   << std::endl;

    }

    void
    Refiner::
    unrefineTheseElements(ElementUnrefineCollection& elements_to_unref)
    {
      bool doPass2=true;
      if (doPass2)
        {
          unrefinePass2(elements_to_unref);
          return;
        }

      if (m_alwaysInitNodeRegistry)
        {
          throw std::logic_error("Refiner::unrefineTheseElements: to use urefinement, you must have setAlwaysInitializeNodeRegistry(false)");
        }
      VERIFY_OP_ON(m_eMesh.getEntityPool().size(), ==, 0, "hmmm");
      
      //m_nodeRegistry->checkDB("unrefine start");

      m_eMesh.get_bulk_data()->modification_begin();

        {
          // mark nodes
          // set<> kept_nodes
          // set<> deleted_nodes

          /* Algorithm Option 1:
             (requires parallel comm)
             1. reset/initialize node registry
             2. loop over elements;
             verify has no children
             get the parent;

             2. loop over elements;
             verify has no children
             get the parent;
             loop over parent sub-dim entities;
             for
          */

          /* Option 2,3: purely local

          0. keep NodeRegistry DB always (could be compressed later with a single int holding element id + subDim +iord)
          1. mark all nodes belonging to non-deleted leaf elements as KEPT
          2. foreach elements_to_unref;
          verify has no children
          get the parent;
          mark nodes of deleted elements as DELETED (but if KEPT is set, then don't delete)
          3. foreach elements_to_unref
          a. delete parent's children (once only of course)
          b. delete any DELETED nodes (if not done automagically by stk_mesh in step 3a)
          c. for sanity, delete DELETED nodes from NodeRegistry DB

          [option 2]:
          4. re-refine using currently marked edges using templates

          [option 3]:
          4. re-refine using currently marked edges using a local triangulation of each parent and its marked edges
          [option 1]: use a local Delaunay method
          [option 2]: project to a master element, use some form of template approach

          5. rebuild when there's a load-balance?
          [option 1]: build a surrogate DB in stk_mesh using face/edge and attributes
          [option 2]: use existing db, add parallel comm

          */

          //std::cout << "tmp elements_to_unref.size() = " << elements_to_unref.size() << std::endl;

          NodeSetType kept_nodes(*m_eMesh.get_bulk_data());
          NodeSetType deleted_nodes(*m_eMesh.get_bulk_data());

          // filter unref set
          filterUnrefSet(elements_to_unref);

          // get kept nodes
          getKeptNodes(kept_nodes, elements_to_unref);

          // get deleted nodes
          getDeletedNodes(deleted_nodes, kept_nodes, elements_to_unref);

          // nothing to be done
          if (0 && deleted_nodes.size() == 0)
            {
              m_eMesh.get_bulk_data()->modification_end();
              std::cout << "Refiner::unrefineTheseElements: deleted_nodes size is 0, nothing to be done, early return." << std::endl;
              return;
            }

          // remove deleted nodes and their associated sub-dim entities from NodeRegistry's db
          //   save kept_nodes_orig on NodeRegistry
          SubDimCellToDataMap to_save;
          NodeSetType kept_nodes_orig(*m_eMesh.get_bulk_data());
          m_nodeRegistry->cleanDeletedNodes(deleted_nodes, kept_nodes_orig, to_save);

          // remove elements to be unrefined
          ElementUnrefineCollection copied_children_to_be_removed = elements_to_unref;

#if DEBUG_UNREF || DEBUG_UNREF_1
          std::cout << "tmp copied_children_to_be_removed.size() [= num elements to be urefined that are children and !ghosts]= " << copied_children_to_be_removed.size() << std::endl;
#endif

          typedef std::set<stk::mesh::Entity, stk::mesh::EntityLess> SetOfEntities;

          SetOfEntities family_trees_to_be_removed(*m_eMesh.get_bulk_data());
          SetOfEntities children_to_be_removed(*m_eMesh.get_bulk_data());
          SetOfEntities children_to_be_removed_with_ghosts(*m_eMesh.get_bulk_data());
          SetOfEntities parent_elements(*m_eMesh.get_bulk_data());

          // set to hold sideset elements to be removed
          SetOfEntities side_elem_set_to_be_removed(*m_eMesh.get_bulk_data());
          SetOfEntities side_elem_family_trees_to_be_removed(*m_eMesh.get_bulk_data());
          SetOfEntities parent_side_elements(*m_eMesh.get_bulk_data());

          // find all elements to be removed from the filtered list and build additional list with ghost elements
          getChildrenToBeRemoved(elements_to_unref,
                                 children_to_be_removed, children_to_be_removed_with_ghosts, copied_children_to_be_removed,
                                 family_trees_to_be_removed,
                                 parent_elements);
          if (DEBUG_UNREF_1) std::cout << "children_to_be_removed.size= " << children_to_be_removed.size() << " elements_to_unref.size= " << elements_to_unref.size() << std::endl;

          // NodeRegistry must be kept in sync
          m_nodeRegistry->clear_element_owner_data(children_to_be_removed_with_ghosts);

          // get the corresponding side elements from the children->side relations
          // FIXME - this is to be deprecated, but leaving it in for now
#if 0
          getSideElemsToBeRemoved(deleted_nodes,
                                  children_to_be_removed, side_elem_set_to_be_removed, side_elem_family_trees_to_be_removed,
                                  parent_side_elements);
#endif

          // first have to delete the family tree (higher ranks have to be deleted first)
          removeFamilyTrees(family_trees_to_be_removed);

          // remove children
          removeChildElements(children_to_be_removed);

          // FIXME - this is to be deprecated, but leaving it in for now
#if 0
          // for sideset elements, remove their family trees first
          removeFamilyTrees(side_elem_family_trees_to_be_removed);

          // remove sideset elements
          removeSideElements(side_elem_set_to_be_removed, children_to_be_removed);
#endif

          // reconnect and remove any dangling side elements
          bool allow_not_found = true;
          fix_side_sets_2(allow_not_found);

          getSideParentsToBeRemeshed(parent_elements, parent_side_elements);

          // remesh the holes left by removing child elems (quad/hex hanging node doesn't need this)
          if (m_needsRemesh)
            {
              remesh(parent_elements);
              remesh(parent_side_elements);
            }

#if CHECK_DEBUG
          check_db("after unrefineTheseElements, b4 mod end");
#endif

          // remove any elements that are empty (these can exist when doing local refinement)
          removeEmptyElements();

          removeDeletedNodes(deleted_nodes);

          set_active_part();

          fix_side_sets_2();

        } // ilevel

      //if (1)  std::cout << "P["<< m_eMesh.get_rank() << "] unrefineTheseElements modification_end start..." << std::endl;
      m_eMesh.get_bulk_data()->modification_end();
      //if (1)  std::cout << "P["<< m_eMesh.get_rank() << "] unrefineTheseElements modification_end ...end" << std::endl;

#if CHECK_DEBUG
      check_db("before clear_element_owner_data_phase_2");
#endif
      //m_nodeRegistry->checkDB("before clear_element_owner_data_phase_2");

      m_nodeRegistry->clear_element_owner_data_phase_2();

      //m_nodeRegistry->checkDB("after clear_element_owner_data_phase_2");
      m_eMesh.get_bulk_data()->modification_begin();
      set_active_part();
      m_eMesh.destroyEntityPool();
      m_eMesh.get_bulk_data()->modification_end();

      //check_sidesets_2(" unrefineTheseElements:: end");

#if CHECK_DEBUG
      check_db("after unrefineTheseElements");
#endif

    }

    // ======================================================================================================================================================
    // ======================================================================================================================================================
    // ======================================================================================================================================================

    // Pass2

    void Refiner::
    remeshRecurse(stk::mesh::Entity element)
    {
      // Note: we rely on the remeshing returning no children if the element isn't marked for refinement
      SetOfEntities element_set_of_one(&element, (&element)+1, *m_eMesh.get_bulk_data());
      remesh(element_set_of_one);

      if (m_eMesh.hasFamilyTree(element))
        {
          std::vector<stk::mesh::Entity> children;
          m_eMesh.getChildren(element, children, true, false);
          for (unsigned ichild=0; ichild < children.size(); ichild++)
            {
              remeshRecurse(children[ichild]);
            }
        }
    }

    void Refiner::
    remeshElements(SetOfEntities& rootElements, stk::mesh::EntityRank rank, int pool_size_hint)
    {
      bool use_pool = true && !use_idServer;
      if (use_pool)
      {
        int pool_size = rootElements.size()*50;
        if (pool_size_hint)
          {
            pool_size = std::max(pool_size, pool_size_hint);
          }
        pool_size *= 2; // safety factor
        m_eMesh.initializeEntityPool(rank, pool_size);
      }
      if (use_idServer)
        {
          m_eMesh.resetIdServer();
          m_eMesh.getNextId(m_eMesh.element_rank()+1u);
          m_eMesh.getNextId(m_eMesh.element_rank());
          m_eMesh.getNextId(m_eMesh.side_rank());
          m_eMesh.getNextId(m_eMesh.node_rank());
        }
      for (SetOfEntities::iterator elIter = rootElements.begin(); elIter != rootElements.end(); ++elIter)
        {
          stk::mesh::Entity element = *elIter;
          //std::cout << "element= " << m_eMesh.identifier(element) << std::endl;
          VERIFY_OP_ON(m_eMesh.numChildren(element), ==, 0, "hmmm");

          remeshRecurse(element);
        }

      m_nodeRegistry->clear_element_owner_data_phase_2(false, false);
      SetOfEntities emptySet(*m_eMesh.get_bulk_data());
      replaceNodeRegistryOwnership(emptySet, rank);
      m_nodeRegistry->clear_element_owner_data_phase_2(true, false);

      set_active_part();
      if (use_pool)
        {
          m_eMesh.destroyEntityPool();
        }
      if (use_idServer)
        {
          m_eMesh.resetIdServer();
        }
    }

    // get a breadth-first list of descendants - if only_leaves is set, then
    //   only insert elements with no children
    //
    // if elements_to_unref is non-null, then return true if all descendants are in the elements_to_unref list, false otherwise
    //   and allow a short-circuit return: Note: this is an efficiency improvement over direct usages of allDescendants
    //

    bool Refiner::allDescendants(stk::mesh::Entity element, SetOfEntities& descendants, unsigned& nlevels, bool only_leaves, ElementUnrefineCollection *elements_to_unref)
    {
      unsigned nlev = nlevels + 1;
      nlevels=0;
      if (m_eMesh.hasFamilyTree(element))
        {
          std::vector<stk::mesh::Entity> children;
          m_eMesh.getChildren(element, children, true, false);
          for (unsigned ichild=0; ichild < children.size(); ichild++)
            {
              if (!only_leaves || m_eMesh.numChildren(children[ichild]) == 0)
                {
                  descendants.insert(children[ichild]);
                  // short-circuit
                  if (elements_to_unref)
                    {
                      if (elements_to_unref->find(children[ichild]) == elements_to_unref->end())
                        return false;
                    }
                }
            }
          for (unsigned ichild=0; ichild < children.size(); ichild++)
            {
              bool all_in_list = allDescendants(children[ichild], descendants, nlev, only_leaves, elements_to_unref);
              if (elements_to_unref && not all_in_list)
                return false;
              nlevels = std::max(nlevels, nlev);
            }
        }
      return true;
    }

    void Refiner::
    filterRecurse(stk::mesh::Entity element, ElementUnrefineCollection& rootElements, ElementUnrefineCollection& elements_to_unref)
    {

      if (1)
        {
          SetOfEntities allD(*m_eMesh.get_bulk_data());
          bool only_leaves=false;
          unsigned nlevels=0;
          //allDescendants(element, allD, nlevels, only_leaves);
          //bool allIn = true;
          bool allIn = (DEBUG_UNREF_3 ? allDescendants(element, allD, nlevels, only_leaves) :
                        allDescendants(element, allD, nlevels, only_leaves, &elements_to_unref) );
          if (DEBUG_UNREF_3)
            {
              for (SetOfEntities::iterator it=allD.begin(); it != allD.end(); ++it)
                {
                  if (elements_to_unref.find(*it) == elements_to_unref.end())
                    {
                      allIn = false;
                      break;
                    }
                }
            }
          if (allD.size() && allIn)
            {
              rootElements.insert(element);
              return;
            }
        }

      if (m_eMesh.hasFamilyTree(element))
        {
          std::vector<stk::mesh::Entity> children;
          m_eMesh.getChildren(element, children, true, false);
          for (unsigned ichild=0; ichild < children.size(); ichild++)
            {
              filterRecurse(children[ichild], rootElements, elements_to_unref);
            }
        }
    }

    void Refiner::
    filterUnrefSetPass2(ElementUnrefineCollection& elements_to_unref,   SetOfEntities& rootElements)
    {
      int print_filter_info = DEBUG_UNREF_1;
      if (print_filter_info)  std::cout << "\n\nP["<< m_eMesh.get_rank() << "] filterUnrefSetPass2: initial elements_to_unref size = " << elements_to_unref.size() << std::endl;

      ElementUnrefineCollection new_set(*m_eMesh.get_bulk_data());
      ElementUnrefineCollection new_root(*m_eMesh.get_bulk_data());

        {
          const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( stk::mesh::MetaData::ELEMENT_RANK );
          for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_entity_in_bucket = bucket.size();
              for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                {
                  stk::mesh::Entity element = bucket[ientity];
                  bool is_root = !m_eMesh.hasFamilyTree(element) || !m_eMesh.isChildElement(element,true);

                  if (is_root)
                    {
                      filterRecurse(element, new_root, elements_to_unref);
                    }
                }
            }

          for (ElementUnrefineCollection::iterator rit=new_root.begin(); rit != new_root.end(); rit++)
            {
              stk::mesh::Entity element = *rit;

              SetOfEntities allD(*m_eMesh.get_bulk_data());
              bool only_leaves=false;
              unsigned nlevels=0;
              //allDescendants(element, allD, nlevels, only_leaves);
              //bool allIn = true;
              bool allIn = (DEBUG_UNREF_3 ? allDescendants(element, allD, nlevels, only_leaves, &elements_to_unref) :
                            allDescendants(element, allD, nlevels, only_leaves, 0) );
              if (DEBUG_UNREF_3)
                {
                  for (SetOfEntities::iterator it=allD.begin(); it != allD.end(); ++it)
                    {
                      if (elements_to_unref.find(*it) == elements_to_unref.end())
                        {
                          allIn = false;
                          break;
                        }
                    }
                }
              VERIFY_OP_ON((allIn && allD.size()), == , true, "hmm");

              if (allIn && allD.size())
                {
                  //if (print_filter_info && (!in_set && !is_root) ) std::cout << "element= " << m_eMesh.identifier(element) << " in_set= " << in_set << " is_root= " << is_root << std::endl;

                  new_root.insert(element);

                  for (SetOfEntities::iterator it=allD.begin(); it != allD.end(); ++it)
                    {
                      new_set.insert(*it);
                    }
                }
            }
          rootElements = new_root;
          elements_to_unref = new_set;
        }

      if (print_filter_info)  std::cout << "\n\nP["<< m_eMesh.get_rank() << "] filterUnrefSetPass2: final elements_to_unref size = " << elements_to_unref.size() << std::endl;

      // check
      if (DEBUG_UNREF)
      {
        for (ElementUnrefineCollection::iterator it=new_root.begin(); it != new_root.end(); ++it)
          {
            stk::mesh::Entity rootElement = *it;
            if (new_set.find(rootElement) != new_set.end())
              throw std::runtime_error("bad set");
          }
      }

    }

    void
    Refiner::
    unrefinePass2(ElementUnrefineCollection& elements_to_unref)
    {
      if (DEBUG_UNREF_1) std::cout << "\n\n\n unrefinePass2:: start \n\nVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV "<< std::endl;

      m_eMesh.get_bulk_data()->modification_begin();
      VERIFY_OP_ON(m_eMesh.getEntityPool().size(), ==, 0, "hmmm");

      ElementUnrefineCollection elements_to_unref_0 = elements_to_unref;

      NodeSetType kept_nodes(*m_eMesh.get_bulk_data());
      NodeSetType kept_nodes_orig(*m_eMesh.get_bulk_data());
      NodeSetType deleted_nodes(*m_eMesh.get_bulk_data());

      // filter unref set pass 2 - only root elements (no parent or not refined)
      SetOfEntities rootElements(*m_eMesh.get_bulk_data());
      SetOfEntities family_trees_to_be_removed(*m_eMesh.get_bulk_data());
      SetOfEntities children_to_be_removed(*m_eMesh.get_bulk_data());

      //unsigned elsizeb4 = elements_to_unref.size();
      filterUnrefSetPass2(elements_to_unref, rootElements);

      if (DEBUG_UNREF_1)
        {
          static int cnt1=0;
          char buf[1000];
          sprintf(buf, "%04d", cnt1);
          if (cnt1==0)
            m_eMesh.save_as("fref.e");
          else
            m_eMesh.save_as("fref.e-s"+std::string(buf));
          ++cnt1;
        }

      // get kept nodes
      // FIXME - check conditions inside getKeptNodes
      //getKeptNodes(kept_nodes, all_kept_nodes_elements);
      if (1)
        {
          const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.element_rank() );
          for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;

              const unsigned num_entity_in_bucket = bucket.size();
              for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                {
                  stk::mesh::Entity element = bucket[ientity];

                  const percept::MyPairIterRelation elem_nodes (m_eMesh, element,  stk::mesh::MetaData::NODE_RANK);

                  bool in_unref_set = elements_to_unref.find( element ) != elements_to_unref.end();
                  if (!in_unref_set)
                    {
                      for (unsigned inode=0; inode < elem_nodes.size(); inode++)
                        {
                          stk::mesh::Entity node = elem_nodes[inode].entity();
                          kept_nodes.insert(node);
                        }
                    }
                }
            }
        }

      // get deleted nodes
      //getDeletedNodes(deleted_nodes, kept_nodes, elements_to_unref);
      {
        for (ElementUnrefineCollection::iterator u_iter = elements_to_unref.begin();
             u_iter != elements_to_unref.end(); ++u_iter)
          {
            stk::mesh::Entity element = *u_iter;

            if (rootElements.find(element) != rootElements.end())
              continue;

            const percept::MyPairIterRelation elem_nodes (m_eMesh, element,  stk::mesh::MetaData::NODE_RANK);

            if (!m_eMesh.isGhostElement(element) && elem_nodes.size() && m_eMesh.hasFamilyTree(element) )
              {
                for (unsigned inode=0; inode < elem_nodes.size(); inode++)
                  {
                    stk::mesh::Entity node = elem_nodes[inode].entity();
                    bool in_kept_nodes_set = kept_nodes.find( node ) != kept_nodes.end();
                    if (!in_kept_nodes_set)
                      {
                        deleted_nodes.insert(node);
                      }
                  }
              }
          }
      }


      if (DEBUG_UNREF_1) std::cout << "deleted_nodes.size= " << deleted_nodes.size() << std::endl;
      // FIXME - additional filter step

      // remove deleted nodes and their associated sub-dim entities from NodeRegistry's db
      //   save kept_nodes_orig on NodeRegistry
      SubDimCellToDataMap to_save;
      m_nodeRegistry->cleanDeletedNodes(deleted_nodes, kept_nodes_orig, to_save);

      // first have to delete the family tree (higher ranks have to be deleted first)

      size_t elements_to_unref_sz = elements_to_unref.size();
      // remove children
      if (1)
        {
          const unsigned FAMILY_TREE_RANK = stk::mesh::MetaData::ELEMENT_RANK + 1u;

          SetOfEntities family_trees(*m_eMesh.get_bulk_data());

          for (ElementUnrefineCollection::iterator it=elements_to_unref.begin(); it != elements_to_unref.end(); ++it)
            {
              stk::mesh::Entity element = *it;

              //unsigned child_ft_level_0 = m_eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, element);
              const percept::MyPairIterRelation child_to_family_tree_relations (m_eMesh, element, FAMILY_TREE_RANK);
              if (child_to_family_tree_relations.size() == 0)
                {
                  throw std::logic_error("Refiner::unrefinePass2 child_to_family_tree_relations.size == 0");
                }

              while (true)
                {
                  percept::MyPairIterRelation rels (m_eMesh, element, FAMILY_TREE_RANK);
                  if (!rels.size())
                    break;
                  stk::mesh::Entity ft_to_rel = rels[0].entity();
                  family_trees.insert(ft_to_rel);
                  stk::mesh::RelationIdentifier ft_to_id = rels[0].relation_ordinal();

                  bool del = m_eMesh.get_bulk_data()->destroy_relation( ft_to_rel, element, ft_to_id);
                  if (!del)
                    throw std::runtime_error("Pass2:: destroy_relation failed 4");
                }

              if (! m_eMesh.get_bulk_data()->destroy_entity( element ) )
                {
                  throw std::runtime_error("pass2 4.2 - couldn't delete element");
                }

            }

          for (SetOfEntities::iterator fiter=family_trees.begin(); fiter != family_trees.end(); fiter++)
            {
              stk::mesh::Entity family_tree = *fiter;
              percept::MyPairIterRelation rels (m_eMesh, family_tree, m_eMesh.element_rank());
              if (rels.size() == 1) // only points to a parent entity
                {
                  if ( ! m_eMesh.get_bulk_data()->destroy_entity( family_tree ) )
                    {
                      throw std::runtime_error("pass2 4.1 - couldn't delete family_tree");
                    }
                }
            }
        }

      // reconnect and remove any dangling side elements
      bool allow_not_found = true;
      fix_side_sets_2(allow_not_found);

      // FIXME
      SetOfEntities parent_side_elements(*m_eMesh.get_bulk_data());
      getSideParentsToBeRemeshed(rootElements, parent_side_elements);

      // remesh the holes left by removing child elems (quad/hex hanging node doesn't need this)
      if (m_needsRemesh)
        {
          remeshElements(rootElements, m_eMesh.element_rank(), elements_to_unref_sz);
          remeshElements(parent_side_elements, m_eMesh.side_rank(), elements_to_unref_sz);
        }
      // FIXME side sets...

      // remove any elements that are empty (these can exist when doing local refinement)
      removeEmptyElements();

      removeDeletedNodes(deleted_nodes);

      set_active_part();

      fix_side_sets_2();

      m_eMesh.get_bulk_data()->modification_end();

      m_nodeRegistry->clear_element_owner_data_phase_2();

      m_eMesh.get_bulk_data()->modification_begin();
      set_active_part();
      m_eMesh.destroyEntityPool();
      m_eMesh.get_bulk_data()->modification_end();

      if (DEBUG_UNREF_1) std::cout << "\n\n ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n\n\n unrefinePass2:: end "<< std::endl;

    }

  } // namespace adapt
} // namespace stk
