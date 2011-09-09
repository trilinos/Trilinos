#include <exception>
#include <fstream>
#include <set>
#include <typeinfo>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <stk_adapt/Refiner.hpp>

#include <stk_percept/MeshUtil.hpp>


namespace stk {
  namespace adapt {
    using namespace std;
    using namespace percept;

    // ====================================================================================================
    // ====================================================================================================
    // ====================================================================================================

    // get set of children with no nieces that are to be removed; put their parents in a list; 
    //   put the children's family trees in a list to be removed
    void Refiner::getChildrenToBeRemoved(ElementUnrefineCollection& elements_to_unref,
                                         SetOfEntities& children_to_be_removed, SetOfEntities& children_to_be_removed_with_ghosts, 
                                         SetOfEntities& copied_children_to_be_removed,
                                         SetOfEntities& family_trees_to_be_removed, 
                                         SetOfEntities& parent_elements)
    {
      const unsigned FAMILY_TREE_RANK = m_eMesh.element_rank() + 1u;

      unsigned nchild_removed = 0;

      for (ElementUnrefineCollection::iterator u_iter = elements_to_unref.begin();
           u_iter != elements_to_unref.end(); ++u_iter)
        {
          stk::mesh::Entity * element_p = *u_iter;
          bool isGhostElement = m_eMesh.isGhostElement(*element_p);
          bool isChild = m_eMesh.isChildElement(*element_p);

          if (!isChild)
            {
              throw std::logic_error(" found parent where child expected");
              continue;
            }

#if DEBUG_UNREF
          //std::cout << "tmp element to be removed id= " << element_p->identifier() << " " << std::endl;
#endif

          bool inCopiedList = copied_children_to_be_removed.find(element_p) != copied_children_to_be_removed.end();

          if (inCopiedList)
            {
              std::vector<stk::mesh::Entity *> siblings;
              stk::mesh::PairIterRelation child_to_family_tree_relations = element_p->relations(FAMILY_TREE_RANK);

              // look for level 0 only - these are children with no children
              unsigned child_ft_level_0 = m_eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, *element_p);

              stk::mesh::Entity *family_tree = child_to_family_tree_relations[child_ft_level_0].entity();
              stk::mesh::PairIterRelation family_tree_relations = family_tree->relations(m_eMesh.element_rank());
              if (family_tree_relations.size() == 0)
                {
                  throw std::logic_error("Refiner::unrefineTheseElements family_tree_relations.size() == 0");
                }

              for (unsigned ichild=1; ichild < family_tree_relations.size(); ichild++)
                {
                  stk::mesh::Entity *child = family_tree_relations[ichild].entity();
                  if (!m_eMesh.isParentElement(*child))
                    {
                      siblings.push_back(child);
                      copied_children_to_be_removed.erase(child);
                    }
                  else
                    {
                      throw std::logic_error("Refiner::unrefineTheseElements found parent where child expected in siblings list");
                    }
                }

#if DEBUG_UNREF
              std::cout << "tmp removing family_tree: " << family_tree->identifier() << std::endl;
              //stk::mesh::EntityId family_tree_id =  family_tree->identifier() ;
#endif

              stk::mesh::Entity *parent = family_tree_relations[FAMILY_TREE_PARENT].entity();
              if (!parent)
                {
                  throw std::logic_error("Refiner::unrefineTheseElements parent == null");
                }
              parent_elements.insert(parent);

              family_trees_to_be_removed.insert(family_tree);

              for (unsigned ichild=0; ichild < siblings.size(); ichild++)
                {
                  stk::mesh::Entity *child = siblings[ichild];

                  bool isSiblingAGhost = m_eMesh.isGhostElement(*child);
                  if (isSiblingAGhost != isGhostElement)
                    {
                      std::cout << "P["<< m_eMesh.getRank() << "] Refiner::buildChildList isSiblingAGhost,  ghost= " << isSiblingAGhost << " " << isGhostElement << std::endl;
                      throw std::logic_error("isSiblingAGhost != isGhostElement");
                    }
#if DEBUG_UNREF
                  //std::cout << "tmp removing child: " << child->identifier() << " " << *child << std::endl;
#endif
                  ++nchild_removed;
                  
                  if (!isGhostElement) 
                    { 
                      children_to_be_removed.insert( child );
                    }
                  children_to_be_removed_with_ghosts.insert( child );
                }

            }
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
          stk::mesh::Entity *family_tree = *family_tree_it;
          if ( ! m_eMesh.getBulkData()->destroy_entity( family_tree ) )
            {
              throw std::logic_error("Refiner::unrefineTheseElements couldn't remove element, destroy_entity returned false for family_tree.");
            }
        }
    }


    void Refiner::getSideElemsToBeRemoved(SetOfEntities& children_to_be_removed, SetOfEntities& side_elem_set_to_be_removed, SetOfEntities& family_trees_to_be_removed, SetOfEntities& parent_side_elements)
    {
      if (getIgnoreSideSets()) return;

      const unsigned FAMILY_TREE_RANK = m_eMesh.element_rank() + 1u;
      for(SetOfEntities::iterator child_it = children_to_be_removed.begin();
          child_it != children_to_be_removed.end(); ++child_it)
        {
          stk::mesh::Entity *child = *child_it;

          if (m_eMesh.hasFamilyTree(*child) && !m_eMesh.isChildWithoutNieces(*child))
            throw std::logic_error("error 34");

          // add sideset elements to list to be removed (and their family tree info)
          mesh::PairIterRelation side_relations = child->relations(m_eMesh.side_rank());
          for (unsigned jside = 0; jside < side_relations.size(); jside++)
            {
              stk::mesh::Entity * side_element = side_relations[jside].entity();
              // FIXME err check
              if (m_eMesh.hasFamilyTree(*side_element) && !m_eMesh.isChildWithoutNieces(*side_element))
                {
                  throw std::logic_error("error 35");
                }
              //if (m_eMesh.hasFamilyTree(*side_element) && m_eMesh.isChildWithoutNieces(*side_element, false))
              if (m_eMesh.hasFamilyTree(*side_element))
              {
                unsigned side_elem_child_ft_level_0 = m_eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, *side_element);

                mesh::PairIterRelation side_element_to_family_tree_relations = side_element->relations(FAMILY_TREE_RANK);
                stk::mesh::Entity *family_tree = side_element_to_family_tree_relations[side_elem_child_ft_level_0].entity();
                stk::mesh::PairIterRelation family_tree_relations = family_tree->relations(side_element->entity_rank());
                if (family_tree_relations.size() == 0)
                  {
                    throw std::logic_error("Refiner::unrefineTheseElements family_tree_relations.size() == 0 [1]");
                  }

                for (unsigned ichild=1; ichild < family_tree_relations.size(); ichild++)
                  {
                    stk::mesh::Entity *side_elem_sibling = family_tree_relations[ichild].entity();
            
                    //std::cout << "tmp unref side element id= " << side_element->identifier() << std::endl;
                    side_elem_set_to_be_removed.insert(side_elem_sibling);
                  }

                stk::mesh::Entity *parent = family_tree_relations[FAMILY_TREE_PARENT].entity();
                parent_side_elements.insert(parent);
                family_trees_to_be_removed.insert(family_tree);
              }
            }
        }
    }



    void Refiner::removeChildElements(SetOfEntities& children_to_be_removed)
    {
      const unsigned FAMILY_TREE_RANK = m_eMesh.element_rank() + 1u;
      for(SetOfEntities::iterator child_it = children_to_be_removed.begin();
          child_it != children_to_be_removed.end(); ++child_it)
        {
          stk::mesh::Entity *child = *child_it;

          if ( ! m_eMesh.getBulkData()->destroy_entity( child ) )
            {
              CellTopology cell_topo(stk::percept::PerceptMesh::get_cell_topology(*child));

              //const mesh::PairIterRelation elem_relations = child->relations(child->entity_rank()+1);
              const mesh::PairIterRelation child_to_ft_relations = child->relations(FAMILY_TREE_RANK);
#if DEBUG_UNREF
              std::cout << "tmp Refiner::unrefineTheseElements couldn't remove element  cell= " << cell_topo.getName() << std::endl;
              std::cout << "tmp child_to_ft_relations.size() = " << child_to_ft_relations.size() << std::endl;
              //std::cout << "tmp ft_id loc, outerloop= " << child_to_family_tree_relations[0].entity()->identifier() << " " << family_tree_id << std::endl;

              m_eMesh.printEntity(std::cout, *child);
#endif

              throw std::logic_error("Refiner::unrefineTheseElements couldn't remove element, destroy_entity returned false.");
            }
        }
    }

    void Refiner::removeSideElements(SetOfEntities& side_elem_set_to_be_removed, SetOfEntities& elements_to_be_deleted)
    {
      for(SetOfEntities::iterator side_elem_it = side_elem_set_to_be_removed.begin();
          side_elem_it != side_elem_set_to_be_removed.end(); ++side_elem_it)
        {
          stk::mesh::Entity *side_elem = *side_elem_it;
          bool del =  m_eMesh.getBulkData()->destroy_entity( side_elem );
          if ( ! del )
            {
              std::cout << "Refiner::unrefineTheseElements couldn't remove side element, destroy_entity returned false." << *side_elem << std::endl;
              del =  m_eMesh.getBulkData()->destroy_entity( side_elem );
              if (1)
                {
                  stk::mesh::PairIterRelation rels = side_elem->relations(m_eMesh.element_rank());
                  for (unsigned irels=0; irels < rels.size(); irels++)
                    {
                      bool in_del = (elements_to_be_deleted.find(rels[irels].entity()) != elements_to_be_deleted.end());
                      std::cout << "Refiner::unrefineTheseElements found element: in del list= " << in_del << std::endl;
                    }
                }
              throw std::logic_error("Refiner::unrefineTheseElements couldn't remove side element, destroy_entity returned false.");
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
      std::cout << "tmp remesh:: parent_elements.size() [elements to be remeshed] = " << parent_elements.size() << std::endl;
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

          vector<stk::mesh::Entity *> new_elements;

          // count num new elements needed on this proc (served by UniformRefinerPattern)
          unsigned num_elem_not_ghost = 0u;

          for (ElementUnrefineCollection::iterator p_iter = parent_elements.begin();
               p_iter != parent_elements.end(); ++p_iter)
            {
              stk::mesh::Entity *parent = *p_iter;

              stk::mesh::Entity& element = *parent;

              const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);
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
          std::cout << "tmp remesh::rank[irank], num_elem_needed= " << m_ranks[irank] << " " << num_elem_needed << std::endl;
#endif

          // create new entities on this proc
          new_elements.resize(0);
          m_eMesh.createEntities( m_ranks[irank], num_elem_needed, new_elements);
          vector<stk::mesh::Entity *>::iterator element_pool_it = new_elements.begin();

          // FIXME - we could directly call this with a refactor to change elementColors passed in here as a generic collection + checking for element Type
          //
          //createElementsAndNodesAndConnectLocal(m_ranks[irank], m_breakPattern[irank], elementColors, needed_entity_ranks, new_elements);

          for (ElementUnrefineCollection::iterator p_iter = parent_elements.begin();
               p_iter != parent_elements.end(); ++p_iter)
            {
              stk::mesh::Entity *parent_p = *p_iter;

              const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(*parent_p);
              CellTopology cell_topo(cell_topo_data);
              unsigned elementType = cell_topo.getKey();
              unsigned bpElementType = m_breakPattern[irank]->getFromTypeKey();
              if (elementType == bpElementType)
                {
                  stk::mesh::Entity& parent = *parent_p;

                  if (!m_eMesh.isGhostElement(parent))
                    {
#if DEBUG_UNREF
                      //std::cout << "P["<< m_eMesh.getRank() << "] parent.owner_rank() = " << parent.owner_rank() << std::endl;
                      std::cout << "tmp Parent to be remeshed = ";
                      m_eMesh.printEntity(std::cout, parent);
#endif
                      if (createNewNeededNodeIds(cell_topo_data, parent, needed_entity_ranks, new_sub_entity_nodes))
                        {
                          //std::cout << "typeid= " << typeid(*breakPattern).name() << std::endl;
                          throw std::logic_error("unrefineTheseElements:: createNewNeededNodeIds failed");
                        }

                      vector<stk::mesh::Entity *>::iterator element_pool_it_b4 = element_pool_it;
                      m_breakPattern[irank]->createNewElements(m_eMesh, *m_nodeRegistry, parent, new_sub_entity_nodes, element_pool_it, m_proc_rank_field);
                      vector<stk::mesh::Entity *>::iterator element_pool_it_af = element_pool_it;
                      num_new_elem_during_remesh += (element_pool_it_af - element_pool_it_b4);
                    }
                }
            }
#if DEBUG_UNREF
          std::cout << "tmp remesh:: nchild elements during remesh for rank[irank] = " << m_ranks[irank] << " " << num_new_elem_during_remesh << std::endl;
#endif

        } // irank

    }


    void
    Refiner::
    unrefineAll()
    {
      ElementUnrefineCollection elements_to_unref;

      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.getBulkData()->buckets( m_eMesh.element_rank() );

      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          //if (removePartSelector(**k))
          {
            stk::mesh::Bucket & bucket = **k ;

            const unsigned num_entity_in_bucket = bucket.size();
            for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
              {
                stk::mesh::Entity& element = bucket[ientity];
                // FIXME
                // skip elements that are already a parent (if there's no family tree yet, it's not a parent, so avoid throwing an error in isParentElement)
                const bool check_for_family_tree = false;  
                bool isParent = m_eMesh.isParentElement(element, check_for_family_tree);
              
                if (isParent)
                  continue;

                const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);

                if (elem_nodes.size() && (m_eMesh.hasFamilyTree(element) && m_eMesh.isChildWithoutNieces(element, false) ) )
                  {
                    bool elementIsGhost = m_eMesh.isGhostElement(element);

                    if (!elementIsGhost)
                      {
                        elements_to_unref.insert(&element);
                      }
                  }
              }
          }
        }
      unrefineTheseElements(elements_to_unref);
    }

#define DEBUG_UNREF 0
#define DEBUG_UNREF_1 0

    void Refiner::
    filterUnrefSet(ElementUnrefineCollection& elements_to_unref)
    {
      int print_filter_info = DEBUG_UNREF_1;
      if (print_filter_info)  std::cout << "P["<< m_eMesh.getRank() << "] filterUnrefSet: initial set size = " << elements_to_unref.size() << std::endl;
      
      const unsigned FAMILY_TREE_RANK = m_eMesh.element_rank() + 1u;
      ElementUnrefineCollection elements_to_unref_copy;

      typedef std::set<stk::mesh::Entity *> SetOfEntities;

      int num_is_parent = 0;
      int num_elem_nodes_0 = 0;
      int num_has_nieces = 0;

      for (ElementUnrefineCollection::iterator u_iter = elements_to_unref.begin();
           u_iter != elements_to_unref.end(); ++u_iter)
        {
          stk::mesh::Entity& element = **u_iter;
          
          const bool check_for_family_tree = false;  
          bool isParent = m_eMesh.isParentElement(element, check_for_family_tree);

          if (isParent)
            {
              ++num_is_parent;
              continue;
            }

          const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);
          int elem_nodes_size = elem_nodes.size();
          bool has_no_nieces = m_eMesh.hasFamilyTree(element) && m_eMesh.isChildWithoutNieces(element, false);

          if (!elem_nodes_size) num_elem_nodes_0++;
          if (!has_no_nieces) num_has_nieces++;
          if (elem_nodes_size && has_no_nieces)
            {
              stk::mesh::PairIterRelation child_to_family_tree_relations = element.relations(FAMILY_TREE_RANK);

              // look for level 0 only - these are children with no children
              unsigned child_ft_level_0 = m_eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, element);

              stk::mesh::Entity *family_tree = child_to_family_tree_relations[child_ft_level_0].entity();
              stk::mesh::PairIterRelation family_tree_relations = family_tree->relations(m_eMesh.element_rank());
              if (family_tree_relations.size() == 0)
                {
                  throw std::logic_error("Refiner::filterUnrefSet family_tree_relations.size() == 0");
                }

              SetOfEntities side_elem_set;

              bool all_siblings_in_unref_set = true;
              bool all_side_sets_ok = true;

              for (unsigned ichild=1; ichild < family_tree_relations.size(); ichild++)
                {
                  stk::mesh::Entity *child = family_tree_relations[ichild].entity();
                  if (m_eMesh.isParentElement(*child))
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
                    mesh::PairIterRelation side_relations = child->relations(m_eMesh.side_rank());
                    for (unsigned jside = 0; jside < side_relations.size(); jside++)
                      {
                        stk::mesh::Entity * side_element = side_relations[jside].entity();
                        side_elem_set.insert(side_element);

                        if (0)
                        {
                            stk::mesh::PairIterRelation side_elem_to_family_tree_relations = side_element->relations(FAMILY_TREE_RANK);
                            stk::mesh::Entity *side_elem_family_tree_0 = side_elem_to_family_tree_relations[0].entity();
                            std::cout << "side_elem_family_tree_0= " << side_elem_family_tree_0 << std::endl;
                        }

                        // FIXME if (!m_eMesh.hasFamilyTree(*side_element) || !m_eMesh.isChildWithoutNieces(*side_element))
                        if (m_eMesh.hasFamilyTree(*side_element) && !m_eMesh.isChildWithoutNieces(*side_element))
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
                      stk::mesh::Entity *child = family_tree_relations[ichild].entity();

                      {
                        mesh::PairIterRelation side_relations = child->relations(m_eMesh.side_rank());
                        for (unsigned jside = 0; jside < side_relations.size(); jside++)
                          {
                            stk::mesh::Entity * side_element = side_relations[jside].entity();
                            //side_elem_set.insert(side_element);

                            if (!m_eMesh.hasFamilyTree(*side_element)) continue;

                            stk::mesh::PairIterRelation side_elem_to_family_tree_relations = side_element->relations(FAMILY_TREE_RANK);

                            // look for level 0 only - these are children with no children
                            unsigned side_elem_ft_level_0 = m_eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, *side_element);

                            stk::mesh::Entity *side_elem_family_tree = side_elem_to_family_tree_relations[side_elem_ft_level_0].entity();
                            stk::mesh::PairIterRelation side_elem_family_tree_relations = side_elem_family_tree->relations(m_eMesh.side_rank());
                            for (unsigned ise_child=1; ise_child < side_elem_family_tree_relations.size(); ise_child++)
                              {
                                stk::mesh::Entity *se_sibling = side_elem_family_tree_relations[ise_child].entity();
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

              if (all_siblings_in_unref_set && all_side_sets_ok)
                {
                  for (unsigned ichild=1; ichild < family_tree_relations.size(); ichild++)
                    {
                      stk::mesh::Entity *child = family_tree_relations[ichild].entity();
                      elements_to_unref_copy.insert(child);
                    }
                }
            }
        }
      if (print_filter_info) std::cout << "tmp filterUnrefSet::elements_to_unref.size = " << elements_to_unref.size() 
                                       << " filtered size= " << elements_to_unref_copy.size() 
                                       << " num_has_nieces= " << num_has_nieces
                                       << " num_elem_nodes_0= " << num_elem_nodes_0
                                       << " num_is_parent= " << num_is_parent
                                       << std::endl;
      elements_to_unref = elements_to_unref_copy;
    }

    void Refiner::
    getKeptNodes(NodeSetType& kept_nodes, ElementUnrefineCollection& elements_to_unref)
    {
      bool doTest=true;

      // mark kept nodes
      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.getBulkData()->buckets( m_eMesh.element_rank() );

      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          //if (removePartSelector(**k))
          {
            stk::mesh::Bucket & bucket = **k ;

            const unsigned num_entity_in_bucket = bucket.size();
            for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
              {
                stk::mesh::Entity& element = bucket[ientity];

                const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);

                //if (m_eMesh.isLeafElement(element) && !m_eMesh.isGhostElement(element))
                if (!doTest || (elem_nodes.size() && m_eMesh.isLeafElement(element)) )
                  {
                    bool in_unref_set = elements_to_unref.find( &element ) != elements_to_unref.end();
                    if (!in_unref_set)
                      {
                        for (unsigned inode=0; inode < elem_nodes.size(); inode++)
                          {
                            stk::mesh::Entity *node = elem_nodes[inode].entity();
                            kept_nodes.insert(node);
#if DEBUG_UNREF
                            std::cout << "tmp kept node: " << *node << " ";
                            m_eMesh.printEntity(std::cout, *node);
#endif

                          }
                      }
                  }
              }
          }
        }


    }

    void Refiner::
    getDeletedNodes(NodeSetType& deleted_nodes, const NodeSetType& kept_nodes, ElementUnrefineCollection& elements_to_unref)
    {
      // mark deleted nodes (nodes in list of elements_to_unref
      for (ElementUnrefineCollection::iterator u_iter = elements_to_unref.begin();
           u_iter != elements_to_unref.end(); ++u_iter)
        {
          stk::mesh::Entity * element = *u_iter;

          //if (m_eMesh.isChildElement(*element) && !m_eMesh.isGhostElement(*element))
          {
            const mesh::PairIterRelation elem_nodes = element->relations(stk::mesh::fem::FEMMetaData::NODE_RANK);

            for (unsigned inode=0; inode < elem_nodes.size(); inode++)
              {
                stk::mesh::Entity *node = elem_nodes[inode].entity();
                bool in_kept_nodes_set = kept_nodes.find( node ) != kept_nodes.end();
                if (!in_kept_nodes_set)
                  {
                    deleted_nodes.insert(node);
#if DEBUG_UNREF
                    std::cout << "tmp deleted node: " << *node << " ";
                    m_eMesh.printEntity(std::cout, *node);
#endif
                  }
              }
          }
        }
      if (DEBUG_UNREF_1) std::cout << "tmp kept_nodes size= " << kept_nodes.size() << " deleted_nodes size= " << deleted_nodes.size() << std::endl;

    }

    void
    Refiner::
    unrefineTheseElements(ElementUnrefineCollection& elements_to_unref)
    {
      if (m_alwaysInitNodeRegistry)
        {
          throw std::logic_error("Refiner::unrefineTheseElements: to use urefinement, you must have setAlwaysInitializeNodeRegistry(false)");
        }

      m_eMesh.getBulkData()->modification_begin();

      //check_sidesets_2(" unrefineTheseElements:: start");

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


      NodeSetType kept_nodes;
      NodeSetType deleted_nodes;

      // filter unref set
      filterUnrefSet(elements_to_unref);

      // get kept nodes
      getKeptNodes(kept_nodes, elements_to_unref);

      // get deleted nodes
      getDeletedNodes(deleted_nodes, kept_nodes, elements_to_unref);

      // nothing to be done
      if (0 && deleted_nodes.size() == 0)
        {
          m_eMesh.getBulkData()->modification_end();
          std::cout << "Refiner::unrefineTheseElements: deleted_nodes size is 0, nothing to be done, early return." << std::endl;
          return;
        }

      // remove deleted nodes and their associated sub-dim entities
      m_nodeRegistry->cleanDeletedNodes(deleted_nodes);

      // remove elements to be unrefined
      ElementUnrefineCollection copied_children_to_be_removed = elements_to_unref;

#if DEBUG_UNREF || DEBUG_UNREF_1
      std::cout << "tmp copied_children_to_be_removed.size() [= num elements to be urefined that are children and !ghosts]= " << copied_children_to_be_removed.size() << std::endl;
#endif

      typedef std::set<stk::mesh::Entity *> SetOfEntities;

      SetOfEntities family_trees_to_be_removed;
      SetOfEntities children_to_be_removed;
      SetOfEntities children_to_be_removed_with_ghosts;
      SetOfEntities parent_elements;

      // set to hold sideset elements to be removed
      SetOfEntities side_elem_set_to_be_removed;
      SetOfEntities side_elem_family_trees_to_be_removed;
      SetOfEntities parent_side_elements;

      // remove elements marked for unrefine (make sure they are children and not ghosts)
      getChildrenToBeRemoved(elements_to_unref,
                             children_to_be_removed, children_to_be_removed_with_ghosts, copied_children_to_be_removed, 
                             family_trees_to_be_removed,
                             parent_elements);

      // NodeRegistry must be kept in sync
      m_nodeRegistry->clear_element_owner_data(children_to_be_removed_with_ghosts);

      // get the corresponding side elements from the children->side relations
      getSideElemsToBeRemoved(children_to_be_removed, side_elem_set_to_be_removed, side_elem_family_trees_to_be_removed,
                              parent_side_elements);

      // first have to delete the family tree (higher ranks have to be deleted first)
      removeFamilyTrees(family_trees_to_be_removed);

      // remove children 
      removeChildElements(children_to_be_removed);

      // for sideset elements, remove their family trees first
      removeFamilyTrees(side_elem_family_trees_to_be_removed);

      // remove sideset elements
      removeSideElements(side_elem_set_to_be_removed, children_to_be_removed);

      // remesh the holes left by removing child elems
      remesh(parent_elements);
      remesh(parent_side_elements);

#if CHECK_DEBUG
      check_db("after unrefineTheseElements, b4 mod end");
#endif
      // remove any elements that are empty (these can exist when doing local refinement)
      removeEmptyElements();

      set_active_part();

      fixElementSides1();

      m_eMesh.getBulkData()->modification_end();

#if CHECK_DEBUG
      check_db("before clear_element_owner_data_phase_2");
#endif

      m_nodeRegistry->clear_element_owner_data_phase_2();

      //check_sidesets_2(" unrefineTheseElements:: end");


#if CHECK_DEBUG
      check_db("after unrefineTheseElements");
#endif

    }

  } // namespace adapt
} // namespace stk
