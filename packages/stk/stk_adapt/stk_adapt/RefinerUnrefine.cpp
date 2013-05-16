#include <exception>
#include <fstream>
#include <set>
#include <typeinfo>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <stk_adapt/Refiner.hpp>

#include <stk_percept/MeshUtil.hpp>

#define DEBUG_UNREF 1
#define DEBUG_UNREF_1 1
#define DEBUG_UNREF_2 0

namespace stk {
  namespace adapt {
    using namespace std;
    using namespace percept;

    // ====================================================================================================
    // ====================================================================================================
    // ====================================================================================================

    void Refiner::remove_dangling_sidesets()
    {
      SetOfEntities sides_to_delete(*m_eMesh.get_bulk_data());
      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.side_rank() );
      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          const unsigned num_entity_in_bucket = bucket.size();
          for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
            {
              stk::mesh::Entity side = bucket[ientity];
              if (!sharesElementFace(side))
                sides_to_delete.insert(side);
            }
        }

      const unsigned FAMILY_TREE_RANK = stk::mesh::MetaData::ELEMENT_RANK + 1u;
      SetOfEntities family_trees(*m_eMesh.get_bulk_data());
      for (SetOfEntities::iterator siter=sides_to_delete.begin(); siter != sides_to_delete.end(); ++siter)
        {
          stk::mesh::Entity side = *siter;

          while (true)
            {
              percept::MyPairIterRelation rels (m_eMesh, side, FAMILY_TREE_RANK);
              if (!rels.size())
                break;
              stk::mesh::Entity to_rel = rels[0].entity();
              family_trees.insert(to_rel);
              stk::mesh::RelationIdentifier to_id = rels[0].relation_ordinal();

              bool del = m_eMesh.get_bulk_data()->destroy_relation( to_rel, side, to_id);
              if (!del)
                throw std::runtime_error("remove_dangling_sidesets:: destroy_relation failed 4");
            }

          if ( ! m_eMesh.get_bulk_data()->destroy_entity( side ) )
            {
              throw std::runtime_error("remove_dangling_sidesets error 4 - couldn't delete");
            }
        }
      for (SetOfEntities::iterator fiter=family_trees.begin(); fiter != family_trees.end(); fiter++)
        {
          stk::mesh::Entity family_tree = *fiter;
          percept::MyPairIterRelation rels (m_eMesh, family_tree, m_eMesh.side_rank());
          if (rels.size() == 1)
            {
              if ( ! m_eMesh.get_bulk_data()->destroy_entity( family_tree ) )
                {
                  throw std::runtime_error("remove_dangling_sidesets error 4.1 - couldn't delete family_tree");
                }
            }
        }

    }

    static void getLevels(PerceptMesh& eMesh, ElementUnrefineCollection& elements_to_unref, int& min_level, int& max_level)
    {
      stk::mesh::FieldBase *refine_level = eMesh.get_field("refine_level");
      if (!refine_level) return;
      min_level=std::numeric_limits<int>::max();
      max_level=0;

      for (ElementUnrefineCollection::iterator it=elements_to_unref.begin(); it != elements_to_unref.end(); ++it)
        {
          stk::mesh::Entity element = *it;
          const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);
          if (!cell_topo_data) continue;
          double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(refine_level) , element );
          if (fdata)
            {
              min_level = std::min(min_level, (int)fdata[0]);
              max_level = std::max(max_level, (int)fdata[0]);
            }
        }
      if (DEBUG_UNREF_1)
        std::cout << "min_level= " << min_level << " max_level= " << max_level << std::endl;
    }

    void getLevels(PerceptMesh& eMesh, unsigned rank, int& min_level, int& max_level)
    {
      const vector<stk::mesh::Bucket*> & buckets = eMesh.get_bulk_data()->buckets( rank );
      stk::mesh::FieldBase *refine_level = eMesh.get_field("refine_level");
      if (!refine_level) return;
      min_level=std::numeric_limits<int>::max();
      max_level=0;

      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(bucket);
          if (!cell_topo_data) continue;
          const unsigned num_entity_in_bucket = bucket.size();
          for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
            {
              stk::mesh::Entity element = bucket[ientity];
              if (refine_level)
                {
                  double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(refine_level) , element );
                  if (fdata)
                    {
                      min_level = std::min(min_level, (int)fdata[0]);
                      max_level = std::max(max_level, (int)fdata[0]);
                    }
                }
            }
        }
      if (DEBUG_UNREF_1)
        std::cout << "min_level= " << min_level << " max_level= " << max_level << std::endl;
    }

    static void filterForLevel(PerceptMesh& eMesh, ElementUnrefineCollection& elements_to_unref, int level)
    {
      stk::mesh::FieldBase *refine_level = eMesh.get_field("refine_level");
      if (!refine_level) return;

      ElementUnrefineCollection new_set(*eMesh.get_bulk_data());
      for (ElementUnrefineCollection::iterator it=elements_to_unref.begin(); it != elements_to_unref.end(); ++it)
        {
          stk::mesh::Entity element = *it;
          double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(refine_level) , element );
          if (fdata && (int)fdata[0] >= level)
            {
              //std::cout << "insert fdata[0] = " << fdata[0] << " level= " << level << std::endl;
              new_set.insert(element);
            }
        }
      elements_to_unref = new_set;
    }

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
                    const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(bucket);

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

    void Refiner::
    remeshDuringRefine(stk::mesh::EntityRank rank)
    {
      const unsigned FAMILY_TREE_RANK = stk::mesh::MetaData::ELEMENT_RANK + 1u;

      SetOfEntities parents(*m_eMesh.get_bulk_data());

      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( rank );

      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
            stk::mesh::Bucket & bucket = **k ;

            const unsigned num_entity_in_bucket = bucket.size();
            for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
              {
                stk::mesh::Entity element = bucket[ientity];

                if (m_eMesh.hasFamilyTree(element))
                  {
                    stk::mesh::Entity parent = m_eMesh.getParent(element, true);
                    if (parent.is_valid())
                      {
                        if (!m_eMesh.hasGrandChildren(parent, true))
                          {
                            parents.insert(parent);
                          }
                      }
                  }
              }
        }

      for (SetOfEntities::iterator pIter = parents.begin(); pIter != parents.end(); ++pIter)
        {
          stk::mesh::Entity parent = *pIter;

          //std::cout << "\n parent= " << parent.identifier() << std::endl;

          std::vector<stk::mesh::Entity> children;
          m_eMesh.getChildren(parent, children, true, false);
          VERIFY_OP_ON(children.size(), !=, 0, "hmmm");
          SetOfEntities familyTrees(*m_eMesh.get_bulk_data());
          SetOfEntities childrenSet(children.begin(), children.end(), *m_eMesh.get_bulk_data());

          for (SetOfEntities::iterator icgp = childrenSet.begin(); icgp != childrenSet.end(); ++icgp)
            {
              stk::mesh::Entity child = *icgp;
              //std::cout << "childOrGrandChild= " << childOrGrandChild.identifier() << " isGranchChild= " << (childrenSet.find(childOrGrandChild) == childrenSet.end()) << std::endl;

              unsigned child_ft_level_0 = m_eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, child);
              const percept::MyPairIterRelation child_to_family_tree_relations (m_eMesh, child,FAMILY_TREE_RANK);
              if (child_to_family_tree_relations.size() == 0)
                {
                  throw std::logic_error("Refiner::remeshDuringRefine child_to_family_tree_relations.size == 0");
                }
              stk::mesh::Entity family_tree = child_to_family_tree_relations[child_ft_level_0].entity();
              percept::MyPairIterRelation family_tree_relations (m_eMesh, family_tree,rank);
              if (family_tree_relations.size() == 0)
                {
                  throw std::logic_error("Refiner::remeshDuringRefine family_tree_relations.size() == 0");
                }

              familyTrees.insert(family_tree);
            }

          removeFamilyTrees(familyTrees);
          removeChildElements(childrenSet);

          SetOfEntities parent1(&parent, (&parent)+1, *m_eMesh.get_bulk_data());
          remesh(parent1);

        }

      m_nodeRegistry->clear_element_owner_data_phase_2(false, false);
      SetOfEntities emptySet(*m_eMesh.get_bulk_data());
      replaceNodeRegistryOwnership(emptySet, rank);
      m_nodeRegistry->clear_element_owner_data_phase_2(true, false);

    }

    void Refiner::
    preUnrefine(ElementUnrefineCollection& elements_to_unref, stk::mesh::EntityRank rank)
    {
      int print_filter_info = 1 || DEBUG_UNREF_1;
      if (print_filter_info)  std::cout << "\n\nP["<< m_eMesh.get_rank() << "] preUnrefine: initial set size = " << elements_to_unref.size() << std::endl;

      const unsigned FAMILY_TREE_RANK = stk::mesh::MetaData::ELEMENT_RANK + 1u;

      SetOfEntities grandParents(*m_eMesh.get_bulk_data());

      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( rank );

      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          {
            stk::mesh::Bucket & bucket = **k ;

            const unsigned num_entity_in_bucket = bucket.size();
            for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
              {
                stk::mesh::Entity element = bucket[ientity];

                if (m_eMesh.hasFamilyTree(element))
                  {
                    stk::mesh::Entity gp = m_eMesh.getGrandParent(element, true);
                    if (gp.is_valid())
                      {
                        if (!m_eMesh.hasGreatGrandChildren(gp, true))
                          {
                            grandParents.insert(gp);
                          }
                      }
                  }
              }
          }
        }

      for (SetOfEntities::iterator gpIter = grandParents.begin(); gpIter != grandParents.end(); ++gpIter)
        {
          stk::mesh::Entity grandParent = *gpIter;

          //std::cout << "\n grandParent= " << grandParent.identifier() << std::endl;

          std::vector<stk::mesh::Entity> children;
          m_eMesh.getChildren(grandParent, children, true, false);
          VERIFY_OP_ON(children.size(), !=, 0, "hmmm");
          SetOfEntities familyTrees(*m_eMesh.get_bulk_data());
          SetOfEntities childrenOrGrandChildrenSet(children.begin(), children.end(), *m_eMesh.get_bulk_data());
          SetOfEntities childrenSet(children.begin(), children.end(), *m_eMesh.get_bulk_data());

          for (unsigned ichild=0; ichild < children.size(); ichild++)
            {
              std::vector<stk::mesh::Entity> grandChildren;
              m_eMesh.getChildren(children[ichild], grandChildren, true, false);
              if (grandChildren.size())
                {
                  childrenOrGrandChildrenSet.insert(grandChildren.begin(), grandChildren.end());
                }
            }

          for (SetOfEntities::iterator icgp = childrenOrGrandChildrenSet.begin(); icgp != childrenOrGrandChildrenSet.end(); ++icgp)
            {
              stk::mesh::Entity childOrGrandChild = *icgp;
              //std::cout << "childOrGrandChild= " << childOrGrandChild.identifier() << " isGranchChild= " << (childrenSet.find(childOrGrandChild) == childrenSet.end()) << std::endl;


              unsigned child_ft_level_0 = m_eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, childOrGrandChild);
              const percept::MyPairIterRelation child_to_family_tree_relations (m_eMesh, childOrGrandChild,FAMILY_TREE_RANK);
              if (child_to_family_tree_relations.size() == 0)
                {
                  throw std::logic_error("Refiner::preUnrefine child_to_family_tree_relations.size == 0");
                }
              stk::mesh::Entity family_tree = child_to_family_tree_relations[child_ft_level_0].entity();
              percept::MyPairIterRelation family_tree_relations (m_eMesh, family_tree,rank);
              if (family_tree_relations.size() == 0)
                {
                  throw std::logic_error("Refiner::preUnrefine family_tree_relations.size() == 0");
                }

              familyTrees.insert(family_tree);
            }

          removeFamilyTrees(familyTrees);
          removeChildElements(childrenOrGrandChildrenSet);

          SetOfEntities grandParent1(&grandParent, (&grandParent)+1, *m_eMesh.get_bulk_data());
          remesh(grandParent1);
          std::vector<stk::mesh::Entity> gpChildren;
          m_eMesh.getChildren(grandParent, gpChildren, false, false);
          for (unsigned igpc=0; igpc < gpChildren.size(); igpc++)
            {
              stk::mesh::Entity gpc = gpChildren[igpc];
              SetOfEntities gpc1(&gpc, (&gpc)+1, *m_eMesh.get_bulk_data());
              remesh(gpc1);
            }

        }

      m_nodeRegistry->clear_element_owner_data_phase_2(false, false);
      SetOfEntities emptySet(*m_eMesh.get_bulk_data());
      replaceNodeRegistryOwnership(emptySet, rank);
      m_nodeRegistry->clear_element_owner_data_phase_2(true, false);

      set_active_part();

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
          if (!parent.is_valid())
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
              throw std::logic_error("Refiner::unrefineTheseElements couldn't remove element, destroy_entity returned false for family_tree.");
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
              //throw std::logic_error("Refiner::unrefineTheseElements couldn't remove node, destroy_entity returned false for node.");
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
                percept::MyPairIterRelation family_tree_relations (m_eMesh, family_tree,side_element.entity_rank());
                if (family_tree_relations.size() == 0)
                  {
                    throw std::logic_error("Refiner::unrefineTheseElements family_tree_relations.size() == 0 [1]");
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
              CellTopology cell_topo(stk::percept::PerceptMesh::get_cell_topology(child));

              //const percept::MyPairIterRelation elem_relations ( child->relations(child->entity_rank(m_eMesh,)+1);
              const percept::MyPairIterRelation child_to_ft_relations (m_eMesh, child,FAMILY_TREE_RANK);
#if DEBUG_UNREF
              std::cout << "tmp Refiner::unrefineTheseElements couldn't remove element  cell= " << cell_topo.getName() << std::endl;
              std::cout << "tmp child_to_ft_relations.size() = " << child_to_ft_relations.size() << std::endl;
              //std::cout << "tmp ft_id loc, outerloop= " << child_to_family_tree_relations[0].entity()->identifier() << " " << family_tree_id << std::endl;

              m_eMesh.print_entity(std::cout, child);
#endif

              throw std::logic_error("Refiner::unrefineTheseElements couldn't remove element, destroy_entity returned false.");
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
              std::cout << "Refiner::unrefineTheseElements couldn't remove side element, destroy_entity returned false." << side_elem << std::endl;
              del =  m_eMesh.get_bulk_data()->destroy_entity( side_elem );
              if (1)
                {
                  percept::MyPairIterRelation rels (m_eMesh, side_elem,stk::mesh::MetaData::ELEMENT_RANK);
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
          //std::cout << "tmp remesh::rank[irank], num_elem_needed= " << m_ranks[irank] << " " << num_elem_needed << std::endl;
#endif


          // create new entities on this proc
          new_elements.resize(0);
          m_eMesh.createEntities( m_ranks[irank], num_elem_needed, new_elements);
          vector<stk::mesh::Entity>::iterator element_pool_it = new_elements.begin();

          // FIXME - we could directly call this with a refactor to change elementColors passed in here as a generic collection + checking for element Type
          //
          //createElementsAndNodesAndConnectLocal(m_ranks[irank], m_breakPattern[irank], elementColors, needed_entity_ranks, new_elements);

          for (ElementUnrefineCollection::iterator p_iter = parent_elements.begin();
               p_iter != parent_elements.end(); ++p_iter)
            {
              stk::mesh::Entity parent_p = *p_iter;
              const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(parent_p);
              CellTopology cell_topo(cell_topo_data);
              unsigned elementType = cell_topo.getKey();
              unsigned bpElementType = m_breakPattern[irank]->getFromTypeKey();
              if (elementType == bpElementType)
                {
                  stk::mesh::Entity parent = parent_p;

                  if (!m_eMesh.isGhostElement(parent))
                    {
#if DEBUG_UNREF_2
                      //std::cout << "P["<< m_eMesh.get_rank() << "] parent.owner_rank() = " << parent.owner_rank() << std::endl;
                      std::cout << "tmp Parent to be remeshed = ";
                      m_eMesh.print_entity(std::cout, parent);
#endif
                      if (createNewNeededNodeIds(cell_topo_data, parent, needed_entity_ranks, new_sub_entity_nodes))
                        {
                          //std::cout << "typeid= " << typeid(*breakPattern).name() << std::endl;
                          throw std::logic_error("unrefineTheseElements:: createNewNeededNodeIds failed");
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

      stk::mesh::FieldBase *refine_field = m_eMesh.get_field("refine_field");
      stk::mesh::FieldBase *refine_field_filtered = m_eMesh.get_field("refine_field_filtered");

      const unsigned FAMILY_TREE_RANK = stk::mesh::MetaData::ELEMENT_RANK + 1u;
      ElementUnrefineCollection elements_to_unref_copy(*m_eMesh.get_bulk_data());

      {
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
                  {
                    double *f_data = PerceptMesh::field_data_entity(refine_field, element);
                    //if (f_data) f_data[0] = -1;
                    if (f_data) f_data[0] = 10;
                  }
                }
            }
          }
        for (ElementUnrefineCollection::iterator u_iter = elements_to_unref.begin();
             u_iter != elements_to_unref.end(); ++u_iter)
          {
            stk::mesh::Entity element = *u_iter;

            {
              double *f_data = PerceptMesh::field_data_entity(refine_field, element);
              //if (f_data) f_data[0] = -1;
              if (f_data) f_data[0] = -30;
            }
          }
      }

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
            double *f_data = PerceptMesh::field_data_entity(refine_field, element);
             if (f_data) f_data[0] = -10;

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
                      //double *f_data = PerceptMesh::field_data_entity(refine_field, child);
                      //   if (f_data) f_data[0] = -20;
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
                            std::cout << "side_elem_family_tree_0= " << side_elem_family_tree_0.identifier() << std::endl;
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

      {
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
                  {
                    double *f_data = PerceptMesh::field_data_entity(refine_field_filtered, element);
                    //if (f_data) f_data[0] = -1;
                    if (f_data) f_data[0] = 10;
                  }
                }
            }
          }

        for (ElementUnrefineCollection::iterator u_iter = elements_to_unref.begin();
             u_iter != elements_to_unref.end(); ++u_iter)
          {
            stk::mesh::Entity element = *u_iter;

            {
              double *f_data = PerceptMesh::field_data_entity(refine_field_filtered, element);
              //if (f_data) f_data[0] = -1;
              if (f_data) f_data[0] = -10;
            }
          }
      }

    }

    void Refiner::
    getKeptNodes(NodeSetType& kept_nodes, ElementUnrefineCollection& elements_to_unref)
    {
      //bool doTest=true;

      stk::mesh::FieldBase *normal_kept_deleted = m_eMesh.get_field("normal_kept_deleted");

      {
        const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.node_rank() );

        for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            //if (removePartSelector(**k))
            {
              stk::mesh::Bucket & bucket = **k ;

              const unsigned num_entity_in_bucket = bucket.size();
              for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                {
                  stk::mesh::Entity node = bucket[ientity];
                  {
                    double *f_data = PerceptMesh::field_data_entity(normal_kept_deleted, node);
                    if (f_data) f_data[0] = 0;  // kept
                  }
                }
            }
          }
      }

      // mark kept nodes
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

                const percept::MyPairIterRelation elem_nodes (m_eMesh, element,  stk::mesh::MetaData::NODE_RANK);

                bool tt = (elem_nodes.size() && m_eMesh.isLeafElement(element));
                tt =  !m_eMesh.hasFamilyTree(element) || m_eMesh.numChildren(element) == 0 || !m_eMesh.hasGrandChildren(element, true);

                if (Util::getFlag(1256))
                  {
                    //tt = !m_eMesh.hasFamilyTree(element) || m_eMesh.numChildren(element) == 0;
                  }

                if (tt)
                  //if (!doTest || elem_nodes.size())
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
                            {
                              double *f_data = PerceptMesh::field_data_entity(normal_kept_deleted, node);
                              if (f_data) f_data[0] = 1;  // kept
                            }
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
      stk::mesh::FieldBase *normal_kept_deleted = m_eMesh.get_field("normal_kept_deleted");
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

          //if (!doTest ||
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
                      {
                        double *f_data = PerceptMesh::field_data_entity(normal_kept_deleted, node);
                        if (f_data) f_data[0] = -1;  // deleted
                      }
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

    static void merge(ElementUnrefineCollection& kept_nodes, ElementUnrefineCollection& kept_nodes_orig)
    {
      for (ElementUnrefineCollection::iterator it=kept_nodes_orig.begin(); it != kept_nodes_orig.end(); ++it)
        {
          kept_nodes.insert(*it);
        }
    }

    static void subtract(ElementUnrefineCollection& kept_nodes, ElementUnrefineCollection& kept_nodes_orig)
    {
      for (ElementUnrefineCollection::iterator it=kept_nodes.begin(); it != kept_nodes.end(); ++it)
        {
          if (kept_nodes_orig.find(*it) != kept_nodes_orig.end())
            kept_nodes_orig.erase(*it);
        }
    }

    void
    Refiner::
    unrefineTheseElements(ElementUnrefineCollection& elements_to_unref)
    {
      if (m_alwaysInitNodeRegistry)
        {
          throw std::logic_error("Refiner::unrefineTheseElements: to use urefinement, you must have setAlwaysInitializeNodeRegistry(false)");
        }

      //m_nodeRegistry->checkDB("unrefine start");

      m_eMesh.get_bulk_data()->modification_begin();

      int min_level=0, max_level=0;
      ElementUnrefineCollection elements_to_unref_0 = elements_to_unref;

      if (getDoLevelBasedUnrefinement())
        getLevels(m_eMesh, elements_to_unref, min_level, max_level);
      for (int ilevel=max_level; ilevel >= min_level; --ilevel)
        {
          if (DEBUG_UNREF_1) std::cout << "ilevel= " << ilevel << std::endl;
          elements_to_unref = elements_to_unref_0;
          std::cout << "tmp srk ilevel= " << ilevel << "  elements_to_unref.size= " << elements_to_unref.size() << std::endl;
          if (getDoLevelBasedUnrefinement())
            {
              std::cout << "tmp srk before filterForLevel elements_to_unref.size= " << elements_to_unref.size() << std::endl;
              filterForLevel(m_eMesh, elements_to_unref, ilevel);
              std::cout << "tmp srk after filterForLevel elements_to_unref.size= " << elements_to_unref.size() << std::endl;
            }

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

          NodeSetType kept_nodes(*m_eMesh.get_bulk_data());
          NodeSetType kept_nodes_orig(*m_eMesh.get_bulk_data());
          NodeSetType kept_nodes_orig_minus_kept_nodes(*m_eMesh.get_bulk_data());
          NodeSetType deleted_nodes(*m_eMesh.get_bulk_data());

          // get kept nodes - original
          getKeptNodes(kept_nodes_orig, elements_to_unref);
          std::cout << "kept_nodes_orig.size= " << kept_nodes_orig.size() << std::endl;
          // filter unref set
          filterUnrefSet(elements_to_unref);

          //           if (Util::getFlag(1256))
          //             {
          //               m_eMesh.save_as("filtered.e");
          //               //exit(1);
          //             }

          //            if (Util::getFlag(1257))
          //              {
          //                m_eMesh.save_as("filtereda.e");
          //                exit(1);
          //              }

          if (0)
            {
              elements_to_unref.clear();
              const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.element_rank() );
              for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
                {
                  stk::mesh::Bucket & bucket = **k ;
                  const unsigned num_entity_in_bucket = bucket.size();
                  for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                    {
                      stk::mesh::Entity element = bucket[ientity];
                      if (m_eMesh.hasFamilyTree(element) && m_eMesh.isChildWithoutNieces(element))
                        elements_to_unref.insert(element);
                    }
                }
            }
          // get kept nodes
          getKeptNodes(kept_nodes, elements_to_unref);
          std::cout << "kept_nodes_orig.size= " << kept_nodes_orig.size() << " kept_nodes.size= " << kept_nodes.size() << std::endl;

          // get deleted nodes
          getDeletedNodes(deleted_nodes, kept_nodes, elements_to_unref);
          unsigned deleted_nodes_size = deleted_nodes.size();
          std::cout << "kept_nodes_orig.size= " << kept_nodes_orig.size() << " before subtract kept_nodes.size= " << kept_nodes.size()
                    << " deleted_nodes.size= " << deleted_nodes_size
                    << std::endl;

          if (Util::getFlag(1256))
            {
              m_eMesh.save_as("filtereda.e");

              if (1)
                {

                  elements_to_unref = elements_to_unref_0;

                  if (1)
                    {
                      preUnrefine(elements_to_unref, m_eMesh.element_rank());
                      elements_to_unref.clear();
                      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.element_rank() );
                      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
                        {
                          stk::mesh::Bucket & bucket = **k ;
                          const unsigned num_entity_in_bucket = bucket.size();
                          for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                            {
                              stk::mesh::Entity element = bucket[ientity];
                              //if (m_eMesh.hasFamilyTree(element) && m_eMesh.isChildWithoutNieces(element))
                              if (m_eMesh.hasFamilyTree(element) && !m_eMesh.isParentElement(element))
                                elements_to_unref.insert(element);
                            }
                        }
                    }

                  filterUnrefSet(elements_to_unref);
                  kept_nodes.clear();
                  deleted_nodes.clear();
                  getKeptNodes(kept_nodes, elements_to_unref);
                  getDeletedNodes(deleted_nodes, kept_nodes, elements_to_unref);

                  m_eMesh.save_as("filteredb.e");
                }

              exit(1);
            }

          if (0)
            {
              if (0) merge(kept_nodes, kept_nodes_orig);
              // kept_nodes_orig -= kept_nodes
              if (0) subtract(kept_nodes, kept_nodes_orig);
              kept_nodes_orig.clear();
              deleted_nodes.clear();
              getDeletedNodes(deleted_nodes, kept_nodes, elements_to_unref);
              std::cout << "kept_nodes_orig.size= " << kept_nodes_orig.size() << " after subtract kept_nodes.size= " << kept_nodes.size()
                        << " deleted_nodes.size= " << deleted_nodes.size()
                        << std::endl;
            }

//           if (Util::getFlag(1257))
//             {
//               m_eMesh.save_as("filtereda.e");
//               exit(1);
//             }

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
          m_nodeRegistry->cleanDeletedNodes(deleted_nodes, kept_nodes_orig, to_save);
          std::cout << "to_save.size= " << to_save.size() << std::endl;

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
          removeChildElements(children_to_be_removed, &elements_to_unref_0);

          // FIXME - this is to be deprecated, but leaving it in for now
#if 0
          // for sideset elements, remove their family trees first
          removeFamilyTrees(side_elem_family_trees_to_be_removed);

          // remove sideset elements
          removeSideElements(side_elem_set_to_be_removed, children_to_be_removed);
#endif

          // reconnect and remove any dangling side elements
          fix_side_sets_2(true);

          getSideParentsToBeRemeshed(parent_elements, parent_side_elements);

          // remesh the holes left by removing child elems (quad/hex hanging node doesn't need this)
          if (m_needsRemesh)
            {
              remesh(parent_elements);
              remesh(parent_side_elements);
            }
//           if (Util::getFlag(1256))
//             {
//               m_eMesh.save_as("filtered.e");
//               exit(1);
//             }

#if CHECK_DEBUG
          check_db("after unrefineTheseElements, b4 mod end");
#endif

          // remove any elements that are empty (these can exist when doing local refinement)
          removeEmptyElements();

          if (0)
          {
            // 1. add kept_nodes_orig back into NodeRegistry
            // 2. remove kept_nodes_orig from deleted_nodes
            // 3. remesh
            // 4.

            SubDimCellToDataMap& map = m_nodeRegistry->getMap();
            for (SubDimCellToDataMap::iterator it=to_save.begin(); it != to_save.end(); ++it)
              {
                map[it->first]=it->second;
              }
            SetOfEntities current_elements(*m_eMesh.get_bulk_data());
            const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.element_rank() );
            for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_entity_in_bucket = bucket.size();
                for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                  {
                    stk::mesh::Entity element = bucket[ientity];
                    if (!m_eMesh.hasFamilyTree(element) || m_eMesh.numChildren(element) == 0)
                      {
                        current_elements.insert(element);
                      }
                  }
              }
            std::cout << "current_elements.size= " << current_elements.size() << " all elem= " << m_eMesh.get_number_elements() << std::endl;
            remesh(current_elements);
            std::cout << "after remesh current_elements.size= " << current_elements.size() << " all elem= " << m_eMesh.get_number_elements() << std::endl;
          }

          removeDeletedNodes(deleted_nodes);

          set_active_part();

          //remove_dangling_sidesets();
          if (Util::getFlag(1256))
            {
              m_eMesh.save_as("filtered.e");
              //exit(1);
            }

          fix_side_sets_2();

          if (deleted_nodes_size != 0)
           break;

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
      m_eMesh.get_bulk_data()->modification_end();

      //check_sidesets_2(" unrefineTheseElements:: end");

#if CHECK_DEBUG
      check_db("after unrefineTheseElements");
#endif

    }

    // ======================================================================================================================================================
    // ======================================================================================================================================================
    // ======================================================================================================================================================

    void Refiner::
    remeshGrandParents(SetOfEntities& grandParents)
    {
      const unsigned FAMILY_TREE_RANK = stk::mesh::MetaData::ELEMENT_RANK + 1u;
      unsigned rank = m_eMesh.element_rank();

      for (SetOfEntities::iterator gpIter = grandParents.begin(); gpIter != grandParents.end(); ++gpIter)
        {
          stk::mesh::Entity grandParent = *gpIter;

          //std::cout << "\n grandParent= " << grandParent.identifier() << std::endl;

          std::vector<stk::mesh::Entity> children;
          m_eMesh.getChildren(grandParent, children, true, false);
          VERIFY_OP_ON(children.size(), !=, 0, "hmmm");
          SetOfEntities familyTrees(*m_eMesh.get_bulk_data());
          SetOfEntities childrenOrGrandChildrenSet(children.begin(), children.end(), *m_eMesh.get_bulk_data());
          SetOfEntities childrenSet(children.begin(), children.end(), *m_eMesh.get_bulk_data());

          for (unsigned ichild=0; ichild < children.size(); ichild++)
            {
              std::vector<stk::mesh::Entity> grandChildren;
              m_eMesh.getChildren(children[ichild], grandChildren, true, false);
              if (grandChildren.size())
                {
                  childrenOrGrandChildrenSet.insert(grandChildren.begin(), grandChildren.end());
                }
            }

          for (SetOfEntities::iterator icgp = childrenOrGrandChildrenSet.begin(); icgp != childrenOrGrandChildrenSet.end(); ++icgp)
            {
              stk::mesh::Entity childOrGrandChild = *icgp;
              //std::cout << "childOrGrandChild= " << childOrGrandChild.identifier() << " isGranchChild= " << (childrenSet.find(childOrGrandChild) == childrenSet.end()) << std::endl;

              unsigned child_ft_level_0 = m_eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, childOrGrandChild);
              const percept::MyPairIterRelation child_to_family_tree_relations (m_eMesh, childOrGrandChild,FAMILY_TREE_RANK);
              if (child_to_family_tree_relations.size() == 0)
                {
                  throw std::logic_error("Refiner::preUnrefine child_to_family_tree_relations.size == 0");
                }
              stk::mesh::Entity family_tree = child_to_family_tree_relations[child_ft_level_0].entity();
              percept::MyPairIterRelation family_tree_relations (m_eMesh, family_tree,rank);
              if (family_tree_relations.size() == 0)
                {
                  throw std::logic_error("Refiner::preUnrefine family_tree_relations.size() == 0");
                }

              familyTrees.insert(family_tree);
            }

          removeFamilyTrees(familyTrees);
          removeChildElements(childrenOrGrandChildrenSet);

          SetOfEntities grandParent_set_of_one(&grandParent, (&grandParent)+1, *m_eMesh.get_bulk_data());
          remesh(grandParent_set_of_one);
          std::vector<stk::mesh::Entity> gpChildren;
          m_eMesh.getChildren(grandParent, gpChildren, false, false);
          for (unsigned igpc=0; igpc < gpChildren.size(); igpc++)
            {
              stk::mesh::Entity gpc = gpChildren[igpc];
              SetOfEntities gpc_set_of_one(&gpc, (&gpc)+1, *m_eMesh.get_bulk_data());
              remesh(gpc_set_of_one);
            }

        }

      m_nodeRegistry->clear_element_owner_data_phase_2(false, false);
      SetOfEntities emptySet(*m_eMesh.get_bulk_data());
      replaceNodeRegistryOwnership(emptySet, rank);
      m_nodeRegistry->clear_element_owner_data_phase_2(true, false);

      set_active_part();

    }

    void Refiner::
    filterUnrefSetPass2(ElementUnrefineCollection& elements_to_unref,   SetOfEntities& grandParents)
    {
      int print_filter_info = 1 || DEBUG_UNREF_1;
      if (print_filter_info)  std::cout << "\n\nP["<< m_eMesh.get_rank() << "] filterUnrefSetPass2: initial elements_to_unref size = " << elements_to_unref.size() << std::endl;

      grandParents.clear();
      for (ElementUnrefineCollection::iterator it = elements_to_unref.begin(); it != elements_to_unref.end(); ++it)
        {
          stk::mesh::Entity element = *it;

          if (m_eMesh.hasFamilyTree(element))
            {
              stk::mesh::Entity gp = m_eMesh.getGrandParent(element, true);
              if (gp.is_valid())
                {
                  if (!m_eMesh.hasGreatGrandChildren(gp, true))
                    {
                      grandParents.insert(gp);
                    }
                }
            }
        }

      if (print_filter_info)
        std::cout << "\n\nP["<< m_eMesh.get_rank() << "] filterUnrefSetPass2: initial grandParents size = "
                  << grandParents.size() << std::endl;

      ElementUnrefineCollection new_set(*m_eMesh.get_bulk_data());
      SetOfEntities new_gp(*m_eMesh.get_bulk_data());
      unsigned count_not_in_set=0;
      for (SetOfEntities::iterator gpIter = grandParents.begin(); gpIter != grandParents.end(); ++gpIter)
        {
          stk::mesh::Entity grandParent = *gpIter;

          std::vector<stk::mesh::Entity> children;
          m_eMesh.getChildren(grandParent, children, true, false);
          VERIFY_OP_ON(children.size(), !=, 0, "hmmm");
          //SetOfEntities childrenOrGrandChildrenSet(children.begin(), children.end(), *m_eMesh.get_bulk_data());
          SetOfEntities childrenOrGrandChildrenSet(*m_eMesh.get_bulk_data());

          for (unsigned ichild=0; ichild < children.size(); ichild++)
            {
              std::vector<stk::mesh::Entity> grandChildren;
              m_eMesh.getChildren(children[ichild], grandChildren, true, false);
              if (grandChildren.size())
                {
                  childrenOrGrandChildrenSet.insert(grandChildren.begin(), grandChildren.end());
                }
              else
                {
                  childrenOrGrandChildrenSet.insert(children[ichild]);
                }
            }

          bool not_in_set = false;
          for (SetOfEntities::iterator icgp = childrenOrGrandChildrenSet.begin(); icgp != childrenOrGrandChildrenSet.end(); ++icgp)
            {
              stk::mesh::Entity childOrGrandChild = *icgp;
              VERIFY_OP_ON(m_eMesh.numChildren(childOrGrandChild), ==, 0, "hmm3");
              if (elements_to_unref.find(childOrGrandChild) == elements_to_unref.end())
                {
                  not_in_set = true;
                  break;
                }
            }
          if (not_in_set)
            {
              ++count_not_in_set;
              continue;
            }
          new_gp.insert(grandParent);

          for (SetOfEntities::iterator icgp = childrenOrGrandChildrenSet.begin(); icgp != childrenOrGrandChildrenSet.end(); ++icgp)
            {
              stk::mesh::Entity childOrGrandChild = *icgp;
              new_set.insert(childOrGrandChild);
            }
        }
      elements_to_unref = new_set;
      grandParents = new_gp;
      if (print_filter_info)
        std::cout << "\n\nP["<< m_eMesh.get_rank() << "] filterUnrefSetPass2: count_not_in_set = " << count_not_in_set << " final grandParents size = "
                  << grandParents.size() << std::endl;
      if (print_filter_info)  std::cout << "\n\nP["<< m_eMesh.get_rank() << "] filterUnrefSetPass2: final elements_to_unref size = " << elements_to_unref.size() << std::endl;

    }

    void
    Refiner::
    unrefinePass2(ElementUnrefineCollection& elements_to_unref)
    {
      std::cout << "\n\n\n unrefinePass2:: start \n\nVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV "<< std::endl;

      m_eMesh.get_bulk_data()->modification_begin();

      ElementUnrefineCollection elements_to_unref_0 = elements_to_unref;

      // all elements candidates for kept nodes
      ElementUnrefineCollection all_kept_nodes_elements(*m_eMesh.get_bulk_data());
      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.element_rank() );
      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          const unsigned num_entity_in_bucket = bucket.size();
          for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
            {
              stk::mesh::Entity element = bucket[ientity];
              if (!m_eMesh.hasFamilyTree(element) || m_eMesh.hasGrandChildren(element) || m_eMesh.hasGreatGrandChildren(element))
                all_kept_nodes_elements.insert(element);
            }
        }

      NodeSetType kept_nodes(*m_eMesh.get_bulk_data());
      NodeSetType kept_nodes_orig(*m_eMesh.get_bulk_data());
      NodeSetType deleted_nodes(*m_eMesh.get_bulk_data());

      // filter unref set pass 2 - only grand parents, their children and grand children
      SetOfEntities grandParents(*m_eMesh.get_bulk_data());
      unsigned elsizeb4 = elements_to_unref.size();
      filterUnrefSetPass2(elements_to_unref, grandParents);
      unsigned elsizeaf = elements_to_unref.size();

      // get kept nodes
      // FIXME - check conditions inside getKeptNodes
      //getKeptNodes(kept_nodes, all_kept_nodes_elements);
      //getKeptNodes(kept_nodes, grandParents);
      //getKeptNodes(kept_nodes, not_elements_to_unref);
      if (0)
        {
          for (SetOfEntities::iterator gpit=grandParents.begin(); gpit != grandParents.end(); ++gpit)
            {
              stk::mesh::Entity gp = *gpit;
              MyPairIterRelation gp_nodes(m_eMesh, gp, m_eMesh.node_rank());
              for (unsigned ii = 0; ii < gp_nodes.size(); ++ii)
                {
                  stk::mesh::Entity node = gp_nodes[ii].entity();
                  kept_nodes.insert(node);
                }
            }
        }

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
                          //   if (0)
                          //    {
                          //      double *f_data = PerceptMesh::field_data_entity(normal_kept_deleted, node);
                          //      if (f_data) f_data[0] = 1;  // kept
                          //    }
                        }
                    }
                }
            }
        }

      // get deleted nodes
      // FIXME - check conditions inside getDeletedNodes
      getDeletedNodes(deleted_nodes, kept_nodes, elements_to_unref);
      unsigned deleted_nodes_size = deleted_nodes.size();
      std::cout << "unrefinePass2:: deleted_nodes_size= " << deleted_nodes_size
                << " kept_nodes.size= " << kept_nodes.size()
                << " grandParents.size= " <<  grandParents.size()
                << " elsizeb4= " << elsizeb4
                << " elsizeaf= " << elsizeaf
                << std::endl;
      //if (deleted_nodes_size) exit(1);

      // remove deleted nodes and their associated sub-dim entities from NodeRegistry's db
      //   save kept_nodes_orig on NodeRegistry
      SubDimCellToDataMap to_save;
      m_nodeRegistry->cleanDeletedNodes(deleted_nodes, kept_nodes_orig, to_save);

      // remesh the holes left by removing child elems (quad/hex hanging node doesn't need this)
      remeshGrandParents(grandParents);

      // remove any elements that are empty (these can exist when doing local refinement)
      removeEmptyElements();

      removeDeletedNodes(deleted_nodes);

      set_active_part();

      fix_side_sets_2();

      m_eMesh.get_bulk_data()->modification_end();

      m_nodeRegistry->clear_element_owner_data_phase_2();

      m_eMesh.get_bulk_data()->modification_begin();
      set_active_part();
      m_eMesh.get_bulk_data()->modification_end();

      std::cout << "\n\n ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n\n\n unrefinePass2:: end "<< std::endl;

    }

  } // namespace adapt
} // namespace stk
