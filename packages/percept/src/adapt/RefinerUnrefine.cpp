// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <exception>
#include <fstream>
#include <set>
#include <typeinfo>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <adapt/Refiner.hpp>
#include <adapt/RefinerUtil.hpp>
#include <adapt/AdaptedMeshVerifier.hpp>
#include <adapt/FixSideSetsSelector.hpp>

#include <percept/MeshUtil.hpp>
#include <percept/RunEnvironment.hpp>
#include <percept/PerceptUtils.hpp>

#define DEBUG_UNREF 0
#define DEBUG_UNREF_1 0
#define DEBUG_UNREF_2 0
#define DEBUG_UNREF_3 0
#define DEBUG_RU 0

#define TIMING(code) code
#define TIMER(name) stk::diag::Timer timer ## name ( #name, rootTimer());  stk::diag::TimeBlock tbTimer ## name (timer ## name)
#define TIMER2(name,parentName) stk::diag::Timer timer ## name ( #name, timer ## parentName);  stk::diag::TimeBlock tbTimer ## name (timer ## name)

#define DO_DETAILED_TIMING 0
#if DO_DETAILED_TIMING
#define DTIMER(name) TIMER(name)
#define DTIMER2(name,parentName) TIMER2(name)
#else
#define DTIMER(name) do {} while(0)
#define DTIMER2(name,parentName) do {} while(0)
#endif

#include <stk_mesh/base/MeshUtils.hpp>

  namespace percept {
    using namespace percept;


    // ====================================================================================================
    // ====================================================================================================
    // ====================================================================================================
    bool use_idServer = false;

    void Refiner::
    update_node_registry()
    {
      SubDimCellToDataMap::iterator iter;
      SubDimCellToDataMap& map = m_nodeRegistry->getMap();
      SubDimCellToDataMap map_new;

      for (unsigned irank=0; irank < m_ranks.size(); irank++)
        {
          if (m_ranks[irank] == m_eMesh.element_rank())
            {
              vector<NeededEntityType> needed_entity_ranks;
              m_breakPattern[irank]->fillNeededEntities(needed_entity_ranks);

              const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( m_ranks[irank] );

              for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
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
                          if (m_eMesh.isGhostElement(element)) continue;

                          {
                            const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);

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

                                for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                                  {
                                    SubDimCell_SDCEntityType subDimEntity(&m_eMesh);
                                    m_nodeRegistry->getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd, cell_topo_data);

                                    SubDimCellData* nodeId_elementOwnerId_ptr = m_nodeRegistry->getFromMapPtr(subDimEntity);
                                    if (nodeId_elementOwnerId_ptr)
                                      {

                                        std::get<SDC_DATA_OWNING_ELEMENT_KEY>(*nodeId_elementOwnerId_ptr) = m_eMesh.entity_key(element);
                                        std::get<SDC_DATA_OWNING_SUBDIM_ORDINAL>(*nodeId_elementOwnerId_ptr) = static_cast<unsigned char>(iSubDimOrd + 1);
                                        std::get<SDC_DATA_OWNING_SUBDIM_RANK>(*nodeId_elementOwnerId_ptr) = static_cast<unsigned char>(needed_entity_rank);

                                        map_new[subDimEntity] = *nodeId_elementOwnerId_ptr;

                                      }
                                  }
                              }
                          }
                        }
                    }
                }
            }
        }
      map = map_new;
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

              const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( rank );

              for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
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
                            // if found an element that is not deleted, make it own the subDimEntity
                            if (elements_to_delete.find(element) == elements_to_delete.end())
                              {
                                if (!m_eMesh.isGhostElement(element))
                                  refineMethodApply(&NodeRegistry::replaceElementOwnership, element, needed_entity_ranks,cell_topo_data);
                              }
                          }
                      }
                  }
                }
            }
        }
    }


    void Refiner::removeFamilyTrees(SetOfEntities& family_trees_to_be_removed)
    {
      for(SetOfEntities::iterator family_tree_it = family_trees_to_be_removed.begin();
          family_tree_it != family_trees_to_be_removed.end(); ++family_tree_it)
        {
          stk::mesh::Entity family_tree = *family_tree_it;
          if ( ! m_eMesh.get_bulk_data()->destroy_entity( family_tree ) )
            {
              throw std::logic_error("Refiner::removeFamilyTrees couldn't remove element, destroy_entity returned false for familytree.");
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



    template<class SetOfEntities>
    static void print_set(PerceptMesh& eMesh, SetOfEntities& side_set, const std::string& msg = "sides_to_remove.size")
    {
      if (DEBUG_RU)
        {
          std::cout << "print_set: " << msg << " side_set.size= " << side_set.size();
          for (typename SetOfEntities::iterator it_side=side_set.begin(); it_side != side_set.end(); ++it_side)
            {
              std::cout << " " << eMesh.id(*it_side)
                // << " E: " << eMesh.print_entity_compact(*it_side)
                        << " nc: " << eMesh.numChildren(*it_side);
              ;
            }
          std::cout << std::endl;
        }
    }

    void Refiner::getSideParentsToBeRemeshed(SetOfEntities& parents_to_be_remeshed, SetOfEntities& parent_side_elements, bool newVersion, SetOfEntities *avoid_sides)
    {
      if (getIgnoreSideSets()) return;

      SetOfEntities neighbors;
      SetOfEntities::iterator p_begin = parents_to_be_remeshed.begin(), p_end = parents_to_be_remeshed.end();

      for(SetOfEntities::iterator parent_it = p_begin;
          parent_it != p_end; ++parent_it)
        {
          stk::mesh::Entity parent = *parent_it;

          if (1)
            {
              neighbors.clear();
              m_eMesh.get_node_neighbors(parent, neighbors, m_eMesh.side_rank());
              SetOfEntities::iterator n_begin = neighbors.begin();
              SetOfEntities::iterator n_end = neighbors.end();
              for (SetOfEntities::iterator neigh_it = n_begin; neigh_it != n_end; ++neigh_it)
                {
                  stk::mesh::Entity side = *neigh_it;
                  if (avoid_sides && avoid_sides->find(side) != avoid_sides->end())
                    continue;
                  int permIndex = -1, permPolarity = 0, iside = -1;
                  bool sc = m_eMesh.should_connect(side, parent, &permIndex, &permPolarity, &iside);
                  //VERIFY_OP_ON(sc, ==, true, "bad sc");
                  if (sc && m_eMesh.numChildren(side) == 0)
                    {
                      parent_side_elements.insert(side);
                    }
                }
            }
        }
    }

    void Refiner::removeChildElements(SetOfEntities& children_to_be_removed, ElementUnrefineCollection *elements_to_unref_0)
    {
      const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(stk::topology::ELEMENT_RANK + 1u);
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

    // for the given set of parent elements, remesh the interior after the children have been removed and any deleted nodes have been removed from the NodeRegistry
    void Refiner::remesh(stk::mesh::Entity parent_element)
    {
      DTIMER(remesh);

      bool debug = false;

      // FIXME for performance
      //static NewSubEntityNodesType s_new_sub_entity_nodes(percept::EntityRankEnd);
      NewSubEntityNodesType s_new_sub_entity_nodes(percept::EntityRankEnd);

      NewSubEntityNodesType& new_sub_entity_nodes = s_new_sub_entity_nodes;

      VERIFY_OP_ON(m_ranks.size(), > , 0, "Logic/configuration error, m_ranks.size=0");

      for (unsigned irank=0; irank < m_ranks.size(); irank++)
        {
          // unsigned num_new_elem_during_remesh = 0;
          vector<NeededEntityType> needed_entity_ranks;
          m_breakPattern[irank]->fillNeededEntities(needed_entity_ranks);
          stk::mesh::Selector& fromPartsSelector = m_fromPartsSelector[irank];

          vector<stk::mesh::Entity> new_elements, ft_new_elements;

          // count num new elements needed on this proc (served by UniformRefinerPattern)
          unsigned num_elem_not_ghost = 0u;

          if (1)
            {
              stk::mesh::Entity element = parent_element;
              const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
              CellTopology cell_topo(cell_topo_data);
              unsigned elementType = cell_topo.getKey();
              unsigned bpElementType = m_breakPattern[irank]->getFromTypeKey();
              if (elementType == bpElementType && fromPartsSelector(m_eMesh.bucket(element)))
                {
                  if (!m_eMesh.isGhostElement(element))
                    {
                      ++num_elem_not_ghost;
                    }
                }
            }

          unsigned num_elem_needed = num_elem_not_ghost * m_breakPattern[irank]->getNumNewElemPerElem();

          if (debug) std::cout << "P[" << m_eMesh.get_rank() << " tmp remesh::rank[" << irank << "] = " << m_ranks[irank] << " , num_elem_needed= " << num_elem_needed << std::endl;

#if DEBUG_UNREF
          //std::cout << "tmp remesh::rank[irank], num_elem_needed= " << m_ranks[irank] << " " << num_elem_needed << std::endl;
#endif


          // create new entities on this proc
          new_elements.clear();
          ft_new_elements.clear();

          const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(stk::topology::ELEMENT_RANK + 1u);
          if (1)
            {
              if (UniformRefinerPatternBase::USE_DECLARE_ELEMENT_SIDE && m_ranks[irank] == m_eMesh.side_rank())
                {
                  new_elements.resize(0);
                }
              else
                {
                  if (!m_eMesh.getEntitiesUsingIdServer(m_ranks[irank], num_elem_needed, new_elements))
                    {
                      std::cout << "P[" << m_eMesh.get_rank() << "] num_elem_needed= " << num_elem_needed << " deplenished " << " rank= " << m_ranks[irank] << std::endl;
                      throw std::logic_error("entity pool deplenished");
                    }
                }
              if (!m_eMesh.getEntitiesUsingIdServer(FAMILY_TREE_RANK, num_elem_needed, ft_new_elements))
                {
                  std::cout << "P[" << m_eMesh.get_rank() << "] FAMILY_TREE_RANK num_elem_needed= " << num_elem_needed << " deplenished " << std::endl;
                  throw std::logic_error("entity pool deplenished");
                }
            }

          vector<stk::mesh::Entity>::iterator element_pool_it = new_elements.begin();
          vector<stk::mesh::Entity>::iterator ft_element_pool_it = ft_new_elements.begin();

          // FIXME - we could directly call this with a refactor to change elementColors passed in here as a generic collection + checking for element Type
          //
          //createElementsAndNodesAndConnectLocal(m_ranks[irank], m_breakPattern[irank], elementColors, needed_entity_ranks, new_elements);

          if (1)
            {
              stk::mesh::Entity parent_p = parent_element;
              const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(parent_p);
              CellTopology cell_topo(cell_topo_data);
              unsigned elementType = cell_topo.getKey();
              unsigned bpElementType = m_breakPattern[irank]->getFromTypeKey();

              if (DEBUG_UNREF)
                {
                  std::string strOld;
                  stk::mesh::PartVector pvOld = m_eMesh.bucket(parent_p).supersets();
                  for (unsigned ii=0; ii < pvOld.size(); ++ii)
                    strOld += pvOld[ii]->name()+" ";

                  std::ostringstream fromP;
                  for (unsigned i_part = 0; i_part < m_breakPattern[irank]->getFromParts().size(); i_part++)
                    {
                      fromP << "i_part = " << i_part << " m_fromParts= " << m_breakPattern[irank]->getFromParts()[i_part]->name()
                            << " is member= " << m_eMesh.bucket(parent_p).member(*m_breakPattern[irank]->getFromParts()[i_part])
                            << std::endl;
                    }


                  auto breakPattern = m_breakPattern[irank];
                  std::cout << "tmp remesh:: irank = " << irank
                            << " parent= " << m_eMesh.identifier(parent_p)
                            << " topo= " << cell_topo_data->name
                            << PerceptMesh::demangle(typeid(*breakPattern).name())
                            << " elementType= " << elementType
                            << " fromTopo= " << m_breakPattern[irank]->getFromTopology()->name
                            << " toTopo= " << m_breakPattern[irank]->getToTopology()->name
                            << " fromPartsSelector(m_eMesh.bucket(parent_p))= " << fromPartsSelector(m_eMesh.bucket(parent_p))
                            << "\n parts= " << strOld
                            << "\n fromP= " << fromP.str()
                            << std::endl;
                }


              if (elementType == bpElementType && fromPartsSelector(m_eMesh.bucket(parent_p)))
                {
                  stk::mesh::Entity parent = parent_p;

                  if (!m_eMesh.isGhostElement(parent))
                    {
#if DEBUG_UNREF_2
                      //std::cout << "P["<< m_eMesh.get_rank() << "] eMesh.owner_rank(parent) = " << eMesh.owner_rank(parent) << std::endl;
                      std::cout << "tmp Parent to be remeshed = ";
                      m_eMesh.print_entity(std::cout, parent);
#endif

                      {
                        //DTIMER2(remesh_createNewNeededNodeIds,remesh);
                        if (createNewNeededNodeIds(cell_topo_data, parent, needed_entity_ranks, new_sub_entity_nodes, m_breakPattern[irank]))
                          {
                            //std::cout << "typeid= " << typeid(*breakPattern).name() << std::endl;
                            throw std::logic_error("remesh:: createNewNeededNodeIds failed");
                          }
                      }
                      vector<stk::mesh::Entity>::iterator element_pool_it_b4 = element_pool_it;
                      {
                        //DTIMER2(remesh_createNewElementsCall,remesh);
                        m_breakPattern[irank]->createNewElements(m_eMesh, *m_nodeRegistry, parent, new_sub_entity_nodes, element_pool_it, ft_element_pool_it, m_proc_rank_field);
                      }
                      vector<stk::mesh::Entity>::iterator element_pool_it_af = element_pool_it;
                      unsigned nnew = (element_pool_it_af - element_pool_it_b4);
                      if (debug)
                        {
                          std::cout << "P[" << m_eMesh.get_rank() << " tmp remesh nnew= " << nnew << std::endl;
                        }
                      // num_new_elem_during_remesh += nnew;
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

      const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( stk::topology::ELEMENT_RANK );

      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
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

                const percept::MyPairIterRelation elem_nodes (m_eMesh, element,stk::topology::NODE_RANK);

                if (elem_nodes.size() && (m_eMesh.hasFamilyTree(element) && m_eMesh.isChildWithoutNieces(element, false) ) )
                  {
                    elements_to_unref.insert(element);
                  }
              }
          }
        }
      unrefineTheseElements(elements_to_unref);
    }

    void
    Refiner::
    unrefineTheseElements(ElementUnrefineCollection& elements_to_unref)
    {
      unrefinePass2(elements_to_unref);
    }

    // ======================================================================================================================================================
    // ======================================================================================================================================================
    // ======================================================================================================================================================

    struct Comp
    {
      PerceptMesh& mesh;
      Comp(PerceptMesh& m) : mesh(m) {}
      bool operator()(const stk::mesh::Entity& a, const stk::mesh::Entity& b) { return mesh.id(a) < mesh.id(b); }
    };

    // Pass2

    void Refiner::
    remeshRecurse(stk::mesh::Entity element, int& s_depth)
    {
      if (s_depth > 20)
        {
          std::cout << "P[" << m_eMesh.get_rank() << "] WARNING: remeshRecurse recursion depth > 20 " << std::endl;
        }
      if (s_depth > 100)
        {
          std::cout << "P[" << m_eMesh.get_rank() << "] WARNING: remeshRecurse recursion depth > 100 " << std::endl;
          VERIFY_MSG("recursion depth too deep");
        }

      // Note: we rely on the remeshing returning no children if the element isn't marked for refinement
      //SetOfEntities element_set_of_one(&element, (&element)+1, *m_eMesh.get_bulk_data());

      if (s_depth > 20)
        {
          std::cout << "remeshing element= " << m_eMesh.identifier(element) << " s_depth= " << s_depth << " " << m_eMesh.print_entity_compact(element) << std::endl;
        }

      remesh(element);

      if (m_eMesh.hasFamilyTree(element))
        {
          std::vector<stk::mesh::Entity> children;
          m_eMesh.getChildren(element, children, true, false);

          if (children.size())
            ++s_depth;

          for (unsigned ichild=0; ichild < children.size(); ichild++)
            {
              remeshRecurse(children[ichild], s_depth);
            }
        }
    }

    void Refiner::
    remeshElements(SetOfEntities& rootElements, stk::mesh::EntityRank rank, int pool_size_hint, SetOfEntities *elemsToBeDeleted)
    {
      DTIMER(remeshElements);

      m_fromPartsSelector.resize(m_ranks.size());
      for (unsigned irank=0; irank < m_ranks.size(); irank++)
        {
          m_fromPartsSelector[irank] = stk::mesh::selectUnion( m_breakPattern[irank]->getFromParts() );
        }

      if (DO_MEMORY) {
        std::string hwm = print_memory_high_water_mark(m_eMesh.parallel());
        if (!m_eMesh.get_rank()) std::cout << "MEM: " << hwm << " unrefiner: remeshElements= " << std::endl;
      }

      {
        DTIMER2(recurse, remeshElements);

        for (SetOfEntities::iterator elIter = rootElements.begin(); elIter != rootElements.end(); ++elIter)
          {
            stk::mesh::Entity element = *elIter;
            VERIFY_OP_ON(m_eMesh.numChildren(element), ==, 0, "hmmm");

            int s_depth=0;
            remeshRecurse(element, s_depth);
          }
      }

      if (m_nodeRegistry->s_use_new_ownership_check)
        {
          SetOfEntities emptySet(*m_eMesh.get_bulk_data());
          if (rank == m_eMesh.element_rank())
            {
              update_node_registry();
            }
        }
      else
      {
        DTIMER2(clear, remeshElements);
        m_nodeRegistry->clear_element_owner_data_phase_2(false, false, elemsToBeDeleted);

        SetOfEntities emptySet(*m_eMesh.get_bulk_data());
        if (rank == m_eMesh.element_rank())
          {
            VERIFY_OP_ON(elemsToBeDeleted, ==, 0, "bad");
            if (elemsToBeDeleted)
              replaceNodeRegistryOwnership(*elemsToBeDeleted, rank);
            else
              replaceNodeRegistryOwnership(emptySet, rank);
          }
        m_nodeRegistry->clear_element_owner_data_phase_2(true, false, elemsToBeDeleted);
      }

      {
        DTIMER2(set_active_part, remeshElements);
        set_active_part();
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
    filterRecurse(stk::mesh::Entity root_element, ElementUnrefineCollection& rootElements, ElementUnrefineCollection& elements_to_unref)
    {

      if (1)
        {
          SetOfEntities allD(*m_eMesh.get_bulk_data());
          bool only_leaves=false;
          unsigned nlevels=0;
          //allDescendants(root_element, allD, nlevels, only_leaves);
          //bool allIn = true;
          bool allIn = (DEBUG_UNREF_3 ? allDescendants(root_element, allD, nlevels, only_leaves) :
                        allDescendants(root_element, allD, nlevels, only_leaves, &elements_to_unref) );
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
              rootElements.insert(root_element);
              return;
            }
        }

      if (m_eMesh.hasFamilyTree(root_element))
        {
          std::vector<stk::mesh::Entity> children;
          m_eMesh.getChildren(root_element, children, true, false);
          for (unsigned ichild=0; ichild < children.size(); ichild++)
            {
              filterRecurse(children[ichild], rootElements, elements_to_unref);
            }
        }
    }

    /** given a list of elements to unrefine, @param elements_to_unref, and an empty list, @param elements_to_be_remeshed,
     *    find a set of elements whose descendants are all in the list of elements_to_unref and who are parents
     *    of at least one generation; put this new list in @param elements_to_be_remeshed, which are then remeshed after
     *    their children are deleted.  Also, delete any elements in elements_to_unref that are in the final
     *    elements_to_be_remeshed set.
     *
     * assertion: after processing, elements_to_be_remeshed should have each element have all descendants in elements_to_unref
     *
     */

    void Refiner::
    filterUnrefSetPass2(ElementUnrefineCollection& elements_to_unref,   SetOfEntities& elements_to_be_remeshed)
    {
      stk::diag::Timer timerFilter_("filterUnrefSetPass2", rootTimer());
      stk::diag::TimeBlock tbTimerFilter_(timerFilter_);

      int print_filter_info = DEBUG_UNREF_1;
      if (print_filter_info)  std::cout << "\n\nP["<< m_eMesh.get_rank() << "] filterUnrefSetPass2: initial elements_to_unref size = " << elements_to_unref.size() << std::endl;

      ElementUnrefineCollection new_elements_to_unref(*m_eMesh.get_bulk_data());
      ElementUnrefineCollection final_root_elements_set(*m_eMesh.get_bulk_data());

        {
          // first, find the true root elements (top of each tree, or have no children), and recurse down to find
          //   the elements that have all descendants in the list to be unref'd
          const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( stk::topology::ELEMENT_RANK );
          for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_entity_in_bucket = bucket.size();
              for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                {
                  stk::mesh::Entity element = bucket[ientity];
                  bool has_family_tree = m_eMesh.hasFamilyTree(element);
                  bool is_child = has_family_tree && m_eMesh.isChildElement(element,true);
                  bool is_root = !has_family_tree || !is_child;

                  if (m_onlyOneLevelUnrefine)
                    {
                      if (has_family_tree && m_eMesh.hasGrandChildren(element, true))
                        {
                          continue;
                        }
                      if (has_family_tree && m_eMesh.isParentElement(element))
                        filterRecurse(element, final_root_elements_set, elements_to_unref);
                    }
                  else if (is_root)
                    {
                      if (m_eMesh.hasFamilyTree(element))
                        VERIFY_OP_ON(0, !=, m_eMesh.numChildren(element),"hmmm");

                      filterRecurse(element, final_root_elements_set, elements_to_unref);
                    }
                }
            }

          for (ElementUnrefineCollection::iterator rit=final_root_elements_set.begin(); rit != final_root_elements_set.end(); ++rit)
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

                  final_root_elements_set.insert(element);

                  for (SetOfEntities::iterator it=allD.begin(); it != allD.end(); ++it)
                    {
                      new_elements_to_unref.insert(*it);
                    }
                }
            }
          elements_to_be_remeshed = final_root_elements_set;
          elements_to_unref = new_elements_to_unref;
        }

      if (print_filter_info)  std::cout << "\n\nP["<< m_eMesh.get_rank() << "] filterUnrefSetPass2: final elements_to_unref size = " << elements_to_unref.size() << std::endl;

      // check
      if (DEBUG_UNREF)
      {
        for (ElementUnrefineCollection::iterator it=final_root_elements_set.begin(); it != final_root_elements_set.end(); ++it)
          {
            stk::mesh::Entity rootElement = *it;
            if (new_elements_to_unref.find(rootElement) != new_elements_to_unref.end())
              throw std::runtime_error("bad set");
          }
      }

    }

    void Refiner::
    get_kept_nodes(SetOfEntities& kept_nodes, ElementUnrefineCollection& elements_to_unref)
    {
      //TIMER2(unref_keptNodes,Unrefine_);

      const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.element_rank() );
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;

          const unsigned num_entity_in_bucket = bucket.size();
          for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
            {
              stk::mesh::Entity element = bucket[ientity];

              bool in_unref_set = elements_to_unref.find( element ) != elements_to_unref.end();
              if (!in_unref_set)
                {
                  const percept::MyPairIterRelation elem_nodes (m_eMesh, element,  stk::topology::NODE_RANK);

                  for (unsigned inode=0; inode < elem_nodes.size(); inode++)
                    {
                      stk::mesh::Entity node = elem_nodes[inode].entity();
                      kept_nodes.insert(node);
                      //std::cout << "P[" << m_eMesh.get_rank() << "] kept node= " << m_eMesh.identifier(node) << std::endl;
                    }
                }
            }
        }

    }

    void Refiner::
    get_deleted_nodes(SetOfEntities& deleted_nodes, SetOfEntities& kept_nodes, ElementUnrefineCollection& elements_to_unref, SetOfEntities& elements_to_be_remeshed)
    {
      if (!m_avoidClearDanglingNodes)
      {
        //TIMER2(getDeletedNodes,Unrefine_);

        for (ElementUnrefineCollection::iterator u_iter = elements_to_unref.begin();
             u_iter != elements_to_unref.end(); ++u_iter)
          {
            stk::mesh::Entity element = *u_iter;

            if (elements_to_be_remeshed.find(element) != elements_to_be_remeshed.end())
              continue;

            const percept::MyPairIterRelation elem_nodes (m_eMesh, element,  stk::topology::NODE_RANK);

            if (0)
              std::cout << "P[" << m_eMesh.get_rank() << "] deleted elem= " << m_eMesh.identifier(element)
                        << " isG = " << m_eMesh.isGhostElement(element) << " nod.siz=" << elem_nodes.size() << " hft= " << m_eMesh.hasFamilyTree(element)
                        << std::endl;

            if ( elem_nodes.size() && m_eMesh.hasFamilyTree(element) )
              {
                for (unsigned inode=0; inode < elem_nodes.size(); inode++)
                  {
                    stk::mesh::Entity node = elem_nodes[inode].entity();
                    bool in_kept_nodes_set = kept_nodes.find( node ) != kept_nodes.end();
                    if (!in_kept_nodes_set)
                      {
                        deleted_nodes.insert(node);
                        // std::cout << "P[" << m_eMesh.get_rank() << "] deleted node= " << m_eMesh.identifier(node) << std::endl;
                      }
                  }
              }
          }
      }
    }

    void Refiner::
    filter_deleted_nodes(SetOfEntities& deleted_nodes)
    {
      if (m_eMesh.get_spatial_dim() == 2)
        return;

      SubDimCellToDataMap& cell_2_data_map = m_nodeRegistry->getMap();
      stk::mesh::BulkData &bulk_data = *m_eMesh.get_bulk_data();
      SetOfEntities keep_nodes(bulk_data);

      for (SubDimCellToDataMap::iterator cell_iter = cell_2_data_map.begin(); cell_iter != cell_2_data_map.end(); ++cell_iter)
        {
          const SubDimCell_SDCEntityType& subDimEntity = (*cell_iter).first;
          SubDimCellData& nodeId_elementOwnerId = (*cell_iter).second;
          NodeIdsOnSubDimEntityType& nodeIds_onSE = std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnerId);
          stk::mesh::EntityId owning_elementId = std::get<SDC_DATA_OWNING_ELEMENT_KEY>(nodeId_elementOwnerId).id();
          stk::mesh::EntityRank owning_elementRank = std::get<SDC_DATA_OWNING_ELEMENT_KEY>(nodeId_elementOwnerId).rank();
          if (!nodeIds_onSE.size())
            continue;

          bool keep_it = subDimEntity.size() == 4;
          if (subDimEntity.size() == 1 && owning_elementId && owning_elementRank == m_eMesh.element_rank())
            {
              stk::mesh::Entity elem = bulk_data.get_entity(m_eMesh.element_rank(), owning_elementId);
              if (m_eMesh.topology(elem) == stk::topology::HEX_8)
                keep_it = true;
            }
          if (keep_it)
            {
              for (unsigned ii=0; ii < nodeIds_onSE.size(); ++ii)
                {
                  if (deleted_nodes.find(nodeIds_onSE[ii]) != deleted_nodes.end())
                    {
                      keep_nodes.insert(nodeIds_onSE[ii]);
                    }
                }
            }
        }

      for (SetOfEntities::iterator it = keep_nodes.begin(); it != keep_nodes.end(); ++it)
        {
          deleted_nodes.erase(*it);
        }
    }

    void Refiner::
    generate_temporary_elements(SetOfEntities& entities_to_delete, stk::mesh::EntityRank rank, std::vector<stk::mesh::Entity>& new_elements)
    {
      // add pseudo-dlements to ensure STK doesn't delete Field values on nodes
      if (!m_do_new_elements)
        return;

      size_t sz = 0;
      std::vector<stk::mesh::Entity> entities;
      for (SetOfEntities::iterator it=entities_to_delete.begin(); it != entities_to_delete.end(); ++it)
        {
          stk::mesh::Entity element = *it;
          const percept::MyPairIterRelation elem_nodes (m_eMesh, element, m_eMesh.node_rank());
          stk::mesh::Entity node0 = stk::mesh::Entity();

          for (unsigned ii=0; ii < elem_nodes.size(); ++ii)
            {
              if (!m_eMesh.aura(elem_nodes[ii].entity()))
                {
                  node0 = elem_nodes[ii].entity();
                  break;
                }
            }
          if (m_eMesh.is_valid(node0))
            {
              ++sz;
              entities.push_back(element);
            }
        }
      if (sz == 0)
        {
          new_elements.resize(0);
          return;
        }

      m_eMesh.getEntitiesUsingIdServer(rank, sz, new_elements);

      size_t ielem = 0;
      VERIFY_OP_ON(new_elements.size(), ==, sz, "bad size");
      for (size_t jel = 0; jel < entities.size(); ++jel)
        {
          stk::mesh::Entity element = entities[jel];
          VERIFY_OP_ON(ielem, <, sz, "bad index");
          stk::mesh::Entity newElement = new_elements[ielem];
          const stk::mesh::PartVector& super = m_eMesh.bucket(element).supersets();
          stk::mesh::PartVector add, rem;
          for (unsigned ipart = 0; ipart < super.size(); ++ipart)
            {
              stk::mesh::Part *  part = super[ipart];

              if ( stk::mesh::is_auto_declared_part(*part) )
                continue;
              bool is_auto_part = part->attribute<AutoPart>() != 0;
              if (is_auto_part)
                continue;
              add.push_back(part);
            }
          const percept::MyPairIterRelation elem_nodes (m_eMesh, element, m_eMesh.node_rank());
          stk::mesh::Entity node0 = stk::mesh::Entity();
          for (unsigned ii=0; ii < elem_nodes.size(); ++ii)
            {
              if (!m_eMesh.aura(elem_nodes[ii].entity()))
                {
                  node0 = elem_nodes[ii].entity();
                  break;
                }
            }
          VERIFY_OP_ON(m_eMesh.is_valid(node0), ==, true, "bad node0");

          // avoid temporary elements pointing to aura nodes, which causes a parallel mesh inconsitency with Parts
          for (unsigned ii=0; ii < elem_nodes.size(); ++ii)
            {
              if (m_eMesh.aura(elem_nodes[ii].entity()))
                m_eMesh.get_bulk_data()->declare_relation(newElement, node0, ii);
              else
                m_eMesh.get_bulk_data()->declare_relation(newElement, elem_nodes[ii].entity(), ii);
            }

          m_eMesh.get_bulk_data()->change_entity_parts( newElement, add, rem );
          ++ielem;
        }
    }

    void Refiner::
    delete_entities(SetOfEntities& entities_to_delete, stk::mesh::EntityRank rank)
    {
      SetOfEntities family_trees(*m_eMesh.get_bulk_data());

      const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(stk::topology::ELEMENT_RANK + 1u);

      for (ElementUnrefineCollection::iterator it=entities_to_delete.begin(); it != entities_to_delete.end(); ++it)
        {
          stk::mesh::Entity element = *it;

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

          if (!m_eMesh.get_bulk_data()->destroy_entity( element ) )
            {
              throw std::runtime_error("pass2 4.2 - couldn't delete element");
            }
        }

      for (SetOfEntities::iterator fiter=family_trees.begin(); fiter != family_trees.end(); ++fiter)
        {
          stk::mesh::Entity family_tree = *fiter;
          percept::MyPairIterRelation rels (m_eMesh, family_tree, rank);
          if (rels.size() == 1) // only points to a parent entity
            {
              if ( ! m_eMesh.get_bulk_data()->destroy_entity( family_tree ) )
                {
                  throw std::runtime_error("pass2 4.1 - couldn't delete family_tree");
                }
            }
        }
    }

    bool Refiner::include_side_bucket(stk::mesh::Bucket& side_bucket, stk::mesh::Selector *excludeSelector)
    {
      if (!bucket_acceptable(side_bucket, side_bucket.entity_rank() ))
        {
          return false;
        }

      if (excludeSelector && (*excludeSelector)(side_bucket))
        {
          return false;
        }
#if  defined(STK_PERCEPT_HAS_GEOMETRY)
      if (m_eMesh.is_in_geometry_parts(m_geomFile, side_bucket))
        {
          return false;
        }
#endif
      return true;
    }

    void Refiner::
    get_deleted_sides(SetOfEntities& sides_to_delete, ElementUnrefineCollection& elements_to_unref, SetOfEntities& elements_to_be_remeshed)
    {
#if 1
      SetOfEntities side_set(*m_eMesh.get_bulk_data());
      //build_side_set(side_set);
      sides_to_delete.clear();
      SetOfEntities neighbors;

      stk::mesh::Selector excludeSelector, *excludeSelectorPtr = 0;
      if (m_excludeParts.size())
        {
          excludeSelector = stk::mesh::selectUnion(m_excludeParts);
          excludeSelectorPtr = &excludeSelector;
        }


      stk::mesh::EntityRank side_rank_iter_begin = m_eMesh.side_rank();
      stk::mesh::EntityRank side_rank_iter_end = m_eMesh.side_rank();
      if (m_eMesh.get_spatial_dim() == 3)
        {
          side_rank_iter_begin = m_eMesh.edge_rank();
        }

      for (stk::mesh::EntityRank side_rank_iter = side_rank_iter_begin; side_rank_iter <= side_rank_iter_end; side_rank_iter++)
        {
          for (auto element : elements_to_unref)
            {
              neighbors.clear();
              m_eMesh.get_node_neighbors(element, neighbors, side_rank_iter);
              for (auto side : neighbors)
                {
                  if (include_side_bucket(m_eMesh.bucket(side), excludeSelectorPtr))
                    {
                      side_set.insert(side);
                    }
                }
            }
        }


      for (auto side : side_set)
        {
          if (m_eMesh.aura(side))
            continue;

          bool keep = false;
          if (1)
            {
              neighbors.clear();
              m_eMesh.get_node_neighbors(side, neighbors, m_eMesh.element_rank());
              SetOfEntities::iterator n_begin = neighbors.begin();
              SetOfEntities::iterator n_end = neighbors.end();
              for (SetOfEntities::iterator neigh_it = n_begin; neigh_it != n_end; ++neigh_it)
                {
                  stk::mesh::Entity element = *neigh_it;
                  // FIXME PERFORMANCE
                  bool notUnref = (elements_to_unref.find(element) == elements_to_unref.end());
                  if (notUnref)
                    {
                      int permIndex = -1, permPolarity = 0, iside = -1;
                      bool sc = m_eMesh.should_connect(side, element, &permIndex, &permPolarity, &iside);
                      if (sc)
                        {
                          keep = true;
                          break;
                        }
                    }
                }
            }
          if (!keep)
            {
              sides_to_delete.insert(side);
            }
        }
#else
      SetOfEntities side_set(*m_eMesh.get_bulk_data());
      build_side_set(side_set);
      sides_to_delete.clear();
      SetOfEntities neighbors;

      for (SetOfEntities::iterator it_side=side_set.begin(); it_side != side_set.end(); ++it_side)
        {
          stk::mesh::Entity side = *it_side;

          if (m_eMesh.aura(side))
            continue;

          bool keep = false;
          if (1)
            {
              neighbors.clear();
              m_eMesh.get_node_neighbors(side, neighbors, m_eMesh.element_rank());
              SetOfEntities::iterator n_begin = neighbors.begin();
              SetOfEntities::iterator n_end = neighbors.end();
              for (SetOfEntities::iterator neigh_it = n_begin; neigh_it != n_end; ++neigh_it)
                {
                  stk::mesh::Entity element = *neigh_it;
                  // FIXME PERFORMANCE
                  bool notUnref = (elements_to_unref.find(element) == elements_to_unref.end());
                  if (notUnref)
                    {
                      int permIndex = -1, permPolarity = 0, iside = -1;
                      bool sc = m_eMesh.should_connect(side, element, &permIndex, &permPolarity, &iside);
                      if (sc)
                        {
                          keep = true;
                          break;
                        }
                    }
                }
            }
          if (!keep)
            {
              sides_to_delete.insert(*it_side);
            }
        }
#endif
    }

    // verify parent and child are on same proc, including the FAMILY_TREE entities
    void Refiner::check_parent_ownership()
    {
      const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(stk::topology::ELEMENT_RANK + 1u);

      for (stk::mesh::EntityRank rank_iter = m_eMesh.side_rank(); rank_iter <= m_eMesh.element_rank(); rank_iter++)
        {
          const stk::mesh::BucketVector & entity_buckets = m_eMesh.get_bulk_data()->buckets( rank_iter );
          for ( stk::mesh::BucketVector::const_iterator it_entity_bucket = entity_buckets.begin() ; it_entity_bucket != entity_buckets.end() ; ++it_entity_bucket )
            {
              stk::mesh::Bucket & entity_bucket = **it_entity_bucket ;
              const unsigned num_elements_in_entity_bucket = entity_bucket.size();
              for (unsigned i_entity = 0; i_entity < num_elements_in_entity_bucket; i_entity++)
                {
                  stk::mesh::Entity entity = entity_bucket[i_entity];
                  percept::MyPairIterRelation entity_ft (m_eMesh, entity, FAMILY_TREE_RANK);
                  for (unsigned ii = 0; ii < entity_ft.size(); ++ii)
                    {
                      VERIFY_OP_ON(m_eMesh.owner_rank(entity_ft[ii].entity()), ==, m_eMesh.owner_rank(entity), "bad entity_ft ownership");
                    }

                  if (m_eMesh.aura(entity))
                    continue;
                  stk::mesh::Entity parent = stk::mesh::Entity();
                  if (m_eMesh.hasFamilyTree(entity))
                    parent = m_eMesh.getParent(entity, true);
                  if (m_eMesh.is_valid(parent))
                    {
                      VERIFY_OP_ON(m_eMesh.owner_rank(parent), ==, m_eMesh.owner_rank(entity), "bad parent ownership");
                    }
                }
            }
        }
    }

    bool Refiner::check_sides_on_same_proc_as_owned_element(const std::string& msg, bool doThrow)
    {
      SetOfEntities side_set(*m_eMesh.get_bulk_data());
      bool gfound = true;
      build_side_set(side_set);
      std::ostringstream str;
      for (SetOfEntities::iterator it_side=side_set.begin(); it_side != side_set.end(); ++it_side)
        {
          stk::mesh::Entity side = *it_side;
          stk::mesh::EntityId cid = 0; //3650;
          bool ldb = m_eMesh.id(side) == cid && m_eMesh.entity_rank(side) == m_eMesh.side_rank();

          if (m_eMesh.aura(side))  // FIXME !owned
            continue;

          // skip edges in 3D
          if (m_eMesh.get_spatial_dim() == 3 && m_eMesh.entity_rank(side) == m_eMesh.edge_rank())
            continue;

          const stk::mesh::Entity *side_elements = m_eMesh.get_bulk_data()->begin_elements(side);
          unsigned side_elements_size = m_eMesh.get_bulk_data()->num_elements(side);
          const stk::mesh::ConnectivityOrdinal *side_element_ords = m_eMesh.get_bulk_data()->begin_element_ordinals(side);

          bool found = false;
          int jj = -1;
          for (unsigned ii = 0; ii < side_elements_size; ++ii)
            {
              stk::mesh::Entity elem = side_elements[ii];
              bool sameOwner = m_eMesh.owner_rank(elem) == m_eMesh.owner_rank(side);
              if (sameOwner)
                {
                  bool isPos = m_eMesh.is_positive_perm(elem, side, side_element_ords[ii]);
                  if (ldb)
                    {
                      str << "P[" << m_eMesh.get_rank() << "] CS TMP srk elem= " << m_eMesh.id(elem) << " ii= " << ii << " side_elements.size= " << side_elements_size
                          << " sameOwner= " << sameOwner << " isPos= " << isPos
                          << std::endl;
                      m_eMesh.print(str,side,true,true);
                      m_eMesh.print(str,elem,true,true);
                    }

                  if (sameOwner && isPos)
                    {
                      jj = ii;
                      found = true;
                      break;
                    }
                }
            }
          if (doThrow && !found)
            {
              std::cerr << m_eMesh.rank() << " CS1 side= " << m_eMesh.id(side) << " msg: " << msg << " jj= " << jj
                        << " side_elements.size= "<< side_elements_size
                        << "\n" << str.str() << std::endl;
              m_eMesh.print(std::cerr, side, true, true);
              VERIFY_MSG(m_eMesh.rank()+" couldn't find an element/side on same proc, msg= " +msg +" side_elements.size= "+toString(side_elements_size)+" side= " +toString(m_eMesh.id(side)));
            }
          if (!found) gfound = false;
        }

      stk::all_reduce( m_eMesh.get_bulk_data()->parallel() , stk::ReduceMin<1>( &gfound ) );

      return !gfound;
    }

    void Refiner::
    require_sides_on_same_proc_as_pos_perm_element()
    {
      stk::diag::Timer timerAdapt_("RefineMesh", rootTimer());
      stk::diag::Timer timerInitRefine_("percept::InitRefine", timerAdapt_);
      stk::diag::Timer timerReqSides_("ReqSidesLoc", timerInitRefine_);
      stk::diag::TimeBlock timerReqSidesBlock_(timerReqSides_);

      const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(stk::topology::ELEMENT_RANK + 1u);
      SetOfEntities side_set(*m_eMesh.get_bulk_data());
      bool only_roots_of_tree = true;

      {
        stk::diag::Timer timerReqSidesBSS_("ReqSidesLocBSS", timerReqSides_);
        stk::diag::TimeBlock timerReqSidesBSSBlock_(timerReqSidesBSS_);
        build_side_set(side_set, only_roots_of_tree);
      }

      SetOfEntities allD(*m_eMesh.get_bulk_data());

      stk::mesh::EntityProcVec change;

      {
        stk::diag::Timer timerReqSidesBCh_("ReqSidesLocBCh", timerReqSides_);
        stk::diag::TimeBlock timerReqSidesBChBlock_(timerReqSidesBCh_);

      for (auto child_side : side_set)
        {
          stk::mesh::Entity side = child_side;
          if (!m_eMesh.owned(side))
            continue;

          stk::topology side_topo = m_eMesh.topology(side);

          stk::mesh::Entity const* side_elements = m_eMesh.get_bulk_data()->begin_elements( side);
          stk::mesh::ConnectivityOrdinal const* side_element_ordinals = m_eMesh.get_bulk_data()->begin_element_ordinals(side);
          unsigned side_elements_size = m_eMesh.get_bulk_data()->num_elements(side);

          bool found = false;
          int element_proc = -1;

          bool is3Dedge = (m_eMesh.get_spatial_dim() == 3 && m_eMesh.entity_rank(side) == m_eMesh.edge_rank());

          for (unsigned ii=0; ii < side_elements_size; ++ii)
            {
              stk::mesh::Entity element = side_elements[ii];

              stk::mesh::ConnectivityOrdinal side_ord = side_element_ordinals[ii];
              stk::mesh::Permutation perm = stk::mesh::INVALID_PERMUTATION;

#if 1
              if (m_eMesh.owner_rank(element) == m_eMesh.owner_rank(side))
                {
                  perm = m_eMesh.find_permutation(element, side, side_ord);
                  if (side_topo.is_positive_polarity(perm))
                    {
                      element_proc = m_eMesh.owner_rank(element);
                      found = true;
                      break;
                    }
                }
              else
                {
                  perm = m_eMesh.find_permutation(element, side, side_ord);
                  if (side_topo.is_positive_polarity(perm))
                    {
                      element_proc = m_eMesh.owner_rank(element);
                    }
                }

#else
              perm = m_eMesh.find_permutation(element, side, side_ord);
              if (side_topo.is_positive_polarity(perm))
                {
                  element_proc = m_eMesh.owner_rank(element);
                  if (m_eMesh.owner_rank(element) == m_eMesh.owner_rank(side))
                    {
                      element_proc = m_eMesh.owner_rank(element);
                      found = true;
                      break;
                    }
                }
#endif
            }
          if (!found && is3Dedge)
            {
              continue;
            }
          if (!found && element_proc < 0)
            {
              VERIFY_MSG("side not pointing to any positive permutation element, side= "+toString(m_eMesh.id(side)));
            }

          if (element_proc != m_eMesh.get_rank())
            {
              bool only_leaves = false;
              allD.clear();
              m_eMesh.allDescendants(side, allD, only_leaves);
              allD.insert(side);

              for (auto aside : allD)
                {
                  VERIFY_OP_ON(m_eMesh.owned(aside), ==, true, "root and child not on same proc");
                  change.push_back(stk::mesh::EntityProc(aside, element_proc));

                  // add FT's
                  stk::mesh::Entity const *entity_ft = m_eMesh.get_bulk_data()->begin(aside, FAMILY_TREE_RANK);
                  unsigned entity_ft_size = m_eMesh.get_bulk_data()->num_connectivity(aside, FAMILY_TREE_RANK);
                  for (unsigned ii = 0; ii < entity_ft_size; ++ii)
                    {
                      if(m_eMesh.owned(entity_ft[ii]))
                      {
                        change.push_back(stk::mesh::EntityProc(entity_ft[ii], element_proc));
                      }
                    }
                }
            }
        }
      }

      if (1)
        {
          stk::diag::Timer timerReqSidesCEO_("ReqSidesLocCEO", timerReqSides_);
          stk::diag::TimeBlock timerReqSidesCEOBlock_(timerReqSidesCEO_);


          size_t csz = change.size();
          stk::all_reduce( m_eMesh.get_bulk_data()->parallel() , stk::ReduceSum<1>( &csz ) );
          if (csz)
            m_eMesh.get_bulk_data()->change_entity_owner(change);
        }
    }

    void
    Refiner::
    unrefinePass2(SetOfEntities& elements_to_unref)
    {
      stk::diag::Timer timerAdapt_("RefineMesh", rootTimer());
      stk::diag::TimeBlock timerAdaptBlock_(timerAdapt_);

      stk::diag::Timer timerDoRefine_("percept::DoRefine", timerAdapt_);
      stk::diag::TimeBlock timerDoRefineBlock_(timerDoRefine_);

      stk::diag::Timer timerUnrefine_("unrefinePass2", timerDoRefine_);
      stk::diag::TimeBlock timerUnrefineBlock_(timerUnrefine_);

      REF_LTRACE_0();

      if (1)
      {
        TIMER2(unref_require_sides,Unrefine_);
        check_sides_on_same_proc_as_owned_element("unrefinePass2", true);
      }

      mod_begin(&timerUnrefine_);

      m_eMesh.initializeIdServer();

      SetOfEntities kept_nodes(*m_eMesh.get_bulk_data());
      SetOfEntities deleted_nodes(*m_eMesh.get_bulk_data());

      // filter unref set pass 2 - only root elements (no parent or not refined)
      SetOfEntities elements_to_be_remeshed(*m_eMesh.get_bulk_data());
      SetOfEntities family_trees_to_be_removed(*m_eMesh.get_bulk_data());
      SetOfEntities children_to_be_removed(*m_eMesh.get_bulk_data());
      SetOfEntities sides_to_be_removed(*m_eMesh.get_bulk_data());

      filterUnrefSetPass2(elements_to_unref, elements_to_be_remeshed);

      FixSideSetsSelectorUnrefine fss_unref(m_eMesh);
      fss_unref.add_elements(elements_to_unref);

      // get kept nodes
      {
        TIMER2(unref_keptNodes,Unrefine_);
        get_kept_nodes(kept_nodes, elements_to_unref);
      }

      // get deleted nodes
      {
        TIMER2(getDeletedNodes,Unrefine_);
        get_deleted_nodes(deleted_nodes, kept_nodes, elements_to_unref, elements_to_be_remeshed);
      }

      // filter deleted nodes - special case for Quad faces
      {
        TIMER2(filterDeletedNodes,Unrefine_);
        filter_deleted_nodes(deleted_nodes);
      }

      if (DEBUG_UNREF_1) std::cout << "deleted_nodes.size= " << deleted_nodes.size() << std::endl;

      // remove deleted nodes and their associated sub-dim entities from NodeRegistry's db
      //   save kept_nodes_orig on NodeRegistry
      SubDimCellToDataMap to_save;
      {
        TIMER2(cleanDeletedNodes,Unrefine_);

        SetOfEntities kept_nodes_orig(*m_eMesh.get_bulk_data());
        m_nodeRegistry->cleanDeletedNodes(deleted_nodes, kept_nodes_orig, to_save);
      }

      print_set(m_eMesh, elements_to_unref, "elements_to_unref");

      // find sides to remove
      if (1)
      {
        TIMER2(get_deleted_sides,Unrefine_);
        get_deleted_sides(sides_to_be_removed, elements_to_unref, elements_to_be_remeshed);
        //std::cout << "sides_to_be_removed.size= " << sides_to_be_removed.size() << std::endl;
        print_set(m_eMesh, sides_to_be_removed, "sides_to_be_removed");
      }

      // remove children
      std::vector<stk::mesh::Entity> new_elements, new_sides;

      SetOfEntities parent_side_elements(*m_eMesh.get_bulk_data());
      SetOfEntities new_element_set(*m_eMesh.get_bulk_data());
      SetOfEntities new_side_set(*m_eMesh.get_bulk_data());

      {
        TIMER2(removeChildren,Unrefine_);
        generate_temporary_elements(elements_to_unref, m_eMesh.element_rank(), new_elements);
        new_element_set.insert(new_elements.begin(), new_elements.end());
        delete_entities(elements_to_unref, m_eMesh.element_rank());
      }

      {
        TIMER2(deleteSides,Unrefine_);
        generate_temporary_elements(sides_to_be_removed, m_eMesh.side_rank(), new_sides);
        new_side_set.insert(new_sides.begin(), new_sides.end());
        delete_entities(sides_to_be_removed, m_eMesh.side_rank());
      }

      // reconnect and remove any dangling side elements
      bool allow_not_found = false;

      // FIXME
      {
        TIMER2(getSideParentsToBeRemeshed,Unrefine_);
        getSideParentsToBeRemeshed(elements_to_be_remeshed, parent_side_elements, true, &new_side_set);
        print_set(m_eMesh, elements_to_be_remeshed, "elements_to_be_remeshed");
        print_set(m_eMesh, parent_side_elements, "parent_side_elements");
      }

      // remesh the holes left by removing child elems (quad/hex hanging node doesn't need this)
      if (m_needsRemesh)
        {
          remeshElements(elements_to_be_remeshed, m_eMesh.element_rank(), 0);

          remeshElements(parent_side_elements, m_eMesh.side_rank(), 0);
        }

      // remove any elements that are empty (these can exist when doing local refinement)
      {
        TIMER2(removeEmptyElements,Unrefine_);
        removeEmptyElements();
      }

      {
        TIMER2(removeDeletedNodes,Unrefine_);
        removeDeletedNodes(deleted_nodes);
      }

      {
        TIMER2(unref_set_active_part,Unrefine_);
        set_active_part();
      }

      if (m_do_new_elements)
        {
          for (SetOfEntities::iterator iter = new_element_set.begin(); iter != new_element_set.end(); ++iter)
            {
              stk::mesh::Entity element = *iter;
              if ( !m_eMesh.get_bulk_data()->destroy_entity( element ) )
                throw std::runtime_error("pass2 4.3.1 - couldn't delete element");
            }
          for (SetOfEntities::iterator iter = new_side_set.begin(); iter != new_side_set.end(); ++iter)
            {
              stk::mesh::Entity element = *iter;
              if ( !m_eMesh.get_bulk_data()->destroy_entity( element ) )
                throw std::runtime_error("pass2 4.3.2 - couldn't delete side");
            }
          m_eMesh.initializeIdServer();
        }

      bool reduced_mod_end = true;
      if (m_eMesh.getProperty("percept_reduced_mod_end") == "false")
        reduced_mod_end = false;

      {
        //m_eMesh.get_bulk_data()->m_autoAuraOption = stk::mesh::BulkData::NO_AUTO_AURA;
        mod_end(&timerUnrefine_,"Unref1");
        //m_eMesh.get_bulk_data()->m_autoAuraOption = stk::mesh::BulkData::AUTO_AURA;
        mod_begin(&timerUnrefine_);
      }

      m_nodeRegistry->cleanInvalidNodes();

      if (1)
      {
        m_eMesh.setProperty("FixSideSets::fix_side_sets_2", "unrefinePass2 second");

        TIMER2(FSS_unref,Unrefine_);
        //fix_side_sets_2(allow_not_found, 0, 0, 0, "Unref1");
        fix_side_sets_2(allow_not_found, 0, 0, &fss_unref, "Unref1");
      }

      {
        TIMER2(ResetFT,Unrefine_);
        reset_family_tree_to_node_relations();
      }

      if (!reduced_mod_end)
      {
        mod_end(&timerUnrefine_,"Unref2");
        mod_begin(&timerUnrefine_);
      }

      if (m_nodeRegistry->s_use_new_ownership_check)
        {
          update_node_registry();
        }
      else
      {
        TIMER2(clear_element_owner_data_phase_2,Unrefine_);
        m_nodeRegistry->clear_element_owner_data_phase_2(true, false);
      }

      {
        TIMER2(unref_set_active_part_2,Unrefine_);
        set_active_part();
      }

      mod_end(&timerUnrefine_,"Unref3");

      REF_LTRACE_0();


      if (1)
      {
        TIMER2(unref_prolong,Unrefine_);
        m_nodeRegistry->prolongate(m_eMesh.get_coordinates_field(), 0, false);
        m_nodeRegistry->prolongateFields();
      }

      m_eMesh.set_parent_element_field();
      if (1)
        {
          RefinerUtil::save_node_registry(m_eMesh, *m_nodeRegistry, "unrefinePass2 end");
          if (m_eMesh.getProperty("NodeRegistry_rebuild_test") == "true")
            {
              static int cnt = 0;
              ++cnt;
              std::string file1 = "node_reg.5."+toString(cnt)+".e";
              m_eMesh.save_as(file1);
              PerceptMesh eMeshNR;
              eMeshNR.set_ioss_read_options("large");
              eMeshNR.open_read_only(file1);
              NodeRegistry nr1(eMeshNR, 0);
              RefinerUtil::rebuild_node_registry(eMeshNR, nr1, true, &m_eMesh, m_nodeRegistry, true);
            }
        }
    }

  } // namespace percept
