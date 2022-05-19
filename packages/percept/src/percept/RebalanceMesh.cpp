// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/RebalanceMesh.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MeshUtils.hpp>

#include <percept/stk_rebalance/Rebalance.hpp>
#include <percept/stk_rebalance/Partition.hpp>
#include <percept/stk_rebalance/ZoltanPartition.hpp>

#include <percept/stk_rebalance_utils/RebalanceUtils.hpp>

#include <stk_util/environment/CPUTime.hpp>

//----------------------------------------------------------------------

namespace percept {

  void RebalanceMesh::fixup_family_trees(PerceptMesh& eMesh, bool debug)
  {
    double ft_time = stk::cpu_time();
    eMesh.get_bulk_data()->modification_begin();
    std::vector<stk::mesh::Entity> vec;
    stk::mesh::Selector on_univ =  eMesh.get_fem_meta_data()->universal_part();
    const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(eMesh.element_rank()+1);
    stk::mesh::get_selected_entities(on_univ, eMesh.get_bulk_data()->buckets(FAMILY_TREE_RANK), vec);

    std::vector<stk::mesh::Entity> new_rels;

    stk::mesh::EntityRank ranks[] = {eMesh.side_rank(), eMesh.element_rank()};
    const unsigned nranks = 2;
    for (unsigned iFT = 0; iFT < vec.size(); iFT++)
      {
        stk::mesh::Entity ft = vec[iFT];
        stk::mesh::Entity const *ft_nodes = eMesh.get_bulk_data()->begin(ft, eMesh.node_rank());
        unsigned ft_nodes_size = eMesh.get_bulk_data()->num_connectivity(ft, eMesh.node_rank());

        new_rels.resize(0);
        for (unsigned irank = 0; irank < nranks; ++irank)
          {
            unsigned ft_elems_size = eMesh.get_bulk_data()->num_connectivity(ft, ranks[irank]);
            if (!ft_elems_size) continue;
            stk::mesh::Entity const* ft_elems = eMesh.get_bulk_data()->begin(ft, ranks[irank]);

            stk::mesh::Entity root = eMesh.rootOfTree(ft_elems[0]);
            VERIFY_OP_ON(eMesh.is_valid(root), ==, true, "bad root");

            // FIXME
            bool found_shared = (eMesh.shared(root) || !eMesh.owned(root))
              || eMesh.get_bulk_data()->in_send_ghost(root);
            if (!found_shared)
              {
                continue;
              }

            stk::mesh::Entity const *root_nodes = eMesh.get_bulk_data()->begin(root, eMesh.node_rank());
            unsigned root_nodes_size = eMesh.get_bulk_data()->num_connectivity(root, eMesh.node_rank());

            found_shared = false;
            for (unsigned mm=0; mm < root_nodes_size; ++mm)
              {
                if (eMesh.shared(root_nodes[mm]))
                  {
                    found_shared = true;
                    break;
                  }
              }
            if (!found_shared)
              continue;

            for (unsigned mm=0; mm < root_nodes_size; ++mm)
              {
                if (!eMesh.shared(root_nodes[mm]))
                  continue;
                bool found = false;
                for (unsigned kk=0; kk < ft_nodes_size; ++kk)
                  {
                    if (ft_nodes[kk] == root_nodes[mm])
                      {
                        found=true;
                        break;
                      }
                  }
                if (!found)
                  {
                    new_rels.push_back(root_nodes[mm]);
                  }
              }
          }

        for (unsigned jj=0; jj < new_rels.size(); ++jj)
          {
            eMesh.get_bulk_data()->declare_relation(ft, new_rels[jj], ft_nodes_size + jj);
          }
      }

    ft_time = stk::cpu_time() - ft_time;

    stk::mesh::fixup_ghosted_to_shared_nodes(*eMesh.get_bulk_data());
    double mod_end_time = stk::cpu_time();
    eMesh.get_bulk_data()->modification_end();
    mod_end_time = stk::cpu_time() - mod_end_time;

    stk::all_reduce(eMesh.get_bulk_data()->parallel(), stk::ReduceSum<1>(&mod_end_time));
    stk::all_reduce(eMesh.get_bulk_data()->parallel(), stk::ReduceSum<1>(&ft_time));
    if (debug && eMesh.get_rank() == 0)
      {
        std::cout << "ft_time = " << ft_time << " mod_end_time= " << mod_end_time << std::endl;
      }
  }

  class RebalanceAdapted : public stk::rebalance::Rebalance
  {
  public:
    RebalanceAdapted(PerceptMesh& eMesh, bool debug = false) : stk::rebalance::Rebalance(), m_eMesh(eMesh), m_debug(debug) {}
  protected:
    PerceptMesh& m_eMesh;
    bool m_debug;

    void check_ownership()
    {
      if (!m_debug) return;
      std::vector<stk::mesh::Entity> vec;
      stk::mesh::Selector on_univ =  m_eMesh.get_fem_meta_data()->universal_part();
      stk::mesh::get_selected_entities(on_univ, m_eMesh.get_bulk_data()->buckets(m_eMesh.element_rank()), vec);

      for (unsigned iElement = 0; iElement < vec.size(); iElement++)
        {
          stk::mesh::Entity element = vec[iElement];
          stk::mesh::Entity root = m_eMesh.rootOfTree(element);
          if (element != root)
            {
              VERIFY_OP_ON(m_eMesh.owner_rank(element), ==, m_eMesh.owner_rank(root), "bad ownership");
            }
        }
    }

    void check_ghosting(const std::string& msg)
    {
      if (!m_debug) return;
      std::vector<stk::mesh::Entity> vec;
      //stk::mesh::Selector on_univ =  m_eMesh.get_fem_meta_data()->universal_part();
      stk::mesh::Selector sel =  m_eMesh.get_fem_meta_data()->locally_owned_part();
      stk::mesh::get_selected_entities(sel, m_eMesh.get_bulk_data()->buckets(m_eMesh.element_rank()), vec);

      std::vector<int>  procs, root_procs;
      for (unsigned iElement = 0; iElement < vec.size(); iElement++)
        {
          stk::mesh::Entity element = vec[iElement];
          stk::mesh::Entity root = m_eMesh.rootOfTree(element);
          if (1)
            {
              bool has_family_tree = m_eMesh.hasFamilyTree(root);
              bool is_child = has_family_tree && m_eMesh.isChildElement(root, true);
              bool is_root = !has_family_tree || !is_child;
              VERIFY_OP_ON(is_root, ==, true, "bad root "+msg);
            }
          if (element != root)
            {
              VERIFY_OP_ON(m_eMesh.owner_rank(element), ==, m_eMesh.owner_rank(root), "bad ownership "+msg);
              m_eMesh.get_bulk_data()->comm_procs(element, procs);
              m_eMesh.get_bulk_data()->comm_procs(root, root_procs);
              std::sort(procs.begin(), procs.end());
              std::sort(root_procs.begin(), root_procs.end());
              VERIFY_OP_ON(procs.size(), ==, root_procs.size(), "bad size "+msg);
              for (unsigned ii=0; ii < procs.size(); ++ii)
                {
                  VERIFY_OP_ON(procs[ii], ==, root_procs[ii], "bad p2 "+msg);
                }
            }
        }
    }

    void fixup_family_trees()
    {
      RebalanceMesh::fixup_family_trees(m_eMesh, m_debug);
    }

    virtual bool full_rebalance(stk::mesh::BulkData  & bulk_data ,
                                stk::rebalance::Partition       * partition,
                                const stk::mesh::EntityRank rank)
    {
      using namespace stk;
      using namespace stk::rebalance;
      mesh::EntityProcVec cs_elem;
      bool rebalancingHasOccurred =  balance_comm_spec_domain( partition, cs_elem );

      if (rebalancingHasOccurred)
        {
          mesh::EntityProcVec cs_elem_new = cs_elem;
          SetOfEntities descendants(*m_eMesh.get_bulk_data());
          for (size_t ii=0; ii < cs_elem.size(); ++ii)
            {
              stk::mesh::Entity element = cs_elem[ii].first;
              descendants.clear();
              bool only_leaves = false;
              m_eMesh.allDescendants(element, descendants, only_leaves);
              for (SetOfEntities::iterator it = descendants.begin(); it != descendants.end(); ++it)
                {
                  VERIFY_OP_ON(bulk_data.parallel_owner_rank(*it), ==, bulk_data.parallel_owner_rank(element), "bad par owner");
                  cs_elem_new.push_back(stk::mesh::EntityProc(*it, cs_elem[ii].second));
                }
            }
          cs_elem = cs_elem_new;
        }

      const stk::mesh::EntityRank node_rank = stk::topology::NODE_RANK;
      const stk::mesh::EntityRank edge_rank = stk::topology::EDGE_RANK;
      const stk::mesh::EntityRank face_rank = stk::topology::FACE_RANK;
      const stk::mesh::EntityRank elem_rank = stk::topology::ELEMENT_RANK;
      // note, this could also be the FAMILY_TREE_RANK...
      const stk::mesh::EntityRank constraint_rank = static_cast<stk::mesh::EntityRank>(elem_rank+1);

      if (1 && rebalancingHasOccurred && partition->partition_dependents_needed() )
        {
          stk::mesh::EntityRank ranks[] = {node_rank, m_eMesh.side_rank(), constraint_rank};
          const unsigned nranks = 3;
          mesh::EntityProcVec cs_elem_new = cs_elem;
          for (unsigned irank=0; irank < nranks; ++irank)
            {
              for (size_t ii=0; ii < cs_elem_new.size(); ++ii)
                {
                  mesh::EntityVector related_entities, side_to_family_tree;
                  mesh::EntityVector elems(1), sides(1);
                  elems[0] = cs_elem_new[ii].first;
                  stk::mesh::get_entities_through_relations(bulk_data, elems, ranks[irank], related_entities);
                  for( size_t j = 0; j < related_entities.size(); ++j )
                    {
                      cs_elem.push_back(stk::mesh::EntityProc(related_entities[j], cs_elem_new[ii].second));
                      if (ranks[irank] == m_eMesh.side_rank())
                        {
                          sides[0] = related_entities[j];
                          stk::mesh::get_entities_through_relations(bulk_data, sides, constraint_rank, side_to_family_tree);
                          for( size_t k = 0; k < side_to_family_tree.size(); ++k )
                            {
                              cs_elem.push_back(stk::mesh::EntityProc(side_to_family_tree[k], cs_elem_new[ii].second));
                            }
                        }
                    }
                }
            }
          // now make it unique
          typedef std::map<stk::mesh::Entity, int, stk::mesh::EntityLess> Emap;
          Emap emap(*m_eMesh.get_bulk_data());
          for (size_t ii=0; ii < cs_elem.size(); ++ii)
            {
              if (bulk_data.parallel_owner_rank(cs_elem[ii].first) == bulk_data.parallel_rank())
                {
                  emap[cs_elem[ii].first] = cs_elem[ii].second;
                }
            }
          if (m_debug) std::cout << "P[" << bulk_data.parallel_rank() << "] cs_elem.size= " << cs_elem.size() << " emap.size= " << emap.size() << std::endl;
          cs_elem.resize(0);
          for (Emap::iterator it=emap.begin(); it != emap.end(); ++it)
            {
              cs_elem.push_back(stk::mesh::EntityProc(it->first, it->second));
            }
          // check
          if (m_debug)
            {
              SetOfEntities descendants(*m_eMesh.get_bulk_data());
              for (size_t ii=0; ii < cs_elem.size(); ++ii)
                {
                  stk::mesh::Entity element = cs_elem[ii].first;
                  int proc = cs_elem[ii].second;
                  if (m_eMesh.entity_rank(element) == elem_rank)
                    {
                      stk::mesh::Entity parent = m_eMesh.getParent(element, false);
                      if (m_eMesh.hasFamilyTree(element))
                        {
                          bool isChild = m_eMesh.isChildElement(element);
                          VERIFY_OP_ON(!isChild, ==, !m_eMesh.is_valid(parent), "bad parent");
                        }
                      if (0)
                        std::cout << "P[" << m_eMesh.get_rank() << "] check e= " << m_eMesh.identifier(element)
                                  << " c= " << m_eMesh.numChildren(element)
                                  << " hft= " << m_eMesh.hasFamilyTree(element)
                                  << " pv= " << m_eMesh.is_valid(parent)
                                  << std::endl;

                      if (!m_eMesh.hasFamilyTree(element) || !m_eMesh.is_valid(parent))
                        {
                          descendants.clear();
                          bool only_leaves = false;
                          m_eMesh.allDescendants(element, descendants, only_leaves);
                          for (SetOfEntities::iterator it = descendants.begin(); it != descendants.end(); ++it)
                            {
                              if(emap.find(*it) == emap.end())
                                throw std::runtime_error("not in map");
                              int new_proc = emap[*it];
                              VERIFY_OP_ON(new_proc, ==, proc, "bad proc");
                              VERIFY_OP_ON(m_eMesh.owner_rank(*it), ==, m_eMesh.owner_rank(element), "bad owner");
                            }
                        }
                    }
                }
            }
          // second check
          if (m_debug)
            {
              stk::mesh::FieldBase *proc_rank_field=0;
              if (!proc_rank_field) proc_rank_field=m_eMesh.get_field(stk::topology::ELEMENT_RANK, "proc_rank");
              if (proc_rank_field)
                {
                  for (size_t ii=0; ii < cs_elem.size(); ++ii)
                    {
                      stk::mesh::Entity element = cs_elem[ii].first;
                      int proc = cs_elem[ii].second;
                      if (m_eMesh.entity_rank(element) == elem_rank)
                        {
                          double *fdata = m_eMesh.field_data(proc_rank_field, element);
                          fdata[0] = double(proc);
                        }
                    }
                }
              //if (m_debug) m_eMesh.save_as("predicted_ranks.e");
            }
        }
      else if(rebalancingHasOccurred && partition->partition_dependents_needed() )
        {
          //stk::mesh::MetaData & fem_meta = stk::mesh::MetaData::get(bulk_data);

          // Don't know the rank of the elements rebalanced, assume all are dependent.
          rebalance_dependent_entities( bulk_data, partition, node_rank, cs_elem, rank );
          if (stk::mesh::InvalidEntityRank != edge_rank && rank != edge_rank)
            rebalance_dependent_entities( bulk_data, partition, edge_rank, cs_elem, rank );
          if (stk::mesh::InvalidEntityRank != face_rank && rank != face_rank)
            rebalance_dependent_entities( bulk_data, partition, face_rank, cs_elem, rank );
          if (stk::mesh::InvalidEntityRank != elem_rank && rank != elem_rank)
            rebalance_dependent_entities( bulk_data, partition, elem_rank, cs_elem, rank );
          if (stk::mesh::InvalidEntityRank != constraint_rank && rank != constraint_rank)
            rebalance_dependent_entities( bulk_data, partition, constraint_rank, cs_elem, rank );
        }

      if ( rebalancingHasOccurred )
        {
          check_ghosting("before");
          double bd_time = stk::cpu_time();
          bulk_data.change_entity_owner( cs_elem );
          bd_time = stk::cpu_time() - bd_time;
          check_ownership();
          double ft_time = stk::cpu_time();
          fixup_family_trees();
          ft_time = stk::cpu_time() - ft_time;
          check_ownership();
          check_ghosting("after");

          stk::all_reduce(m_eMesh.get_bulk_data()->parallel(), stk::ReduceSum<1>(&bd_time));
          stk::all_reduce(m_eMesh.get_bulk_data()->parallel(), stk::ReduceSum<1>(&ft_time));
          if (m_debug && m_eMesh.get_rank() == 0)
            {
              std::cout << "fixup_family_trees time = " << ft_time << " change_entity_owner time= " << bd_time << std::endl;
            }
        }

      //: Finished
      return rebalancingHasOccurred;
    }

  };

  RebalanceMesh::RebalanceMesh(PerceptMesh& mesh, ScalarFieldType *weight_field, bool deb, bool doAdaptedMesh)
    : m_eMesh(mesh), m_debug(deb), m_weight_field(weight_field), m_rank_to_rebalance(m_eMesh.element_rank()), m_doAdaptedMesh(doAdaptedMesh)
  {
  }

  void RebalanceMesh::build_entities_to_rebalance_list(stk::mesh::EntityVector& vec)
  {
    vec.resize(0);
    const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.element_rank() );
    for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
      {
        if ((**k).owned())
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_elements_in_bucket = bucket.size();
            for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
              {
                stk::mesh::Entity element = bucket[iElement];
                stk::mesh::Entity parent = m_eMesh.getParent(element, false);
                if (0)
                  std::cout << "P[" << m_eMesh.get_rank() << "] e= " << m_eMesh.identifier(element)
                            << " c= " << m_eMesh.numChildren(element)
                            << " hft= " << m_eMesh.hasFamilyTree(element)
                            << " pv= " << m_eMesh.is_valid(parent)
                            << std::endl;
                int nc = m_eMesh.numChildren(element);
                VERIFY_OP_ON(nc, !=, 1, "bad nc");

                if (!m_eMesh.hasFamilyTree(element) || !m_eMesh.is_valid(parent))
                  {
                    vec.push_back(element);
                  }
              }
          }
      }
  }

  double RebalanceMesh::rebalance()
  {
    const unsigned p_rank = m_eMesh.get_parallel_rank();

    stk::mesh::MetaData & fem_meta  = *m_eMesh.get_fem_meta_data();
    stk::mesh::BulkData & bulk  = *m_eMesh.get_bulk_data();

    // Put weights field on all elements
    //const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(m_eMesh.element_rank() + 1u);
    m_rank_to_rebalance = m_eMesh.element_rank();

    // Use Zoltan to determine new partition
    Teuchos::ParameterList emptyList;
    stk::rebalance::Zoltan zoltan_partition(bulk, bulk.parallel(), m_eMesh.get_spatial_dim(), emptyList);

    stk::mesh::Selector selector(fem_meta.universal_part());
    stk::mesh::Selector owned_selector = fem_meta.locally_owned_part();

    double imbalance_threshold_before = 0;
    double imbalance_threshold_after = 0;

    std::vector<double> load_factor;
    m_eMesh.get_load_factor(load_factor, m_debug, "before rebalance");

    if (m_doAdaptedMesh)
      {
        if (m_rank_to_rebalance == m_eMesh.element_rank())
          {
            std::vector<stk::mesh::Entity> entities_to_rebalance;

            build_entities_to_rebalance_list(entities_to_rebalance);
            set_weights(entities_to_rebalance);

            imbalance_threshold_before = stk::rebalance::check_balance(bulk, m_weight_field, m_rank_to_rebalance, 0, &entities_to_rebalance);

            RebalanceAdapted rb(m_eMesh, m_debug);
            rb.rebalance(bulk, selector, m_eMesh.get_coordinates_field(), m_weight_field, zoltan_partition, m_rank_to_rebalance, &entities_to_rebalance, true, m_debug);

            build_entities_to_rebalance_list(entities_to_rebalance);
            set_weights(entities_to_rebalance);

            // do it a second time?
            if (0)
              {
                stk::rebalance::Zoltan zoltan_partition1(bulk, bulk.parallel(), m_eMesh.get_spatial_dim(), emptyList);
                RebalanceAdapted rb1(m_eMesh, m_debug);
                rb1.rebalance(bulk, selector, m_eMesh.get_coordinates_field(), m_weight_field, zoltan_partition1, m_rank_to_rebalance, &entities_to_rebalance, false, m_debug);
                build_entities_to_rebalance_list(entities_to_rebalance);
                set_weights(entities_to_rebalance);
              }

            imbalance_threshold_after = stk::rebalance::check_balance(bulk, m_weight_field, m_rank_to_rebalance, 0, &entities_to_rebalance);
          }
        else
          {
            throw std::runtime_error("RebalanceMesh:: not implemented");
          }
      }
    else
      {
        imbalance_threshold_before = compute_imbalance();
        stk::rebalance::Rebalance rb;
        rb.rebalance(bulk, selector, m_eMesh.get_coordinates_field(), 0, zoltan_partition, m_rank_to_rebalance);
        imbalance_threshold_after = compute_imbalance();
      }

    if(m_debug && 0 == p_rank )
      std::cout << std::endl
                << "RebalanceMesh:: imbalance_threshold before rebalance = " << imbalance_threshold_before << std::endl;

    m_eMesh.get_load_factor(load_factor, m_debug, "after rebalance");

    size_t num_local_elems = stk::mesh::count_selected_entities(owned_selector, bulk.buckets(m_rank_to_rebalance));

    if(m_debug && 0 == p_rank )
      std::cout << std::endl
                << "RebalanceMesh:: imbalance_threshold before rebalance = " << imbalance_threshold_before
                << ", after = " << imbalance_threshold_after
                << " num_local_elems= " << num_local_elems
                << std::endl;

    m_eMesh.set_proc_rank_field();

    return imbalance_threshold_after;
  }

  void
  RebalanceMesh::set_weights(const std::vector<stk::mesh::Entity>& entities_to_rebalance)
  {
    stk::mesh::BulkData & bulk  = *m_eMesh.get_bulk_data();

    if (m_weight_field) {

      SetOfEntities descendants(bulk);

      for (size_t iElement=0; iElement < entities_to_rebalance.size(); ++iElement) {

        stk::mesh::Entity element = entities_to_rebalance[iElement];

        double w = 0.0;
        descendants.clear();
        bool only_leaves = false;
        m_eMesh.allDescendants(element, descendants, only_leaves);
        for (SetOfEntities::iterator it = descendants.begin(); it != descendants.end(); ++it) {
          if (m_eMesh.numChildren(*it) == 0)
            w += 1.0;
        }
        if (m_eMesh.numChildren(element)) {
          VERIFY_OP_ON(descendants.size(), !=, 0, "bad");
        }
        if (w == 0.0)
          w = 1.0;
        double *f_data = m_eMesh.field_data(*m_weight_field, element);
        f_data[0] = w;
      }

      //if (m_debug) m_eMesh.save_as("weights.e");
    }
  }

  double
  RebalanceMesh::compute_imbalance()
  {
    stk::mesh::BulkData & bulk  = *m_eMesh.get_bulk_data();

    double imbalance = 0.0;
    if (m_doAdaptedMesh)
      {
        std::vector<stk::mesh::Entity> eb;
        build_entities_to_rebalance_list(eb);
        set_weights(eb);

        imbalance = stk::rebalance::check_balance(bulk, m_weight_field, m_rank_to_rebalance, 0, &eb);
      }
    else
      {
        imbalance = stk::rebalance::check_balance(bulk, m_weight_field, m_rank_to_rebalance);
      }

    return imbalance;
  }
}
