// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_TransitionElementAdapter_hpp
#define adapt_TransitionElementAdapter_hpp

#include <functional>
#include <queue>
#include <cmath>
#include <typeinfo>

#include <adapt/AdaptedMeshVerifier.hpp>
#include <stk_util/environment/CPUTime.hpp>
//#include <percept/GeometryVerifier.hpp>
#include <percept/Histograms.hpp>

#include <adapt/HangingNodeAdapter.hpp>

#define USE_ADAPT_PERFORMANCE_TESTING_CALLGRIND 0
#if USE_ADAPT_PERFORMANCE_TESTING_CALLGRIND
#include "/usr/netpub/valgrind-3.8.1/include/valgrind/callgrind.h"
#endif

#define DO_ALT_TIMER 0

#define TIMING(code) code
#define TIMER(name) stk::diag::Timer timer ## name ( #name, Base::rootTimer());  stk::diag::TimeBlock tbTimer ## name (timer ## name)
#define TIMER2(name,parentName) stk::diag::Timer timer ## name ( #name, timer ## parentName);  stk::diag::TimeBlock tbTimer ## name (timer ## name)
#define TIMER3(name,parentName) \
   stk::diag::Timer timer ## name ( #name, timer ## parentName); \
   stk::diag::Timer * altTimer ## name = Base::getAlternateRootTimer(); \
   if (altTimer ## name && DO_ALT_TIMER) \
     Base::setAlternateRootTimer(& timer ## name ); \
   stk::diag::TimeBlock tbTimer ## name (timer ## name)

#define START_TIMER(name) do { timersLocal_[#name] = stk::cpu_time(); } while(0)
#define STOP_TIMER(name)  do { timersLocal_[#name] = stk::cpu_time() - timersLocal_[#name]; timers_[#name] += timersLocal_[#name]; } while(0)

#include <stk_mesh/base/MeshUtils.hpp>

  namespace percept {

    extern bool s_do_transition_break;
    static std::map<std::string, double> timers_, timersLocal_;
    static bool s_allow_skip_edge_neighbor = true;

    template<class RefinePredicate>
    class TransitionElementAdapter : public HangingNodeAdapter<RefinePredicate>
    {

    protected:
      RefineFieldType                    *m_refine_field;
      bool                                m_refine_field_set;
      TransitionElementType              *m_transition_element_field;
      bool                                m_transition_element_field_set;
      int                                 m_fileid;
      bool                                m_debug;
      bool                                m_debug_print;
      AdaptedMeshVerifier                *m_adaptedMeshVerifier;

    public:
      typedef HangingNodeAdapter<RefinePredicate> Base;

      ~TransitionElementAdapter() override {
        if( m_adaptedMeshVerifier ) {
          delete m_adaptedMeshVerifier;
        }
      }

      TransitionElementAdapter(RefinePredicate& predicate_refine,
                               percept::PerceptMesh& eMesh,
                               UniformRefinerPatternBase & bp,
                               stk::mesh::FieldBase *proc_rank_field=0,
                               bool debug = false) :
        Base(predicate_refine, eMesh, bp, proc_rank_field),
        m_refine_field(0), m_refine_field_set(false),
        m_transition_element_field(0), m_transition_element_field_set(false),
        m_fileid(0), m_debug(debug), m_debug_print(false), m_adaptedMeshVerifier(0)
      {
        const std::type_info& bpt = typeid(bp);
        if (bpt == typeid(Local_Tri3_Tri3_N_HangingNode)
            || bpt == typeid(Local_Tet4_Tet4_N_HangingNode)
            || bpt == typeid(Local_Quad4_Quad4_N_Transition)
            || bpt == typeid(Local_Quad4_Tri3_Hybrid_Transition)
            || bpt == typeid(Local_Hex8_Hex8_N_Transition)
            || bpt == typeid(Local_Pyr5_Pyr5_N_Transition)
            || bpt == typeid(Local_Wedge6_Wedge6_N_Transition)
            || bpt == typeid(Local_Hybrid_3D)
            )
          {
          }
        else
          {
            throw std::logic_error("TransitionElementAdapter can only use "
                                   "\nLocal_Tet4_Tet4_N_HangingNode,"
                                   "\nLocal_Tri3_Tri3_N_HangingNode,"
                                   "\nLocal_Quad4_Quad4_N_Transition,"
                                   "\nLocal_Pyr5_Pyr5_N_Transition,"
                                   "\nLocal_Wedge6_Wedge6_N_Transition,"
                                   "\nLocal_Hybrid_3D,"
                                   "\nLocal_Hex8_Hex8_N_Transition refine patterns");
          }
        Base:: m_predicate_refine.setRefineStage(-1);

        if (bp.m_mark_centroid_always)
          {
            predicate_refine.m_mark_centroid_always = true;
          }
      }

      void print_counts(const std::string& msg)
      {
        PerceptMesh& eMesh = Base::m_eMesh;
        size_t nmarked_refine=0, nmarked_unrefine=0;
        size_t nte=0;
        size_t nele=0, nele_parents=0;
        stk::mesh::Selector on_locally_owned_part =  ( eMesh.get_fem_meta_data()->locally_owned_part() );
        const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( eMesh.element_rank() );
        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            stk::mesh::Bucket & bucket = **k ;
            if (!on_locally_owned_part(bucket)) continue;
            const unsigned num_entity_in_bucket = bucket.size();
            for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
              {
                stk::mesh::Entity element = bucket[ientity];
                const bool check_for_family_tree = false;
                ++nele_parents;
                bool isParent = eMesh.isParentElement(element, check_for_family_tree);
                if (isParent)
                  continue;
                ++nele;
                int *refine_field_elem = stk::mesh::field_data( *m_refine_field , element );
                int *transition_element = stk::mesh::field_data( *m_transition_element_field , element );
                if (transition_element[0]) ++nte;
                if (refine_field_elem[0] > 0)
                  nmarked_refine++;
                else if (refine_field_elem[0] < 0)
                  nmarked_unrefine++;

                if (eMesh.hasFamilyTree(element))
                  {
                    stk::mesh::Entity parent = eMesh.getParent(element, true);
                    if (eMesh.is_valid(parent))
                      {
                        int *refine_field_parent = stk::mesh::field_data( *m_refine_field , parent );
                        if (refine_field_parent[0] > 0)
                          nmarked_refine++;
                        else if (refine_field_parent[0] < 0)
                          nmarked_unrefine++;
                      }
                  }

              }
          }
        stk::ParallelMachine pm = Base::m_eMesh.get_bulk_data()->parallel();
        stk::all_reduce( pm, stk::ReduceSum<1>( &nmarked_refine ));
        stk::all_reduce( pm, stk::ReduceSum<1>( &nmarked_unrefine ));
        stk::all_reduce( pm, stk::ReduceSum<1>( &nte ));
        stk::all_reduce( pm, stk::ReduceSum<1>( &nele ));
        stk::all_reduce( pm, stk::ReduceSum<1>( &nele_parents ));
        const unsigned numRefinedElemPerElem = eMesh.get_spatial_dim() == 3 ? 8 : 4;
        size_t neleNew = nele + numRefinedElemPerElem*nmarked_refine;

        if (eMesh.get_rank() == 0)
          {
            std::cout << "TransitionElementAdapter:: " << msg
                      << " nte= " << nte
                      << " marked ref= "
                      << nmarked_refine << " marked unref= " << nmarked_unrefine
                      << " nele active= " << nele << " nele total= " << nele_parents
                      << " nele new (est)= " << neleNew
                      << std::endl;
          }
      }


      void save(std::string msg)
      {
        if (m_debug)
          {
            std::ostringstream fileid_ss;
            fileid_ss << std::setfill('0') << std::setw(4) << (m_fileid);

            std::string oname = "debug-adapt-seq.e";
            if (m_fileid > 0) oname += "-s" + fileid_ss.str();
            if (Base::m_eMesh.get_rank() == 0)
              std::cout << "oname= " << oname << " " << msg << std::endl;
            Base::m_eMesh.save_as(oname);
            ++m_fileid;
          }
      }

      void save1(std::string msg, std::string onameIn)
      {
        if (m_debug)
          {
            static std::map<std::string, int> fid;
            if (fid.find(onameIn) == fid.end())
              fid[onameIn] = 0;
            std::ostringstream fileid_ss;
            fileid_ss << std::setfill('0') << std::setw(4) << (fid[onameIn]);

            std::string oname = "debug-adapt-seq.e";
            if (onameIn.length()) oname = onameIn;
            if (fid[onameIn] > 0) oname += "-s" + fileid_ss.str();
            if (Base::m_eMesh.get_rank() == 0)
              std::cout << "oname= " << oname << " " << msg << std::endl;
            bool sv = Base::m_eMesh.get_output_active_children_only();
            Base::m_eMesh.output_active_children_only(false);
            Base::m_eMesh.save_as(oname);
            Base::m_eMesh.output_active_children_only(sv);
            ++fid[onameIn];
          }
      }

      bool is_edge_neighbor_across_given_edge(stk::mesh::Entity element_0, stk::mesh::Entity element_1, int edge_0, int *edge_1,
                                              const CellTopologyData * const bucket_topo_data_element_0
                                              )
      {
        PerceptMesh& eMesh = Base::m_eMesh;
        const MyPairIterRelation element_0_nodes(*eMesh.get_bulk_data(), element_0, stk::topology::NODE_RANK );
        const MyPairIterRelation element_1_nodes(*eMesh.get_bulk_data(), element_1, stk::topology::NODE_RANK );

        const CellTopologyData * const element_0_topo_data = (bucket_topo_data_element_0 ? bucket_topo_data_element_0 : eMesh.get_cell_topology(element_0));
        shards::CellTopology element_0_topo(element_0_topo_data);

        const CellTopologyData * const element_1_topo_data = (element_0_nodes.size() == element_1_nodes.size() ? element_0_topo_data : eMesh.get_cell_topology(element_1));
        shards::CellTopology element_1_topo(element_1_topo_data);

        int iedge_0 = edge_0;
        {
          unsigned num_nodes_on_edge_0 = element_0_topo_data->edge[iedge_0].topology->vertex_count;
          for (unsigned iedge_1 = 0; iedge_1 <  element_1_topo_data->edge_count; iedge_1++)
            {
              unsigned num_nodes_on_edge_1 = element_1_topo_data->edge[iedge_1].topology->vertex_count;
              if (num_nodes_on_edge_0 != num_nodes_on_edge_1)
                continue;
              bool edges_match = true;
              for (unsigned jnode_0 = 0; jnode_0 < num_nodes_on_edge_0; jnode_0++)
                {
                  stk::mesh::EntityId edge_0_id = eMesh.identifier(element_0_nodes[ element_0_topo_data->edge[iedge_0].node[jnode_0] ].entity());
                  bool found = false;
                  for (unsigned jnode_1 = 0; jnode_1 < num_nodes_on_edge_1; jnode_1++)
                    {
                      stk::mesh::EntityId edge_1_id = eMesh.identifier(element_1_nodes[ element_1_topo_data->edge[iedge_1].node[jnode_1] ].entity());
                      if (edge_1_id == edge_0_id)
                        {
                          found = true;
                          break;
                        }
                    }
                  if (!found)
                    {
                      edges_match = false;
                      break;
                    }
                }
              if (edges_match)
                {
                  if (edge_1) *edge_1 = iedge_1;
                  return true;
                }
            }
        }
        return false;
      }

      bool process_element(stk::mesh::Entity element, bool& did_change, int& nchange_1, int& nchange_2,
                           LocalSetType& selected_neighbors, bool check_only=false)
      {
        PerceptMesh& eMesh = Base::m_eMesh;

        bool ret_val = false;

        int *transition_element = stk::mesh::field_data( *m_transition_element_field , element );

        if (!transition_element[0])
          return false;

        bool isElemParent = eMesh.hasFamilyTree(element) && eMesh.isParentElement(element, true);
        VERIFY_OP_ON(isElemParent, ==, false, "hmmm");

        selected_neighbors.clear();
        Base::get_node_neighbors(element, selected_neighbors);

        // for all node neighbors that are leaves
        LocalSetType::iterator neighbor_begin = selected_neighbors.begin();
        LocalSetType::iterator neighbor_end = selected_neighbors.end();

        for (LocalSetType::iterator neighbor = neighbor_begin;
             neighbor != neighbor_end; ++neighbor)
          {
            stk::mesh::Entity neigh = *neighbor;

            bool isNParent = eMesh.hasFamilyTree(neigh) && eMesh.isParentElement(neigh, true);
            int *refine_field_neigh = stk::mesh::field_data( *m_refine_field , neigh );
            int *refine_field_elem = stk::mesh::field_data( *m_refine_field , element );

            if (refine_field_neigh[0] >= 1 )
              {
                int edge_0 = -1, edge_1 = -1;
                bool isEdgeN = eMesh.is_edge_neighbor(element, neigh, &edge_0, &edge_1, 0);
                if (isEdgeN)
                  {
                    stk::mesh::Entity parent = eMesh.getParent(element, true);
                    VERIFY_OP_ON(eMesh.is_valid(parent), ==, true, "hmmm");
                    VERIFY_OP_ON(eMesh.isGhostElement(parent), ==, false, "bad parent 2");

                    if (isEdgeN)
                      {
                        int jedge_1 = -1;
                        if (s_allow_skip_edge_neighbor)
                          {
                            bool ien = is_edge_neighbor_across_given_edge(element, parent, edge_0, &jedge_1, 0);
                            if (ien)
                              continue;
                          }
                      }

                    int *refine_field_parent = stk::mesh::field_data( *m_refine_field , parent );
                    if (refine_field_parent[0] < 2)
                      {
                        refine_field_parent[0] = 2;

                        if (check_only)
                          {
                            std::cout << "check_only: upgrading parent,eleme,neigh= " << eMesh.identifier(parent) << " "
                                      << eMesh.identifier(element) << " "
                                      << eMesh.identifier(neigh)
                                      << " refine_field_neigh= " << refine_field_neigh[0]
                                      << " isNParent= " << isNParent
                                      << std::endl;
                            return true;
                          }

                        static std::vector<stk::mesh::Entity> children;
                        eMesh.getChildren(parent, children, true, false);
                        for (unsigned ii=0; ii < children.size(); ++ii)
                          {
                            int *refine_field_child = stk::mesh::field_data( *m_refine_field , children[ii] );
                            refine_field_child[0] = 0;
                            if (0)
                               std::cout << "upgrading parent,child,neigh= " << eMesh.identifier(parent) << " "
                                        << eMesh.identifier(children[ii]) << " "
                                        << eMesh.identifier(neigh)
                                        << " refine_field_neigh= " << refine_field_neigh[0]
                                        << " isNParent= " << isNParent
                                        << std::endl;
                          }

                        VERIFY_OP_ON(refine_field_elem[0], ==, 0, "hmmm");
                        ++nchange_2;
                        did_change = true;
                        ret_val = true;
                        return true;
                      }
                  }
              }
          }
        return ret_val;
      }

      void communicate_field_data(std::vector< const stk::mesh::FieldBase *> fields)
      {
          stk::mesh::communicate_field_data(Base::m_eMesh.get_bulk_data()->aura_ghosting(), fields);
      }

      void downgrade_marked_transition_elements_and_upgrade_parents(std::vector<stk::mesh::Entity>& vec, int &nchange_1, bool &did_change)
      {
          PerceptMesh& eMesh = Base::m_eMesh;
          std::vector<stk::mesh::Entity> children;

          // First loop over all locally owned elements
          // downgrade any TE's marked for refine (but upgrade parents to be refined isotropically)
          for (unsigned iElement = 0; iElement < vec.size(); iElement++) {
              stk::mesh::Entity element = vec[iElement];

              int *transition_element = stk::mesh::field_data(
                      *m_transition_element_field, element);
              int *refine_field_elem = stk::mesh::field_data(*m_refine_field,
                      element);
              if (transition_element[0] && refine_field_elem[0] >= 1) {
                  refine_field_elem[0] = 0;

                  stk::mesh::Entity parent = stk::mesh::Entity();
                  if (eMesh.hasFamilyTree(element))
                      parent = eMesh.getParent(element, true);

                  VERIFY_OP_ON(eMesh.is_valid(parent), ==, true, "valid");
                  VERIFY_OP_ON(eMesh.isGhostElement(parent), ==, false,
                          "bad parent 1");

                  int *refine_field_parent = stk::mesh::field_data(
                          *m_refine_field, parent);

                  if (0 <= refine_field_parent[0] && refine_field_parent[0] < 2) {
                      refine_field_parent[0] = 2;
                      refine_field_elem[0] = 0;

                      eMesh.getChildren(parent, children, true, false);
                      for (unsigned ii = 0; ii < children.size(); ++ii) {
                          int *refine_field_child = stk::mesh::field_data(
                                  *m_refine_field, children[ii]);
                          refine_field_child[0] = 0;
                      }

                      ++nchange_1;
                      did_change = true;
                  }
              }
          }
      }

      void gather_transition_element_leaves(std::vector<stk::mesh::Entity>& vec, std::queue<stk::mesh::Entity> &elementQ)
      {
          PerceptMesh& eMesh = Base::m_eMesh;
          for (unsigned iElement = 0; iElement < vec.size(); iElement++)
          {
              stk::mesh::Entity element = vec[iElement];

              // only leaves
              bool isParent = eMesh.hasFamilyTree(element) && eMesh.isParentElement(element, true);
              if (isParent) continue;

              int *transition_element = stk::mesh::field_data( *m_transition_element_field , element );
              if (transition_element[0])
                  elementQ.push(element);
          }
      }


      bool enforce_transition_element_consistency(int iter)
      {
        // todo: sum values of fields before/after comm field data, check they are the same for locally owned...
        // double-check locally owned

    	PerceptMesh& eMesh = Base::m_eMesh;
        bool did_change = false;
        int nchange_1 = 0, nchange_2 = 0;

        communicate_field_data({m_refine_field, m_transition_element_field});

        std::vector<stk::mesh::Entity> vec;
        stk::mesh::Selector on_locally_owned_part =  ( eMesh.get_fem_meta_data()->locally_owned_part() );
        stk::mesh::get_selected_entities(on_locally_owned_part, eMesh.get_bulk_data()->buckets(eMesh.element_rank()), vec);

        downgrade_marked_transition_elements_and_upgrade_parents(vec, nchange_1, did_change);

        communicate_field_data({m_refine_field, m_transition_element_field});

        std::queue<stk::mesh::Entity> elementQ;
        gather_transition_element_leaves(vec, elementQ);

        // now enforce no neighbors of TEs are marked without marking TE parents
        SetOfEntities selected_neighbors;

        while(elementQ.size())
          {
            stk::mesh::Entity element = elementQ.front();
            elementQ.pop();

            int *refine_field_elem = stk::mesh::field_data( *m_refine_field , element );
            if (refine_field_elem[0] > 1)  // bcarnes: what would it mean for a TE leaf element to have a value > 1?
              continue;

            bool changed = process_element(element, did_change, nchange_1, nchange_2, selected_neighbors);
            if (changed)
              {
                // FIXME for edge neigh...

                // for all node neighbors that are leaves and TE's and not yet marked, put on queue
                for (SetOfEntities::iterator neighbor  = selected_neighbors.begin();
                                             neighbor != selected_neighbors.end(); ++neighbor)
                  {
                    stk::mesh::Entity neigh = *neighbor;
                    if (neigh == element) continue;
                    if (eMesh.isGhostElement(neigh)) continue;
                    bool isNParent = eMesh.hasFamilyTree(neigh) && eMesh.isParentElement(neigh, true);
                    if (isNParent) continue;
                    int *transition_element_neigh = stk::mesh::field_data( *m_transition_element_field , neigh );
                    if (!transition_element_neigh[0]) continue;
                    int *refine_field_neigh = stk::mesh::field_data( *m_refine_field , neigh );
                    if (refine_field_neigh[0] < 2 )
                      {
                        elementQ.push(neigh);
                      }
                  }
              }
          }

        communicate_field_data({m_refine_field, m_transition_element_field});

        stk::all_reduce( eMesh.get_bulk_data()->parallel(), stk::ReduceMax<1>( &did_change ) );

        void *data = &did_change;
        this->special_processing("tea_after_single_pass_enforce_te", data);

        return did_change;
      }

      void clearMarkedParentRefineField()
      {
        PerceptMesh& eMesh = Base::m_eMesh;

        const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( eMesh.element_rank() );
        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_entity_in_bucket = bucket.size();
            for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
              {
                stk::mesh::Entity element = bucket[ientity];
                int *refine_field = stk::mesh::field_data( *m_refine_field , element );

                if (refine_field[0] == 2)
                  refine_field[0] = 0;
              }
          }
      }

      void printGeometryStats(const std::string& msg)
      {
        PerceptMesh& eMesh = Base::m_eMesh;
        if (m_debug)
          {
            if ( eMesh.get_rank()==0)
              {
                std::cout << "printGeometryStats:: checking geometry stats for " << msg << std::endl;
              }
            //std::string histogram_options = "{mesh: [edge_length, quality_edge, quality_vol_edge_ratio, volume] }";
            std::string histogram_options = "{mesh: [edge_length, quality_edge, quality_vol_edge_ratio, volume] }";

            Histograms<double> histograms;
            HistogramsParser<double> hparser(histogram_options);
            hparser.create(histograms);

            eMesh.mesh_field_stats(&histograms);
            histograms.compute_uniform_bins(10);

            if (!eMesh.get_rank()) {
              std::cout << "user-requested histograms= " << std::endl;
              histograms.print(true);
            }

          }
      }

      virtual void refine()
      {
#if USE_ADAPT_PERFORMANCE_TESTING_CALLGRIND
  CALLGRIND_START_INSTRUMENTATION;
  CALLGRIND_TOGGLE_COLLECT;
#endif

        stk::diag::Timer timerRefine_("TE_Refine", Base::rootTimer());
        stk::diag::TimeBlock timerRefineBlock_(timerRefine_);
        stk::diag::Timer *altTimer = Base::getAlternateRootTimer();
        if (altTimer && DO_ALT_TIMER)
          Base::setAlternateRootTimer(&timerRefine_);

        PerceptMesh& eMesh = Base::m_eMesh;

          {
            {
              TIMER2(TER_ReqSides_CPO,Refine_);
              Base::check_parent_ownership();
            }
            bool csop = false;
            {
              TIMER2(TER_ReqSides_CSOP,Refine_);
              csop = Base::check_sides_on_same_proc_as_owned_element("TEA:start refine", false);
            }
            if (csop)
              {
                if (eMesh.get_rank()==0) std::cout << eMesh.rank() << " check_sides_on_same_proc_as_owned_element: "
                                                   << " failed, doing fix_side_sets_2" << std::endl;
                {
                  TIMER2(FSS_TER_ReqSides,Refine_);
                  Base::fix_side_sets_2(false,0,0,0,"TEARef");
                }
                Base::mod_end(&timerRefine_,"TEARef");
                {
                  TIMER2(TER_ReqSides_RSO,Refine_);
                  Base::require_sides_on_same_proc_as_pos_perm_element();
                }
              }
            {
              TIMER2(TER_ReqSides2,Refine_);
              Base::check_sides_on_same_proc_as_owned_element("TEA:start refine after fix_side_sets_2", true);
            }
          }

        Base::m_removeFromNewNodesPart = false;
        this->special_processing("tea_pre_refine");

        START_TIMER(Refine);
        if (m_debug)
          {
            std::vector<double> load_factor;
            eMesh.get_load_factor(load_factor, true, "before refine");
          }

        save("before-refine");

        printGeometryStats("refine start");

        if (!m_transition_element_field_set)
          {
            m_transition_element_field = eMesh.get_transition_element_field();
            VERIFY_OP_ON(m_transition_element_field, !=, 0, "must have transition element field");
            m_transition_element_field_set = true;
          }

        if (!m_refine_field_set)
          {
            m_refine_field_set = true;
            m_refine_field = eMesh.get_fem_meta_data()->template get_field<RefineFieldType::value_type>(stk::topology::ELEMENT_RANK, "refine_field");
          }

        {
          START_TIMER(0TER_Verifier);
          if (m_debug) {

            if (NULL == m_adaptedMeshVerifier) {
              TIMER2(TER_Verifier,Refine_);
              m_adaptedMeshVerifier = new percept::AdaptedMeshVerifier(m_debug);
              if (!m_adaptedMeshVerifier->isValid(Base::m_eMesh, true))
                throw std::runtime_error("Invalid initial mesh");
            }
          }
          STOP_TIMER(0TER_Verifier);
        }

        if (m_debug)
          {
            START_TIMER(0Refine_getOrCheckTransitionElementSet);
            getOrCheckTransitionElementSet("before refine");
            STOP_TIMER(0Refine_getOrCheckTransitionElementSet);
          }

        if (m_debug)
          print_counts("before enforce te consistency");


        // ensure no remaining marks from prior run
        {
          TIMER2(TER_clearTEField,Refine_);
          clearMarkedParentRefineField();
        }

        this->special_processing("tea_before_enforce_te_consistency");

          {
            {
              TIMER2(TER_EnforceTE,Refine_);
              START_TIMER(0TER_EnforceTE);
              int max_iter=100;
              int iter=0;
              bool did_change=false;
              while ((iter++ < max_iter) && (did_change = this->enforce_transition_element_consistency(iter)) )
                {
                  if (m_debug && eMesh.get_rank()==0)
                    std::cout << "P[" << Base::m_eMesh.get_rank() << " iter= " << iter << " did_change= " << did_change
                              << std::endl;
                  if (iter > max_iter-5)
                    {
                      throw std::runtime_error("too many iterations");
                    }
                }
              if (m_debug)
                print_counts("after ute-globalIter");
              if (m_debug)
                if (Base::m_eMesh.get_rank()==0)
                  std::cout << "P[" << Base::m_eMesh.get_rank() << "] TransitionElementAdapter::refine took niter= " << iter << std::endl;

              if (m_debug)
                {
                  TIMER2(TER_checkTE,Refine_);
                  bool vv = checkTransitionElementsAfterEnforce();
                  if (eMesh.get_rank() == 0)
                    std::cout << "checkTransitionElements after enforce = " << vv << std::endl;
                  VERIFY_OP_ON(vv,==,true,"checkTransitionElements after enforce");
                }

              save("after-ref-enforce");
              STOP_TIMER(0TER_EnforceTE);
            }
            this->special_processing("tea_after_enforce_te_consistency");

            {
              TIMER3(TER_doMark0,Refine_);
              START_TIMER(0TER_doMark0);
              Base:: m_predicate_refine.setRefineStage(-10);
              Base::m_predicate_refine.setMarkNone(false);

              Base:: initializeRefine();

              Base::doMark(1);
              Base:: m_predicate_refine.setRefineStage(0);

              Base::setRefinerSelector(0);
              STOP_TIMER(0TER_doMark0);
            }

            if (m_debug)
              {
                START_TIMER(0Refine_getOrCheckTransitionElementSet);
                getOrCheckTransitionElementSet("before ute");
                STOP_TIMER(0Refine_getOrCheckTransitionElementSet);
              }

            save1("after-mark","0after-mark.e");

            {
              {
                TIMER3(TER_RefUTE,Refine_);
                START_TIMER(0TER_RefUTE);
                Base::setNeedsRemesh(true);
                s_do_transition_break = true;
                Base::setAvoidFixSideSets(false);
                Base::setAvoidFixSideSetChecks(true);
                Base::setAvoidClearDanglingNodes(true);
                unrefineTransitionElements();
                Base::setAvoidFixSideSetChecks(false);
                Base::setAvoidClearDanglingNodes(false);
                save("after-ute-before-ref");
                STOP_TIMER(0TER_RefUTE);
              }

              {
                TIMER3(TEA_Refine1,Refine_);
                START_TIMER(0TER_Refine1);
                clearMarkedParentRefineField();
                Base::m_removeFromNewNodesPart = true;
                Base::doRefine();
                Base::setRefinerSelector(0);
                if (0 && Base::m_eMesh.get_rank() == 0)
                  Base::getRefinementInfo().printTable(std::cerr, 0, false, "TEA_Refine1");
                STOP_TIMER(0TER_Refine1);
              }

              Base:: m_predicate_refine.setRefineStage(-1);

              s_do_transition_break = false;
              Base::setAvoidClearDanglingNodes(false);
              Base::setNeedsRemesh(false);
              Base::setAvoidFixSideSets(false);
              save("after-ute");
            }
            save1("after-ute","0after-ute.e");

            if (m_debug)
              {
                TIMER2(TER_checkTE,Refine_);
                START_TIMER(0TER_checkTE);
                bool vv = checkTransitionElements();
                if (eMesh.get_rank() == 0)
                  std::cout << "checkTransitionElements after ute = " << vv << std::endl;
                VERIFY_OP_ON(vv,==,true,"checkTransitionElements after ute");
                STOP_TIMER(0TER_checkTE);
              }

            {
              TIMER2(TER_Verifier,Refine_);
              START_TIMER(0TER_Verifier);
              if (m_adaptedMeshVerifier && !m_adaptedMeshVerifier->isValid(Base::m_eMesh, false))
                {
                  throw std::runtime_error("Invalid mesh after UTE/doRefine");
                }
              STOP_TIMER(0TER_Verifier);
            }

            if (m_debug)
              {
                START_TIMER(0Refine_getOrCheckTransitionElementSet);
                getOrCheckTransitionElementSet("after ute ref 1");
                STOP_TIMER(0Refine_getOrCheckTransitionElementSet);
              }

          }

        if (m_debug)
          {
            TIMER2(TER_checkTE,Refine_);
            START_TIMER(0TER_checkTE);
            bool vv = checkTransitionElements(true);
            if (m_debug && eMesh.get_rank() == 0)
              std::cout << "checkTransitionElements after ute converged = " << vv << std::endl;
            if (!vv)
              {
                if (eMesh.get_rank() == 0)
                  std::cout << "bad checkTransitionElements after ute converged = " << vv << std::endl;
                throw std::runtime_error("bad checkTransitionElements after ute converged");
              }
            STOP_TIMER(0TER_checkTE);
          }

        s_do_transition_break = true;
        Base:: m_predicate_refine.setMarkNone(false);
        Base:: m_predicate_refine.setRefineStage(0);

        save1("after-ref-doRefine2","0after-ref-doRefine2.e");

        {
          TIMER2(TER_Verifier, Refine_);
          START_TIMER(0TER_Verifier);
          if (m_adaptedMeshVerifier && !m_adaptedMeshVerifier->isValid(Base::m_eMesh, false))
            {
              throw std::runtime_error("Invalid mesh after refine doRefine2");
            }
          STOP_TIMER(0TER_Verifier);
        }

        {
          START_TIMER(0Refine_getOrCheckTransitionElementSet);
          if (m_debug)
            getOrCheckTransitionElementSet("after ref 2");
          STOP_TIMER(0Refine_getOrCheckTransitionElementSet);
        }

        printGeometryStats("refine end");

        STOP_TIMER(Refine);
        if (m_debug)
          {
            std::vector<double> load_factor;
            eMesh.get_load_factor(load_factor, true, "after refine");
          }

        this->special_processing("tea_post_refine");

        if (1)
          {
            TIMER2(TER_CheckSides, Refine_);
            //Base::check_parent_ownership();
            Base::check_sides_on_same_proc_as_owned_element("TEA:end refine", true);
          }

        if (m_debug)
          {
            stk::ParallelMachine pm = Base::m_eMesh.get_bulk_data()->parallel();

            for (std::map<std::string, double>::iterator it = timers_.begin();
                 it != timers_.end(); ++it)
              {
                stk::all_reduce( pm, stk::ReduceSum<1>( &it->second ));
              }

            if (eMesh.get_rank() == 0)
              {
                std::cout << "\nTransitionElementAdapter::refine timing info: " << std::endl;
                double time_main = timers_["Refine"];
                for (std::map<std::string, double>::iterator it = timers_.begin();
                     it != timers_.end(); ++it)
                  {
                    std::string a1 = it->first;
                    std::string a0 = a1;
                    if (a1[0] == '0') a1[0] = '\t';
                    std::cout << (a0[0] == '0' ? '\t' : ' ') << (it->second/time_main)*100 << "% " << it->second << "\t"  << a1 << std::endl;
                  }
              }
          }

        if (altTimer && DO_ALT_TIMER)
          Base::setAlternateRootTimer(altTimer);

#if USE_ADAPT_PERFORMANCE_TESTING_CALLGRIND
  CALLGRIND_TOGGLE_COLLECT;
  CALLGRIND_STOP_INSTRUMENTATION;
#endif

      }

      size_t upgrade_element_siblings(ElementUnrefineCollection& elements_to_unref)
      {
        PerceptMesh& eMesh = Base::m_eMesh;
        size_t nmod_elements = 0;
        ElementUnrefineCollection enew = elements_to_unref;
        for (ElementUnrefineCollection::iterator it = elements_to_unref.begin(); it != elements_to_unref.end(); ++it)
          {
            stk::mesh::Entity element = *it;
            static std::vector<stk::mesh::Entity> children;
            stk::mesh::Entity parent = eMesh.getParent(element, true);
            VERIFY_OP_ON(eMesh.is_valid(parent), ==, true, "found unrefined element with no parent");
            bool hasChildren = eMesh.getChildren(parent, children, true, false);
            (void)hasChildren;
            for (unsigned ii=0; ii < children.size(); ++ii)
              {
                if (element != children[ii])
                  {
                    if (enew.find(children[ii]) == enew.end())
                      ++nmod_elements;
                    enew.insert(children[ii]);
                  }
              }
          }
        elements_to_unref = enew;
        return nmod_elements;
      }

      virtual void unrefine() {

        stk::diag::Timer timerUnrefine_("TE_UnRef", Base::rootTimer());
        stk::diag::TimeBlock timerUnrefineBlock_(timerUnrefine_);

        stk::diag::Timer *altTimer = Base::getAlternateRootTimer();
        if (altTimer && DO_ALT_TIMER)
          Base::setAlternateRootTimer(&timerUnrefine_);

        save1("before-unrefine","0before-unref.e");

        this->special_processing("tea_pre_unrefine");

        printGeometryStats("unrefine start");

        IAdapter& breaker = *this;
        Base::setNeedsRemesh(true);
        ElementUnrefineCollection elements_to_unref(*Base::m_eMesh.get_bulk_data());
        breaker.buildUnrefineList(elements_to_unref);

        if (1)
          {
            size_t nmod_elements = upgrade_element_siblings(elements_to_unref);
            if (m_debug && 0 == Base::m_eMesh.get_rank())
              std::cout << "upgrade_element_siblings::nmod_elements= " << nmod_elements << std::endl;
          }
        ElementUnrefineCollection enew(*Base::m_eMesh.get_bulk_data());
        getOrCheckTransitionElementSet("after unref ref 1", &enew);
        elements_to_unref.insert(enew.begin(), enew.end());
        Base::setAvoidClearDanglingNodes(false);
        Base::m_onlyOneLevelUnrefine = true;
        breaker.unrefineTheseElements(elements_to_unref);
        Base::m_onlyOneLevelUnrefine = false;


        if (1)
          {
            Base::mod_begin(&timerUnrefine_);
            Base::removeDanglingNodes();
            Base::mod_end(&timerUnrefine_,"TEARemoveDN");
          }

        if (m_adaptedMeshVerifier && !m_adaptedMeshVerifier->isValid(Base::m_eMesh, false))
          {
            throw std::runtime_error("Invalid mesh after unrefine");
          }
        printGeometryStats("unrefine end");
        save1("after-unrefine","0after-unref.e");

        this->special_processing("tea_post_unrefine");

        if (altTimer && DO_ALT_TIMER)
          Base::setAlternateRootTimer(altTimer);
      }

      virtual void
      doBreak(int num_registration_loops=1)
      {
        refine();
      }

      virtual void
      unrefineTheseElements(ElementUnrefineCollection& elements_to_unref)
      {
        Base::unrefineTheseElements(elements_to_unref);
      }

      void reportParentTransitionElementErrors(stk::mesh::Entity element, const std::string& msg)
      {
          PerceptMesh& eMesh = Base::m_eMesh;

          const bool check_for_family_tree = false;
          bool isParent = eMesh.isParentElement(element, check_for_family_tree);

          if (!isParent) return;

          std::cout << "P[" << eMesh.get_rank() << "] <" << msg << "> found bad element - transition_element found that is a parent "
                    << eMesh.identifier(element)
                    << " isParent= " << isParent
                    << " isGhostElement = " << eMesh.isGhostElement(element)
                    << "\n Note: This may be due to not marking all elements including ghosts - check stk::mesh::Selector usage.\n"
                    << std::endl;
          eMesh.print(element);
          {
            std::string file = "te.err.dat."+toString(eMesh.get_bulk_data()->parallel_size()) + "."+toString(eMesh.get_rank())+".vtk";
            std::set<stk::mesh::Entity> elem_set_debug;
            static std::vector<stk::mesh::Entity> children;
            elem_set_debug.insert(element);
            stk::mesh::Entity parentOfElem = eMesh.getParent(element, true);
            if (eMesh.is_valid(parentOfElem))
              {
                std::cout << "P[" << eMesh.get_rank() << "] parentOfElem is:  " << std::endl;
                eMesh.print(parentOfElem);
                std::cout << "P[" << eMesh.get_rank() << "] my siblings= " << std::endl;
                eMesh.getChildren(parentOfElem, children, true, false);
                for (unsigned ich=0; ich < children.size(); ++ich)
                  {
                    eMesh.print(children[ich]);
                    elem_set_debug.insert(children[ich]);
                  }

              }
            else
              std::cout << "P[" << eMesh.get_rank() << "] parentOfElem is null " << std::endl;
            std::cout << "P[" << eMesh.get_rank() << "] children= " << std::endl;
            eMesh.getChildren(element, children, true, false);
            for (unsigned ich=0; ich < children.size(); ++ich)
              {
                eMesh.print(children[ich]);
                elem_set_debug.insert(children[ich]);
              }
            eMesh.dump_vtk(file, false, &elem_set_debug);
          }

          throw std::runtime_error("bad element - transition_element found that is a parent");
      }

      void getOrCheckTransitionElementSet(const std::string& msg, ElementUnrefineCollection * elements_to_unref=0)
      {
        PerceptMesh& eMesh = Base::m_eMesh;
        TransitionElementType *transition_element_field = m_transition_element_field;

        const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( eMesh.element_rank() );
        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_entity_in_bucket = bucket.size();
            for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
              {
                stk::mesh::Entity element = bucket[ientity];

                int *transition_element = 0;
                transition_element = stk::mesh::field_data( *transition_element_field , element );

                if (transition_element[0])
                  {
                	reportParentTransitionElementErrors(element, msg);

                     if (elements_to_unref) {
                        elements_to_unref->insert(element);
                    }
                  }
              }
          }
      }

      void checkNoTransitionElementsRemain(const std::string& msg)
      {
        PerceptMesh& eMesh = Base::m_eMesh;
        TransitionElementType *transition_element_field = m_transition_element_field;

        const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( eMesh.element_rank() );
        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_entity_in_bucket = bucket.size();
            for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
              {
                stk::mesh::Entity element = bucket[ientity];
                const bool check_for_family_tree = false;
                bool isParent = eMesh.isParentElement(element, check_for_family_tree);
                (void)isParent;

                int *transition_element = 0;
                transition_element = stk::mesh::field_data( *transition_element_field , element );
                //% ghost

                if (transition_element[0])
                  {
                    throw std::runtime_error("Found transition element "+msg);
                  }
              }
          }
      }

      void unrefineTransitionElements()
      {
        PerceptMesh& eMesh = Base::m_eMesh;
        ElementUnrefineCollection elements_to_unref(*eMesh.get_bulk_data());
        getOrCheckTransitionElementSet("in unrefineTransitionElements", &elements_to_unref);
        Base::m_onlyOneLevelUnrefine = true;
        Refiner::unrefineTheseElements(elements_to_unref);
        Base::m_onlyOneLevelUnrefine = false;
        if (m_debug)
          getOrCheckTransitionElementSet("end of unrefineTransitionElements");
      }

      bool checkTransitionElements(bool checkMarks=false)
      {
        bool val = checkTransitionElements1(checkMarks);

        stk::ParallelMachine pm = Base::m_eMesh.get_bulk_data()->parallel();
        stk::all_reduce( pm, stk::ReduceMin<1>( &val ));

        return val;
      }

      bool checkTransitionElements1(bool checkMarks=false)
      {
        PerceptMesh& eMesh = Base::m_eMesh;
        bool val = true;
        SubDimCell_SDCEntityType subDimEntity(&eMesh);

        TransitionElementType *transition_element_field = m_transition_element_field;

        {
          std::vector< const stk::mesh::FieldBase *> fields;
          fields.push_back(m_refine_field);
          fields.push_back(transition_element_field);
          stk::mesh::communicate_field_data(eMesh.get_bulk_data()->aura_ghosting(), fields);
        }

        LocalSetType neighbors;
        const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( eMesh.element_rank() );
        stk::mesh::Selector on_locally_owned_part =  ( eMesh.get_fem_meta_data()->locally_owned_part() );

        size_t nbad_elements=0;
        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (on_locally_owned_part(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              const CellTopologyData * const bucket_topo_data = eMesh.get_cell_topology(bucket);
              const unsigned num_elements_in_bucket = bucket.size();
              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  stk::mesh::Entity element = bucket[iElement];

                  bool isParent = (eMesh.hasFamilyTree(element) && eMesh.isParentElement(element));
                  if (isParent)
                    continue;

                  int *transition_element = 0;
                  transition_element = stk::mesh::field_data( *transition_element_field , element );
                  int *refine_field_elem = stk::mesh::field_data( *m_refine_field , element );
                  bool isTransition = transition_element[0] > 0;
                  if (isTransition)
                    {
                      if (refine_field_elem[0] >= 1)
                        {
                          throw std::runtime_error("checkTransitionElements1: bad transition element mark");
                        }
                      bool local_val = true;

                      neighbors.clear();
                      Base::get_node_neighbors(element, neighbors);

                      static std::vector<stk::mesh::Entity> children;
                      stk::mesh::Entity parent = eMesh.getParent(element, true);
                      bool hasChildren = eMesh.getChildren(parent, children, true, false);
                      (void)hasChildren;

                      LocalSetType::iterator neighbor_begin = neighbors.begin();
                      LocalSetType::iterator neighbor_end = neighbors.end();
                      for (LocalSetType::iterator it = neighbor_begin; it != neighbor_end; ++it)
                        {
                          stk::mesh::Entity neigh = *it;

                          bool isNParent = (eMesh.hasFamilyTree(neigh) && eMesh.isParentElement(neigh));

                          int *refine_field_neigh = stk::mesh::field_data( *m_refine_field , neigh );

                          if (refine_field_neigh[0] >= 1)
                            {

                              bool isEdgeN = false; // FIXME
                              int edge_0=-1, edge_1=-1;
                              isEdgeN = eMesh.is_edge_neighbor(element, neigh, &edge_0, &edge_1, 0);

                              if (isEdgeN)
                                {

                                  if (isEdgeN)  // FIXME
                                    {
                                      int jedge_1 = -1;
                                      if (s_allow_skip_edge_neighbor)
                                        {
                                          bool ien = is_edge_neighbor_across_given_edge(element, parent, edge_0, &jedge_1, 0);
                                          if (ien)
                                            continue;
                                        }
                                    }

                                  ++nbad_elements;
                                  val = false;
                                  local_val = false;
                                  transition_element = stk::mesh::field_data( *transition_element_field , neigh );
                                  if (1)
                                    std::cout << "bad element = " << eMesh.print_entity_compact(element) << "\n neigh = " << eMesh.print_entity_compact(neigh)
                                              << " transition_element= " << transition_element[0]
                                              << " refine_field_neigh= " << refine_field_neigh[0]
                                              << " isNParent= " << isNParent
                                              << std::endl;
                                }
                            }
                        }
                      if (local_val)
                        {
                          unsigned nedges = (unsigned)bucket_topo_data->edge_count;
                          for (unsigned iSubDimOrd = 0; iSubDimOrd < nedges; ++iSubDimOrd)
                            {
                              Base::getNodeRegistry().getSubDimEntity(subDimEntity, element, eMesh.edge_rank(), iSubDimOrd, bucket_topo_data);

                              static SubDimCellData empty_SubDimCellData;

                              SubDimCellData* nodeId_elementOwnderId_ptr = Base::getNodeRegistry().getFromMapPtr(subDimEntity);
                              SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
                              bool is_empty = nodeId_elementOwnderId_ptr == 0;
                              bool is_marked = (!is_empty && std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnderId).size() != 0);
                              if (is_marked)
                                throw std::runtime_error("bad mark in checkTransitionElements1");
                            }
                        }
                    }
                }
            }
          }
        if (m_debug)
          {
            stk::ParallelMachine pm = Base::m_eMesh.get_bulk_data()->parallel();
            stk::all_reduce( pm, stk::ReduceSum<1>( &nbad_elements ));
            if (eMesh.get_rank()==0) std::cout << "checkTransitionElements1: nbad_elements = " << nbad_elements << std::endl;
          }

        return val;
      }

      bool checkTransitionElementsAfterEnforce(bool checkMarks=false)
      {
        bool val = checkTransitionElementsAfterEnforce1(checkMarks);

        stk::ParallelMachine pm = Base::m_eMesh.get_bulk_data()->parallel();
        stk::all_reduce( pm, stk::ReduceSum<1>( &val ));

        return val;
      }

      bool checkTransitionElementsAfterEnforce1(bool checkMarks=false)
      {
        PerceptMesh& eMesh = Base::m_eMesh;
        bool val = true;
        SubDimCell_SDCEntityType subDimEntity(&eMesh);

        TransitionElementType *transition_element_field = m_transition_element_field;

        LocalSetType neighbors;

        {
          std::vector< const stk::mesh::FieldBase *> fields;
          fields.push_back(m_refine_field);
          fields.push_back(transition_element_field);
          stk::mesh::communicate_field_data(eMesh.get_bulk_data()->aura_ghosting(), fields);
        }

        const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( eMesh.element_rank() );
        stk::mesh::Selector on_locally_owned_part =  ( eMesh.get_fem_meta_data()->locally_owned_part() );

        size_t nbad_elements=0;
        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (on_locally_owned_part(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              //const CellTopologyData * const bucket_topo_data = eMesh.get_cell_topology(bucket);
              const unsigned num_elements_in_bucket = bucket.size();
              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  stk::mesh::Entity element = bucket[iElement];

                  bool isParent = (eMesh.hasFamilyTree(element) && eMesh.isParentElement(element));
                  if (isParent)
                    continue;

                  int *transition_element = 0;
                  transition_element = stk::mesh::field_data( *transition_element_field , element );
                  int *refine_field_elem = stk::mesh::field_data( *m_refine_field , element );
                  bool isTransition = transition_element[0] > 0;
                  if (isTransition)
                    {
                      if (refine_field_elem[0] >= 1)
                        {
                          throw std::runtime_error("checkTransitionElements1: bad transition element mark");
                        }

                      neighbors.clear();
                      Base::get_node_neighbors(element, neighbors);

                      static std::vector<stk::mesh::Entity> children;
                      stk::mesh::Entity parent = eMesh.getParent(element, true);
                      bool hasChildren = eMesh.getChildren(parent, children, true, false);
                      (void)hasChildren;

                      LocalSetType::iterator neighbor_begin = neighbors.begin();
                      LocalSetType::iterator neighbor_end = neighbors.end();
                      for (LocalSetType::iterator it = neighbor_begin; it != neighbor_end; ++it)
                        {
                          stk::mesh::Entity neigh = *it;

                          bool isNParent = (eMesh.hasFamilyTree(neigh) && eMesh.isParentElement(neigh));

                          int *refine_field_neigh = stk::mesh::field_data( *m_refine_field , neigh );

                          if (refine_field_neigh[0] >= 1)
                            {
                              bool isEdgeN = false;
                              int edge_0=-1, edge_1=-1;
                              isEdgeN = eMesh.is_edge_neighbor(element, neigh, &edge_0, &edge_1, 0);

                              if (isEdgeN)
                                {

                                    {
                                      int jedge_1 = -1;
                                      if (s_allow_skip_edge_neighbor)
                                        {
                                          bool ien = is_edge_neighbor_across_given_edge(element, parent, edge_0, &jedge_1, 0);
                                          if (ien)
                                            continue;
                                        }
                                    }

                                  int *refine_field_parent = stk::mesh::field_data( *m_refine_field , parent );

                                  if (refine_field_parent[0] <= 1)
                                    {
                                      ++nbad_elements;
                                      val = false;
                                      transition_element = stk::mesh::field_data( *transition_element_field , neigh );
                                      LocalSetType selected_neighbors_1;
                                      bool did_change_1 = false;
                                      int nchange_1=0, nchange_2=0;
                                      bool pe = process_element(element, did_change_1, nchange_1, nchange_2, selected_neighbors_1, true);
                                      if (0)
                                        std::cout << "bad element = " << eMesh.identifier(element) << " neigh = " << eMesh.identifier(neigh)
                                                  << " transition_element= " << transition_element[0]
                                                  << " refine_field_neigh= " << refine_field_neigh[0]
                                                  << " isNParent= " << isNParent
                                                  << std::endl;
                                      VERIFY_OP_ON(pe, ==, false, "Bad checkTransitionElementsAfterEnforce1");
                                    }
                                }
                            }
                        }
                    }
                }
            }
          }
        if (m_debug)
          {
            stk::ParallelMachine pm = Base::m_eMesh.get_bulk_data()->parallel();
            stk::all_reduce( pm, stk::ReduceSum<1>( &nbad_elements ));
            if (eMesh.get_rank()==0) std::cout << "checkTransitionElements1: nbad_elements = " << nbad_elements << std::endl;
          }

        return val;
      }


    };


  }

#endif
