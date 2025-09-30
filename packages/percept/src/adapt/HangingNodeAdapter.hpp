// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_HangingNodeAdapter_hpp
#define adapt_HangingNodeAdapter_hpp

#include <functional>
#include <cmath>
#include <unordered_set>

#include <adapt/PredicateBasedElementAdapter.hpp>
#include <percept/SetOfEntities.hpp>
#include <stk_mesh/base/HashEntityAndEntityKey.hpp>

#include <percept/pooled_alloc.h>

#include <stk_mesh/base/FieldParallel.hpp>

#define HNA_TIMING(code) do { code } while(0)
#define HNA_TIMER(name) stk::diag::Timer timer ## name ( #name, Base::rootTimer());  stk::diag::TimeBlock tbTimer ## name (timer ## name)
#define HNA_TIMER2(name,parentName) stk::diag::Timer timer ## name ( #name, timer ## parentName);  stk::diag::TimeBlock tbTimer ## name (timer ## name)
#define HNA_TIMER2A(name,parentName) stk::diag::Timer timer ## name ( #name, parentName);  stk::diag::TimeBlock tbTimer ## name (timer ## name)

  namespace percept {

    template<typename Value, typename Less = std::less<Value>, typename Allocator = std::allocator<Value> >
    struct PooledStdSet
    {
      typedef std::set<Value, Less, Allocator> Type;
      typedef Less less;
    };

    typedef SetOfEntities LocalSetType;

    class HangingNodeAdapterBase
    {

    public:

      HangingNodeAdapterBase(percept::PerceptMesh& eMesh) : m_pMesh(eMesh), m_hn_transition_element_field(m_pMesh.get_transition_element_field())
      {
        for (int i=0; i < 3; i++)
          m_will_refine[i] = false;
        m_will_refine[eMesh.side_rank()] = true;
      }

      /// optionally mark sidesets/edges/faces for refinement
      /// @param will_refine: {node-, edge-, face-} neighbors
      void set_sub_dim_refine_ranks(const bool will_refine[3]);
      void get_sub_dim_refine_ranks(bool will_refine[3]);

      // low-level interface:

    protected:

      /// wrappers around refine and unrefine
      void refine(IAdapter& breaker);
      void unrefine(IAdapter& breaker);


      /** after marking using field "refine_field", revisit all elements and upgrade
       * any that need it to enforce two to one
       * @return true if there was a change to refine_field; useful in a while(did_change) loop
       */
      bool enforce_two_to_one_refine();

      /** after marking using field "refine_field", revisit all elements and downgrade
       * any that need it to enforce two to one during unrefinement
       * @return true if there was a change to refine_field; useful in a while(did_change) loop
       */
      bool enforce_two_to_one_unrefine();

      bool hasChildren(stk::mesh::Entity element);

    public:
      void get_node_neighbors(stk::mesh::Entity element, LocalSetType& neighbors);

    protected:

      bool is_face_neighbor(stk::mesh::Entity element, int element_level, stk::mesh::Entity neighbor, int neighbor_level,
                            bool check_level_early_return = true,
                            const CellTopologyData * const bucket_topo_data_element_0 =0
                            );
      bool is_edge_neighbor(stk::mesh::Entity element, int element_level, stk::mesh::Entity neighbor, int neighbor_level,
                            bool check_level_early_return = true,
                            const CellTopologyData * const bucket_topo_data_element_0 =0
                            );
      bool is_node_neighbor(stk::mesh::Entity element, int element_level, stk::mesh::Entity neighbor, int neighbor_level,
                            bool check_level_early_return = true
                            );
      //bool min_max_neighbors_level(stk::mesh::Entity element, int min_max[2], TransitionElementType *refine_level, bool check_what[3] );
      void get_neighbors(stk::mesh::Entity element, RefineLevelType *refine_level, bool get_what[3],
                         std::set<stk::mesh::Entity>& selected_neighbors,
                         const CellTopologyData * const bucket_topo_data =0
                         );
      bool should_check(stk::mesh::Entity element, int refine_level_elem,
                        stk::mesh::Entity neighbor, int refine_level_neighbor,
                        const CellTopologyData * const bucket_topo_data
                        );

    public:
      percept::PerceptMesh& m_pMesh;
    protected:
      TransitionElementType *m_hn_transition_element_field;

    protected:
      bool m_will_refine[3];
    };

    template<class RefinePredicate>
    class HangingNodeAdapter : public PredicateBasedElementAdapter<RefinePredicate>, public HangingNodeAdapterBase
    {

      typedef PredicateBasedElementAdapter<RefinePredicate> Base;

      bool m_debug_print;
      bool m_debug;
    protected:
      bool m_only_enforce_unrefine;
      bool m_transition_element_break;
    public:

      HangingNodeAdapter(RefinePredicate& predicate_refine,
                           percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0) :
        Base(predicate_refine, eMesh, bp, proc_rank_field), HangingNodeAdapterBase(eMesh), m_debug_print(false), m_debug(false),
        m_only_enforce_unrefine(false), m_transition_element_break(false)
      {
        Base::setNeedsRemesh(false);
        //m_debug = true;

      }

      /// wrappers around refine and unrefine
      virtual void refine()
      {
        IAdapter& breaker = *this;
        {
          HNA_TIMER2A(HNA_RefEnf, breaker.rootTimer());

          int max_iter=100;
          int iter=0;
          bool did_change=false;
          while ((iter++ < max_iter) && (did_change = this->enforce_two_to_one_refine()) )
            {
              if (m_debug_print)
                std::cout << "P[" << m_pMesh.get_rank() << " iter= " << iter << " did_change= " << did_change
                          << std::endl;
            }
        }

        {
          HNA_TIMER2A(HNA_RefDoBreak, breaker.rootTimer());
          breaker.doBreak();
        }

      }

      virtual void unrefine()
      {
        IAdapter& breaker = *this;
        {
          HNA_TIMER2A(HNA_UnRefEnf, breaker.rootTimer());
          int max_iter=100;
          int iter=0;
          bool did_change=false;
          while ((iter++ < max_iter) && (did_change = this->enforce_two_to_one_unrefine()) )
            {
              if (m_debug_print)
                std::cout << "P[" << m_pMesh.get_rank() << "]  iter= " << iter
                          << " did_change= " << did_change
                          << std::endl;
            }

        }

        if (m_only_enforce_unrefine)
          return;

        ElementUnrefineCollection elements_to_unref(*m_pMesh.get_bulk_data());
        {
          HNA_TIMER2A(HNA_UnRefBuild, breaker.rootTimer());
          breaker.buildUnrefineList(elements_to_unref);
        }
        {
          HNA_TIMER2A(HNA_UnRefUnRef, breaker.rootTimer());
          breaker.unrefineTheseElements(elements_to_unref);
        }
      }

      static const bool ters_active = true;

      bool enforce_two_to_one_refine()
      {
        stk::ParallelMachine pm = m_pMesh.get_bulk_data()->parallel();
        bool did_change = false;
        RefineLevelType *refine_level = m_pMesh.get_refine_level_field();
        if (!refine_level)
          {
            throw std::logic_error("must have refine_level field for hanging-node refinement");
          }
        RefineFieldType *refine_field = m_pMesh.get_fem_meta_data()-> template get_field<RefineFieldType::value_type>(stk::topology::ELEMENT_RANK, "refine_field");
        if (!refine_field)
          {
            throw std::logic_error("must have refine_field field for hanging-node refinement");
          }

        {
          std::vector< const stk::mesh::FieldBase *> fields;
          fields.push_back(refine_level);
          fields.push_back(refine_field);
          //stk::mesh::copy_owned_to_shared( *m_pMesh.get_bulk_data(), fields);
          stk::mesh::communicate_field_data(m_pMesh.get_bulk_data()->aura_ghosting(), fields);

        }

        if (0 && m_debug_print)
          {
            std::ostringstream ostr;
            stk::mesh::EntityId id[]={12,13};
            for (unsigned i=0; i < 2; i++)
              {
                stk::mesh::Entity element = m_pMesh.get_entity(m_pMesh.element_rank(), id[i]);
                if (m_pMesh.is_valid(element))
                  {
                    int *refine_level_elem = stk::mesh::field_data( *refine_level , element );
                    int *refine_field_elem = stk::mesh::field_data( *refine_field , element );
                    ostr << "P[" << m_pMesh.get_rank() << "] elem " << id[i] << " lev= " << refine_level_elem[0]
                         << " ref_field= " << refine_field_elem[0] << "\n";
                  }
              }
            std::cout << ostr.str() << std::endl;
          }

        SelectIfRefined sir(*this);
        HNRefinerSelector<SelectIfRefined> ters(*this, sir);

        stk::mesh::Selector on_locally_owned_part =  ( m_pMesh.get_fem_meta_data()->locally_owned_part() );
        const stk::mesh::BucketVector & buckets = m_pMesh.get_bulk_data()->buckets( m_pMesh.element_rank() );

        // count and print
        if (0)
          {
            unsigned nref=0;
            unsigned nele=0;
            for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                if (on_locally_owned_part(**k))
                  {
                    stk::mesh::Bucket & bucket = **k ;
                    const unsigned num_elements_in_bucket = bucket.size();
                    for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                      {
                        stk::mesh::Entity element = bucket[iElement];

                        if (hasChildren(element))
                          continue;
                        ++nele;

                        //int *refine_level_elem = stk::mesh::field_data( *refine_level , element );
                        int *refine_field_elem = stk::mesh::field_data( *refine_field , element );
                        if (refine_field_elem[0] >= 1)
                          {
                            ++nref;
                          }
                      }
                  }
              }
            //if (m_pMesh.get_rank()==0)
            std::cout << "P[" << m_pMesh.get_rank() << "] enforce_two_to_one_refine:: nref0 for this iter= " << nref << " nele= " << nele << std::endl;
            stk::all_reduce( pm, stk::ReduceSum<1>( &nref ) );
            stk::all_reduce( pm, stk::ReduceSum<1>( &nele ) );
            //if (m_pMesh.get_rank()==0)
            std::cout << "P[" << m_pMesh.get_rank() << "] enforce_two_to_one_refine:: nref for this iter= " << nref << " nele= " << nele << std::endl;
          }

        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (1 || on_locally_owned_part(**k))
              {
                stk::mesh::Bucket & bucket = **k ;
                const CellTopologyData * const bucket_topo_data = m_pMesh.get_cell_topology(bucket);
                const unsigned num_elements_in_bucket = bucket.size();
                for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                  {
                    stk::mesh::Entity element = bucket[iElement];
                    if (ters_active && !ters(element))
                      continue;

                    if (hasChildren(element))
                      continue;

                    int *refine_level_elem = stk::mesh::field_data( *refine_level , element );
                    int *refine_field_elem = stk::mesh::field_data( *refine_field , element );
                    if (refine_field_elem[0] <= 0)
                      {
                        LocalSetType selected_neighbors;

                        get_node_neighbors(element, selected_neighbors);

                        for (LocalSetType::iterator neighbor = selected_neighbors.begin();
                             neighbor != selected_neighbors.end(); ++neighbor)
                          {
                            if (hasChildren(*neighbor))
                              continue;

                            int *refine_level_neigh = stk::mesh::field_data( *refine_level , *neighbor );
                            int *refine_field_neigh = stk::mesh::field_data( *refine_field , *neighbor );

                            // if any neighbor is my level + 1 and marked for refine,
                            //   and I am not marked for refine, mark me for refine
                            if ( (refine_level_neigh[0] > refine_level_elem[0])
                                 && refine_field_neigh[0] > 0)
                              {
                                if (should_check(element, refine_level_elem[0], *neighbor, refine_level_neigh[0],
                                                 bucket_topo_data))
                                  {
                                    refine_field_elem[0] = 1;

                                    if (m_debug_print)
                                      std::cout << "P[" << m_pMesh.get_rank() << "] enforce_two_to_one_refine: upgrading element ("
                                                << m_pMesh.identifier(element) << "," << m_pMesh.isGhostElement(element) << ")"
                                                << " due to neighbor: "
                                                << m_pMesh.identifier(*neighbor) << "," << m_pMesh.isGhostElement(*neighbor) << ")"
                                                << std::endl;
                                    did_change = true;
                                    break;
                                  }
                              }
                          }
                      }
                  }
              }
          }

        if (0 && m_debug_print)
          {
            std::ostringstream ostr;
            stk::mesh::EntityId id[]={12,13};
            unsigned nid=1;
            for (unsigned i=0; i < nid; i++)
              {
                stk::mesh::Entity element = m_pMesh.get_entity(m_pMesh.element_rank(), id[i]);
                if (m_pMesh.is_valid(element))
                  {
                    int *refine_level_elem = stk::mesh::field_data( *refine_level , element );
                    int *refine_field_elem = stk::mesh::field_data( *refine_field , element );
                    ostr << "P[" << m_pMesh.get_rank() << "] t4 b4 copy_owned_to_shared elem " << id[i] << " lev= " << refine_level_elem[0]
                         << " ref_field= " << refine_field_elem[0] << "\n";
                  }
              }
            std::cout << ostr.str() << std::endl;
          }

        {
          std::vector< const stk::mesh::FieldBase *> fields;
          fields.push_back(refine_level);
          fields.push_back(refine_field);
          //stk::mesh::copy_owned_to_shared( *m_pMesh.get_bulk_data(), fields);
          stk::mesh::communicate_field_data( m_pMesh.get_bulk_data()->aura_ghosting(), fields);

        }

        if (0 && m_debug_print)
          {
            MPI_Barrier( MPI_COMM_WORLD );
            std::ostringstream ostr;
            stk::mesh::EntityId id[]={12,13};
            unsigned nid=1;
            for (unsigned i=0; i < nid; i++)
              {
                stk::mesh::Entity element = m_pMesh.get_entity(m_pMesh.element_rank(), id[i]);
                if (m_pMesh.is_valid(element))
                  {
                    int *refine_level_elem = stk::mesh::field_data( *refine_level , element );
                    int *refine_field_elem = stk::mesh::field_data( *refine_field , element );
                    ostr << "P[" << m_pMesh.get_rank() << "] t4 after copy_owned_to_shared elem " << id[i] << " lev= " << refine_level_elem[0]
                         << " ref_field= " << refine_field_elem[0] << "\n";
                  }
              }
            std::cout << ostr.str() << std::endl;
          }

        stk::all_reduce( pm, stk::ReduceMax<1>( &did_change ) );

        // count and print
        if (0 && !did_change)
          {
            unsigned nref=0;
            for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                if (on_locally_owned_part(**k))
                  {
                    stk::mesh::Bucket & bucket = **k ;
                    const unsigned num_elements_in_bucket = bucket.size();
                    for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                      {
                        stk::mesh::Entity element = bucket[iElement];

                        if (hasChildren(element))
                          continue;

                        //int *refine_level_elem = stk::mesh::field_data( *refine_level , element );
                        int *refine_field_elem = stk::mesh::field_data( *refine_field , element );
                        if (refine_field_elem[0] >= 1)
                          {
                            ++nref;
                          }
                      }
                  }
              }
            stk::all_reduce( pm, stk::ReduceSum<1>( &nref ) );
            if (m_pMesh.get_rank()==0)
              std::cout << "P[" << m_pMesh.get_rank() << "] enforce_two_to_one_refine:: nref1 for this iter= " << nref << std::endl;
          }

        return did_change;
      }

      bool enforce_two_to_one_unrefine()
      {
        bool did_change = false;
        RefineLevelType *refine_level = m_pMesh.get_refine_level_field();
        if (!refine_level)
          {
            throw std::logic_error("must have refine_level field for hanging-node refinement");
          }
        RefineFieldType *refine_field = m_pMesh.get_fem_meta_data()-> template get_field<RefineFieldType::value_type>(stk::topology::ELEMENT_RANK, "refine_field");
        if (!refine_level)
          {
            throw std::logic_error("must have refine_field field for hanging-node refinement");
          }

        {
          std::vector< const stk::mesh::FieldBase *> fields;
          fields.push_back(refine_level);
          fields.push_back(refine_field);
          //stk::mesh::copy_owned_to_shared( *m_pMesh.get_bulk_data(), fields);
          stk::mesh::communicate_field_data(m_pMesh.get_bulk_data()->aura_ghosting(), fields);
        }

        SelectIfUnrefined sir;
        HNRefinerSelector<SelectIfUnrefined> ters(*this, sir);

        stk::mesh::Selector on_locally_owned_part =  ( m_pMesh.get_fem_meta_data()->locally_owned_part() );
        const stk::mesh::BucketVector & buckets = m_pMesh.get_bulk_data()->buckets( m_pMesh.element_rank() );
        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (1 || on_locally_owned_part(**k))
              {
                stk::mesh::Bucket & bucket = **k ;
                const CellTopologyData * const bucket_topo_data = m_pMesh.get_cell_topology(bucket);
                const unsigned num_elements_in_bucket = bucket.size();
                for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                  {
                    stk::mesh::Entity element = bucket[iElement];
                    if (ters_active && !ters(element))
                      continue;

                    if (hasChildren(element))
                      continue;


                    int *refine_level_elem = stk::mesh::field_data( *refine_level , element );
                    int *refine_field_elem = stk::mesh::field_data( *refine_field , element );
                    if (refine_field_elem[0] < 0)
                      {
                        LocalSetType selected_neighbors;

                        get_node_neighbors(element, selected_neighbors);

                        for (LocalSetType::iterator neighbor = selected_neighbors.begin();
                             neighbor != selected_neighbors.end(); ++neighbor)
                          {
                            if (hasChildren(*neighbor))
                              continue;

                            int *refine_level_neigh = stk::mesh::field_data( *refine_level , *neighbor );
                            int *refine_field_neigh = stk::mesh::field_data( *refine_field , *neighbor );

                            // if any neighbor is more refined (level is higher)
                            //   and I am marked for unrefine, unmark me for unrefine
                            if ( (refine_level_neigh[0] > refine_level_elem[0])
                                 && refine_field_elem[0] < 0)
                              {
                                if (should_check(element, refine_level_elem[0], *neighbor, refine_level_neigh[0],
                                                 bucket_topo_data))
                                  {
                                    refine_field_elem[0] = 0;
                                    if (m_debug_print)
                                      std::cout << "P[" << m_pMesh.get_rank() << "] enforce_two_to_one_unrefine: downgrading element " << m_pMesh.identifier(element)
                                                << " due to neighbor: " << m_pMesh.identifier(*neighbor)
                                                << " with refine_field_neigh= " << refine_field_neigh[0]
                                                << std::endl;
                                    did_change = true;
                                    break;
                                  }
                              }
                          }
                      }
                  }
              }
          }
        {
          std::vector< const stk::mesh::FieldBase *> fields;
          fields.push_back(refine_level);
          fields.push_back(refine_field);
          //stk::mesh::copy_owned_to_shared( *m_pMesh.get_bulk_data(), fields);
          stk::mesh::communicate_field_data(m_pMesh.get_bulk_data()->aura_ghosting(), fields);
        }

        stk::ParallelMachine pm = m_pMesh.get_bulk_data()->parallel();
        stk::all_reduce( pm, stk::ReduceMax<1>( &did_change ) );

        return did_change;
      }

      class SelectIfRefined {
        SubDimCell_SDCEntityType m_subDimEntity;
        vector<NeededEntityType> m_needed_entity_ranks;

      public:
        SelectIfRefined(HangingNodeAdapter& hna) : m_subDimEntity(&hna.m_pMesh) {
          hna.m_breakPattern[0]->fillNeededEntities(m_needed_entity_ranks);
        }
        bool operator()(HangingNodeAdapter& hna, stk::mesh::Entity element) {
          bool isRefined = (hna.m_predicate_refine(element) & DO_REFINE);
          if (hna.m_predicate_refine.m_mark_centroid_always)
            {
              const CellTopologyData * const cell_topo_data = hna.m_pMesh.get_cell_topology(element);
              CellTopology cell_topo(cell_topo_data);


              bool isMarked = false;
              for (unsigned ineed_ent=0; ineed_ent < m_needed_entity_ranks.size(); ineed_ent++)
                {
                  unsigned numSubDimNeededEntities = 0;
                  stk::mesh::EntityRank needed_entity_rank = m_needed_entity_ranks[ineed_ent].first;

                  if (needed_entity_rank == hna.m_pMesh.edge_rank())
                    {
                      numSubDimNeededEntities = cell_topo_data->edge_count;
                    }
                  else if (needed_entity_rank == hna.m_pMesh.face_rank())
                    {
                      numSubDimNeededEntities = cell_topo_data->side_count;
                    }
                  else if (needed_entity_rank == stk::topology::ELEMENT_RANK)
                    {
                      //numSubDimNeededEntities = 1;
                      numSubDimNeededEntities = 0; // all elements are marked... FIXME
                    }

                  for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                    {

                      hna.getNodeRegistry().getSubDimEntity(m_subDimEntity, element, stk::topology::ELEMENT_RANK, iSubDimOrd, cell_topo_data);
                      static SubDimCellData empty_SubDimCellData;
                      SubDimCellData* nodeId_elementOwnderId_ptr = hna.getNodeRegistry().getFromMapPtr(m_subDimEntity);
                      //SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
                      bool is_empty = nodeId_elementOwnderId_ptr == 0;
                      if (!is_empty)
                        {
                          isMarked = true;
                          break;
                        }
                    }
                }
              return isRefined || isMarked;
            }
          return isRefined;
        }
      };

      class SelectIfUnrefined {
      public:
        bool operator()(HangingNodeAdapter& hna, stk::mesh::Entity element) {
          return (hna.m_predicate_refine(element) & DO_UNREFINE);
        }
      };

      class SelectIfInUnrefineTESet {
        ElementUnrefineCollection& m_to_unref;
      public:
        SelectIfInUnrefineTESet(ElementUnrefineCollection& elements_to_unref) : m_to_unref(elements_to_unref) {}
        bool operator()(HangingNodeAdapter& hna, stk::mesh::Entity element) {
          return (m_to_unref.find(element) != m_to_unref.end());
        }
      };

      template<class LocalSelector>
      class HNRefinerSelector : public RefinerSelector {
        HangingNodeAdapter& m_hna;
        LocalSelector& m_localSelector;
        std::unordered_set<stk::mesh::Entity,std::hash<stk::mesh::Entity>> m_elements;
        std::unordered_set<stk::mesh::Entity,std::hash<stk::mesh::Entity>> m_nodes;
        typedef HangingNodeAdapter<RefinePredicate> Base;

      public:
        HNRefinerSelector(HangingNodeAdapter& hna, LocalSelector& ls) : m_hna(hna), m_localSelector(ls)
        {
          init();
        }

        void init()
        {
          m_elements.clear();
          m_nodes.clear();
          PerceptMesh& eMesh = m_hna.m_eMesh;
          const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( eMesh.element_rank() );
          for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_entity_in_bucket = bucket.size();
              for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                {
                  stk::mesh::Entity element = bucket[ientity];
                  //if (m_hna.m_predicate_refine(element) & DO_REFINE)
                  if (m_localSelector(m_hna, element))
                    {
                      unsigned nnode= eMesh.get_bulk_data()->num_nodes(element);
                      stk::mesh::Entity const *elem_nodes = eMesh.get_bulk_data()->begin_nodes(element);
                      for (unsigned ii=0; ii < nnode; ii++)
                        {
                          m_nodes.insert(elem_nodes[ii]);
                        }
                    }
                }
            }
          for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_entity_in_bucket = bucket.size();
              for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                {
                  stk::mesh::Entity element = bucket[ientity];
                  unsigned nnode= eMesh.get_bulk_data()->num_nodes(element);
                  stk::mesh::Entity const *elem_nodes = eMesh.get_bulk_data()->begin_nodes(element);
                  for (unsigned ii=0; ii < nnode; ii++)
                    {
                      if (m_nodes.find(elem_nodes[ii]) != m_nodes.end())
                        {
                          m_elements.insert(element);
                          break;
                        }
                    }
                }
            }
        }
        virtual bool operator()(stk::mesh::Entity element) override {
          if (m_elements.find(element) != m_elements.end())
            return true;
          return false;
        }
      };

    };

  }

#endif
