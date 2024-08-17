// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_Refiner_hpp
#define adapt_Refiner_hpp

#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <cmath>
#include <utility>
#include <math.h>
#include <map>
#include <set>
#include <vector>

#include <stk_util/diag/PrintTimer.hpp>
#include <stk_util/environment/Env.hpp>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

#include <percept/stk_mesh.hpp>
#include <percept/ProgressMeter.hpp>
#include <adapt/UniformRefinerPattern.hpp>

#include <adapt/NodeRegistry.hpp>

#include <adapt/SubDimCell.hpp>

#include <adapt/RefinementInfoByType.hpp>

#if defined( STK_PERCEPT_HAS_GEOMETRY )
#include <percept/mesh/geometry/kernel/MeshGeometry.hpp>
#endif

#include <percept/PerceptMesh.hpp>

#define DO_REF_LTRACE 0
#define REF_LTRACE_PROC 1
#define REF_LTRACE_TEST (m_eMesh.get_rank() < REF_LTRACE_PROC)
#if DO_REF_LTRACE
#define REF_LTRACE(msg) do { if (REF_LTRACE_TEST) std::cout << msg << " " << __FILE__ << ":" << __LINE__ << std::endl; } while(0)
#define REF_LTRACE_0() REF_LTRACE("")
#else
#define REF_LTRACE(msg) do { } while(0)
#define REF_LTRACE_0() REF_LTRACE("")
#endif

  namespace percept {

    // free functions used in NodeRegistry, FixSideSets, etc.
    void mod_begin_timer(stk::mesh::BulkData& bulk_data, stk::diag::Timer& parent_timer);
    void mod_end_timer(stk::mesh::BulkData& bulk_data, stk::diag::Timer& parent_timer, const std::string& msg);

    typedef std::map<stk::mesh::Part*, stk::mesh::Part*> SideElementPartMap;
    typedef std::map<stk::mesh::Part*, stk::mesh::PartVector> SidePartMap;

    using std::vector;
    using std::map;
    using std::set;

    typedef percept::SetOfEntities SetOfEntities;
    typedef percept::SetOfEntities elements_to_be_destroyed_type;
    typedef elements_to_be_destroyed_type ElementUnrefineCollection;

    class RefinerSelector {
    public:
      // if participates in adaptivity, return true - e.g., for elements
      // far away from the "action", this can return false - be aware
      // that even though an element may not be refined, if it's possible
      // that its neighbor is refined, this must return true
      virtual bool operator()(stk::mesh::Entity element) { throw std::runtime_error("not impl"); }
      virtual bool use_batch_filter() { return false; }
      virtual void batch_filter(stk::mesh::EntityRank rank, std::vector<stk::mesh::Entity>& elements) { throw std::runtime_error("not impl"); }
      virtual void batch_filter(stk::mesh::EntityRank rank, SetOfEntities& elements) { throw std::runtime_error("not impl"); }
    };

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    //template<class UniformRefinerPattern>
    class Refiner
#ifndef SWIG //NLM
  : public percept::Observable<ProgressMeterData>
#endif
    {
    public:
      Refiner(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0);

      virtual ~Refiner();

      void doProgressPrint(const std::string& str);
      virtual void
      doBreak(int num_registration_loops=1);

      void doRebalance();

      void
      setRemoveOldElements(bool do_remove);
      bool
      getRemoveOldElements();

      void setAddChildrenToParts(bool doit) { m_doAddChildrenToParts = doit; }
      bool getAddChildrenToParts() { return m_doAddChildrenToParts; }

      // if set, this will avoid refining this part vecto
      void setExcludeParts(const stk::mesh::PartVector& parts) { m_excludeParts = parts; }
      const stk::mesh::PartVector&  getExcludeParts() { return m_excludeParts; }

      void setAvoidFixSideSets(bool val) { m_avoidFixSideSets = val; }
      bool getAvoidFixSideSets() { return m_avoidFixSideSets; }

      void setAvoidFixSideSetChecks(bool val) { m_avoidFixSideSetChecks = val; }
      bool getAvoidFixSideSetChecks() { return m_avoidFixSideSetChecks; }

      void setAvoidClearDanglingNodes(bool val) { m_avoidClearDanglingNodes = val; }
      bool getAvoidClearDanglingNodes() { return m_avoidClearDanglingNodes; }

      /* for future
      void
      setIOSaveInactiveElements(bool do_save) { m_doIOSaveInactiveElements = do_save; }
      bool
      getIOSaveInactiveElements() { return m_doIOSaveInactiveElements; }
      */

      void
      setGeometryFile(std::string file_name);

      void
      setSmoothGeometry(bool do_smooth) { m_doSmoothGeometry = do_smooth; }
      bool
      getSmoothGeometry() { return m_doSmoothGeometry; }

      void
      setRemoveGeometryBlocks(bool do_remove) { m_removeGeometryBlocks = do_remove; }
      bool
      getRemoveGeometryBlocks() { return m_removeGeometryBlocks; }

      void
      setIgnoreSideSets(bool ignore_sidesets) ;

      bool
      getIgnoreSideSets();

      void
      setDoRebalance(bool doRebal, double threshold=1.0) { m_doRebalance = doRebal; m_rebalThreshold = threshold; }
      bool
      getDoRebalance(double *threshold = 0) { if (threshold) *threshold = m_rebalThreshold; return m_doRebalance; }

      RefinementInfo&
      getRefinementInfo();

      void
      setDoProgressMeter(bool do_progress);
      bool
      getDoProgressMeter();

      void setFixAllBlockBoundaries(bool val) { m_fixAllBlockBoundaries=val; }
      bool getFixAllBlockBoundaries() { return m_fixAllBlockBoundaries; }

      void setRefinerSelector(RefinerSelector *sel) { m_refinerSelector = sel; }
      RefinerSelector *getRefinerSelector() { return m_refinerSelector; }

#if defined(STK_BUILT_FOR_SIERRA)
      void set_rbar_special_treatment(BlockNamesType& rbar_names) { m_rbar_names = rbar_names; }
#endif

      // ================================ unrefine

      typedef SetOfEntities NodeSetType;

      virtual void
      unrefineTheseElements(ElementUnrefineCollection& elements_to_unref);

      void
      unrefineAll();

      void
      setAlwaysInitializeNodeRegistry(bool do_init) { m_alwaysInitNodeRegistry = do_init; }

      bool
      getAlwaysInitializeNodeRegistry() { return m_alwaysInitNodeRegistry; }

      stk::diag::Timer *getAlternateRootTimer() { return m_alternateRootTimer; }
      void setAlternateRootTimer(stk::diag::Timer *timer) { m_alternateRootTimer = timer; }

      stk::diag::Timer *getModBegEndRootTimer() { return m_modBegEndRootTimer; }
      void setModBegEndRootTimer(stk::diag::Timer *timer) { m_modBegEndRootTimer = timer; }

      void reset_family_tree_to_node_relations();

      /// for quad/hex hangin node topology, we don't need to remesh during unrefinement
      ///   default is on (for tri/tet local refinement)
      void
      setNeedsRemesh(bool needsRemesh) { m_needsRemesh = needsRemesh; }

      bool
      getNeedsRemesh() { return m_needsRemesh; }

      /// for tri/tet local refinement, set this to true for much better
      /// unrefinement behavior (elements are unrefined in order they were
      /// refined which leads to less locking behavior)
      void
      setDoLevelBasedUnrefinement(bool doLevelBasedUnrefinement) { m_doLevelBasedUnrefinement = doLevelBasedUnrefinement; }

      bool
      getDoLevelBasedUnrefinement() { return m_doLevelBasedUnrefinement; }

      void mod_begin(stk::diag::Timer *timer=0);
      void mod_end(stk::diag::Timer *timer=0, const std::string& msg="");


#if  defined(STK_PERCEPT_HAS_GEOMETRY)

      enum SMOOTHING_OPTIONS {
        // snaps and discards original coord field, tries to smooth
        SNAP_PLUS_SMOOTH,
        // keeps original and snapped states; does line search between; tries to keep always-valid mesh
        USE_LINE_SEARCH_WITH_MULTIPLE_STATES
        //,END_OPTIONS
      };

      void snapAndSmooth(bool geomSnap, std::string geomFile, bool use_ref_mesh=true);

      void smoothGeometry(MeshGeometry* mesh_geometry, stk::mesh::Selector *selector, SMOOTHING_OPTIONS option, bool use_ref_mesh=true);
#endif

      void deleteParentElements();

      void fix_side_sets_2(bool allow_not_found=false, SetOfEntities *avoid_elems = 0, SetOfEntities *avoid_sides = 0, RefinerSelector *sel=0, const std::string& msg="");

      /// determine side part to elem part relations
      static
      void get_side_part_relations(PerceptMesh& eMesh, bool checkParentChild, SidePartMap& side_part_map, bool debug = false);

      bool connect(stk::mesh::Entity side, bool& valid_side_part_map, SetOfEntities *avoid=0);
      bool connectSidesForced(stk::mesh::Entity element, stk::mesh::Entity side_elem, bool& valid_side_part_map, bool use_coordinate_compare=false);

      NodeRegistry& getNodeRegistry() { return *m_nodeRegistry; }
      percept::PerceptMesh& getMesh() { return m_eMesh; }

      //============= unrefine

    public:

      void
      replaceNodeRegistryOwnership(ElementUnrefineCollection& elements_to_delete, stk::mesh::EntityRank rank);

      void
      initializeDB(bool use_rebuild_node_registry=true);

      void
      unrefinePass2(ElementUnrefineCollection& elements_to_unref);

      stk::diag::Timer &rootTimer();

      void
      set_active_part();

      static void
      set_active_part(PerceptMesh& eMesh);

      virtual void
      initializeRefine();

      void check_parent_ownership();
      void require_sides_on_same_proc_as_pos_perm_element();
      bool check_sides_on_same_proc_as_owned_element(const std::string& msg, bool doThrow = true);
      void build_side_set(SetOfEntities& side_set, bool only_roots = false);
      bool bucket_acceptable(stk::mesh::Bucket& bucket, stk::mesh::EntityRank rank);
      bool include_side_bucket(stk::mesh::Bucket& side_bucket, stk::mesh::Selector *excludeSelector);

    protected:
      void collectElemsToRefine(const unsigned irank, stk::mesh::EntityRank rank, const unsigned elementType,
                  std::vector<stk::mesh::Entity>& elems, int& jele);

      void fillElementRankTypeInfo(std::vector<stk::mesh::EntityRank>& ranks);

      void getRefinementInfo(std::vector<stk::mesh::EntityRank>& ranks);

      void filterUsingRefinerSelector(stk::mesh::EntityRank rank, std::vector<stk::mesh::Entity>& elements);

      typedef  std::pair<stk::mesh::EntityRank, unsigned > ElementRankTypeInfo;


      // these are called by doBreak
      virtual void
      doMark(int num_registration_loops=1);
      virtual void
      doRefine();
      virtual void
      finalizeRefine();

      // called by doMark
      virtual void
      preMark(int iter, int num_registration_loops);
      //virtual void
      //mark(int iter, int num_registration_loops);

      // return whether to break out of the registration loop
      virtual bool
      postMark(int iter, int num_registration_loops);


      // support methods
      void
      remeshRecurse(stk::mesh::Entity element, int& s_depth);

      void
      remeshElements(SetOfEntities& rootElements, stk::mesh::EntityRank rank, int pool_size_hint=0, SetOfEntities* elemsToBeDeleted=0);

      bool allDescendants(stk::mesh::Entity element, SetOfEntities& descendants, unsigned& nlevels, bool only_leaves=false,
                          ElementUnrefineCollection *elements_to_unref=0);

      void
      filterUnrefSetPass2(ElementUnrefineCollection& elements_to_unref,   SetOfEntities& grandParents);
      void
      filterRecurse(stk::mesh::Entity element, ElementUnrefineCollection& rootElements, ElementUnrefineCollection& elements_to_unref);


      void
      removeDeletedNodes(NodeSetType& deleted_nodes);

      void removeFamilyTrees(SetOfEntities& family_trees_to_be_removed);


      void getSideParentsToBeRemeshed(SetOfEntities& children_to_be_removed, SetOfEntities& parent_side_elements, bool newVersion = false, SetOfEntities *avoid_sides = 0);

      void removeChildElements(SetOfEntities& children_to_be_removed, ElementUnrefineCollection* elements_to_unref_0=0);

      void get_kept_nodes(SetOfEntities& kept_nodes, ElementUnrefineCollection& elements_to_unref);
      void get_deleted_nodes(SetOfEntities& deleted_nodes, SetOfEntities& kept_nodes, ElementUnrefineCollection& elements_to_unref, SetOfEntities& elements_to_be_remeshed);
      void filter_deleted_nodes(SetOfEntities& deleted_nodes);

      void delete_entities(ElementUnrefineCollection& elements_to_unref, stk::mesh::EntityRank rank);
      void generate_temporary_elements(SetOfEntities& entities_to_delete, stk::mesh::EntityRank rank, std::vector<stk::mesh::Entity>& new_elements);

      void get_deleted_sides(SetOfEntities& sides_to_delete, ElementUnrefineCollection& elements_to_unref, SetOfEntities& elements_to_be_remeshed);

      void remesh(stk::mesh::Entity parent_element);
      //============= unrefine end

      /**  Overrides start =======>
       */

      unsigned
	  countAndGatherAllElements(unsigned irank, stk::mesh::EntityRank rank, unsigned elementType, std::vector<stk::mesh::Entity> &elements);

      /** Overrides
       *   m_nodeRegistry data member.  The same loop should be executed every time this method is called (i.e. the same elements should be visited).  Also, there
       *   is policy associated with the @param only_count and @param doAllElements inputs.  If doAllElements==true, then ghost elements should not be skipped.
       *   If only_count==true, then *do not* call the supplied @param function, rather just count the number of elements to be visited and return that number.
       *   @return Always return the number of elements visited (modulo the doAllElements parameter).
       *
       *   Note: client code can just use the supplied implementation in this base class - it is not necessary to override
       *
       *   @see UniformRefiner implementation of this method for an example.
       */

      virtual unsigned
      doForAllElements(unsigned irank, std::string function_info,
                       stk::mesh::EntityRank rank, NodeRegistry::ElementFunctionPrototype function,
                       unsigned elementType,
                       vector<NeededEntityType>& needed_entity_ranks,
                       bool doAllElements=true) ;

      /** Create a list of nodes from the new nodes that can be easily deciphered by the UniformRefinerPattern.
       *
       *  This is a helper function that gets all the nodes ready for use in creating new elements.  For each rank of those needed by
       *  the algorithm (supplied information from the break pattern), it fills the 3D array new_sub_entity_nodes with the new nodes.
       *
       *  Returns the 3D array new_sub_entity_nodes[entity_rank_of_entity_or_sub_dim_entity][ordinal_of_sub_dim_entity][ordinal_of_node_on_sub_dim_entity]
       *
       */
      virtual bool
      createNewNeededNodeIds(const CellTopologyData * const cell_topo_data,
                             const stk::mesh::Entity element, vector<NeededEntityType>& needed_entity_ranks, NewSubEntityNodesType& nodes, UniformRefinerPatternBase *breakPattern) ;

      /** Method that actually creates new elements by first calling createNewNeededNodeIds then calls the break pattern's createNewElements method.
       *
       *  A sample implementation is shown in @see UniformRefiner
       */
      virtual size_t
      createElementsAndNodesAndConnectLocal(unsigned irank,  stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern,
                                            unsigned elementType,
                                            vector<NeededEntityType>& needed_entity_ranks,
                                            vector<stk::mesh::Entity>& new_elements_pool,
                                            vector<stk::mesh::Entity>& ft_new_elements_pool,
                                            vector<stk::mesh::Entity>::iterator * new_elements_pool_end_iter,
                                            vector<stk::mesh::Entity>::iterator * ft_new_elements_pool_end_iter
                                            );

      /** This is a helper method that loops over all sub-dimensional entities whose rank matches on of those in @param needed_entity_ranks
       *    and registers that sub-dimensional entity as needing a new node, or whatever other function NodeRegistry requires (getFromRemote(), etc)
       *  Override it to only apply the @param function to the desired sub-entities (e.g. for non-uniform/local refinement)
       *
       *  It is a copy of NodeRegistry's doForAllSubEntities method.  Provided here so it can be overridden.
       *
       *  Note: this is the minimal function that needs to be overridden to get different marking/refining behavior
       */

      virtual void
      refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity element,
                        vector<NeededEntityType>& needed_entity_ranks, const CellTopologyData * const bucket_topo_data);

      /// allow for insertion of calls throughout derived classes algorithms that can callback to other derived classes overloads
      ///  For example, for wedge special processing we have special patterns, special marks (using modify_marks above),
      ///    and special refine_field setting for wedges in boundary layers - see TEA_SpecialWedgeRefinement class
      /// @param step is some arbitrary step info that can be used in a switch statement;
      /// @param data is arbitrary data needed by specializations of this method
      virtual void special_processing(const std::string& step, void *data = 0)
      {
      }

      /// =========>  Overrides  end

      void update_node_registry();

      void
      removeFamilyTrees();

      void
      removeOldElements(unsigned irank, stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern );

      void
      removeElements( elements_to_be_destroyed_type& elements_to_be_destroyed, unsigned irank=0);

      void
      removeEmptyElements();

      void
      removeEmptyFamilyTrees();

      void
      removeFromNewNodesPart();

      // empty nodes (nodes not referred to by any elements) are possibly created during refine, this method removes them
    public:
      void
      removeDanglingNodes();

    protected:

      void
      addOldElementsToPart(stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern, unsigned *elementType = 0u);

      void
      renameNewParts(stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern);

      void
      fixSurfaceAndEdgeSetNames(stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern);

      void
      buildElementSideDB(SubDimCellToDataMap& cell_2_data_map);

      void
      checkBreakPatternValidityAndBuildRanks(std::vector<stk::mesh::EntityRank>& ranks, stk::diag::Timer *timer);

      void
      add_children_to_parts();

      std::string get_parent_element_topology(const std::string& surfaceName);

    protected:
      percept::PerceptMesh& m_eMesh;

      std::vector<UniformRefinerPatternBase *> m_breakPattern;

      NodeRegistry* m_nodeRegistry;
      stk::mesh::FieldBase *m_proc_rank_field;
      bool m_doRemove;

      // for future:
      // bool m_doIOSaveInactiveElements; // default false

      std::vector<stk::mesh::EntityRank> m_ranks;
      bool m_ignoreSideSets;
      std::string m_geomFile;
      bool m_geomSnap;

      RefinementInfo m_refinementInfo;

      int m_progress_meter_frequency;
      bool m_doProgress;

      bool m_alwaysInitNodeRegistry;
      bool m_doSmoothGeometry;
      bool m_allocated;

      bool m_removeGeometryBlocks;
      SidePartMap m_side_part_map;
      bool m_fixAllBlockBoundaries;

#if defined(STK_BUILT_FOR_SIERRA)
      BlockNamesType m_rbar_names;
#endif
      bool m_needsRemesh;
      bool m_doLevelBasedUnrefinement;

      vector< ElementRankTypeInfo > m_elementRankTypeInfo;

      stk::diag::Timer *m_alternateRootTimer;
      stk::diag::Timer *m_modBegEndRootTimer;
      RefinerSelector *m_refinerSelector;
      RefinerSelector *m_fixSideSetsSelector;
      bool m_doAddChildrenToParts;
      stk::mesh::PartVector m_excludeParts;
      bool m_avoidFixSideSets;
      bool m_avoidFixSideSetChecks;
      bool m_avoidClearDanglingNodes;
      std::vector<stk::mesh::Entity> m_element_ft_pool;
      bool m_onlyOneLevelUnrefine;
      bool m_doRebalance;
      double m_rebalThreshold;
      bool m_removeFromNewNodesPart;
      bool m_do_new_elements;
      stk::diag::TimerSet m_timerSet;
      stk::diag::Timer m_timer;
      std::vector<stk::mesh::Selector> m_fromPartsSelector;
    };



  }

#endif
