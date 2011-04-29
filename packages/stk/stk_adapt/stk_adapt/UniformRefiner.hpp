#ifndef stk_adapt_UniformRefiner_hpp
#define stk_adapt_UniformRefiner_hpp

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

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

#include <boost/tuple/tuple_io.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include <stk_percept/stk_mesh.hpp>
#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/ProgressMeter.hpp>
#include <stk_adapt/UniformRefinerPattern.hpp>
#include <stk_adapt/Colorer.hpp>

#include <stk_adapt/NodeRegistry.hpp>

#include <stk_adapt/SubDimCell.hpp>

#include <stk_adapt/RefinementInfoByType.hpp>

#define UNIFORM_REF_REMOVE_OLD_STD_SET 1
#define UNIFORM_REF_REMOVE_OLD_STD_VECTOR 0
#define UNIFORM_REF_REMOVE_OLD_BOOST_SET 0

#if UNIFORM_REF_REMOVE_OLD_BOOST_SET
#include <boost/unordered_set.hpp>
#endif


namespace stk {
  namespace adapt {

    using std::vector;
    using std::map;
    using std::set;


#if UNIFORM_REF_REMOVE_OLD_STD_SET
    typedef std::set<stk::mesh::Entity *> elements_to_be_destroyed_type;
#endif
#if UNIFORM_REF_REMOVE_OLD_STD_VECTOR
    typedef std::vector<stk::mesh::Entity *> elements_to_be_destroyed_type;
#endif
#if UNIFORM_REF_REMOVE_OLD_BOOST_SET
    typedef boost::unordered_set<stk::mesh::Entity *> elements_to_be_destroyed_type;
#endif


    /// e.g. UniformRefiner<shards::Hex<8>, shards::Tet<4> >
    //template<typename FromTopology, typename ToTopology>
#if 0
    class ParallelMeshModAlgorithm
    {
    public:
      virtual void planActions()=0;
      virtual void performActions()=0;
    protected:
      void helperFunction1();
      void helperFunction2();
      //...
    };
#endif

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    //template<class UniformRefinerPattern>
  class UniformRefiner : public stk::percept::Observable<ProgressMeterData> 
    {
    public:
      UniformRefiner(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0);
      UniformRefiner(percept::PerceptMesh& eMesh, std::vector<UniformRefinerPatternBase *>&  bp, stk::mesh::FieldBase *proc_rank_field=0);


      //UniformRefiner(percept::PerceptMesh& eMesh, UniformRefinerPattern<void, void, 0>& bp);
  
      void 
      doBreak();

      void 
      setRemoveOldElements(bool do_remove);

      void
      setGeometryFile(std::string file_name);

      bool 
      getRemoveOldElements();
      
      static BlockNamesType 
      getBlockNames(std::string& block_name, unsigned proc_rank, percept::PerceptMesh& eMesh);

      static BlockNamesType 
      correctBlockNamesForPartPartConsistency(percept::PerceptMesh& eMesh, BlockNamesType& blocks);

      void 
      setIgnoreSideSets(bool ignore_sidesets) ;

      bool 
      getIgnoreSideSets();

      std::vector< RefinementInfoByType >& 
      getRefinementInfoByType();

      void 
      setQueryPassOnly(bool doQueryOnly);

    protected:
  
      //void checkParallelConsitency();

      unsigned
      doForAllElements(stk::mesh::EntityRank rank, NodeRegistry::ElementFunctionPrototype function, vector< ColorerSetType >& elementColors, vector<NeededEntityType>& needed_entity_ranks,
                       bool only_count=false, bool doAllElements=true);

      void 
      createElementsAndNodesAndConnectLocal(unsigned irank,  UniformRefinerPatternBase* breakPattern, 
                   vector< ColorerSetType >& elementColors,   vector<NeededEntityType>& needed_entity_ranks,  vector<stk::mesh::Entity *>& new_elements_pool);

      bool
      createNewNeededNodeIds(const CellTopologyData * const cell_topo_data, 
                             const stk::mesh::Entity& element, vector<NeededEntityType>& needed_entity_ranks, NewSubEntityNodesType& nodes);

      void 
      removeOldElements(stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern );

      void 
      removeOldElements( elements_to_be_destroyed_type& elements_to_be_destroyed);

      void 
      addOldElementsToPart(stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern, unsigned *elementType = 0u);

      void 
      renameNewParts(stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern);

      void 
      fixSurfaceAndEdgeSetNames(stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern);

      void 
      fixElementSides();

      void 
      fixElementSides1();

      void 
      fixElementSides(stk::mesh::EntityRank side_rank);

      void 
      fixElementSides1(stk::mesh::EntityRank side_rank);

      void 
      checkFixElementSides(stk::mesh::EntityRank side_rank, stk::mesh::EntityRank elem_rank);

      void 
      buildElementSideDB(SubDimCellToDataMap& cell_2_data_map);

      void 
      trace_print(std::string msg);

      void 
      checkBreakPatternValidityAndBuildRanks(std::vector<stk::mesh::EntityRank>& ranks);



    private:
      percept::PerceptMesh& m_eMesh;

      //UniformRefinerPatternBase & m_breakPattern;
      std::vector<UniformRefinerPatternBase *> m_breakPattern;

      NodeRegistry* m_nodeRegistry;
      stk::mesh::FieldBase *m_proc_rank_field;
      bool m_doRemove;

      std::vector<stk::mesh::EntityRank> m_ranks;
      bool m_ignoreSideSets;
      std::string m_geomFile;
      bool m_geomSnap;

      std::vector< RefinementInfoByType > m_refinementInfoByType;
      bool m_doQueryOnly;

      int m_progress_meter_frequency;
      bool m_doProgress;
    };



  }
}
#endif
