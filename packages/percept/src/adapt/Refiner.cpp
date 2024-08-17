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
#include <sys/stat.h>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <adapt/Refiner.hpp>
#include <adapt/RefinerUtil.hpp>
#include <adapt/FixSideSets.hpp>
#include <percept/PerceptMesh.hpp>
#include <percept/PerceptUtils.hpp>

#include <percept/MeshUtil.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/util/human_bytes.hpp>

#include <adapt/NodeRegistryDef.hpp>
#include <percept/RebalanceMesh.hpp>

#if defined( STK_PERCEPT_HAS_GEOMETRY )
#include <percept/mesh/geometry/kernel/GeometryKernelGregoryPatch.hpp>
#if HAVE_OPENNURBS
#include <percept/mesh/geometry/kernel/GeometryKernelOpenNURBS.hpp>
#endif
#if HAVE_CUBIT
#include <percept/mesh/geometry/kernel/GeometryKernelPGEOM.hpp>
#endif
#include <percept/mesh/geometry/kernel/MeshGeometry.hpp>
#include <percept/mesh/geometry/kernel/GeometryFactory.hpp>

#include <percept/mesh/mod/smoother/MeshSmoother.hpp>
#include <percept/mesh/mod/smoother/ReferenceMeshSmootherAlgebraic.hpp>
#endif

#include <adapt/FixSideSetsSelector.hpp>
#include <adapt/AdaptedMeshVerifier.hpp>

#include <percept/mesh/geometry/volume/VolumeUtil.hpp>

#define TIMER2(name,parentName) \
  stk::diag::Timer       timer ## name (        #name, timer ## parentName);  \
  stk::diag::TimeBlock tbTimer ## name (timer ## name)

#define USE_PERCEPT_PERFORMANCE_TESTING_CALLGRIND 0
#if USE_PERCEPT_PERFORMANCE_TESTING_CALLGRIND
#include "/usr/netpub/valgrind-3.8.1/include/valgrind/callgrind.h"
#endif

#include <stk_mesh/base/MeshUtils.hpp>

  namespace percept {
    extern bool s_do_transition_break;

    void mod_begin_timer(stk::mesh::BulkData& bulk_data, stk::diag::Timer& parent_timer)
    {
      stk::diag::Timer timerModBeg_("ModBeg", parent_timer);
      stk::diag::TimeBlock timerModBegBlock_(timerModBeg_);
      
      bulk_data.modification_begin();
    }
    
    void mod_end_timer(stk::mesh::BulkData& bulk_data, stk::diag::Timer& parent_timer, const std::string& msg)
    {
      stk::diag::Timer timerModEndAll_("ModEndAll",    parent_timer);
      stk::diag::Timer timerModEnd_(   "ModEnd: "+msg, timerModEndAll_);
      stk::diag::TimeBlock timerModEndBlock_(timerModEnd_);
      stk::diag::TimeBlock timerModEndBlockAll_(timerModEndAll_);
      
      bulk_data.modification_end();
    }

    using namespace percept;

    NodeRegistry *s_nodeRegistry = 0;

#define CHK(eMesh) do { } while (0)

#if 0
      static void check(NodeRegistry* nr, const std::string& where)
      {

        PerceptMesh& eMesh = nr->getMesh();
        int id = 2;
        stk::mesh::Entity entity = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), id);
        std::ostringstream str;
        //str << where;
        if(eMesh.is_valid(entity))
          {
            str << "P[" << eMesh.get_rank() << "] CHK: " << where << " elem= " << eMesh.print_entity_compact(entity);
            str << "\n";
            nr->query(str, id, (int)eMesh.side_rank(), 4);
          }
        else
          {
            str << " not valid";
          }
        std::cout << str.str() << std::endl;
      }
#endif


    stk::diag::Timer& Refiner::rootTimer()
    {
      if (m_alternateRootTimer)
        {
          return *m_alternateRootTimer;
        }
      else
        {
          return m_timer;
        }
    }

    Refiner::Refiner(percept::PerceptMesh& eMesh, UniformRefinerPatternBase &  bp, stk::mesh::FieldBase *proc_rank_field) :
      m_eMesh(eMesh), m_breakPattern(),
      m_nodeRegistry(0),
      m_proc_rank_field(proc_rank_field), m_doRemove(true), m_ranks(), m_ignoreSideSets(false),
      m_geomFile(""), m_geomSnap(false),
      m_refinementInfo(this),
      m_progress_meter_frequency(20),
      m_doProgress(false),
      m_alwaysInitNodeRegistry(true),
      m_doSmoothGeometry(false),
      m_removeGeometryBlocks(false),
      m_fixAllBlockBoundaries(false),
      m_needsRemesh(true),
      m_doLevelBasedUnrefinement(false)
      ,m_alternateRootTimer(0)
      ,m_modBegEndRootTimer(0)
      ,m_refinerSelector(0)
      ,m_doAddChildrenToParts(true)
      ,m_avoidFixSideSets(false)
      ,m_avoidFixSideSetChecks(false)
      ,m_avoidClearDanglingNodes(false)
      ,m_onlyOneLevelUnrefine(false)
      ,m_doRebalance(false)
      ,m_rebalThreshold(1.0)
      ,m_removeFromNewNodesPart(true)
      ,m_do_new_elements(true)
      ,m_timerSet(sierra::Diag::TIMER_ALL)
      ,m_timer(stk::diag::createRootTimer("Refiner", m_timerSet))

    {
      bp.setSubPatterns(m_breakPattern, eMesh);
      m_nodeRegistry = new NodeRegistry (m_eMesh, this);
      s_nodeRegistry = m_nodeRegistry;
      m_nodeRegistry->initialize();
      m_nodeRegistry->init_comm_all();
      m_allocated = false;
      m_timer.start();
    }

    Refiner::~Refiner()
    {
      if (m_nodeRegistry)
        delete m_nodeRegistry;
      if (m_allocated)
      {
        //for (unsigned int i = 0; i < m_breakPattern.size(); i++)
        // SRK only delete the first one - others are handled by the parent 0'th break pattern
        for (unsigned int i = 0; i < 1; i++)
        {
          delete m_breakPattern[i];
        }
      }
      stk::diag::deleteRootTimer(m_timer);
    }

    void Refiner::doProgressPrint(const std::string& msg)
    {
      if (m_doProgress && m_eMesh.get_rank() == 0)
        {
          size_t now=0, hwm=0;
          stk::get_memory_usage(now, hwm);
          std::cout << std::left << std::setw(50) << msg.c_str();
          std::cout << " cpu: " << m_eMesh.cpu_time() << " [sec] mem= " << stk::human_bytes(now) << " [now_proc0] " << stk::human_bytes(hwm) << " [hwm_proc0]" << std::endl;
        }
    }

#define EXTRA_PRINT_UR_GETBLOCKS 0

    void Refiner::
    setGeometryFile(std::string file_name) { m_geomFile = file_name;
      m_geomSnap = true; }


    void Refiner::
    setRemoveOldElements(bool do_remove) { m_doRemove = do_remove; }

    bool Refiner::
    getRemoveOldElements() { return m_doRemove; }

    void Refiner::
    setDoProgressMeter(bool do_progress) { m_doProgress = do_progress; }

    bool Refiner::
    getDoProgressMeter() { return m_doProgress; }

    void Refiner::
    setIgnoreSideSets(bool ignore_ss)
    {
      m_ignoreSideSets= ignore_ss;
    }

    bool Refiner::
    getIgnoreSideSets()
    {
      return m_ignoreSideSets;
    }

    void Refiner::
    addOldElementsToPart(stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern, unsigned *elementType)
    {
      EXCEPTWATCH;
      std::string oldPartName = breakPattern->getOldElementsPartName()+toString(rank);
      stk::mesh::Part *oldPart = m_eMesh.get_fem_meta_data()->get_part(oldPartName);

      if (!oldPart)
        {
          std::cout << "oldPartName= " << oldPartName << std::endl;
          throw std::runtime_error("oldpart is null");
        }

      stk::mesh::PartVector add_parts(1, oldPart);
      stk::mesh::PartVector remove_parts;
      stk::mesh::Selector on_locally_owned_part =  ( m_eMesh.get_fem_meta_data()->locally_owned_part() );

      // The list of Parts that this break pattern will refine.  Only remove elements belonging to these parts.
      stk::mesh::Selector fromPartsSelector = stk::mesh::selectUnion( breakPattern->getFromParts() );

      std::vector<stk::mesh::Entity> elems;
      const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( rank );

      //unsigned nele=0;
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (on_locally_owned_part(**k) && fromPartsSelector(**k) )
            {
              stk::mesh::Bucket & bucket = **k ;
              const CellTopologyData * const bucket_cell_topo_data = m_eMesh.get_cell_topology(bucket);
              shards::CellTopology topo(bucket_cell_topo_data);

              const unsigned num_elements_in_bucket = bucket.size();

              for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
                {
                  EXCEPTWATCH;
                  stk::mesh::Entity element = bucket[i_element];
                  if (!m_eMesh.is_valid(element))
                    {
                      std::cout << "element = 0" << std::endl;
                      throw std::runtime_error("element = 0");
                    }

                  if (elementType && (topo.getKey() != *elementType))
                    {
                    }
                  else
                    {
                      elems.push_back(element);
                      //++nele;
                    }
                }
            }
        }


      for (unsigned ielem=0; ielem < elems.size(); ielem++)
        {
          m_eMesh.get_bulk_data()->change_entity_parts( elems[ielem], add_parts, remove_parts );
        }

    }

    struct myVec
    {
      double *data;
      int len;
      int res;
    };

    void Refiner::
    checkBreakPatternValidityAndBuildRanks(std::vector<stk::mesh::EntityRank>& ranks, stk::diag::Timer *timer)
    {
      ranks.resize(0);

      if (m_doRemove)
        mod_begin(timer);

      for (unsigned ibp = 0; ibp < m_breakPattern.size(); ibp++)
        {
          if (m_breakPattern[ibp])
            {
              stk::mesh::EntityRank irank = m_breakPattern[ibp]->getPrimaryEntityRank();
              stk::mesh::EntityRank irank_prev = static_cast<stk::mesh::EntityRank>(percept::EntityRankEnd);

              if (ibp > 0) irank_prev = m_breakPattern[ibp-1]->getPrimaryEntityRank();
              if (irank > irank_prev)
                {
                  for (unsigned jbp = 0; jbp < m_breakPattern.size(); jbp++)
                    {
                      if (m_breakPattern[jbp])
                        {
                          stk::mesh::EntityRank jrank = m_breakPattern[jbp]->getPrimaryEntityRank();
                          std::cout << "tmp jbp= " << jbp << " jrank= " << jrank
                                    << " fromTopo= " << m_breakPattern[jbp]->getFromTopology()->name
                                    << " toTopo= " << m_breakPattern[jbp]->getToTopology()->name
                                    << std::endl;
                        }
                    }

                  std::cout << "tmp irank= " << irank << " irank_prev= " << irank_prev << std::endl;
                  throw std::logic_error("m_breakPattern: must be in decreasing order of rank");
                }
              ranks.push_back(irank);
              if (m_doRemove)
                {
                  unsigned elementType = m_breakPattern[ibp]->getFromTypeKey();
                  addOldElementsToPart(irank, m_breakPattern[ibp], &elementType);
                }
            }
          else
            {
              std::cout << "ibp = " << ibp << std::endl;
              throw std::logic_error("m_breakPattern is null");
            }
        }
      if (m_doRemove)
        mod_end_timer(*m_eMesh.get_bulk_data(), *timer, "RefinerCBR");

    }

    void Refiner::
    refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity element, std::vector<NeededEntityType>& needed_entity_ranks, const CellTopologyData * const bucket_topo_data)
    {
      const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);

      shards::CellTopology cell_topo(cell_topo_data);

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
              /// note: at this level of granularity we can do single edge refinement, hanging nodes, etc.
              bool doMark = true;
              if (needed_entity_ranks[ineed_ent].third.size())
                {
                  VERIFY_OP_ON(needed_entity_ranks[ineed_ent].third.size(), ==, numSubDimNeededEntities, "bad size");
                  if (!needed_entity_ranks[ineed_ent].third[iSubDimOrd])
                    doMark = false;
                }
              if (doMark)
                (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, true, bucket_topo_data);

            } // iSubDimOrd
        } // ineed_ent
    }

    // called by doMark
    void Refiner::
    preMark(int iter, int num_registration_loops) {}

    bool Refiner::
    postMark(int iter, int num_registration_loops) { return false; }

    void Refiner::
    doBreak(int num_registration_loops)
    {
      EXCEPTWATCH;

      m_eMesh.get_bulk_data()->delete_face_adjacent_element_graph();

      initializeRefine();

      doMark(num_registration_loops);

      doRefine();
    }

    void Refiner::
    initializeRefine()
    {
      stk::diag::Timer timerAdapt_("RefineMesh", rootTimer());
      stk::diag::TimeBlock timerAdaptBlock_(timerAdapt_);

      stk::diag::Timer timerInitRefine_("percept::InitRefine", timerAdapt_);
      stk::diag::TimeBlock timerInitRefineBlock_(timerInitRefine_);

      REF_LTRACE("initializeRefine: start");

      {
        TIMER2(IR_dumpDB,InitRefine_);
        m_nodeRegistry->dumpDB("start of doBreak");
      }

      {
        TIMER2(IR_req_same_proc,InitRefine_);
        if (m_eMesh.getProperty("percept_Refiner_require_sides") != "false")
          require_sides_on_same_proc_as_pos_perm_element();
      }

      {
        TIMER2(IR_side_part,InitRefine_);

        if (m_eMesh.getProperty("use_side_map") == "true")
          get_side_part_relations(m_eMesh, false, m_side_part_map);
      }

      {
        TIMER2(IR_checkPolarity,InitRefine_);
        m_eMesh.setProperty("AdaptedMeshVerifier::checkPolarity","initializeRefine");
        AdaptedMeshVerifier::checkPolarity(m_eMesh);
      }

      REF_LTRACE("initializeRefine: after checkPolarity");

      {
        TIMER2(IR_checkBPVali,InitRefine_);
        std::vector<stk::mesh::EntityRank>& ranks = m_ranks;
        // check logic of break pattern setup and also build ranks used vector
        checkBreakPatternValidityAndBuildRanks(ranks, &timerInitRefine_);
      }

      {
        TIMER2(IR_getRInf, InitRefine_);
        getRefinementInfo(m_ranks);
      }
    }

    void Refiner::
    getRefinementInfo(std::vector<stk::mesh::EntityRank>& ranks)
    {
      ///////////////////////////////////////////////////////////
      ///// Get info on refinements that will be done
      ///////////////////////////////////////////////////////////
      if (1)
        {

          unsigned num_orig_nodes=0;
          {
            std::vector<size_t> count1 ;

            stk::mesh::count_entities(stk::mesh::Selector(m_eMesh.get_fem_meta_data()->universal_part()) , *m_eMesh.get_bulk_data(), count1 );
            if (count1.size() < 3)
              {
                throw std::logic_error("logic error in Refiner m_refinementInfoaByType");
              }
            num_orig_nodes = count1[0];
          }

          m_refinementInfo.m_refinementInfoByType.resize(ranks.size());

          stk::mesh::PartVector fromPartsAll;

          for (unsigned irank = 0; irank < ranks.size(); irank++)
            {
              stk::mesh::PartVector * fromParts = &(m_breakPattern[irank]->getFromParts());
              if (fromParts)
                {
                  for (unsigned ipart = 0; ipart < fromParts->size(); ipart++)
                    {
                      // print message about beams && shells
                      if (1 && m_eMesh.get_rank()==0)
                        {
                          stk::mesh::Part *part = (*fromParts)[ipart];
                          const CellTopologyData * part_cell_topo_data = m_eMesh.get_cell_topology(*part);
                          if (part_cell_topo_data)
                            {
                              shards::CellTopology part_cell_topo(part_cell_topo_data);
                              std::string ptopo_name = part_cell_topo.getName();
                              bool printAll = false;
                              if (printAll)
                                {
                                  std::cout << "P[0] INFO: will refine block: " << part->name() << " with topology= " << part_cell_topo.getName() << std::endl;
                                }
                            }
                        }

                      fromPartsAll.push_back((*fromParts)[ipart]);
                    }
                }
            }

          for (unsigned irank = 0; irank < ranks.size(); irank++)
            {
              EXCEPTWATCH;
              unsigned elementType = m_breakPattern[irank]->getFromTypeKey();
              shards::CellTopology cell_topo(m_breakPattern[irank]->getFromTopology());
              //std::cout << "a1 toposrk cell_topo= " << cell_topo.getName() << std::endl;

              stk::mesh::Selector selector(m_eMesh.get_fem_meta_data()->locally_owned_part());
              if (fromPartsAll.size())
                {
                  selector = stk::mesh::Selector();
                  for (unsigned ipart = 0; ipart < fromPartsAll.size(); ipart++)
                    {
                      stk::mesh::Part *part = fromPartsAll[ipart];
                      const CellTopologyData * part_cell_topo_data = m_eMesh.get_cell_topology(*part);
                      if (part_cell_topo_data)
                        {
                          shards::CellTopology part_cell_topo(part_cell_topo_data);
                          if (part_cell_topo.getKey() == elementType)
                            {
                              selector = selector | *part;
                            }
                        }
                    }
                  selector = selector & (stk::mesh::Selector(m_eMesh.get_fem_meta_data()->locally_owned_part()));
                }
              std::vector<size_t> count ;
              stk::mesh::count_entities( selector, *m_eMesh.get_bulk_data(), count );
              if (count.size() < 3)
                {
                  throw std::logic_error("logic error in Refiner m_refinementInfoByType");
                }
              unsigned n_ele = count[ ranks[irank] ];

              m_refinementInfo.m_refinementInfoByType[irank].m_rank = ranks[irank];
              m_refinementInfo.m_refinementInfoByType[irank].m_numOrigElems = n_ele;

              m_refinementInfo.m_refinementInfoByType[irank].m_numNewElems = n_ele * m_breakPattern[irank]->getNumNewElemPerElem();
              m_refinementInfo.m_refinementInfoByType[irank].m_topology = cell_topo;
              m_refinementInfo.m_refinementInfoByType[irank].m_numOrigNodes = num_orig_nodes;
              m_refinementInfo.m_refinementInfoByType[irank].m_numNewNodes = 0; // can't predict this
            }

          // sum info from all procs
          {
            stk::ParallelMachine pm = m_eMesh.get_bulk_data()->parallel();

            for (unsigned irank = 0; irank < ranks.size(); irank++)
              {
                stk::all_reduce( pm, stk::ReduceSum<1>( &m_refinementInfo.m_refinementInfoByType[irank].m_numOrigElems ) );
                stk::all_reduce( pm, stk::ReduceSum<1>( &m_refinementInfo.m_refinementInfoByType[irank].m_numNewElems ) );
                stk::all_reduce( pm, stk::ReduceSum<1>( &m_refinementInfo.m_refinementInfoByType[irank].m_numOrigNodes ) );
              }
          }

        }
      m_elementRankTypeInfo.resize(ranks.size());

      fillElementRankTypeInfo(ranks);
    }

    void Refiner::
    doMark(int num_registration_loops)
    {
#if USE_PERCEPT_PERFORMANCE_TESTING_CALLGRIND
  CALLGRIND_START_INSTRUMENTATION;
  CALLGRIND_TOGGLE_COLLECT;
#endif

      stk::diag::Timer timerAdapt_("RefineMesh", rootTimer());
      stk::diag::TimeBlock timerAdaptBlock_(timerAdapt_);

      stk::diag::Timer timerDoMark_("percept::DoMark", timerAdapt_);
      stk::diag::TimeBlock timerDoMarkBlock_(timerDoMark_);

      getRefinementInfo().full_stats_before_mark();

      REF_LTRACE("doMark: start");

      static SubDimCellData empty_SubDimCellData;
      std::vector<stk::mesh::EntityRank>& ranks = m_ranks;

      // do elements first, then any faces or edge elements

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // do the top level, all elements of this rank operation
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      fillElementRankTypeInfo(ranks);

      // FIXME warn if a topology shows up without a break pattern

      ///////////////////////////////////////////////////////////
      /////  // start top-level ranks
      ///////////////////////////////////////////////////////////

      for (int ireg=0; ireg < num_registration_loops; ++ireg) {
        // start top-level ranks

        {
          TIMER2(preMark,DoMark_);
          preMark(ireg, num_registration_loops);
        }

        ///////////////////////////////////////////////////////////
        /////  // node registration step
        ///////////////////////////////////////////////////////////

        {
          // node registration step
          EXCEPTWATCH;
          TIMER2(nodeReg,DoMark_);

          /**/  TRACE_PRINT("Refiner: beginRegistration (top-level rank)... ");
          m_nodeRegistry->dumpDB("before init");
          if (m_alwaysInitNodeRegistry)
            {
              m_nodeRegistry->initialize();
            }

          m_nodeRegistry->init_comm_all();

          // register non-ghosted elements needs for new nodes, parallel create new nodes
          m_nodeRegistry->beginRegistration(ireg,num_registration_loops);
          {
            TIMER2(LB_NodeReg,DoMark_);

            for (unsigned irank = 0; irank < ranks.size(); irank++)
            {
              unsigned elementType = m_breakPattern[irank]->getFromTypeKey();
              {
                EXCEPTWATCH;

                ElementRankTypeInfo& e_info = m_elementRankTypeInfo[irank];
                VERIFY_OP_ON(ranks[irank], ==, e_info.first,"er1");
                VERIFY_OP_ON(elementType, ==, e_info.second,"er2");

                std::vector<NeededEntityType> needed_entity_ranks;
                m_breakPattern[irank]->fillNeededEntities(needed_entity_ranks);

                bool doAllElements = true;

                doForAllElements(irank, "[0/16] Register New Nodes",
                		ranks[irank], &NodeRegistry::registerNeedNewNode,
						elementType, needed_entity_ranks,
						doAllElements);
              }
            }
          }

          if (DO_MEMORY) {
            std::string hwm = print_memory_both(m_eMesh.parallel());
            if (!m_eMesh.get_rank()) std::cout << "MEM: " << hwm << " before createNewNodesInParallel= " << std::endl;
          }

          {
            TIMER2(NR_EndReg,DoMark_);
            m_nodeRegistry->endRegistration(ireg,num_registration_loops);
          }

          if (DO_MEMORY) {
            std::string hwm = print_memory_both(m_eMesh.parallel());
            if (!m_eMesh.get_rank()) std::cout << "MEM: " << hwm << " after  createNewNodesInParallel= " << std::endl;
          }

        }
        m_nodeRegistry->dumpDB("after registration");

        ///////////////////////////////////////////////////////////
        /////  Check for remote
        ///////////////////////////////////////////////////////////

        {   // beginCheckForRemote()
          EXCEPTWATCH;
          TIMER2(checkForRemote,DoMark_);

          /**/ TRACE_PRINT("Refiner: beginCheckForRemote (top-level rank)... ");

          // now register ghosted elements needs for new nodes (this does a pack operation)
          m_nodeRegistry->beginCheckForRemote(ireg,num_registration_loops);
          {
            TIMER2(LB_CheckForRem,DoMark_);
            for (unsigned irank = 0; irank < ranks.size(); irank++)
            {
              unsigned elementType = m_breakPattern[irank]->getFromTypeKey();
              {
                EXCEPTWATCH;

                std::vector<NeededEntityType> needed_entity_ranks;
                m_breakPattern[irank]->fillNeededEntities(needed_entity_ranks);

                bool doAllElements = false;  // only do ghost elements if false
                if (!m_nodeRegistry->s_use_new_ownership_check)
                  {
                    doAllElements = false;
                  }
                else
                  {
                    doAllElements = true;
                  }
                doForAllElements(irank, "[1/16] Check For Comm",
                		ranks[irank], &NodeRegistry::checkForRemote,
						elementType, needed_entity_ranks, doAllElements);
              }
            }
          }
          m_nodeRegistry->endCheckForRemote(ireg,num_registration_loops);                /**/   TRACE_PRINT("Refiner: endCheckForRemote (top-level rank)... ");

        }

        ///////////////////////////////////////////////////////////
        /////  Get from remote
        ///////////////////////////////////////////////////////////
        /// communicate all-to-all the new node creation information which also updates the node registry so it can
        /// be queried locally now for any ghost or non-ghost element

        { // get from remote

          {
            TIMER2(beginAndGetFromRemote,DoMark_);

            m_nodeRegistry->beginGetFromRemote(ireg,num_registration_loops);
            {
              TIMER2(LB_GetFromRem,DoMark_);
              for (unsigned irank = 0; irank < ranks.size(); irank++)
                {
                  unsigned elementType = m_breakPattern[irank]->getFromTypeKey();
                  {
                    EXCEPTWATCH;

                    std::vector<NeededEntityType> needed_entity_ranks;
                    m_breakPattern[irank]->fillNeededEntities(needed_entity_ranks);

                    bool doAllElements = false;   // if false, ghost elements only
                    if (!m_nodeRegistry->s_use_new_ownership_check)
                      {
                        doAllElements = false;
                      }
                    else
                      {
                        doAllElements = true;   // if false, ghost elements only
                      }
                    doForAllElements(irank, "[2/16] Get From Remote",
                    		ranks[irank], &NodeRegistry::getFromRemote,
							elementType, needed_entity_ranks, doAllElements);
                  }
                }
            }
          }

          {
            TIMER2(endGetFromRemote, DoMark_);
            m_nodeRegistry->endGetFromRemote(ireg,num_registration_loops);
          }

          m_nodeRegistry->dumpDB("after endGetFromRemote");

        }  // get from remote

        bool break_loop = false;
        {
          TIMER2(postMark, DoMark_);
          break_loop = postMark(ireg, num_registration_loops);
        }

        {
          mod_begin_timer(*m_eMesh.get_bulk_data(), timerDoMark_);
          mod_end_timer(  *m_eMesh.get_bulk_data(), timerDoMark_, "Mark");
        }

        if (break_loop)
          break;

      } // start top-level ranks - num_registration_loops end loop

      getRefinementInfo().full_stats_after_mark();

      REF_LTRACE("doMark: end");
    }

    void Refiner::
    fillElementRankTypeInfo(std::vector<stk::mesh::EntityRank>& ranks)
    {
      m_elementRankTypeInfo.resize(ranks.size());

      for (unsigned irank = 0; irank < ranks.size(); irank++)
        {
          EXCEPTWATCH;
          unsigned elementType = m_breakPattern[irank]->getFromTypeKey();
          shards::CellTopology cell_topo(m_breakPattern[irank]->getFromTopology());

          if (TRACE_STAGE_PRINT) std::cout << "tmp Refiner:: irank = " << irank << " ranks[irank] = " << ranks[irank]
                                           << " elementType= " << elementType
                                           << " cell_topo= " << cell_topo.getName()
                                           << std::endl;

          m_elementRankTypeInfo[irank] = ElementRankTypeInfo(ranks[irank], elementType);
        }
    }

    void Refiner::
    doRefine()
    {
#if USE_PERCEPT_PERFORMANCE_TESTING_CALLGRIND
  CALLGRIND_START_INSTRUMENTATION;
  CALLGRIND_TOGGLE_COLLECT;
#endif
      stk::diag::Timer timerAdapt_("RefineMesh", rootTimer());
      stk::diag::TimeBlock timerAdaptBlock_(timerAdapt_);

      stk::diag::Timer timerDoRefine_("percept::DoRefine", timerAdapt_);
      stk::diag::TimeBlock timerDoRefineBlock_(timerDoRefine_);

      getRefinementInfo().full_stats_before_refine();

      m_eMesh.m_nodeRegistry = (void *) &getNodeRegistry();
      REF_LTRACE("doRefine: start");

      static SubDimCellData empty_SubDimCellData;
      std::vector<stk::mesh::EntityRank>& ranks = m_ranks;

      m_eMesh.initializeIdServer();

      FixSideSetsSelectorRefine fss_ref(m_eMesh);

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // for each element type, in top-down rank order, do the rest of the refinement operations
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      mod_begin_timer(*m_eMesh.get_bulk_data(), timerDoRefine_);

      size_t num_new_elements[2] = {0,0};
      for (unsigned irank = 0; irank < ranks.size(); irank++)
        {
          EXCEPTWATCH;

          unsigned elementType = m_breakPattern[irank]->getFromTypeKey();

          ElementRankTypeInfo& e_info = m_elementRankTypeInfo[irank];
          VERIFY_OP_ON(elementType, ==, e_info.second, "mismatch");

          // loop over elements, build faces, edges in threaded mode (guaranteed no mem conflicts)
          // (note: invoke UniformRefinerPattern: what entities are needed)
          std::vector<NeededEntityType> needed_entity_ranks;
          m_breakPattern[irank]->fillNeededEntities(needed_entity_ranks);
          std::vector<stk::mesh::Entity> new_elements, ft_new_elements;

          {
            EXCEPTWATCH;

            // count num new elements needed on this proc (served by UniformRefinerPattern)
            unsigned num_elem_not_ghost = 0;
            {
              TIMER2(LB_RegNewNodes,DoRefine_);
              TIMER2(RegNewNodes,DoRefine_);
              std::vector<stk::mesh::Entity> elements;
              num_elem_not_ghost = countAndGatherAllElements(irank, ranks[irank], elementType, elements);
            }

            unsigned num_elem_needed = 0;
            {
              TIMER2(EstimateNumNewElems,DoRefine_);
              getRefinementInfo().full_stats_before_estimate(ranks[irank], num_elem_not_ghost);
              num_elem_needed = m_breakPattern[irank]->estimateNumberOfNewElements(m_eMesh, ranks[irank], *m_nodeRegistry, num_elem_not_ghost);
              // FIXME use this simple formula when doing UMR
              //num_elem_needed = num_elem_not_ghost * m_breakPattern[irank]->getNumNewElemPerElem();
              unsigned nelNeedMarked = num_elem_needed / m_breakPattern[irank]->getNumNewElemPerElem();
              getRefinementInfo().full_stats_after_estimate(ranks[irank], nelNeedMarked);
              num_new_elements[0] += num_elem_needed;
            }

            // create new entities on this proc
            {
              TIMER2(CreateEntities,DoRefine_);
              new_elements.resize(0);                                                /**/ TRACE_PRINT("Refiner: createEntities... ranks[irank]==ranks[0] ");
              ft_new_elements.resize(0);
              auto breakPattern = m_breakPattern[irank];
              std::string bpName = PerceptMesh::demangle(typeid(*breakPattern).name());

              if (UniformRefinerPatternBase::USE_DECLARE_ELEMENT_SIDE && ranks[irank] == m_eMesh.side_rank())
                {
                  new_elements.resize(0);
                }
              else
                {
#if USE_CREATE_ENTITIES
                  m_eMesh.createEntities( ranks[irank], num_elem_needed, new_elements);
#else
                  m_eMesh.getEntitiesUsingIdServer( ranks[irank], num_elem_needed, new_elements);
#endif
                }

              const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(stk::topology::ELEMENT_RANK + 1u);
#if USE_CREATE_ENTITIES
              m_eMesh.createEntities( FAMILY_TREE_RANK, num_elem_needed, ft_new_elements);  /**/ TRACE_PRINT("Refiner: createEntities... ranks[irank]==ranks[0] done ");
#else
              m_eMesh.getEntitiesUsingIdServer( FAMILY_TREE_RANK, num_elem_needed, ft_new_elements);
#endif
            }
          }

          /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          ///  Global element ops: here's where we e.g. connect the new elements by declaring new relations
          /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          /**/                                                TRACE_PRINT("Refiner: createElementsAndNodesAndConnectLocal... ");
          /**/                                                TRACE_CPU_TIME_AND_MEM_0(CONNECT_LOCAL);

          size_t num_actual_new_elems = 0;
          {
            TIMER2(LB_CreateElements,DoRefine_);
            std::vector<stk::mesh::Entity>::iterator new_elements_pool_end_iter;
            std::vector<stk::mesh::Entity>::iterator ft_new_elements_pool_end_iter;

            num_actual_new_elems =
              createElementsAndNodesAndConnectLocal(irank, ranks[irank], m_breakPattern[irank], e_info.second, needed_entity_ranks, new_elements, ft_new_elements,
                                                    &new_elements_pool_end_iter,
                                                    &ft_new_elements_pool_end_iter);

            if (ranks[irank] == m_eMesh.element_rank())
              {
                fss_ref.add_elements(new_elements.begin(), new_elements_pool_end_iter);
              }
            num_new_elements[1] += num_actual_new_elems;
          }

          /**/                                                TRACE_CPU_TIME_AND_MEM_1(CONNECT_LOCAL);
          /**/                                                TRACE_PRINT("Refiner: createElementsAndNodesAndConnectLocal...done ");

          /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          ///  Global node loop operations:  this is where we perform ops like adding new nodes to the right parts, interpolating fields, etc.
          /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        } // irank

      /**/                                                TRACE_PRINT("Refiner: addToExistingParts [etc.]... ");
#if !STK_ADAPT_URP_LOCAL_NODE_COMPS
      if (1)
        {
          EXCEPTWATCH;
          // only need to do this once: the map is fully built and we loop over the map's faces/edges, which are fixed after the getFromRemote step

          doProgressPrint("Stage: [5/16] Add to existing parts...");

          {
            TIMER2(AddToParts,DoRefine_);
            m_nodeRegistry->addToExistingPartsNew();
          }

          doProgressPrint("Stage: [6/16] Prolongate fields...");

          {
            TIMER2(ProlongCoord,DoRefine_);
            m_nodeRegistry->prolongate(m_eMesh.get_coordinates_field());
          }

#if defined(STK_BUILT_FOR_SIERRA)
          if (m_rbar_names.size())
            m_nodeRegistry->add_rbars(m_rbar_names);
#endif
        }
#endif

      /**/                                                TRACE_PRINT("Refiner: addToExistingParts [etc.] ...done ");

      REF_LTRACE("doRefine: fix_side_sets...");

      doProgressPrint("Stage: [7/16] Remove empty elements...");

      if (1)
      {
        TIMER2(RemEmptyElem,DoRefine_);
        removeEmptyElements();
      }

      bool reduced_mod_end = true;
      if (m_eMesh.getProperty("percept_reduced_mod_end") == "false")
        reduced_mod_end = false;

      /***********************/                           TRACE_PRINT("Refiner: fix_side_sets_2 ");
      {
        doProgressPrint("Stage: [8/16] Fix side sets...");

        mod_end_timer(  *m_eMesh.get_bulk_data(), timerDoRefine_, "Refine0");
        mod_begin_timer(*m_eMesh.get_bulk_data(), timerDoRefine_);

        if (!getIgnoreSideSets())
          {
            TIMER2(FSS_refine,DoRefine_);
            m_eMesh.setProperty("FixSideSets::fix_side_sets_2","refiner");

            fix_side_sets_2(false,0,0, &fss_ref, "Ref1");
          }
      }
      /***********************/                           TRACE_PRINT("Refiner: fix_side_sets_2...done ");

      {
        TIMER2(ProlongFields,DoRefine_);
        m_nodeRegistry->prolongateFields(); // Needed after mod-end so that parts are parallel consistent so that fields are parallel consistent
      }

      REF_LTRACE("doRefine: fix_side_sets...done");

      if (m_doRemove)
        {
          doProgressPrint("Stage: [9/16] Remove old elements...");

          EXCEPTWATCH;

          TIMER2(DoRemove,DoRefine_);

          for (unsigned irank = 0; irank < ranks.size(); irank++)
            {
#if PERCEPT_USE_FAMILY_TREE
              if (irank == 0)
                removeFamilyTrees();
#endif
              removeOldElements(irank, ranks[irank], m_breakPattern[irank]);
              renameNewParts(ranks[irank], m_breakPattern[irank]);
//              fixSurfaceAndEdgeSetNames(ranks[irank], m_breakPattern[irank]);
            }
        }


      // remove any elements that are empty (these can exist when doing local refinement)
      {
        doProgressPrint("Stage: [10/16] Remove empty elements part 2...");

        TIMER2(RemEmptyElem,DoRefine_);
        removeEmptyElements();
      }

      if (m_removeFromNewNodesPart)
        {
          TIMER2(RemoveFromNewNodesPart,DoRefine_);
          removeFromNewNodesPart();
        }

      doProgressPrint("Stage: [11/16] Modification_end...");
      {
        mod_end_timer(  *m_eMesh.get_bulk_data(), timerDoRefine_, "Refine1");
        mod_begin_timer(*m_eMesh.get_bulk_data(), timerDoRefine_);
      }

      doProgressPrint("Stage: [12/16] Remove unattached nodes...");

      // remove nodes not referred to by elements
      {
        TIMER2(RemDangling,DoRefine_);
        removeDanglingNodes();
      }

      REF_LTRACE("doRefine: removeDanglingNodes...done");

      doProgressPrint("Stage: [13/16] Remove empty family trees...");

      {
        TIMER2(RemEmptyFT,DoRefine_);
        removeEmptyFamilyTrees();
      }

      set_active_part();
      if (m_doAddChildrenToParts)
        add_children_to_parts();

      reset_family_tree_to_node_relations();

      doProgressPrint("Stage: [14/16] Modification end # 2...");

      /**/                                                TRACE_PRINT("Refiner: mod_end...start... ");
      if (!reduced_mod_end)
        {
          // force a flush of all pending deletes, etc
          mod_end_timer(  *m_eMesh.get_bulk_data(), timerDoRefine_, "Flush1");
          mod_begin_timer(*m_eMesh.get_bulk_data(), timerDoRefine_);
          mod_end_timer(  *m_eMesh.get_bulk_data(), timerDoRefine_, "Flush2");
        }
      else
        {
          mod_end_timer(*m_eMesh.get_bulk_data(), timerDoRefine_, "Refine3");
        }

      /**/                                                TRACE_PRINT("Refiner: mod_end...done ");

      doProgressPrint("Stage: [15/16] Checks on ownership...");

      if (1)
      {
        TIMER2(ref_check_own,DoRefine_);
        check_parent_ownership();
        check_sides_on_same_proc_as_owned_element("doRefine:b4setPEF", true);
      }

      m_eMesh.set_parent_element_field();
      if (1)
        {
          RefinerUtil::save_node_registry(m_eMesh, *m_nodeRegistry, "Refiner: end");
        }

      m_eMesh.setProperty("AdaptedMeshVerifier::checkPolarity", "doRefine end");
      AdaptedMeshVerifier::checkPolarity(m_eMesh);

#if  defined(STK_PERCEPT_HAS_GEOMETRY)
      bool use_ref_mesh = true;
      if (m_eMesh.getProperty("smooth_use_reference_mesh") == "0")
        use_ref_mesh = false;
      snapAndSmooth(m_geomSnap, m_geomFile, use_ref_mesh);
#endif

      REF_LTRACE("doRefine: doRebalance...");

      doProgressPrint("Stage: [16/16] Rebalance...");

      /**/                                                TRACE_PRINT( "Refiner:doBreak ... done");
      {
        TIMER2(Rebalance,DoRefine_);
        doRebalance();
      }
      getRefinementInfo().countCurrentNodes(m_eMesh);
      getRefinementInfo().full_stats_after_refine();

      m_nodeRegistry->dumpDB("after doBreak");

      doProgressPrint("Stage: [16/16] Rebalance...done, Refinement done.");

#if USE_PERCEPT_PERFORMANCE_TESTING_CALLGRIND
  CALLGRIND_TOGGLE_COLLECT;
  CALLGRIND_STOP_INSTRUMENTATION;
#endif

    }

    void Refiner::doRebalance()
    {
      if (m_doRebalance)
        {
          if (m_eMesh.get_rank() == 0)
            {
              std::cout << "Refiner:: rebalancing... weights field= " << m_eMesh.m_weights_field << std::endl;
            }
          if (1)
            {
              RefinerUtil::save_node_registry(m_eMesh, *m_nodeRegistry, "doRebalance");
            }

          RebalanceMesh rb(m_eMesh, m_eMesh.m_weights_field, false);
          const double imb_before = rb.compute_imbalance();
          if (imb_before > m_rebalThreshold)
            {
              const double imb_after = rb.rebalance();
              if (m_eMesh.get_rank() == 0)
                {
                  std::cout << "Refiner:: imbalance before= " << imb_before << " imbalance after= " << imb_after << std::endl;
                }
              m_eMesh.m_markNone = true;
              m_nodeRegistry->setCheckForGhostedNodes(true);
              initializeDB();
              m_nodeRegistry->setCheckForGhostedNodes(false);
              m_eMesh.m_markNone = false;

              if (1)
                {
                  std::vector< const stk::mesh::FieldBase *> fields;
                  fields.push_back(m_eMesh.m_refine_field);
                  fields.push_back(m_eMesh.m_transition_element_field);
                  stk::mesh::communicate_field_data(m_eMesh.get_bulk_data()->aura_ghosting(), fields);
                }
            }
          else
            {
              if (m_eMesh.get_rank() == 0)
                {
                  std::cout << "Refiner:: no rebalance done, imbalance before= " << imb_before << " which is < threshold =  " << m_rebalThreshold << std::endl;
                }
            }
        }
    }

    void Refiner::
    finalizeRefine()
    {
    }

    /// Delete all elements that aren't child elements
    void Refiner::deleteParentElements()
    {
      stk::diag::Timer timerAdapt_("RefineMesh", rootTimer());
      stk::diag::Timer timerDoRefine_("percept::DoRefine", timerAdapt_);

      doProgressPrint("Stage: Deleting parent elements...");

      //check_sidesets_2(" deleteParentElements:: start");
      //check_sidesets(" deleteParentElements:: start");
      //check_sidesets_1(" deleteParentElements:: start");

      std::vector<stk::mesh::EntityRank> ranks_to_be_deleted;
      ranks_to_be_deleted.push_back(stk::topology::ELEMENT_RANK);
      ranks_to_be_deleted.push_back(m_eMesh.side_rank());
      if (m_eMesh.get_spatial_dim() == 3)
        ranks_to_be_deleted.push_back(m_eMesh.edge_rank());

      //std::cout << "tmp srk ranks_to_be_deleted= " << ranks_to_be_deleted << std::endl;

      elements_to_be_destroyed_type parents(*m_eMesh.get_bulk_data());
      for (unsigned irank=0; irank < ranks_to_be_deleted.size(); irank++)
        {

          const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( ranks_to_be_deleted[irank] );
          int nchild=0;
          for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;

              // only do "old" elements
              //if (!oldPartSelector(bucket))
              //  continue;

              const unsigned num_elements_in_bucket = bucket.size();
              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  stk::mesh::Entity element = bucket[iElement];
                  if (!m_eMesh.isParentElement(element, false))
                  //if (!m_eMesh.hasFamilyTree(element) || m_eMesh.isChildElement(element, true))
                  //if (!m_eMesh.hasFamilyTree(element) || m_eMesh.isChildElementLeaf(element, true))
                    {
                      // it has no family tree, so it's a leaf, or it has no children
                      ++nchild;
                    }
                  else
                    {
                      parents.insert(element);
                    }
                }
            }
          (void)nchild;
        }

      mod_begin_timer(*m_eMesh.get_bulk_data(), timerDoRefine_);

      removeFamilyTrees();

      removeElements(parents);

      fix_side_sets_2(false,0,0,0,"Ref2");

      mod_end_timer(*m_eMesh.get_bulk_data(), timerDoRefine_, "RefinerDelPar");

      doProgressPrint("Stage: Deleting parent elements, number of parent elements = ");

    }

    void Refiner::removeEmptyElements()
    {
      for (unsigned irank=0; irank < m_ranks.size(); irank++)
        {
          elements_to_be_destroyed_type list(*m_eMesh.get_bulk_data());

          const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( m_ranks[irank] );

          for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;

              const unsigned num_elements_in_bucket = bucket.size();
              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  stk::mesh::Entity element = bucket[iElement];
                  if (0 == m_eMesh.get_bulk_data()->num_connectivity(element, stk::topology::NODE_RANK))
                    {
                      list.insert(element);
                    }
                }
            }

          removeElements(list);
        }
    }

    void Refiner::removeEmptyFamilyTrees()
    {
      const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(stk::topology::ELEMENT_RANK + 1u);

      elements_to_be_destroyed_type list(*m_eMesh.get_bulk_data());

      const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( FAMILY_TREE_RANK );

      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;

          const unsigned num_elements_in_bucket = bucket.size();
          for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
            {
              stk::mesh::Entity element = bucket[iElement];
              if (m_eMesh.get_bulk_data()->has_no_relations(element))
                {
                  list.insert(element);
                }
            }
        }

      for (elements_to_be_destroyed_type::iterator it = list.begin(); it != list.end(); ++it)
        {
          if ( ! m_eMesh.get_bulk_data()->destroy_entity( *it ) )
            {
              throw std::runtime_error("couldn't destroy family tree entity");
            }
        }
    }

    void Refiner::removeDanglingNodes()
    {
      if (m_avoidClearDanglingNodes)
        return;

      SetOfEntities node_list(*m_eMesh.get_bulk_data());
      SetOfEntities pseudos(*m_eMesh.get_bulk_data());

      const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.node_rank() );

      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          stk::mesh::BulkData &mesh = bucket.mesh();

          const unsigned num_nodes_in_bucket = bucket.size();
          for (unsigned iElement = 0; iElement < num_nodes_in_bucket; iElement++)
            {
              stk::mesh::Entity node = bucket[iElement];
              size_t num_rels = mesh.count_relations(node);
              if (0) std::cout << "node= " << m_eMesh.identifier(node) << " delete= " << (0==num_rels) << std::endl;
              if (0 == num_rels)
                {
                  node_list.insert(node);
                }
            }
        }

      //if (1 && !m_eMesh.get_rank()) std::cout << "P[" << m_eMesh.get_rank() << "] tmp number of dangling nodes = " << node_list.size() << " and pseudos= " << pseudos.size() << std::endl;
      //!srk
      getNodeRegistry().clear_dangling_nodes(&node_list);

      if (0)
        {
          for (SetOfEntities::iterator itbd = pseudos.begin(); itbd != pseudos.end();  ++itbd)
            {
              stk::mesh::Entity pseudo_p = *itbd;

              if ( ! m_eMesh.get_bulk_data()->destroy_entity( pseudo_p ) )
                {
                  throw std::logic_error("Refiner::removeDanglingNodes couldn't remove pseudo");

                }
            }
        }

      for (SetOfEntities::iterator itbd = node_list.begin(); itbd != node_list.end();  ++itbd)
        {
          stk::mesh::Entity node_p = *itbd;

          if ( ! m_eMesh.get_bulk_data()->destroy_entity( node_p ) )
            {
              throw std::logic_error("Refiner::removeDanglingNodes couldn't remove node");

            }
        }

      // check for any null entities
      //std::cout << "check for any null entities..." << std::endl;
      const stk::mesh::BucketVector & elem_buckets = m_eMesh.get_bulk_data()->buckets( stk::topology::ELEMENT_RANK );

      for ( stk::mesh::BucketVector::const_iterator k = elem_buckets.begin() ; k != elem_buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;

          const unsigned num_nodes_in_bucket = bucket.size();
          for (unsigned iElement = 0; iElement < num_nodes_in_bucket; iElement++)
            {
              stk::mesh::Entity elem = bucket[iElement];
              const percept::MyPairIterRelation rels (m_eMesh, elem, m_eMesh.node_rank());
              for (unsigned j=0; j < rels.size(); j++)
                {
                  if (!m_eMesh.is_valid(rels[j].entity())) throw std::runtime_error("bad node in an element");
                }
            }
        }


    }

#if defined( STK_PERCEPT_HAS_GEOMETRY )
#if 0
    static stk::mesh::Selector getNodeWasSnappedSelector(MeshGeometry& mesh_geometry)
    {
      stk::mesh::Selector selector;
      const std::vector<GeometryEvaluator*>& geomEvaluators = mesh_geometry.getGeomEvaluators();

      for (unsigned s=0; s < geomEvaluators.size(); s++)
        {
          selector |= geomEvaluators[s]->mMesh;
        }
      return selector;
    }
#endif
#endif


#if defined( STK_PERCEPT_HAS_GEOMETRY )
    void Refiner::snapAndSmooth(bool geomSnap, std::string geomFile, bool use_ref_mesh)
    {
      //SMOOTHING_OPTIONS option = SNAP_PLUS_SMOOTH;
      SMOOTHING_OPTIONS option = USE_LINE_SEARCH_WITH_MULTIPLE_STATES;

      GeometryKernel *geomKernel = 0;
      if (geomFile.length() != 0.0)
        {
        if (geomFile.find(".3dm") != std::string::npos)
          {
            std::string m2gFile = geomFile.substr(0,geomFile.length()-3) + "m2g";

            struct stat s;
            if (0 == stat(m2gFile.c_str(), &s))
              {
#if HAVE_CUBIT
                geomKernel = new GeometryKernelPGEOM();
#else
                throw std::runtime_error("CUBIT not supported on this platform");
#endif
              }
            else
              {
#if HAVE_OPENNURBS
                geomKernel = new GeometryKernelOpenNURBS();
#else
                throw std::runtime_error("OPENNURBS not supported on this platform");
#endif
              }
          }
        else if (geomFile.find(".e") != std::string::npos || geomFile.find(".g") != std::string::npos || geomFile.find(".exo") != std::string::npos)
          {
            geomKernel = new GeometryKernelGregoryPatch(m_eMesh, false);
            if (m_eMesh.get_ioss_read_options().find("auto-decomp:yes") != std::string::npos) {
              geomKernel->set_property("auto_decomp", "true");
            }
            if (m_eMesh.get_ioss_read_options().find("large") != std::string::npos) {
              geomKernel->set_property("exo_large", "true");
            }
          }
        else if(geomFile.find(".sat") != std::string::npos)
          {
#ifdef HAVE_ACIS
    	    geomKernel = new GeometryKernelPGEOM();
#else
    	    throw std::runtime_error("ACIS not supported on this platform");
#endif
          }
        else
          {
            VERIFY_MSG("invalid file extension on --input_geometry file \n   "
                       "-- valid extensions are .3dm (OpenNURBS) or .e,.g,.exo \n"
                       "for GregoryPatch Exodus files or .sat for ACIS (assumes \n"
                       "there is also a file with the same name ending in  .m2g) - file= " + geomFile);
          }
        }

      // set to 0.0 for no checks, > 0.0 for a fixed check delta, < 0.0 (e.g. -0.5) to check against local edge length average times this |value|
      double doCheckMovement = 0.0;
      //double doCheckMovement = -1.0;

      // anything exceeding a value > 0.0 will be printed
      double doCheckCPUTime = 0.0;
      //double doCheckCPUTime = 0.1;

      MeshGeometry mesh_geometry(m_eMesh, geomKernel, doCheckMovement, doCheckCPUTime);
      GeometryFactory factory(geomKernel, &mesh_geometry);
      if (geomFile != "") {
        factory.read_file(geomFile, &m_eMesh);
      }

      switch(option) {
      case SNAP_PLUS_SMOOTH:
        {
          mesh_geometry.snap_points_to_geometry(&m_eMesh);
          if (doCheckMovement != 0.0)
            mesh_geometry.print_node_movement_summary();

          if (m_doSmoothGeometry)
            {
              smoothGeometry(&mesh_geometry, 0, option, use_ref_mesh);
              mesh_geometry.snap_points_to_geometry(&m_eMesh);
            }
        }
        break;
      case USE_LINE_SEARCH_WITH_MULTIPLE_STATES:
        {
          //VERIFY_OP_ON(m_eMesh.get_coordinates_field()->number_of_states(), ==, 3, "Must use PerceptMesh::set_num_coordinate_field_states(3) to use new smoothing.");
          stk::mesh::FieldBase *nm1_field = m_eMesh.get_field(stk::topology::NODE_RANK, "coordinates_NM1");

          if (m_doSmoothGeometry && nm1_field)
            {
              // make a copy of current non-snapped state (dst,src)
              m_eMesh.copy_field(m_eMesh.get_field(stk::topology::NODE_RANK, "coordinates_NM1"), m_eMesh.get_coordinates_field() );
            }

          // do the snap
          if (geomSnap)
            {
              mesh_geometry.snap_points_to_geometry(&m_eMesh);
            }
          if (doCheckMovement != 0.0)
            mesh_geometry.print_node_movement_summary();

          if (m_doSmoothGeometry && nm1_field)
            {
              if (!geomSnap)
                mesh_geometry.pre_process(&m_eMesh);

              // make a copy of current snapped state
              m_eMesh.copy_field(m_eMesh.get_field(stk::topology::NODE_RANK, "coordinates_N"), m_eMesh.get_coordinates_field() );

              // reset current state to non-snapped state
              m_eMesh.copy_field(m_eMesh.get_coordinates_field(), m_eMesh.get_field(stk::topology::NODE_RANK, "coordinates_NM1") );

              // option to smooth without geometry
              if (geomFile == "")
                {
                  bool option_fix_all_internal_and_outer_boundary_nodes=getFixAllBlockBoundaries();
                  stk::mesh::Selector boundarySelector;
                  if (option_fix_all_internal_and_outer_boundary_nodes)
                    {
                      stk::mesh::Part *skin_part = m_eMesh.get_skin_part("inner_skin_part", true);
                      boundarySelector = boundarySelector | *skin_part;
                    }
                  else
                    {
                      // build a selector from all surface parts
                      const stk::mesh::PartVector parts = m_eMesh.get_fem_meta_data()->get_parts();
                      for (unsigned ip=0; ip < parts.size(); ip++)
                        {
                          bool stk_auto= stk::mesh::is_auto_declared_part(*parts[ip]);

                          if (stk_auto) continue;
                          stk::mesh::EntityRank per = parts[ip]->primary_entity_rank();

                          if (per == m_eMesh.side_rank())
                            {
                              std::cout << "INFO::smoothing: freezing points on boundary: " << parts[ip]->name() << std::endl;
                              boundarySelector = boundarySelector | *parts[ip];
                            }
                        }
                    }
		  stk::mesh::Selector owned = m_eMesh.get_bulk_data()->mesh_meta_data().locally_owned_part();
                  boundarySelector = boundarySelector & owned;
                  smoothGeometry(0, &boundarySelector, option, use_ref_mesh);
                }
              else
                smoothGeometry(&mesh_geometry, 0, option, use_ref_mesh);

            }

        }
        break;
      }

      if (geomKernel) delete geomKernel;
    }
#endif

#if  defined(STK_PERCEPT_HAS_GEOMETRY)
    void Refiner::smoothGeometry(MeshGeometry* mesh_geometry, stk::mesh::Selector* selector, SMOOTHING_OPTIONS option, bool use_ref_mesh)
    {
      bool do_smoothing = true;
      if (do_smoothing)
        {
          switch(option) {
          case SNAP_PLUS_SMOOTH:
            {
              /// deprecated
            }
            break;
          case USE_LINE_SEARCH_WITH_MULTIPLE_STATES:
            {
              // geometry used for classification of fixed/non-fixed nodes
              int niter = 1001;
              double tol = 1.e-4;
              std::string sni = m_eMesh.getProperty("smoother_niter");
              std::string snt = m_eMesh.getProperty("smoother_tol");
              if (sni.length()) niter = std::stoi(sni);
              if (snt.length()) tol = std::stod(snt);

              if (m_eMesh.getProperty("smoother_type") == "algebraic")
                {
                  double drop_off_coeffs[3] = {1,1,1};  // FIXME
                  int nlayers_drop_off = 40;
                  niter = 1;
                  percept::ReferenceMeshSmootherAlgebraic pmmpsi(&m_eMesh, selector, mesh_geometry, niter, 1.e-4, drop_off_coeffs, nlayers_drop_off);
                  pmmpsi.m_do_animation = 0;
                  pmmpsi.m_use_ref_mesh = use_ref_mesh;
                  pmmpsi.run();
                }
              else
                {
                  percept::ReferenceMeshSmootherConjugateGradientImpl<STKMesh> pmmpsi(&m_eMesh, selector, mesh_geometry, niter, tol);
                  pmmpsi.m_use_ref_mesh = use_ref_mesh;
                  pmmpsi.run();
                }
            }
            break;
          }
          return;
        }
    }
#endif

    unsigned Refiner::
    countAndGatherAllElements(unsigned irank, stk::mesh::EntityRank rank, unsigned elementType, std::vector<stk::mesh::Entity> &elements)
    {
    	unsigned num_elem = 0;
    	elements.resize(0);
        stk::mesh::Selector selector(m_eMesh.get_fem_meta_data()->universal_part());
        stk::mesh::PartVector * fromParts = &(m_breakPattern[irank]->getFromParts());
        if (fromParts)
            selector = stk::mesh::selectUnion(*fromParts);

        size_t nn_ele=0;
        const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( rank );
        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            stk::mesh::Bucket & bucket = **k ;
            if (selector(bucket))
              {
                const CellTopologyData * const bucket_cell_topo_data = m_eMesh.get_cell_topology(bucket);
                shards::CellTopology topo(bucket_cell_topo_data);
                if (topo.getKey() == elementType)
                  {
                    unsigned num_elements_in_bucket = bucket.size();
                    elements.resize(nn_ele + num_elements_in_bucket);
                    for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                      {
                        stk::mesh::Entity element = bucket[iElement];
                        elements[nn_ele++] = element;

                        bool elementIsGhost = m_eMesh.isGhostElement(element);
                        if (!elementIsGhost)
                          ++num_elem;
                      }
                  }
              }
          }

        // bcarnes: not sure if this code is part of counting or applying the function to the elements?
        if (m_refinerSelector)
          {
            const bool stats = false;
            if (stats)
              {
                getRefinementInfo().full_stats_before_filter(rank, elements.size());
              }
            filterUsingRefinerSelector(rank, elements);
            if (stats)
              {
                getRefinementInfo().full_stats_after_filter(rank, elements.size());
              }
          }

        return num_elem;
    }
   unsigned Refiner::
    doForAllElements(unsigned irank, std::string function_info,
                     stk::mesh::EntityRank rank, NodeRegistry::ElementFunctionPrototype function,
                     unsigned elementType,
                     std::vector<NeededEntityType>& needed_entity_ranks,
                     bool doAllElements)
    {
      EXCEPTWATCH;

      std::vector<stk::mesh::Entity> elements;
      unsigned num_elem = countAndGatherAllElements(irank, rank, elementType, elements);

      int progress_meter_num_total = 0;
      if (m_doProgress && m_eMesh.get_rank() == 0)
        {
          progress_meter_num_total = num_elem;
          if (progress_meter_num_total)
            {
              std::ostringstream oss; oss << function_info << " [" << 100.0*((double)irank)/((double)m_ranks.size()) << " %]" << " cpu: " << m_eMesh.cpu_time() << " [sec]";
              ProgressMeterData pd(ProgressMeterData::INIT, 0.0, oss.str());
              notifyObservers(&pd);
            }
        }
      int progress_meter_when_to_post = progress_meter_num_total / m_progress_meter_frequency;
      if (0 == progress_meter_when_to_post)
        progress_meter_when_to_post = 1;
      double d_progress_meter_num_total = progress_meter_num_total;

      for (size_t iElement=0; iElement < elements.size(); ++iElement)
        {
          stk::mesh::Entity element = elements[iElement];
          const CellTopologyData * const bucket_cell_topo_data = m_eMesh.get_cell_topology(element);
          shards::CellTopology topo(bucket_cell_topo_data);

          VERIFY_OP_ON(m_eMesh.is_valid(element), ==, true, "doForAllElements bad element");
          bool elementIsGhost = m_eMesh.isGhostElement(element);

          if (doAllElements || elementIsGhost)
            {
              refineMethodApply(function, element, needed_entity_ranks, bucket_cell_topo_data);
            }

          if (m_doProgress && m_eMesh.get_rank() == 0 && iElement && (iElement % progress_meter_when_to_post == 0) )
            {
              double progress_meter_percent = 100.0*((double)num_elem)/std::max(d_progress_meter_num_total,1.0);
              std::ostringstream oss; oss << function_info << " [" << 100.0*((double)irank)/((double)m_ranks.size()) << " %]" << " cpu: " << m_eMesh.cpu_time() << " [sec]";
              ProgressMeterData pd(ProgressMeterData::RUNNING, progress_meter_percent, oss.str());
              notifyObservers(&pd);
              if (0) std::cout << "progress_meter_percent = " << progress_meter_percent << std::endl;
            }
        } // elements

      if (m_doProgress && m_eMesh.get_rank() == 0 && progress_meter_num_total)
        {
          std::ostringstream oss; oss << function_info << " [" << 100.0*((double)irank)/((double)m_ranks.size()) << " %]" << " cpu: " << m_eMesh.cpu_time() << " [sec]";
          ProgressMeterData pd(ProgressMeterData::FINI, 0.0, oss.str());
          notifyObservers(&pd);
        }

      return num_elem;
    }

    void Refiner::filterUsingRefinerSelector(stk::mesh::EntityRank rank, std::vector<stk::mesh::Entity>& elements)
    {
      if (m_refinerSelector)
        {
          if (m_refinerSelector->use_batch_filter())
            {
              m_refinerSelector->batch_filter(rank, elements);
            }
          else
            {
              size_t nele=0;
              for (size_t ii=0; ii < elements.size(); ++ii)
                {
                  if ((*m_refinerSelector)(elements[ii]) || m_eMesh.isGhostElement(elements[ii]))
                    {
                      ++nele;
                    }
                }
              //std::cout << "P[" << m_eMesh.get_rank() << "] tmp srk filterUsingRefinerSelector:: orig size= " << elements.size()
              //          << " filtered size= " << nele << std::endl;
              std::vector<stk::mesh::Entity> elements_new(nele);
              nele=0;
              for (size_t ii=0; ii < elements.size(); ++ii)
                {
                  if ((*m_refinerSelector)(elements[ii]) || m_eMesh.isGhostElement(elements[ii]))
                    {
                      elements_new[nele++] = elements[ii];
                    }
                }
              //std::cout << "P[" << m_eMesh.get_rank() << "] tmp srk 1 filterUsingRefinerSelector:: orig size= " << elements.size()
              //          << " filtered size= " << nele << std::endl;
              elements = elements_new;
            }
        }
    }

    void Refiner::
    collectElemsToRefine(const unsigned irank, stk::mesh::EntityRank rank, const unsigned elementType,
            std::vector<stk::mesh::Entity>& elems, int& jele)
    {
        stk::mesh::Selector selector(m_eMesh.get_fem_meta_data()->universal_part());
        stk::mesh::PartVector * fromParts = &(m_breakPattern[irank]->getFromParts());
        if (fromParts)
        {
            selector = stk::mesh::selectUnion(*fromParts);
        }
        if (m_excludeParts.size())
        {
            selector = selector & !stk::mesh::selectUnion(m_excludeParts);
            if (!m_eMesh.get_rank() && m_eMesh.getProperty("MeshAdapt.debug") == "true")
                std::cout << "Refiner:createElementsAndNodesAndConnectLocal with excluded parts, selector= " << selector << std::endl;
        }

        // create new elements and connect them up

        const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( rank );
        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
            stk::mesh::Bucket & bucket = **k ;
            if (selector(bucket))
            {
                const CellTopologyData * const bucket_cell_topo_data = m_eMesh.get_cell_topology(bucket);
                shards::CellTopology topo(bucket_cell_topo_data);
                if (topo.getKey() == elementType)
                {
                    unsigned num_elements_in_bucket = bucket.size();
                    jele += num_elements_in_bucket;
                    for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                    {
                        stk::mesh::Entity element = bucket[iElement];
                        if (!m_eMesh.is_valid(element))
                          {
                            throw std::runtime_error("invalid element in Refiner::collectElemsToRefine");
                          }

                        elems.push_back(element);
                    }
                }
            }
        }
    }

    size_t Refiner::
    createElementsAndNodesAndConnectLocal(unsigned irank, stk::mesh::EntityRank rank, UniformRefinerPatternBase *breakPattern,
                                          unsigned elementType,   std::vector<NeededEntityType>& needed_entity_ranks,
                                          std::vector<stk::mesh::Entity>& new_elements_pool,
                                          std::vector<stk::mesh::Entity>& ft_element_pool,
                                          std::vector<stk::mesh::Entity>::iterator * new_elements_pool_end_iter,
                                          std::vector<stk::mesh::Entity>::iterator * ft_new_elements_pool_end_iter
                                          )
    {
      EXCEPTWATCH;
      static NewSubEntityNodesType s_new_sub_entity_nodes(percept::EntityRankEnd);

      NewSubEntityNodesType& new_sub_entity_nodes = s_new_sub_entity_nodes;

      std::vector<stk::mesh::Entity>::iterator element_pool_it = new_elements_pool.begin();
      std::vector<stk::mesh::Entity>::iterator ft_element_pool_it = ft_element_pool.begin();

      int jele = 0;

      std::vector<stk::mesh::Entity> elems;
      collectElemsToRefine(irank, rank, elementType, elems, jele);

      const int nele = jele;
      jele = 0;

      if (nele == 0)
      {
          *new_elements_pool_end_iter = element_pool_it;
          *ft_new_elements_pool_end_iter = ft_element_pool_it;
          return 0;
      }

      stk::mesh::Entity first_element = elems[0];
      const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(first_element);
      shards::CellTopology cell_topo(cell_topo_data);

      for (int iElement = 0; iElement < nele; iElement++)
        {
          stk::mesh::Entity element = elems[iElement];

          if (m_proc_rank_field && rank == stk::topology::ELEMENT_RANK)
            {
              double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(m_proc_rank_field) , element );
              fdata[0] = double(m_eMesh.owner_rank(element));
            }
          // FIXME

          // skip elements that are already a parent (if there's no family tree yet, it's not a parent, so avoid throwing an error if isParentElement)
          const bool check_for_family_tree = false;
          bool isParent = m_eMesh.isParentElement(element, check_for_family_tree);
          if (isParent)
            continue;

          if (m_eMesh.isGhostElement(element))
              continue;

          if (createNewNeededNodeIds(cell_topo_data, element, needed_entity_ranks,
                  new_sub_entity_nodes, breakPattern))
          {
              throw std::logic_error("needed_entity_ranks[ineed_ent].second");
          }

          breakPattern->createNewElements(m_eMesh, *m_nodeRegistry, element, new_sub_entity_nodes,
                  element_pool_it, ft_element_pool_it, m_proc_rank_field);

          ++jele;
        }

      size_t num_new_elems = element_pool_it - new_elements_pool.begin();

      *new_elements_pool_end_iter = element_pool_it;
      *ft_new_elements_pool_end_iter = ft_element_pool_it;

      return num_new_elems;
    }

    /// create a list of nodes from the new nodes that can be easily deciphered by the UniformRefinerPattern
    /// Returns the 3D array new_sub_entity_nodes[entity_rank][ordinal_of_sub_dim_entity][ordinal_of_node_on_sub_dim_entity]

    bool Refiner::
    createNewNeededNodeIds(const CellTopologyData * const cell_topo_data,
                           const stk::mesh::Entity element, std::vector<NeededEntityType>& needed_entity_ranks,
                           NewSubEntityNodesType& new_sub_entity_nodes, UniformRefinerPatternBase *breakPattern)
    {
      EXCEPTWATCH;

      NodeRegistry& nodeRegistry = *m_nodeRegistry;

      const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);

      // CHECK - cache this
      for (unsigned ineed_ent=0; ineed_ent < needed_entity_ranks.size(); ineed_ent++)
        {
          unsigned numSubDimNeededEntities = 0;

          // special case of face in 3d or edge in 2d
          if (needed_entity_ranks[ineed_ent].first == m_eMesh.entity_rank(element))
            {
              numSubDimNeededEntities = 1;
            }
          else if (needed_entity_ranks[ineed_ent].first == m_eMesh.edge_rank())
            {
              numSubDimNeededEntities = cell_topo_data->edge_count;
            }
          else if (needed_entity_ranks[ineed_ent].first == m_eMesh.face_rank())
            {
              numSubDimNeededEntities = cell_topo_data->side_count;
            }
          else if (needed_entity_ranks[ineed_ent].first == stk::topology::ELEMENT_RANK)
            {
              numSubDimNeededEntities = 1;
            }

          if (static_cast<unsigned>(needed_entity_ranks[ineed_ent].first) >= new_sub_entity_nodes.size())
            {
              throw std::logic_error("Refiner::createNewNeededNodeIds logic err #1");
            }
          new_sub_entity_nodes[needed_entity_ranks[ineed_ent].first].resize(numSubDimNeededEntities);

          // ensure nodes don't get inadvertently reused
          for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
            {
              new_sub_entity_nodes[needed_entity_ranks[ineed_ent].first][iSubDimOrd].resize(0);
            }

          for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
            {
              // CHECK
              NodeIdsOnSubDimEntityType* nodeIds_onSE_ptr = nodeRegistry.getNewNodesOnSubDimEntity(element, needed_entity_ranks[ineed_ent].first, iSubDimOrd);
              if (nodeIds_onSE_ptr == 0)
                {
                  new_sub_entity_nodes[needed_entity_ranks[ineed_ent].first][iSubDimOrd].resize(0);
                  continue;
                }
              NodeIdsOnSubDimEntityType& nodeIds_onSE = *nodeIds_onSE_ptr;

              if (nodeIds_onSE.size() == 0)
                {
                  new_sub_entity_nodes[needed_entity_ranks[ineed_ent].first][iSubDimOrd].resize(0);
                  continue;
                }

              if (needed_entity_ranks[ineed_ent].third.size())
                {
                  VERIFY_OP_ON(needed_entity_ranks[ineed_ent].third.size(), ==, numSubDimNeededEntities, "bad size");
                  if (!needed_entity_ranks[ineed_ent].third[iSubDimOrd])
                    {
                      new_sub_entity_nodes[needed_entity_ranks[ineed_ent].first][iSubDimOrd].resize(0);
                      continue;
                    }
                }

              unsigned num_new_nodes_needed = needed_entity_ranks[ineed_ent].second;

              if (num_new_nodes_needed == 0)
                  return true;

              if (iSubDimOrd >= new_sub_entity_nodes[needed_entity_ranks[ineed_ent].first].size())
                  throw std::logic_error("Refiner::createNewNeededNodeIds logic err #2");

              new_sub_entity_nodes[needed_entity_ranks[ineed_ent].first][iSubDimOrd].resize(num_new_nodes_needed);

              for (unsigned i_new_node = 0; i_new_node < num_new_nodes_needed; i_new_node++)
                {
                  if (!m_eMesh.is_valid(nodeIds_onSE[i_new_node]))
                    {
                      stk::mesh::Entity node1 = m_eMesh.get_bulk_data()->get_entity(stk::topology::NODE_RANK, nodeIds_onSE.m_entity_id_vector[i_new_node]);

                      if (!m_eMesh.is_valid(node1))
                        {
                          if (!m_nodeRegistry->getUseCustomGhosting())
                          {
                            static stk::mesh::PartVector empty_parts;
                            node1 = m_eMesh.get_bulk_data()->declare_node(nodeIds_onSE.m_entity_id_vector[i_new_node], empty_parts);
                          }
                          if (!m_eMesh.is_valid(node1))
                          {
                            throw std::logic_error("Refiner::createNewNeededNodeIds logic err #4");
                          }
                        }
                      nodeIds_onSE[i_new_node] = node1;
                      VERIFY_OP_ON(m_eMesh.identifier(node1), ==, nodeIds_onSE.m_entity_id_vector[i_new_node], "Refiner::createNewNeededNodeIds logic err #4.1");
                    }

                  new_sub_entity_nodes[needed_entity_ranks[ineed_ent].first][iSubDimOrd][i_new_node] = m_eMesh.identifier(nodeIds_onSE[i_new_node]);
                }
            }
        }
      return false;
    }

    /** Creates a map of element sides to their higher-dimensional base elements
     */

#define EXTRA_PRINT_UR_BESDB 0

    void Refiner::
    buildElementSideDB(SubDimCellToDataMap& cell_2_data_map)
    {

    }


    /** Sets orientations and associativity of elements to sub-dimensional faces/edges after refinement.
     */
#define EXTRA_PRINT_UR_FES 0

    struct FindPair {
      const std::string& m_side_part_name;
      const std::string& m_elem_part_name;
      FindPair(const std::string& side_part_name, const std::string& elem_part_name) : m_side_part_name(side_part_name), m_elem_part_name(elem_part_name) {}

      bool operator()(std::pair<stk::mesh::Part* const, stk::mesh::PartVector >& iter)
      {
        bool side_found =  iter.first->name() == m_side_part_name;
        bool elem_found = false;
        for (unsigned i=0; i < iter.second.size(); i++)
          {
            if (m_elem_part_name == iter.second[i]->name())
              {
                elem_found = true;
                break;
              }
          }
        return side_found && elem_found;
      }
    };

#define DEBUG_GSPR 0

    static void add_if_not_present(stk::mesh::Part *side_part, stk::mesh::Part *elem_part, SidePartMap& side_part_map)
    {
      SidePartMap::iterator found = std::find_if(side_part_map.begin(), side_part_map.end(), FindPair(side_part->name(), elem_part->name()));
      if (found != side_part_map.end())
        return;

      stk::mesh::PartVector& epv_add = side_part_map[side_part];
      epv_add.push_back(elem_part);
      if (DEBUG_GSPR) std::cout << "Refiner::get_side_part_relations: add_to_map = " << side_part->name() << " elem_part= " << elem_part->name() << std::endl;
    }


    // determine side part to elem part relations
    //static
    void Refiner::
    get_side_part_relations(PerceptMesh& eMesh, bool checkParentChild, SidePartMap& side_part_map, bool debug)
    {
      EXCEPTWATCH;

      //stk::mesh::EntityRank node_rank = eMesh.node_rank();
      stk::mesh::EntityRank edge_rank = eMesh.edge_rank();
      stk::mesh::EntityRank side_rank = eMesh.side_rank();
      stk::mesh::EntityRank element_rank = stk::topology::ELEMENT_RANK;

      int spatialDim = eMesh.get_spatial_dim();

      stk::mesh::EntityRank side_rank_iter_begin = side_rank;
      stk::mesh::EntityRank side_rank_iter_end = side_rank;
      if (spatialDim == 3)
        {
          side_rank_iter_begin = edge_rank;
        }

      if (debug || DEBUG_GSPR)
        {
          bool doAll = true;
          const stk::mesh::PartVector parts = eMesh.get_fem_meta_data()->get_parts();
          for (unsigned ip=0; ip < parts.size(); ip++)
            {
              bool stk_auto= stk::mesh::is_auto_declared_part(*parts[ip]);
              //const CellTopologyData *const topology = eMesh.get_cell_topology(*parts[ip]);
              if (stk_auto && !doAll)
                continue;
              if (eMesh.get_rank() == 0) std::cout << "tmp srk get_side_part_relations:: parts[ip]-> == " << parts[ip]->name() << std::endl;
            }
        }

      // get super-relations (side_part.name() --> elem_part.name())
      for (stk::mesh::EntityRank side_rank_iter = side_rank_iter_begin; side_rank_iter <= side_rank_iter_end; side_rank_iter++)
        {
          const stk::mesh::BucketVector & side_buckets = eMesh.get_bulk_data()->buckets( side_rank_iter );
          for ( stk::mesh::BucketVector::const_iterator it_side_bucket = side_buckets.begin() ; it_side_bucket != side_buckets.end() ; ++it_side_bucket )
            {
              stk::mesh::Bucket & side_bucket = **it_side_bucket ;
              stk::mesh::PartVector const& side_parts = side_bucket.supersets();

              if (DEBUG_GSPR)
                {
                  std::cout << "\nside_bucket.supersets() =  " << std::endl;
                  for (unsigned isp=0; isp < side_parts.size(); isp++)
                    {
                      bool stk_auto= stk::mesh::is_auto_declared_part(*side_parts[isp]);
                      const CellTopologyData *const topology = eMesh.get_cell_topology(*side_parts[isp]);
                      unsigned per = side_parts[isp]->primary_entity_rank();
                      std::cout << "superset= " << side_parts[isp]->name() << " stk_auto= " << stk_auto << " topology= " << topology << " primary_entity_rank= " << per << std::endl;
                    }
                }

              stk::mesh::PartVector elem_parts = side_parts; // since high-rank parts are in side_parts already
              for (unsigned isp=0; isp < side_parts.size(); isp++)
                {
                  if ( stk::mesh::is_auto_declared_part(*side_parts[isp]) )
                    continue;

                  const AutoPart *side_auto_part = side_parts[isp]->attribute<AutoPart>();
                  if (side_auto_part)
                    continue;

                  stk::mesh::EntityRank per = side_parts[isp]->primary_entity_rank();
                  if (per != side_rank_iter)
                    continue;

                  for (unsigned iep=0; iep < elem_parts.size(); iep++)
                    {
                      if ( stk::mesh::is_auto_declared_part(*elem_parts[iep]) )
                        continue;

                      const AutoPart *auto_part = elem_parts[iep]->attribute<AutoPart>();
                      if (elem_parts[iep]->name().find(UniformRefinerPatternBase::getOldElementsPartName()) != std::string::npos)
                        {
                          if (!auto_part) throw std::runtime_error("Refiner::get_side_part_relations: bad old part attribute for auto");
                        }

                      if (auto_part)
                        {
                          continue;
                        }

                      if (elem_parts[iep]->primary_entity_rank() == element_rank)
                        {
                          if (elem_parts[iep] != side_parts[isp])
                            {
                              SidePartMap::iterator found = side_part_map.find(side_parts[isp]);
                              if (found == side_part_map.end())
                                {
                                  side_part_map[side_parts[isp]] = stk::mesh::PartVector(1, elem_parts[iep]);
                                }
                              else
                                {
                                  stk::mesh::PartVector& epv = found->second;
                                  stk::mesh::PartVector::iterator epv_found = std::find(epv.begin(), epv.end(), elem_parts[iep]);
                                  if (epv_found == epv.end())
                                    {
                                      epv.push_back(elem_parts[iep]);
                                    }
                                }
                            }
                        }
                    }
                }
            }

          // add parts created by UnrefinerPattern
          const stk::mesh::PartVector parts = eMesh.get_fem_meta_data()->get_parts();
          SidePartMap side_part_map_copy = side_part_map;

          SidePartMap::iterator iter;
          for (iter = side_part_map_copy.begin(); iter != side_part_map_copy.end(); iter++)
            {
              stk::mesh::Part *side_part = iter->first;
              std::string side_part_name = side_part->name();
              const stk::mesh::PartVector *side_pv  = side_part->attribute<stk::mesh::PartVector>();
              const stk::mesh::PartVector& epv = iter->second;
              for (unsigned iepv=0; iepv < epv.size(); iepv++)
                {
                  stk::mesh::Part *elem_part = epv[iepv];
                  std::string elem_part_name = elem_part->name();
                  if (DEBUG_GSPR) std::cout << "looking for other parts: side_part = " << side_part->name() << " elem_part= " << elem_part_name << std::endl;

                  const stk::mesh::PartVector *elem_pv  = elem_part->attribute<stk::mesh::PartVector>();

                  stk::mesh::PartVector side_pv1;
                  if (side_pv) side_pv1 = *side_pv;
                  stk::mesh::PartVector elem_pv1;
                  if (elem_pv) elem_pv1 = *elem_pv;
                  side_pv1.push_back(side_part);
                  elem_pv1.push_back(elem_part);
                  for (unsigned iside_pv = 0; iside_pv < side_pv1.size(); iside_pv++)
                    {
                      for (unsigned ielem_pv = 0; ielem_pv < elem_pv1.size(); ielem_pv++)
                        {
                          stk::mesh::Part *sidep = side_pv1[iside_pv];
                          stk::mesh::Part *elemp = elem_pv1[ielem_pv];
                          add_if_not_present(sidep, elemp, side_part_map);
                        }
                    }
                }
            }
        }

      if ((debug || DEBUG_GSPR) && !eMesh.get_rank())
        {
          SidePartMap::iterator iter;
          for (iter = side_part_map.begin(); iter != side_part_map.end(); iter++)
            {
              stk::mesh::PartVector& epv = iter->second;
              for (unsigned iepv=0; iepv < epv.size(); iepv++)
                {
                  std::cout << "Refiner::get_side_part_relations: side_part = " << std::setw(50) << iter->first->name() << " elem_part= " << std::setw(50) << epv[iepv]->name() << std::endl;
                }
            }
        }
    }

#undef EXTRA_PRINT_UR_FES

    void Refiner::removeFamilyTrees()
    {
      EXCEPTWATCH;

      doProgressPrint("Stage: Removing family_trees... ");

      elements_to_be_destroyed_type elements_to_be_destroyed(*m_eMesh.get_bulk_data());

      const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(stk::topology::ELEMENT_RANK + 1u);
      const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( FAMILY_TREE_RANK );

      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          //if (removePartSelector(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_elements_in_bucket = bucket.size();

              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  stk::mesh::Entity element   = bucket[iElement];
                  stk::mesh::Entity element_p = element;

                  elements_to_be_destroyed.insert(element_p);
                }
            }
        }
      removeElements(elements_to_be_destroyed);

      doProgressPrint("Stage: Removing family_trees, size() = ");
    }

    void Refiner::
    removeOldElements(unsigned irank, stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern)
    {
      EXCEPTWATCH;

      const stk::mesh::Part *oldPart = m_eMesh.getPart(breakPattern->getOldElementsPartName()+toString(rank));

      if (1 && oldPart)
        {
          const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(*oldPart);
          std::string ct_name = (cell_topo_data ? cell_topo_data->name : "");
          //std::cout << "tmp removeOldElements::name= " << oldPart->name() << " for rank= " << rank << " topology= " << ct_name << std::endl;
        }

      if (!oldPart)
        {
          std::cout << "name= " << breakPattern->getOldElementsPartName()+toString(rank) << std::endl;
          throw std::runtime_error("oldPart is null");
        }

      stk::mesh::Selector removePartSelector (*oldPart);

      const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( rank );

      elements_to_be_destroyed_type elements_to_be_destroyed(*m_eMesh.get_bulk_data());

      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (removePartSelector(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_elements_in_bucket = bucket.size();

              if (0)
                {
                  std::string str;
                  stk::mesh::PartVector const& pv = bucket.supersets();
                  for (unsigned ip = 0; ip < pv.size(); ip++)
                    {
                      str += " "+pv[ip]->name();
                    }
                  std::cout << "P[" << m_eMesh.get_rank() << "] removing elements in bucket of parts: " << str << std::endl;
                }

              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  stk::mesh::Entity element = bucket[iElement];
                  stk::mesh::Entity element_p = element;

                  if (!m_eMesh.isGhostElement(element))
                    {
                      elements_to_be_destroyed.insert(element_p);
                    }
                }
            }
        }

      removeElements(elements_to_be_destroyed, irank);
    }

    void Refiner::removeElements(elements_to_be_destroyed_type& elements_to_be_destroyed, unsigned irank)
    {
      elements_to_be_destroyed_type elements_to_be_destroyed_pass2(*m_eMesh.get_bulk_data());

      int progress_meter_num_total = elements_to_be_destroyed.size();
      if (m_doProgress && m_eMesh.get_rank() == 0 && progress_meter_num_total)
        {
          std::ostringstream oss; oss << "Delete Elements " << " pass [" << 100.0*((double)irank)/((double)m_ranks.size()) << " %]" << " cpu: " << m_eMesh.cpu_time() << " [sec]";
          ProgressMeterData pd(ProgressMeterData::INIT, 0.0, oss.str());
          notifyObservers(&pd);
        }
      unsigned num_elem = 0;
      int progress_meter_when_to_post = progress_meter_num_total / m_progress_meter_frequency;
      if (0 == progress_meter_when_to_post)
        progress_meter_when_to_post = 1;
      double d_progress_meter_num_total = progress_meter_num_total;

      for (elements_to_be_destroyed_type::iterator itbd = elements_to_be_destroyed.begin(); itbd != elements_to_be_destroyed.end();  ++itbd)
        {
          stk::mesh::Entity element_p = *itbd;

          if (m_doProgress && m_eMesh.get_rank() == 0 && (num_elem % progress_meter_when_to_post == 0) )
            {
              double progress_meter_percent = 100.0*((double)num_elem)/std::max(d_progress_meter_num_total,1.0);
              std::ostringstream oss; oss << "Delete Elements " << " pass [" << 100.0*((double)irank)/((double)m_ranks.size()) << " %]" << " cpu: " << m_eMesh.cpu_time() << " [sec]";
              ProgressMeterData pd(ProgressMeterData::RUNNING, progress_meter_percent, oss.str());
              notifyObservers(&pd);
            }

          ++num_elem;

          if ( ! m_eMesh.get_bulk_data()->destroy_entity( element_p ) )
            {
              elements_to_be_destroyed_pass2.insert(element_p);
              //throw std::logic_error("Refiner::removeElements couldn't remove element");

            }
        }

      if (m_doProgress && m_eMesh.get_rank() == 0 && progress_meter_num_total)
        {
          std::ostringstream oss; oss << "Delete Elements " << " pass [" << 100.0*((double)irank)/((double)m_ranks.size()) << " %]" << " cpu: " << m_eMesh.cpu_time() << " [sec]";
          ProgressMeterData pd(ProgressMeterData::FINI, 0.0, oss.str());
          notifyObservers(&pd);
        }

      //std::cout << "tmp Refiner::removeElements pass2 size = " << elements_to_be_destroyed_pass2.size() << std::endl;
      for (elements_to_be_destroyed_type::iterator itbd = elements_to_be_destroyed_pass2.begin();
           itbd != elements_to_be_destroyed_pass2.end();  ++itbd)
        {
          stk::mesh::Entity element_p = *itbd;

          if (1)
            {
              const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(stk::topology::ELEMENT_RANK + 1u);
              stk::mesh::EntityRank rank= m_eMesh.entity_rank(element_p);
              for (stk::mesh::EntityRank higher_rank = static_cast<stk::mesh::EntityRank>(rank+1); higher_rank <= FAMILY_TREE_RANK; ++higher_rank)
                {
                  while (true)
                    {
                      percept::MyPairIterRelation rels (m_eMesh, element_p, higher_rank);
                      if (!rels.size())
                        break;
                      stk::mesh::Entity to_rel = rels[0].entity();
                      stk::mesh::RelationIdentifier to_id = rels[0].relation_ordinal();

                      bool del = m_eMesh.get_bulk_data()->destroy_relation( to_rel, element_p, to_id);
                      if (!del)
                        throw std::runtime_error("removeOldElements:: destroy_relation failed 1");
                    }
                }
            }

          if ( ! m_eMesh.get_bulk_data()->destroy_entity( element_p ) )
            {
              shards::CellTopology cell_topo(m_eMesh.get_cell_topology(element_p));
              std::cout << "tmp Refiner::removeElements couldn't remove element in pass2,...\n tmp destroy_entity returned false: cell= " << cell_topo.getName() << std::endl;
              const percept::MyPairIterRelation elem_relations (m_eMesh, element_p, static_cast<stk::mesh::EntityRank>(m_eMesh.entity_rank(element_p)+1));
              std::cout << "tmp elem_relations[rank+1].size() = " << elem_relations.size() << std::endl;
              const percept::MyPairIterRelation elem_relations_elem (m_eMesh, element_p, m_eMesh.element_rank());
              std::cout << "tmp elem_relations[3].size() = " << elem_relations_elem.size() << std::endl;

              throw std::logic_error("Refiner::removeElements couldn't remove element, destroy_entity returned false.");
            }
        }
    }

    /// fix names of surfaces (changing for example surface_hex8_quad4 to surface_tet4_tri3)
    void Refiner::
    fixSurfaceAndEdgeSetNames(stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern)
    {
      EXCEPTWATCH;
      stk::mesh::PartVector toParts = breakPattern->getToParts();

      bool debug = true;

      //std::cout << "toParts.size()= " << toParts.size() << " typeid= " << typeid(*breakPattern).name()  << std::endl;

      for (unsigned i_part = 0; i_part < toParts.size(); i_part++)
        {
          //const std::string & partName = toParts[i_part]->name();
          std::string * toPartName_p = const_cast<std::string *> (&toParts[i_part]->name());

          std::string toPartName = toParts[i_part]->name();
          if ( toPartName.find("surface_", 0) == std::string::npos)
            {
              if (debug) std::cout << "tmp fixSurfaceAndEdgeSetNames:: skipping toPartName= " << toPartName << " typeid= " << typeid(*breakPattern).name()  << std::endl;
              continue;
            }

          std::string newToPartName = toPartName;

          StringStringMap::iterator map_it;
          StringStringMap str_map =  breakPattern->fixSurfaceAndEdgeSetNamesMap();
          if (debug) std::cout << "tmp fixSurfaceAndEdgeSetNamesMap:: str_map.size()= " << str_map.size()
            //<< " " << breakPattern->getFromTopoPartName() << "__" << breakPattern->getToTopoPartName()
                           << " typeid= " << typeid(*breakPattern).name()
                           << std::endl;

          for (map_it = str_map.begin(); map_it != str_map.end(); map_it++)
            {
              std::string from_str = map_it->first;
              std::string to_str = map_it->second;
              Util::replace(newToPartName, from_str, to_str);
              if (debug)
                std::cout << "tmp fixSurfaceAndEdgeSetNamesMap: old= " << toPartName << " new= " << newToPartName << std::endl;
            }

          *toPartName_p = newToPartName;

          if (debug)
            std::cout << "tmp fixSurfaceAndEdgeSetNamesMap:: P[" << m_eMesh.get_rank() << "] new part name= " << toParts[i_part]->name()
                      << " old part name = " << toPartName
                      << std::endl;
        }
    }

    // FIXME this is a hack to rename parts
    /// Renames as follows:
    ///   originalPartName -> originalPartName_uo_1000    The original part holds the elements to be converted, and is renamed to be the "old" part
    ///   originalPartName_urpconv -> originalPartName    The new part has the same name as the original part with urpconv appended, which
    ///                                                      is then changed back to the original part name
    ///   fromPartName -> fromPartName+"_uo_1000"
    ///   toPartName_urpconv -> toPartName
    ///
    /// So, after the renaming, the original part name holds the new elements, and the original elements are
    ///   in the part with the original name appended with _uo_1000.  These parts are ignored on subsequent input.
    ///
#define DEBUG_RENAME_NEW_PARTS 0
#define DEBUG_RENAME_NEW_PARTS_1 0
    static std::string strip_hashes(std::string in, stk::mesh::PartVector& parts, const std::string& sep, bool check=true)
    {
      std::string out=in;
      size_t pos=0;
      pos = in.find(sep);
      if (pos != std::string::npos)
        {
          std::string o1 = in.substr(pos+1);
          size_t pos2 = o1.find(sep);
          out = in.substr(0,pos)+o1.substr(pos2+1);
        }

      if (DEBUG_RENAME_NEW_PARTS_1) std::cout << "tmp srk in= " << std::setw(50) << in << " out= " << std::setw(50) << out << std::endl;
      if (check && out != in)
        {
          for (unsigned ii=0; ii < parts.size(); ++ii)
            {
              if (parts[ii]->name() == out)
                {
                  std::cout << "bad part[" << ii << "] = " << parts[ii]->name() << std::endl;
                  for (unsigned jj=0; jj < parts.size(); ++jj)
                    {
                      std::cout << "part[" << jj << "] = " << parts[jj]->name() << std::endl;
                    }
                  throw std::runtime_error("bad name change");
                }
            }
        }
      return out;
    }

    std::string Refiner::get_parent_element_topology(const std::string& surfaceName)
    {
      // If the sideset has a "canonical" name as in "surface_{id}",
      // Then the sideblock name will be of the form:
      //  * "surface_eltopo_sidetopo_id" or
      //  * "surface_block_id_sidetopo_id"
      // If the sideset does *not* have a canonical name, then
      // the sideblock name will be of the form:
      //  * "{sideset_name}_eltopo_sidetopo" or
      //  * "{sideset_name}_block_id_sidetopo"
      // Generated mesh will create sidesets of the form
      //  * "surface_id_sidetopo

      const stk::mesh::BulkData& bulk = *m_eMesh.get_bulk_data();
      const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
      std::vector<std::string> tokens;
      stk::util::tokenize(surfaceName, "_", tokens);

      size_t tokenSize = tokens.size();

      std::string parent_element_topology;

      if(tokenSize >= 4) {
        parent_element_topology = tokens[1];
      } else if(tokenSize == 3) {
        const bool allDigits = tokens[1].find_first_not_of("0123456789") == std::string::npos;

        std::string parentSurfaceName;

        if (allDigits) {
          // Generated mesh format
          parentSurfaceName = tokens[0] + "_" + tokens[1];
          stk::mesh::Part* parentSurface = meta.get_part(parentSurfaceName);

          if(nullptr != parentSurface) {
            std::vector<const stk::mesh::Part*> touchingBlocks = meta.get_blocks_touching_surface(parentSurface);
            if(touchingBlocks.size() == 1) {
              convert_stk_topology_to_ioss_name(touchingBlocks[0]->topology(), parent_element_topology);
            }
          }
        } else {
          // non-canonical format
          parent_element_topology = tokens[1];
        }
      }

      return parent_element_topology;
    }

    void Refiner::
    renameNewParts(stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern)
    {
      EXCEPTWATCH;

      stk::mesh::PartVector toParts = breakPattern->getToParts();
      stk::mesh::PartVector fromParts = breakPattern->getFromParts();
      bool do_strip_hashes = breakPattern->m_do_strip_hashes;
      bool do_strip_hashes_from = false;

      stk::mesh::MetaData* meta = m_eMesh.get_fem_meta_data();
      stk::mesh::PartVector all_parts = meta->get_parts();

      if (DEBUG_RENAME_NEW_PARTS)
        {
          std::cout << "\n\ntmp srk getFromTopoPartName()= " << breakPattern->getFromTopoPartName() << " getToTopoPartName()= " << breakPattern->getToTopoPartName() << std::endl;
          std::cout << "fromParts.size() = " << fromParts.size() << " toParts.size() = " << toParts.size() << std::endl;
          for (unsigned i_part = 0; i_part < toParts.size(); i_part++)
            {
              std::cout << "tmp toParts[i_part]->name() = " << toParts[i_part]->name() << std::endl;
            }
          for (unsigned i_part = 0; i_part < fromParts.size(); i_part++)
            {
              std::cout << " fromParts[i_part]->name() = " << fromParts[i_part]->name()  << std::endl;
            }
        }

      VERIFY_OP_ON(fromParts.size(), ==, toParts.size(), "sizes wrong");

      if (fromParts.size() == toParts.size())
        {
          for (unsigned i_part = 0; i_part < toParts.size(); i_part++)
            {
              if (DEBUG_RENAME_NEW_PARTS) std::cout << "tmp before: fromPartName= " << fromParts[i_part]->name()
                                                    << " toPartName= " << toParts[i_part]->name() << std::endl;

              std::string toPartName = toParts[i_part]->name();
              if (toParts[i_part]->name() == fromParts[i_part]->name())
                {
                  continue;
                }
              std::string fromPartName = toPartName;
              int len = fromPartName.length();
              int clen = breakPattern->getAppendConvertString().length();
              fromPartName.erase(len - clen, clen);
              //mesh::Part *fromPart = m_eMesh.get_non_const_part(fromPartName);
              stk::mesh::Part *fromPart = fromParts[i_part];
              VERIFY_OP_ON(fromPart, !=, 0, std::string("Refiner::renameNewParts null fromPart found, fromPart= ")+fromPartName);

              if (1)
                {
                  std::string newToPartName = fromPartName;
                  std::string newFromPartName = fromPartName + breakPattern->getAppendOriginalString();
                  if (do_strip_hashes) newToPartName = strip_hashes(newToPartName, all_parts, breakPattern->getConvertSeparatorString(), false);
                  if (do_strip_hashes_from) newFromPartName = strip_hashes(newFromPartName, all_parts, breakPattern->getConvertSeparatorString(), false);

                  meta->delete_part_alias_case_insensitive(*fromParts[i_part], newToPartName);

                  meta->add_part_alias(*toParts[i_part], newToPartName);
                  meta->add_part_alias(*fromParts[i_part], newFromPartName);

                  if (DEBUG_RENAME_NEW_PARTS) {
                    std::cout << "tmp renameNewParts: to alias for " << toParts[i_part]->name() << " = " << newToPartName
                              << " parent topo = " << get_parent_element_topology(toParts[i_part]->name()) << std::endl;
                    std::cout << "tmp renameNewParts: from alias for " << fromParts[i_part]->name() << " = " << newFromPartName
                              << " parent topo = " << get_parent_element_topology(fromParts[i_part]->name()) << std::endl;
                  }

                  StringStringMap::iterator map_it;
                  StringStringMap str_map =  breakPattern->fixSurfaceAndEdgeSetNamesMap();
                  if(rank == meta->side_rank()) {
                    if (DEBUG_RENAME_NEW_PARTS) std::cout << "tmp renameNewParts:: str_map.size()= " << str_map.size()
                                       << " typeid= " << typeid(*breakPattern).name()
                                       << std::endl;

                    for (map_it = str_map.begin(); map_it != str_map.end(); map_it++)
                    {
                      std::string from_str = map_it->first;
                      std::string to_str = map_it->second;
                      Util::replace(newToPartName, from_str, to_str);
                      Util::replace(newFromPartName, from_str, to_str);
                      if (DEBUG_RENAME_NEW_PARTS) {
                        std::cout << "tmp renameNewParts: old toPartNane= " << toPartName << " new toPartName= " << newToPartName << std::endl;
                        std::cout << "tmp renameNewParts: old fromPartNane= " << fromPartName << " new fromPartName= " << newFromPartName << std::endl;
                      }
                    }

                    // This is to prevent the collapse of refined subset parts onto the same name
                    newToPartName += ("." + get_parent_element_topology(toParts[i_part]->name()));
                    newFromPartName += ("." + get_parent_element_topology(fromParts[i_part]->name()));
                  }

                  stk::io::set_alternate_part_name(*toParts[i_part], newToPartName);
                  stk::io::set_alternate_part_name(*fromParts[i_part], newFromPartName);

                  if (DEBUG_RENAME_NEW_PARTS) {
                    std::cout << "tmp  after: fromPartName= " << fromParts[i_part]->name() << " (" << newFromPartName << ") "
                              <<              " toPartName= " << toParts[i_part]->name()  << " (" << newToPartName << ") " << std::endl;
                    std::cout << "tmp P[" << m_eMesh.get_rank() << "] fromPartName: " << fromPartName << " part= " << toParts[i_part]->name()
                              << " old part name = " << fromPart->name()
                              << std::endl;
                  }
                }
            }
        }
    }


    RefinementInfo&
    Refiner::
    getRefinementInfo()
    {
      return m_refinementInfo;
    }

    void Refiner::
    add_children_to_parts()
    {
      static std::vector<stk::mesh::Part*> add_parts(1, static_cast<stk::mesh::Part*>(0));
      static std::vector<stk::mesh::Part*> remove_parts;

      std::vector<stk::mesh::Entity> children;
      const std::vector< stk::mesh::Part * > & parts = m_eMesh.get_fem_meta_data()->get_parts();
      stk::mesh::Selector on_locally_owned_part =  ( m_eMesh.get_fem_meta_data()->locally_owned_part() );

      std::vector<std::pair<stk::mesh::Entity,stk::mesh::Part*> > toChange;

      unsigned nparts = parts.size();
      for (unsigned ipart=0; ipart < nparts; ipart++)
        {
          stk::mesh::Part& part = *parts[ipart];

          if (stk::mesh::is_auto_declared_part(part))
            continue;

          bool percept_auto_part = part.attribute<percept::AutoPart>() != 0;
          if (percept_auto_part)
            continue;

          const stk::mesh::EntityRank part_rank = part.primary_entity_rank();

          if (part_rank == stk::topology::NODE_RANK || part_rank == stk::topology::INVALID_RANK)
            continue;

          stk::mesh::Selector selector(part);
          //std::cout << "part= " << part.name() << std::endl;
          const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( part_rank );

          for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;

              if (selector(bucket))
                {
                  const unsigned num_entity_in_bucket = bucket.size();
                  for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                    {
                      stk::mesh::Entity element = bucket[ientity];
                      if (m_eMesh.hasFamilyTree(element))
                        {
                          m_eMesh.getChildren(element, children, false);
                          for (unsigned ich = 0; ich < children.size(); ich++)
                            {
                              if (on_locally_owned_part(m_eMesh.bucket(children[ich])) && !selector(m_eMesh.bucket(children[ich]))
                                  && m_eMesh.bucket(children[ich]).topology() == m_eMesh.bucket(element).topology())
                                {
                                  toChange.push_back(std::make_pair(children[ich],&part));
                                }
                            }
                        }
                    }
                }
            }
        }
      for (unsigned ii=0; ii < toChange.size(); ++ii)
        {
          stk::mesh::Entity element = toChange[ii].first;
          add_parts[0] = toChange[ii].second;
          m_eMesh.get_bulk_data()->change_entity_parts( element, add_parts, remove_parts );
        }
    }

#if 1
    void Refiner::set_active_part()
    {
      set_active_part(m_eMesh);
    }

    // static
    void Refiner::set_active_part(PerceptMesh& eMesh)
    {
      // deal with parts
      stk::mesh::EntityRank part_ranks[] = {eMesh.element_rank(), eMesh.side_rank()};
      for (unsigned irank=0; irank < 2; irank++)
        {
          std::string active_part_name = "refine_active_elements_part_"+toString(part_ranks[irank]);
          std::string inactive_part_name = "refine_inactive_elements_part_"+toString(part_ranks[irank]);
          stk::mesh::Part* child_elements_part = eMesh.get_non_const_part(active_part_name);
          stk::mesh::Part* parent_elements_part = eMesh.get_non_const_part(inactive_part_name);
          stk::mesh::Selector in_child_part(*child_elements_part);
          stk::mesh::Selector in_parent_part(*parent_elements_part);

          if (child_elements_part && parent_elements_part)
            {
              std::vector<stk::mesh::Part*> child_parts(1, child_elements_part);
              std::vector<stk::mesh::Part*> parent_parts(1, parent_elements_part);
              stk::mesh::Selector on_locally_owned_part =  ( eMesh.get_fem_meta_data()->locally_owned_part() );

              std::vector<stk::mesh::Entity> child_entities;
              std::vector<stk::mesh::Entity> parent_entities;

              const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( part_ranks[irank] );

              for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
                {
                  stk::mesh::Bucket & bucket = **k ;

                  if (on_locally_owned_part(bucket))
                    {
                      bool in_child_part_element = in_child_part(bucket);
                      bool in_parent_part_element = in_parent_part(bucket);
                      const unsigned num_entity_in_bucket = bucket.size();
                      for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                        {
                          stk::mesh::Entity element = bucket[ientity];
                          if (eMesh.hasFamilyTree(element) && eMesh.isParentElement(element, true))
                            {
                              if (in_child_part_element || !in_parent_part_element)
                              {
                                parent_entities.push_back(element);
                              }
                            }
                          else
                            {
                              if (!in_child_part_element || in_parent_part_element)
                              {
                                child_entities.push_back(element);
                              }
                            }
                        }
                    }
                }

              //std::cout << "tmp Refiner::set_active_part: child_entities= " << child_entities.size() << " parent_entities= " << parent_entities.size() << std::endl;
              eMesh.get_bulk_data()->change_entity_parts(child_entities, child_parts, parent_parts);
              eMesh.get_bulk_data()->change_entity_parts(parent_entities,  parent_parts, child_parts);
            }
        }

    }

#else
    void Refiner::set_active_part()
    {
      // deal with parts
      stk::mesh::EntityRank part_ranks[] = {m_eMesh.element_rank(), m_eMesh.side_rank()};
      for (unsigned irank=0; irank < 2; irank++)
        {
          std::string active_part_name = "refine_active_elements_part_"+toString(part_ranks[irank]);
          std::string inactive_part_name = "refine_inactive_elements_part_"+toString(part_ranks[irank]);
          stk::mesh::Part* child_elements_part = m_eMesh.get_non_const_part(active_part_name);
          stk::mesh::Part* parent_elements_part = m_eMesh.get_non_const_part(inactive_part_name);
          stk::mesh::Selector in_child_part(*child_elements_part);
          stk::mesh::Selector in_parent_part(*parent_elements_part);

          if (child_elements_part && parent_elements_part)
            {
              std::vector<stk::mesh::Part*> empty;
              std::vector<stk::mesh::Part*> child_parts(1, child_elements_part);
              std::vector<stk::mesh::Part*> parent_parts(1, parent_elements_part);
              stk::mesh::Selector on_locally_owned_part =  ( m_eMesh.get_fem_meta_data()->locally_owned_part() );

              std::vector<stk::mesh::Entity> child_entities_add_to_child;
              std::vector<stk::mesh::Entity> child_entities_remove_from_parent;
              std::vector<stk::mesh::Entity> parent_entities_add_to_parent;
              std::vector<stk::mesh::Entity> parent_entities_remove_from_child;

              const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( part_ranks[irank] );

              for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
                {
                  stk::mesh::Bucket & bucket = **k ;

                  if (on_locally_owned_part(bucket))
                    {
                      bool in_child_part_element = in_child_part(bucket);
                      bool in_parent_part_element = in_parent_part(bucket);
                      const unsigned num_entity_in_bucket = bucket.size();
                      for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                        {
                          stk::mesh::Entity element = bucket[ientity];
                          if (m_eMesh.hasFamilyTree(element) && m_eMesh.isParentElement(element, true))
                            {
                              if (in_child_part_element)
                                {
                                  parent_entities_remove_from_child.push_back(element);
                                }
                              if (!in_parent_part_element)
                                {
                                  parent_entities_add_to_parent.push_back(element);
                                }
                            }
                          else
                            {
                              if (!in_child_part_element)
                              {
                                child_entities_add_to_child.push_back(element);
                              }
                              if (in_parent_part_element)
                              {
                                child_entities_remove_from_parent.push_back(element);
                              }
                            }
                        }
                    }
                }

              //std::cout << "tmp Refiner::set_active_part: child_entities= " << child_entities.size() << " parent_entities= " << parent_entities.size() << std::endl;
              for (unsigned iv=0; iv < child_entities_add_to_child.size(); iv++)
                {
                  m_eMesh.get_bulk_data()->change_entity_parts( child_entities_add_to_child[iv],   child_parts, empty);
                }
              for (unsigned iv=0; iv < parent_entities_add_to_parent.size(); iv++)
                {
                  m_eMesh.get_bulk_data()->change_entity_parts( parent_entities_add_to_parent[iv],  parent_parts, empty);
                }
              for (unsigned iv=0; iv < child_entities_remove_from_parent.size(); iv++)
                {
                  m_eMesh.get_bulk_data()->change_entity_parts( child_entities_remove_from_parent[iv],   empty, parent_parts );
                }
              for (unsigned iv=0; iv < parent_entities_remove_from_child.size(); iv++)
                {
                  m_eMesh.get_bulk_data()->change_entity_parts( parent_entities_remove_from_child[iv], empty, child_parts );
                }
            }
        }

    }
#endif
    //    ========================================================================================================================
    //    ========================================================================================================================
    //    ========================================================================================================================

    void Refiner::initializeDB(bool use_rebuild_node_registry)
    {
      if (use_rebuild_node_registry)
        {
          RefinerUtil::rebuild_node_registry(m_eMesh, *m_nodeRegistry, true);
          return;
        }
      this->special_processing("refiner_pre_initializeDB");
      bool debug=false;
      m_nodeRegistry->initialize();
      m_nodeRegistry->init_comm_all();

      for (unsigned irank = 0; irank < m_ranks.size(); irank++)
        {
          unsigned elementType = m_breakPattern[irank]->getFromTypeKey();
          {
            EXCEPTWATCH;

            ElementRankTypeInfo& e_info = m_elementRankTypeInfo[irank];
            VERIFY_OP_ON(m_ranks[irank], ==, e_info.first,"er1");
            VERIFY_OP_ON(elementType, ==, e_info.second,"er2");

            std::vector<NeededEntityType> needed_entity_ranks;
            m_breakPattern[irank]->fillNeededEntities(needed_entity_ranks);

            bool doAllElements = true;

            unsigned num_elem_not_ghost_0_incr = doForAllElements(irank, "initializeEmpty",
                                                                  m_ranks[irank], &NodeRegistry::initializeEmpty,
                                                                  elementType, needed_entity_ranks,
                                                                  doAllElements);

            if (debug) std::cout << "P[" << m_eMesh.get_rank() << "] tmp srk1 irank= " << irank << " ranks[irank]= " << m_ranks[irank]
                             << " nodeRegistry size= " << m_nodeRegistry->getMap().size()
                             << " num_elem_not_ghost_0_incr= " << num_elem_not_ghost_0_incr
                             << std::endl;
          }
        }

      if (use_rebuild_node_registry)
        {
          RefinerUtil::rebuild_node_registry(m_eMesh, *m_nodeRegistry, false);
        }
      else {

      SubDimCellToDataMap& cell_2_data_map = m_nodeRegistry->getMap();
      if (debug) std::cout << "tmp srk1 cell_2_data_map size= " << cell_2_data_map.size() << std::endl;

      for (SubDimCellToDataMap::iterator cell_iter = cell_2_data_map.begin(); cell_iter != cell_2_data_map.end(); ++cell_iter)
        {
          const SubDimCell_SDCEntityType& subDimEntity = (*cell_iter).first;
          SubDimCellData& nodeId_elementOwnderId = (*cell_iter).second;
          stk::mesh::EntityId owning_elementId = std::get<SDC_DATA_OWNING_ELEMENT_KEY>(nodeId_elementOwnderId).id();
          NodeIdsOnSubDimEntityType& nodeIds_onSE = std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnderId);
          stk::mesh::EntityRank owning_elementRank = std::get<SDC_DATA_OWNING_ELEMENT_KEY>(nodeId_elementOwnderId).rank();

          if (0 && debug) std::cout << "P[" << m_eMesh.get_rank() << "] tmp srk1 subDimEntity.size= " << subDimEntity.size() << " owning_elementId= " << owning_elementId << " nodeIds_onSE.size= " << nodeIds_onSE.size() << std::endl;
          // if (subDimEntity.size() == 1)
          //   continue;

          if (1)
            {
              stk::mesh::Entity owning_element = m_eMesh.get_bulk_data()->get_entity(owning_elementRank, owning_elementId);
              VERIFY_OP_ON(m_eMesh.is_valid(owning_element), ==, true, "hmmm");
              //double ela = m_eMesh.edge_length_ave(owning_element);
              double centroid[3] = {0,0,0};
              int nDim = m_eMesh.get_spatial_dim();
              if (subDimEntity.size() == 1)
                {
                  computeCentroid(subDimEntity[0], centroid, *(m_eMesh.get_coordinates_field()));
                }
              else
                {
                  for (unsigned ii=0; ii < subDimEntity.size(); ++ii)
                    {
                      stk::mesh::Entity node = subDimEntity[ii];
                      double *coords = m_eMesh.field_data_inlined(m_eMesh.get_coordinates_field(), node);
                      for (int jj=0; jj < nDim; ++jj)
                        {
                          centroid[jj] += coords[jj]/double(subDimEntity.size());
                        }
                    }
                }

              std::vector<stk::mesh::Entity> children;
              double min_dist = std::numeric_limits<double>::max();
              stk::mesh::Entity min_node = stk::mesh::Entity();
              const percept::MyPairIterRelation parent_nodes(m_eMesh, owning_element, m_eMesh.node_rank());

              m_eMesh.getChildren(owning_element, children, false);
              if (children.size() == 0)
                continue;

              double emin = std::numeric_limits<double>::max();
              for (unsigned ich = 0; ich < children.size(); ich++)
                {
                  double min_edge_len = 0, max_edge_len = 0;
                  m_eMesh.edge_length_ave(owning_element, m_eMesh.get_coordinates_field(), &min_edge_len, &max_edge_len);
                  if (min_edge_len < emin)
                    emin = min_edge_len;

                  const percept::MyPairIterRelation child_nodes(m_eMesh, children[ich], m_eMesh.node_rank());
                  for (unsigned jj=0; jj < child_nodes.size(); ++jj)
                    {
                      stk::mesh::Entity cnode = child_nodes[jj].entity();
                      bool fnd=false;
                      for (unsigned kk=0; kk < parent_nodes.size(); ++kk)
                        {
                          if (cnode == parent_nodes[kk].entity())
                            {
                              fnd=true;
                              break;
                            }
                        }
                      if (fnd) continue;
                      double *c_coords = m_eMesh.field_data_inlined(m_eMesh.get_coordinates_field(), cnode);
                      double dist=0;
                      for (int k=0; k < nDim; ++k)
                        {
                          dist += (centroid[k] - c_coords[k])*(centroid[k] - c_coords[k]);
                        }
                      dist = std::sqrt(dist);
                      if (dist < min_dist)
                        {
                          min_dist = dist;
                          min_node = cnode;
                        }
                    }
                }
              bool found = min_dist < 1.e-3*emin;

              // FIXME - what about quad faces?
              if (0 && !found)
                {
                  double ela_local = 0.0;
                  double nela = 0.0;
                  for (unsigned ii=0; ii < subDimEntity.size()-1; ++ii)
                    {
                      double *coords_ii = m_eMesh.field_data_inlined(m_eMesh.get_coordinates_field(), subDimEntity[ii]);
                      for (unsigned jj=ii+1; jj < subDimEntity.size(); ++jj)
                        {
                          double *coords_jj = m_eMesh.field_data_inlined(m_eMesh.get_coordinates_field(), subDimEntity[jj]);
                          nela += 1.0;
                          double dist=0.0;
                          for (int k=0; k < nDim; ++k)
                            {
                              dist += (coords_ii[k] - coords_jj[k])*(coords_ii[k] - coords_jj[k]);
                            }
                          dist = std::sqrt(dist);
                          ela_local += dist;
                        }
                    }
                  ela_local /= nela;
                  found = min_dist < 1.e-4*ela_local;
                }
              if (found)
                {
                  if (0)
                    {
                      std::cout << "P[" << m_eMesh.get_rank() << "] found node= " << m_eMesh.identifier(min_node)
                                << " " << m_eMesh.print_entity_compact(min_node)
                                << std::endl;
                    }

                  nodeIds_onSE.resize(1);
                  nodeIds_onSE.m_entity_id_vector.resize(1);
                  nodeIds_onSE[0] = min_node;
                  nodeIds_onSE.m_entity_id_vector[0] = m_eMesh.identifier(min_node);
                  nodeIds_onSE.m_mark = NodeRegistry::NR_MARK;
                }
              else
                {
                  nodeIds_onSE.resize(0);
                  nodeIds_onSE.m_entity_id_vector.resize(0);
                }
            }
        }
      }

      if (1)
        {
          mod_begin();
          this->removeDanglingNodes();
          mod_end(0,"RefinerInitDB");
        }

      this->special_processing("refiner_post_initializeDB");
    }

    void Refiner::removeFromNewNodesPart()
    {
      stk::mesh::Part* new_nodes_part = m_eMesh.get_non_const_part("refine_new_nodes_part");
      VERIFY_OP_ON(new_nodes_part, !=, 0, "new_nodes_part");

      std::vector<stk::mesh::Part*> remove_parts(1, new_nodes_part);
      std::vector<stk::mesh::Part*> add_parts;
      std::vector<stk::mesh::Entity> node_vec;

      stk::mesh::Selector removePartSelector(*new_nodes_part & m_eMesh.get_fem_meta_data()->locally_owned_part() );
      const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.node_rank() );
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          if (removePartSelector(bucket))
            {
              const unsigned num_entity_in_bucket = bucket.size();
              for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                {
                  stk::mesh::Entity node = bucket[ientity];
                  node_vec.push_back(node);
                  if (m_eMesh.m_new_nodes_field)
                    {
                      NewNodesType::value_type *ndata = stk::mesh::field_data(*m_eMesh.m_new_nodes_field, node);
                      if (ndata)
                        {
                          ndata[0] = static_cast<NewNodesType::value_type>(0);
                        }
                    }
                }
            }
        }
        m_eMesh.get_bulk_data()->change_entity_parts( node_vec, add_parts, remove_parts );
    }

  void Refiner::
  fix_side_sets_2(bool allow_not_found, SetOfEntities *avoid_elems, SetOfEntities *avoid_sides, RefinerSelector *sel, const std::string& msg)
  {
    if (m_avoidFixSideSets)
      return;

    FixSideSets fss(this, m_eMesh, m_excludeParts, m_side_part_map, m_geomFile, m_avoidFixSideSetChecks, sel, m_doProgress);
    fss.fix_side_sets_2(allow_not_found, avoid_elems, avoid_sides, msg);
  }

  void Refiner::
  build_side_set(SetOfEntities& side_set, bool only_roots)
  {
    FixSideSets fss(this, m_eMesh, m_excludeParts, m_side_part_map, m_geomFile, m_avoidFixSideSetChecks, 0, m_doProgress); // FIXME - RefinerSelector?
    fss.build_side_set(side_set, only_roots);
  }

  bool Refiner::bucket_acceptable(stk::mesh::Bucket& bucket, stk::mesh::EntityRank rank)
  {
    stk::mesh::PartVector const& side_parts = bucket.supersets();
    for (unsigned isp=0; isp < side_parts.size(); ++isp)
      {
        stk::mesh::Part& part = *side_parts[isp];
        bool is_auto = stk::mesh::is_auto_declared_part(part);
        const AutoPart *side_auto_part = part.attribute<AutoPart>();
        bool is_percept_auto_part = side_auto_part != 0;
        if (!is_percept_auto_part && !is_auto && part.primary_entity_rank() == rank)
          {
            return true;
          }
      }
    return false;
  }

  void Refiner::mod_begin(stk::diag::Timer *timer)
  {
    stk::diag::Timer& timerRoot_ = (getModBegEndRootTimer() ? *getModBegEndRootTimer() : (timer ? *timer :  rootTimer()));
    stk::diag::Timer timerModBeg_("percept::ModBeg", timerRoot_);
    stk::diag::TimeBlock timerModBegBlock_(timerModBeg_);

    m_eMesh.get_bulk_data()->modification_begin();

  }

  void Refiner::mod_end(stk::diag::Timer *timer, const std::string& msg)
  {
    stk::diag::Timer& timerRoot_ = (getModBegEndRootTimer() ? *getModBegEndRootTimer() : (timer ? *timer :  rootTimer()));
    stk::diag::Timer timerModEndAll_("ModEndAll", timerRoot_);
    stk::diag::Timer timerModEnd_("ModEnd"+msg, timerModEndAll_);
    stk::diag::TimeBlock timerModEndBlock_(timerModEnd_);
    stk::diag::TimeBlock timerModEndBlockAll_(timerModEndAll_);

    //stk::mesh::fixup_ghosted_to_shared_nodes(*m_eMesh.get_bulk_data());
    m_eMesh.get_bulk_data()->modification_end();
  }

  void add_ft_nodes(PerceptMesh& eMesh, stk::mesh::Entity parent_elem, stk::mesh::Entity child_elem)
  {
    const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(eMesh.element_rank() + 1u);
    SetOfEntities to_add(*eMesh.get_bulk_data());

    MyPairIterRelation parent_ft(eMesh, parent_elem, FAMILY_TREE_RANK);
    unsigned parent_elem_ft_level_1 = 0;
    if (parent_ft.size() == 2)
      parent_elem_ft_level_1 = eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_1, parent_elem);

    stk::mesh::Entity family_tree = parent_ft[parent_elem_ft_level_1].entity();

    bool checkInShared = true;

    percept::MyPairIterRelation parent_elem_nodes (eMesh, parent_elem,  stk::topology::NODE_RANK );
    for (unsigned i = 0; i < parent_elem_nodes.size(); i++)
      {
        if (checkInShared && !eMesh.get_bulk_data()->in_shared(parent_elem_nodes[i].entity())) continue;

        bool found = false;
        percept::MyPairIterRelation ft_nodes (eMesh, family_tree,  stk::topology::NODE_RANK );
        for (unsigned j = 0; j < ft_nodes.size(); j++)
          {
            if (ft_nodes[j].entity() == parent_elem_nodes[i].entity())
              {
                found = true;
                break;
              }
          }
        if (!found)
          {
            to_add.insert(parent_elem_nodes[i].entity());
            VERIFY_OP_ON(eMesh.get_bulk_data()->in_index_range(parent_elem_nodes[i].entity()), ==, true, "parent_elem_nodes bad");
            //eMesh.get_bulk_data()->declare_relation(family_tree, parent_elem_nodes[i].entity(), ft_nodes.size());
          }
      }

    percept::MyPairIterRelation child_elem_nodes (eMesh, child_elem,  stk::topology::NODE_RANK );
    if (child_elem_nodes.size() == 0)
      {
        throw std::runtime_error("child_elem has no nodes");
      }
    for (unsigned i = 0; i < child_elem_nodes.size(); i++)
      {
        if (checkInShared && !eMesh.get_bulk_data()->in_shared(child_elem_nodes[i].entity())) continue;

        bool found = false;
        percept::MyPairIterRelation ft_nodes (eMesh, family_tree,  stk::topology::NODE_RANK );
        for (unsigned j = 0; j < ft_nodes.size(); j++)
          {
            if (ft_nodes[j].entity() == child_elem_nodes[i].entity())
              {
                found = true;
                break;
              }
          }
        if (!found)
          {
            to_add.insert(child_elem_nodes[i].entity());
            VERIFY_OP_ON(eMesh.get_bulk_data()->in_index_range(child_elem_nodes[i].entity()), ==, true, "child_elem_nodes bad");
            //eMesh.get_bulk_data()->declare_relation(family_tree, child_elem_nodes[i].entity(), ft_nodes.size());

          }
      }

    // check for second level and subsequent refinement
    if (parent_ft.size() == 2)
      {
        unsigned parent_elem_ft_level_0 = eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, parent_elem);
        stk::mesh::Entity family_tree_level_0 = parent_ft[parent_elem_ft_level_0].entity();

        percept::MyPairIterRelation ft_level_0_nodes (eMesh, family_tree_level_0,  stk::topology::NODE_RANK );
        for (unsigned i = 0; i < ft_level_0_nodes.size(); i++)
          {
            if (checkInShared && !eMesh.get_bulk_data()->in_shared(ft_level_0_nodes[i].entity())) continue;

            bool found = false;
            percept::MyPairIterRelation ft_nodes (eMesh, family_tree,  stk::topology::NODE_RANK );
            for (unsigned j = 0; j < ft_nodes.size(); j++)
              {
                if (ft_nodes[j].entity() == ft_level_0_nodes[i].entity())
                  {
                    found = true;
                    break;
                  }
              }
            if (!found)
              {
                VERIFY_OP_ON(eMesh.get_bulk_data()->in_index_range(ft_level_0_nodes[i].entity()), ==, true, "ft_level_0_nodes bad 0");
                //eMesh.get_bulk_data()->declare_relation(family_tree, ft_level_0_nodes[i].entity(), ft_nodes.size());
                to_add.insert(ft_level_0_nodes[i].entity());
              }
          }
      }

    // add nodes to family_tree
    {
      percept::MyPairIterRelation ft_nodes (eMesh, family_tree,  stk::topology::NODE_RANK );
      unsigned ftns=ft_nodes.size();

      std::vector<stk::mesh::Entity> to_add_vec(to_add.begin(), to_add.end());

      for (unsigned ita=0; ita < to_add_vec.size(); ita++)
        {
          eMesh.get_bulk_data()->declare_relation(family_tree, to_add_vec[ita], ftns+ita);
        }
    }

  }

  void add_ft_nodes(PerceptMesh& eMesh, stk::mesh::Entity elem)
  {
    std::vector<stk::mesh::Entity> children;
    if (eMesh.hasFamilyTree(elem))
      {
        eMesh.getChildren(elem, children, true, false);
        if (children.size() == 0)
          {
            return;
          }
      }

    for (auto child : children)
      {
        add_ft_nodes(eMesh, elem, child);
      }
  }

  void delete_ft_nodes(PerceptMesh& eMesh, stk::mesh::Entity ft)
  {
    std::vector<stk::mesh::Entity> nodes(eMesh.get_bulk_data()->begin_nodes(ft), eMesh.get_bulk_data()->end_nodes(ft));
    std::vector<stk::mesh::ConnectivityOrdinal> nords( eMesh.get_bulk_data()->begin_node_ordinals(ft),  eMesh.get_bulk_data()->end_node_ordinals(ft));

    for (unsigned jj=0; jj < nodes.size(); ++jj)
      {
        bool del = eMesh.get_bulk_data()->destroy_relation( ft, nodes[jj], nords[jj]);

        VERIFY_OP_ON(del, ==, true, "bad del");
      }
  }

  void Refiner::reset_family_tree_to_node_relations()
  {

    const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(m_eMesh.element_rank() + 1u);

    std::vector<stk::mesh::Entity> ftvec;
    stk::mesh::get_selected_entities(m_eMesh.get_fem_meta_data()->locally_owned_part(), m_eMesh.get_bulk_data()->buckets(FAMILY_TREE_RANK), ftvec, false/*don't sort*/);
    for (auto ft : ftvec)
      {
        delete_ft_nodes(m_eMesh, ft);
      }

    for (stk::mesh::EntityRank rank_iter = m_eMesh.side_rank(); rank_iter <= m_eMesh.element_rank(); ++rank_iter)
      {
        SetOfEntities eset(*m_eMesh.get_bulk_data());
        std::vector<stk::mesh::Entity> evec;
        stk::mesh::get_selected_entities(m_eMesh.get_fem_meta_data()->locally_owned_part() , m_eMesh.get_bulk_data()->buckets(rank_iter), evec, false/*don't sort*/);
        for (auto elem : evec)
          {
            //eset.insert(elem);

            if (m_eMesh.numChildren(elem) == 0)
              {
                eset.insert(elem);

                stk::mesh::Entity parent = m_eMesh.getParent(elem, false);
                if (m_eMesh.is_valid(parent))
                  {
                    eset.insert(parent);
                    if (0)
                      {
                        stk::mesh::Entity grand_parent = m_eMesh.getParent(parent, false);
                        if (m_eMesh.is_valid(grand_parent))
                          {
                            eset.insert(grand_parent);
                          }
                      }
                  }
              }
          }
        for (auto elem : eset)
          {
            add_ft_nodes(m_eMesh, elem);
          }
      }
  }

} // namespace percept
