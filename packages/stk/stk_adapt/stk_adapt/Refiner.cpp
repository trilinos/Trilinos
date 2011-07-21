#include <exception>
#include <fstream>
#include <set>
#include <typeinfo>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <stk_adapt/Refiner.hpp>

#include <stk_percept/MeshUtil.hpp>

#if defined( STK_ADAPT_HAS_GEOMETRY )
#include <stk_adapt/geometry/GeometryKernelOpenNURBS.hpp>
#include <stk_adapt/geometry/MeshGeometry.hpp>
#include <stk_adapt/geometry/GeometryFactory.hpp>
#endif

#include <stk_adapt/NodeRegistryDef.hpp>

namespace stk {
  namespace adapt {
    using namespace std;
    using namespace percept;

    Refiner::Refiner(percept::PerceptMesh& eMesh, UniformRefinerPatternBase &  bp, stk::mesh::FieldBase *proc_rank_field) :
      m_eMesh(eMesh), m_breakPattern(),
      m_nodeRegistry(0),
      m_proc_rank_field(proc_rank_field), m_doRemove(true), m_ranks(), m_ignoreSideSets(false),
      m_geomFile(""), m_geomSnap(false),
      m_doQueryOnly(false),
      m_progress_meter_frequency(20),
      m_doProgress(true && (0 == eMesh.getRank()) ),
      m_alwaysInitNodeRegistry(true),
      m_doSmoothGeometry(true)
    {
      bp.setSubPatterns(m_breakPattern, eMesh);
      m_nodeRegistry = new NodeRegistry (m_eMesh);
      m_nodeRegistry->initialize();
      m_nodeRegistry->init_comm_all();
    }

    Refiner::~Refiner()
    {
      if (m_nodeRegistry)
        delete m_nodeRegistry;
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
    removeFromOldPart(stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern)
    {
      EXCEPTWATCH;

      std::string oldPartName = breakPattern->getOldElementsPartName()+toString(rank);
      mesh::Part *oldPart = m_eMesh.getFEM_meta_data()->get_part(oldPartName);
      //std::cout << "tmp removeFromOldPart:: oldPartName= " << oldPartName << std::endl;
      if (!oldPart)
        {
          std::cout << "oldPartName= " << oldPartName << std::endl;
          throw std::runtime_error("oldpart is null");
        }

      mesh::PartVector remove_parts(1, oldPart);
      mesh::PartVector add_parts;
      mesh::Selector on_locally_owned_part =  ( m_eMesh.getFEM_meta_data()->locally_owned_part() );

      std::vector<stk::mesh::Entity*> elems;
      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.getBulkData()->buckets( rank );

      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (on_locally_owned_part(**k)) // && fromPartsSelector(**k) )
            {
              stk::mesh::Bucket & bucket = **k ;

              const unsigned num_elements_in_bucket = bucket.size();

              for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
                {
                  stk::mesh::Entity& element = bucket[i_element];
                  elems.push_back(&element);
                }
            }
        }

      for (unsigned ielem=0; ielem < elems.size(); ielem++)
        {
          //std::cout << "tmp removeFromOldPart:: element = " << *elems[ielem] << std::endl;
          m_eMesh.getBulkData()->change_entity_parts( *elems[ielem], add_parts, remove_parts );
        }
    }

    void Refiner::
    addOldElementsToPart(stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern, unsigned *elementType)
    {
      EXCEPTWATCH;
      //m_eMesh.getBulkData()->modification_begin();
      std::string oldPartName = breakPattern->getOldElementsPartName()+toString(rank);
      mesh::Part *oldPart = m_eMesh.getFEM_meta_data()->get_part(oldPartName);
#define DEBUG_REMOVE_OLD_PARTS 0

      if (DEBUG_REMOVE_OLD_PARTS) std::cout << "tmp addOldElementsToPart:: oldPartName= " << oldPartName << std::endl;

      if (!oldPart)
        {
          std::cout << "oldPartName= " << oldPartName << std::endl;
          throw std::runtime_error("oldpart is null");
        }

      mesh::PartVector add_parts(1, oldPart);
      mesh::PartVector remove_parts;
      mesh::Selector on_locally_owned_part =  ( m_eMesh.getFEM_meta_data()->locally_owned_part() );

      // The list of Parts that this break pattern will refine.  Only remove elements belonging to these parts.
      mesh::Selector fromPartsSelector = mesh::selectUnion( breakPattern->getFromParts() );

      std::vector<stk::mesh::Entity*> elems;
      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.getBulkData()->buckets( rank );

      unsigned nele=0;
      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (on_locally_owned_part(**k) && fromPartsSelector(**k) )
            //if ( on_locally_owned_part(**k) )
            {
              stk::mesh::Bucket & bucket = **k ;
              const CellTopologyData * const bucket_cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(bucket);
              shards::CellTopology topo(bucket_cell_topo_data);

              const unsigned num_elements_in_bucket = bucket.size();

              for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
                {
                  EXCEPTWATCH;
                  stk::mesh::Entity& element = bucket[i_element];
                  if (&element == 0)
                    {
                      std::cout << "element = 0" << std::endl;
                      throw std::runtime_error("element = 0");
                      //exit(1);
                    }

                  if (elementType && (topo.getKey() != *elementType))
                    {
                    }
                  else
                    {
                      elems.push_back(&element);
                      ++nele;
                      //std::cout << "tmp adding to oldParts = " << element << std::endl;
                    }
                }
            }
        }


      for (unsigned ielem=0; ielem < elems.size(); ielem++)
        {
          //std::cout << "tmp addOldElementsToPart element = " << *elems[ielem] << std::endl;
          m_eMesh.getBulkData()->change_entity_parts( *elems[ielem], add_parts, remove_parts );
        }

      //m_eMesh.getBulkData()->modification_end();
    }

    void Refiner::
    trace_print(std::string msg)
    {
      size_t heap_in_Mb = 0;
      size_t memory_in_Mb = 0;
      double cpu = 0.0;

      if (TRACE_STAGE_PRINT)
        {
          memory_in_Mb = Util::memory(heap_in_Mb);
          memory_in_Mb = memory_in_Mb / (1024*1024);
          heap_in_Mb = heap_in_Mb / (1024*1024);
        }

      cpu = Util::cpu_time();
      std::cout
        << msg
        << " cpu_time= " << cpu/(60.) << " [min] "
        << " mem= " << memory_in_Mb << " [Mb] "
        //<< " heap= " << heap_in_Mb << " [Mb] "
        << std::endl;

    }


    struct myVec
    {
      double *data;
      int len;
      int res;
    };

    static void doPrintSizes()
    {
      if (0)
        {
          std::cout
            << "sizeof(myVec) = " << sizeof(myVec) << " "
            << "sizeof(Relation) = " << sizeof(stk::mesh::Relation) << " "
            << "sizeof(Entity) = " << sizeof(stk::mesh::Entity) << " "
            << "sizeof(EntityImpl) = " << sizeof(stk::mesh::impl::EntityImpl) << " "
            << "\nsizeof(EntityKey) = " << sizeof(stk::mesh::EntityKey) << " "
            << "\nsizeof(RelationVector) = " << sizeof(stk::mesh::RelationVector) << " "
            << "\nsizeof(EntityCommInfoVector) = " << sizeof(stk::mesh::EntityCommInfoVector) << " "
            << "\nsizeof(Bucket *) = " << sizeof(stk::mesh::Bucket *) << " "
            << "\nsizeof(unsigned) = " << sizeof(unsigned) << " "
            << "\nsizeof(size_t) = " << sizeof(size_t) << " "
            << "\nsizeof(EntityModificationLog) = " << sizeof(stk::mesh::EntityModificationLog) << std::endl;

        }
    }

    void Refiner::
    checkBreakPatternValidityAndBuildRanks(std::vector<stk::mesh::EntityRank>& ranks)
    {
      ranks.resize(0);

      m_eMesh.getBulkData()->modification_begin();
      for (unsigned ibp = 0; ibp < m_breakPattern.size(); ibp++)
        {
          if (m_breakPattern[ibp])
            {
              stk::mesh::EntityRank irank = m_breakPattern[ibp]->getPrimaryEntityRank();
              stk::mesh::EntityRank irank_prev = stk::percept::EntityRankEnd;
              if (ibp > 0) irank_prev = m_breakPattern[ibp-1]->getPrimaryEntityRank();
              if (irank > irank_prev)
                {
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
      m_eMesh.getBulkData()->modification_end();

    }

    void Refiner::
    applyNodeRegistryFunctionForSubEntities(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity& element, vector<NeededEntityType>& needed_entity_ranks)
    {
      const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);

      CellTopology cell_topo(cell_topo_data);
      const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);

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
          else if (needed_entity_rank == m_eMesh.element_rank())
            {
              numSubDimNeededEntities = 1;
            }

          for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
            {
              /// note: at this level of granularity we can do single edge refinement, hanging nodes, etc.
              //SubDimCell_SDSEntityType subDimEntity;
              //getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);

              (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd);

            } // iSubDimOrd
        } // ineed_ent
    }

    void Refiner::
    doBreak()
    {
      EXCEPTWATCH;

      /**/                                                TRACE_PRINT( "Refiner:doBreak start...");

      m_nodeRegistry->dumpDB("start of doBreak");

      if (0) doPrintSizes();

      //NodeRegistry nr (m_eMesh);
      //m_nodeRegistry = &nr;

//       if (m_nodeRegistry)
//         {
//           delete m_nodeRegistry;
//         }

//        if (!m_nodeRegistry)
//          {
//            m_nodeRegistry = new NodeRegistry (m_eMesh);
//          }

      CommDataType buffer_entry;

      stk::mesh::BulkData& bulkData = *m_eMesh.getBulkData();
      static SubDimCellData empty_SubDimCellData;

      // color elements
#if 0
      struct EntityExcluder
      {
        virtual bool exclude(stk::mesh::Entity& element) = 0;
      };
#endif

      std::vector<stk::mesh::EntityRank>& ranks = m_ranks;
      // check logic of break pattern setup and also build ranks used vector
      checkBreakPatternValidityAndBuildRanks(ranks);

      ///////////////////////////////////////////////////////////
      ///// Get info on refinements that will be done
      ///////////////////////////////////////////////////////////
      if (1)
        {

          m_refinementInfoByType.resize(ranks.size());

          stk::mesh::PartVector fromPartsAll;

          for (unsigned irank = 0; irank < ranks.size(); irank++)
            {
              stk::mesh::PartVector * fromParts = &(m_breakPattern[irank]->getFromParts());
              if (fromParts)
                {
                  for (unsigned ipart = 0; ipart < fromParts->size(); ipart++)
                    {
                      fromPartsAll.push_back((*fromParts)[ipart]);
                    }
                }
            }

          for (unsigned irank = 0; irank < ranks.size(); irank++)
            {
              EXCEPTWATCH;
              unsigned elementType = m_breakPattern[irank]->getFromTypeKey();
              shards::CellTopology cell_topo(m_breakPattern[irank]->getFromTopology());

              mesh::Selector selector(m_eMesh.getFEM_meta_data()->locally_owned_part());
              if (fromPartsAll.size())
                {
                  selector = mesh::Selector();
                  for (unsigned ipart = 0; ipart < fromPartsAll.size(); ipart++)
                    {
                      mesh::Part *part = fromPartsAll[ipart];
                      const CellTopologyData * part_cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(*part);
                      if (part_cell_topo_data)
                        {
                          CellTopology part_cell_topo(part_cell_topo_data);
                          if (part_cell_topo.getKey() == elementType)
                            {
                              selector = selector | *part;
                            }
                        }
                    }
                  selector = selector & (mesh::Selector(m_eMesh.getFEM_meta_data()->locally_owned_part()));
                }
              std::vector<unsigned> count ;
              stk::mesh::count_entities( selector, *m_eMesh.getBulkData(), count );
              if (count.size() < 3)
                {
                  throw std::logic_error("logic error in Refiner m_refinementInfoByType");
                }
              unsigned n_ele = count[ ranks[irank] ];

              m_refinementInfoByType[irank].m_numOrigElems = n_ele;

              m_refinementInfoByType[irank].m_numNewElems = n_ele * m_breakPattern[irank]->getNumNewElemPerElem();
              m_refinementInfoByType[irank].m_topology = cell_topo;
            }

          // sum info from all procs
          {
            stk::ParallelMachine pm = m_eMesh.getBulkData()->parallel();

            for (unsigned irank = 0; irank < ranks.size(); irank++)
              {
                stk::all_reduce( pm, stk::ReduceSum<1>( &m_refinementInfoByType[irank].m_numOrigElems ) );
                stk::all_reduce( pm, stk::ReduceSum<1>( &m_refinementInfoByType[irank].m_numNewElems ) );
              }
          }

        }

      if (m_doQueryOnly)
        {
          return;
        }

      // do elements first, then any faces or edge elements

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // do the top level, all elements of this rank operation
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      unsigned num_elem_not_ghost_0 = 0;

      ///////////////////////////////////////////////////////////
      ///// Do the mesh coloring step for each type of element
      ///////////////////////////////////////////////////////////

      vector< vector< ColorerSetType > > elementColorsByType = vector < vector< ColorerSetType > > (ranks.size());

      for (unsigned irank = 0; irank < ranks.size(); irank++)
        {
          EXCEPTWATCH;
          unsigned elementType = m_breakPattern[irank]->getFromTypeKey();
          shards::CellTopology cell_topo(m_breakPattern[irank]->getFromTopology());

          if (TRACE_STAGE_PRINT) std::cout << "tmp Refiner:: irank = " << irank << " ranks[irank] = " << ranks[irank]
                                           << " elementType= " << elementType
                                           << " cell_topo= " << cell_topo.getName()
                                           << std::endl;

          std::vector<stk::mesh::EntityRank> ranks_one(1, ranks[irank]);

          // this gives a list of colored elements for this element type only
          stk::mesh::PartVector * fromParts = 0;
          fromParts = &(m_breakPattern[irank]->getFromParts());

          //!FIXME add part info
          Colorer meshColorerThisTypeOnly(elementColorsByType[irank], ranks_one);   TRACE_PRINT("Refiner: Color mesh (all top level rank elements)... ");
          meshColorerThisTypeOnly.color(m_eMesh, &elementType, fromParts);          TRACE_PRINT("Refiner: Color mesh (all top level rank elements)...done ");

          if (0 && elementColorsByType[irank].size() == 0)
            {
              std::cout << "WARNING: no elements found of this type: " << cell_topo.getName() << " key= " << elementType << std::endl;
            }
        }

      // FIXME warn if a topology shows up without a break pattern

      ///////////////////////////////////////////////////////////
      /////  // start top-level ranks
      ///////////////////////////////////////////////////////////

      {   // start top-level ranks


        ///////////////////////////////////////////////////////////
        /////  // node registration step
        ///////////////////////////////////////////////////////////

        {  // node registration step
          EXCEPTWATCH;

          /**/  TRACE_PRINT("Refiner: beginRegistration (top-level rank)... ");
          m_nodeRegistry->dumpDB("before init");
          if (m_alwaysInitNodeRegistry)
            {
              m_nodeRegistry->initialize();
            }
          else
            {
              //m_nodeRegistry->clear_element_owner_data();
              m_nodeRegistry->init_entity_repo();
            }

          m_nodeRegistry->init_comm_all();                           

          m_eMesh.adapt_parent_to_child_relations().clear();

          // register non-ghosted elements needs for new nodes, parallel create new nodes
          m_nodeRegistry->beginRegistration();

          for (unsigned irank = 0; irank < ranks.size(); irank++)
            {
              unsigned elementType = m_breakPattern[irank]->getFromTypeKey();
              {
                EXCEPTWATCH;

                vector< ColorerSetType >& elementColors = elementColorsByType[irank];

                vector<NeededEntityType> needed_entity_ranks;
                m_breakPattern[irank]->fillNeededEntities(needed_entity_ranks);

                bool count_only = false;
                bool doAllElements = true;

                unsigned num_elem_not_ghost_0_incr = doForAllElements(ranks[irank], &NodeRegistry::registerNeedNewNode, elementColors, elementType, needed_entity_ranks,
                                                                      count_only, doAllElements);

                num_elem_not_ghost_0 += num_elem_not_ghost_0_incr;
              }
            }

          m_nodeRegistry->endRegistration();                    /**/  TRACE_PRINT("Refiner: endRegistration (top-level rank)... ");
        }
        m_nodeRegistry->dumpDB("after registration");

#define CHECK_DEBUG 0
        if (CHECK_DEBUG)
          {
            MPI_Barrier( MPI_COMM_WORLD );
            std::cout << "P["<< m_eMesh.getRank()
                      <<"] ========================================================================================================================" << std::endl;
            m_nodeRegistry->checkDB();
            check_db("after registerNeedNewNode");
            MPI_Barrier( MPI_COMM_WORLD );
            std::cout << "P["<< m_eMesh.getRank()
                      <<"] ========================================================================================================================" << std::endl;
          }

        ///////////////////////////////////////////////////////////
        /////  Check for remote
        ///////////////////////////////////////////////////////////

        {   // beginCheckForRemote()
          EXCEPTWATCH;

          /**/                                                        TRACE_PRINT("Refiner: beginCheckForRemote (top-level rank)... ");

          // now register ghosted elements needs for new nodes (this does a pack operation)
          m_nodeRegistry->beginCheckForRemote();
          unsigned num_elem = 0;
          for (unsigned irank = 0; irank < ranks.size(); irank++)
            {
              unsigned elementType = m_breakPattern[irank]->getFromTypeKey();
              //if (ranks[irank] >= m_eMesh.face_rank())
              {
                EXCEPTWATCH;

                vector< ColorerSetType >& elementColors = elementColorsByType[irank];

                vector<NeededEntityType> needed_entity_ranks;
                m_breakPattern[irank]->fillNeededEntities(needed_entity_ranks);

                bool count_only = false;
                bool doAllElements = false;  // only do ghost elements
                num_elem = doForAllElements(ranks[irank], &NodeRegistry::checkForRemote, elementColors, elementType, needed_entity_ranks, count_only, doAllElements);
              }
            }
          m_nodeRegistry->endCheckForRemote();                /**/   TRACE_PRINT("Refiner: endCheckForRemote (top-level rank)... ");

          if (1 && CHECK_DEBUG)
            {
              std::cout << "num_elem= " << num_elem << std::endl;
              MPI_Barrier( MPI_COMM_WORLD );
              std::cout << "P["<< m_eMesh.getRank()
                        <<"] ========================================================================================================================" << std::endl;
              m_nodeRegistry->checkDB();
              check_db("after checkForRemote");
              MPI_Barrier( MPI_COMM_WORLD );
              std::cout << "P["<< m_eMesh.getRank()
                        <<"] ========================================================================================================================" << std::endl;
            }

        }

        ///////////////////////////////////////////////////////////
        /////  Get from remote
        ///////////////////////////////////////////////////////////
        /// communicate all-to-all the new node creation information which also updates the node registry so it can
        /// be queried locally now for any ghost or non-ghost element

        { // get from remote

          EXCEPTWATCH;

          /**/                                                        TRACE_PRINT("Refiner: beginGetFromRemote (top-level rank)... ");
          m_nodeRegistry->beginGetFromRemote();
          unsigned num_elem = 0;
          for (unsigned irank = 0; irank < ranks.size(); irank++)
            {
              unsigned elementType = m_breakPattern[irank]->getFromTypeKey();
              {
                EXCEPTWATCH;

                vector< ColorerSetType >& elementColors = elementColorsByType[irank];

                vector<NeededEntityType> needed_entity_ranks;
                m_breakPattern[irank]->fillNeededEntities(needed_entity_ranks);

                bool count_only = false;
                bool doAllElements = false;   // ghost elements only
                num_elem = doForAllElements(ranks[irank], &NodeRegistry::getFromRemote, elementColors, elementType, needed_entity_ranks, count_only, doAllElements);
              }
            }

          m_nodeRegistry->endGetFromRemote();                    /**/  TRACE_PRINT("Refiner: endGetFromRemote (top-level rank)... ");
          m_nodeRegistry->dumpDB("after endGetFromRemote");

          //stk::diag::printTimersTable(std::cout, perceptTimer(), stk::diag::METRICS_ALL, false);

          if (CHECK_DEBUG)
            {
              std::cout << "num_elem= " << num_elem << std::endl;
              MPI_Barrier( MPI_COMM_WORLD );
              std::cout << "P["<< m_eMesh.getRank()
                        <<"] ========================================================================================================================" << std::endl;
              m_nodeRegistry->checkDB();
              check_db("end getFromRemote");
              MPI_Barrier( MPI_COMM_WORLD );
              
              std::cout << "P["<< m_eMesh.getRank()
                        <<"] ========================================================================================================================" << std::endl;
            }
        }  // get from remote
      } // start top-level ranks

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // for each element type, in top-down rank order, do the rest of the refinement operations
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      bulkData.modification_begin();
      for (unsigned irank = 0; irank < ranks.size(); irank++)
        {
          EXCEPTWATCH;

          unsigned elementType = m_breakPattern[irank]->getFromTypeKey();
          if (TRACE_STAGE_PRINT)
            std::cout << "tmp Refiner:: irank = " << irank
                      << " ranks[irank] = " << ranks[irank] << " elementType= " << elementType << std::endl;

          std::vector<stk::mesh::EntityRank> ranks_one(1, ranks[irank]);

          vector< ColorerSetType >& elementColors = elementColorsByType[irank];

          // loop over elements, build faces, edges in threaded mode (guaranteed no mem conflicts)
          // (note: invoke UniformRefinerPattern: what entities are needed)
          vector<NeededEntityType> needed_entity_ranks;
          m_breakPattern[irank]->fillNeededEntities(needed_entity_ranks);

          vector<stk::mesh::Entity *> new_elements;

          //bulkData.modification_begin();

          {
            EXCEPTWATCH;

            // count num new elements needed on this proc (served by UniformRefinerPattern)
            bool count_only = true;
            bool doAllElements = true;
            /**/                                                TRACE_PRINT("Refiner: registerNeedNewNode count_only(true) ranks[irank]==ranks[0]... ");
            unsigned num_elem_not_ghost = doForAllElements(ranks[irank], &NodeRegistry::registerNeedNewNode, elementColors, elementType, needed_entity_ranks, count_only, doAllElements);
            /**/                                                TRACE_PRINT("Refiner: registerNeedNewNode count_only(true) ranks[irank]==ranks[0]... done ");

            unsigned num_elem_needed = num_elem_not_ghost * m_breakPattern[irank]->getNumNewElemPerElem();

            if (0 && num_elem_not_ghost != num_elem_not_ghost_0)
              {
                std::cout << "num_elem_not_ghost_0 = " << num_elem_not_ghost_0 << " num_elem_not_ghost= " << num_elem_not_ghost << std::endl;
                throw std::runtime_error("num_elem_not_ghost_0 != num_elem_not_ghost");
              }

            // create new entities on this proc
            m_nodeRegistry->beginLocalMeshMods();
            new_elements.resize(0);                                                /**/ TRACE_PRINT("Refiner: createEntities... ranks[irank]==ranks[0] ");
            m_eMesh.createEntities( ranks[irank], num_elem_needed, new_elements);  /**/ TRACE_PRINT("Refiner: createEntities... ranks[irank]==ranks[0] done ");
            m_nodeRegistry->endLocalMeshMods();

          }
          m_nodeRegistry->dumpDB("after endLocalMeshMods");

          /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          ///  Global element ops: here's where we e.g. connect the new elements by declaring new relations
          /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          /**/                                                TRACE_PRINT("Refiner: createElementsAndNodesAndConnectLocal... ");
          /**/                                                TRACE_CPU_TIME_AND_MEM_0(CONNECT_LOCAL);

          createElementsAndNodesAndConnectLocal(ranks[irank], m_breakPattern[irank], elementColors, needed_entity_ranks, new_elements);

          /**/                                                TRACE_CPU_TIME_AND_MEM_1(CONNECT_LOCAL);
          /**/                                                TRACE_PRINT("Refiner: createElementsAndNodesAndConnectLocal...done ");

          /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          ///  Global node loop operations:  this is where we perform ops like adding new nodes to the right parts, interpolating fields, etc.
          /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

          if (TRACE_STAGE_PRINT && !m_eMesh.getRank()) {
            Util::trace_cpu_time_and_mem_print(CONNECT_LOCAL, "CONNECT_LOCAL");
            Util::trace_cpu_time_and_mem_print(CONNECT_LOCAL_createNewNeededNodes, "CONNECT_LOCAL_createNewNeededNodes");
            Util::trace_cpu_time_and_mem_print(CONNECT_LOCAL_createNewElements, "CONNECT_LOCAL_createNewElements");
            Util::trace_cpu_time_and_mem_print(CONNECT_LOCAL_URP_createOrGetNode, "CONNECT_LOCAL_URP_createOrGetNode");
            Util::trace_cpu_time_and_mem_print(CONNECT_LOCAL_URP_declare_relation, "CONNECT_LOCAL_URP_declare_relation");
          }

        } // irank

      /**/                                                TRACE_PRINT("Refiner: addToExistingParts [etc.]... ");
#if !STK_ADAPT_URP_LOCAL_NODE_COMPS
      if (1)
        {
          EXCEPTWATCH;
          // only need to do this once: the map is fully built and we loop over the map's faces/edges, which are fixed after the getFromRemote step

          m_nodeRegistry->addToExistingPartsNew();
          //std::cout << "tmp makeCentroid... " << std::endl;
          m_nodeRegistry->makeCentroid(m_eMesh.getCoordinatesField());
          //std::cout << "tmp makeCentroid...done " << std::endl;
          //std::cout << "tmp interpolateFields... " << std::endl;
          //FIXME m_nodeRegistry->interpolateFields();
          //std::cout << "tmp interpolateFields...done " << std::endl;
        }
      //std::cout << "tmp dumpElements 1" << std::endl;
      // m_eMesh.dumpElements();
#endif
      /**/                                                TRACE_PRINT("Refiner: addToExistingParts [etc.] ...done ");

      /***********************/                           TRACE_PRINT("Refiner: fixElementSides1 ");
      fixElementSides1();
      m_eMesh.adapt_parent_to_child_relations().clear();
      /***********************/                           TRACE_PRINT("Refiner: fixElementSides1...done ");

      //std::cout << "tmp dumpElements 2" << std::endl;
      //m_eMesh.dumpElements();

#if CHECK_DEBUG
      std::cout << "m_doRemove= " << m_doRemove << std::endl;
      check_db("b4 remove");
#endif

      if (m_doRemove)
        {
          EXCEPTWATCH;

          //bulkData.modification_begin();

          /***********************/                           TRACE_PRINT("Refiner: fixElementSides1 ");
          //           fixElementSides1();
          //           m_eMesh.adapt_parent_to_child_relations().clear();
          /***********************/                           TRACE_PRINT("Refiner: fixElementSides1...done ");

          if (m_doProgress)
            {
              ProgressMeterData pd(ProgressMeterData::INIT, 0.0, "removeOldElements");
              notifyObservers(&pd);
            }

          for (unsigned irank = 0; irank < ranks.size(); irank++)
            {
              if (m_doProgress)
                {
                  ProgressMeterData pd(ProgressMeterData::RUNNING, 100.0*((double)irank)/((double)ranks.size()), "removeOldElements" );
                  notifyObservers(&pd);
                }

#if PERCEPT_USE_FAMILY_TREE
              if (irank == 0)
                removeFamilyTrees();
#endif

              removeOldElements(ranks[irank], m_breakPattern[irank]);
              renameNewParts(ranks[irank], m_breakPattern[irank]);
              fixSurfaceAndEdgeSetNames(ranks[irank], m_breakPattern[irank]);
            }
          if (m_doProgress)
            {
              ProgressMeterData pd(ProgressMeterData::FINI, 0.0, "removeOldElements");
              notifyObservers(&pd);
            }

        }
      else  // m_doRemove
        {
          if (0)
            {
              for (unsigned irank = 0; irank < ranks.size(); irank++)
                removeFromOldPart(ranks[irank], m_breakPattern[irank]);
            }
        }

      /**/                                                TRACE_PRINT("Refiner: modification_end...start... ");
      bulkData.modification_end();
      /**/                                                TRACE_PRINT("Refiner: modification_end...done ");

      //std::cout << "tmp dumpElements 3" << std::endl;
      //m_eMesh.dumpElements();

#if defined( STK_ADAPT_HAS_GEOMETRY )
      if (m_geomSnap)
        {

            GeometryKernelOpenNURBS gk;
            MeshGeometry mesh_geometry(&gk);
            GeometryFactory factory(&gk, &mesh_geometry);
            factory.read_file(m_geomFile, &m_eMesh);
            mesh_geometry.snap_points_to_geometry(&m_eMesh);

            if (m_doSmoothGeometry)
              {
                smoothGeometry(mesh_geometry);
                mesh_geometry.snap_points_to_geometry(&m_eMesh);
              }



        }
#endif

      /**/                                                TRACE_PRINT( "Refiner:doBreak ... done");

      m_nodeRegistry->dumpDB("after doBreak");
#if CHECK_DEBUG
      check_db("after doBreak");
#endif
      //std::cout << "tmp m_nodeRegistry.m_gee_cnt= " << m_nodeRegistry->m_gee_cnt << std::endl;
      //std::cout << "tmp m_nodeRegistry.m_gen_cnt= " << m_nodeRegistry->m_gen_cnt << std::endl;

    } // doBreak

    
    // FIXME - temp until we figure out what to do with parent/child, persistence, etc.
    // FIXME - just deletes elements, not family trees for now

    /// Delete all elements that aren't child elements
    void Refiner::deleteParentElements()
    {

      elements_to_be_destroyed_type parents;

      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.getBulkData()->buckets( m_eMesh.element_rank() );

      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
        {
          stk::mesh::Bucket & bucket = **k ;

          // only do "old" elements
          //if (!oldPartSelector(bucket))
          //  continue;

          const unsigned num_elements_in_bucket = bucket.size();
          for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
            {
              stk::mesh::Entity& element = bucket[iElement];
              if (m_eMesh.isParentElement(element, false))
                {
#if UNIFORM_REF_REMOVE_OLD_STD_VECTOR
                      parents.push_back(&element);
#else
                      parents.insert(&element);
#endif
                }
            }
        }

      m_eMesh.getBulkData()->modification_begin();
#if PERCEPT_USE_FAMILY_TREE
      removeFamilyTrees();
#endif
      //std::cout << "tmp removeOldElements(parents) " << std::endl;
      removeOldElements(parents);
      m_eMesh.getBulkData()->modification_end();

    }

#if defined( STK_ADAPT_HAS_GEOMETRY )
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

    // orthogonal corrections
    static void move_quad_centroid_node_method1(int spatialDim, double *centroid_node_orig_pos, double *centroid_node_new_pos, double *edge_centroids[4])
    {
      double *edges[2][2] = {{edge_centroids[0], edge_centroids[2]}, {edge_centroids[1], edge_centroids[3]} };
      double proposed_delta[3] = {0,0,0};

      for (int idim=0; idim < spatialDim; idim++)
        {
          centroid_node_new_pos[idim] = centroid_node_orig_pos[idim];
        }

      for (int iedge=0; iedge < 2; iedge++)
        {
          double dxyz[3]={0,0,0};
          double d=0;
          for (int idim=0; idim < spatialDim; idim++)
            {
              //proposed_delta[idim] = (centroid_node_new_pos[idim] - centroid_node_orig_pos[idim]);
              double midpoint= 0.5*(edges[iedge][1][idim] + edges[iedge][0][idim]);
              proposed_delta[idim] = (midpoint - centroid_node_orig_pos[idim]);
              
              dxyz[idim] = edges[iedge][1][idim] - edges[iedge][0][idim];
              d += dxyz[idim]*dxyz[idim];
            }
          double proj=0;
          if (d<1.e-10) d=1;
          d = std::sqrt(d);
          for (int idim=0; idim < spatialDim; idim++)
            {
              dxyz[idim] /= d;
              proj += dxyz[idim]*proposed_delta[idim];
            }
          for (int idim=0; idim < spatialDim; idim++)
            {
              centroid_node_new_pos[idim] += proj*dxyz[idim];
            }
        }
    }
#endif

#if defined( STK_ADAPT_HAS_GEOMETRY )
    void Refiner::smoothGeometry(MeshGeometry& mesh_geometry)
    {
      /**
       *  0. cache set of nodes involved in moving/snapping to geometry
       *  1. get parent elements touching the boundary - just the set of elements with a node in the set from (0)
       *  2. for each elem, loop over its faces, create subDimEntity, find new node from NodeRegistry
       *  3. if new node is not in set (0), compute new location as average of the face's edge mid-nodes
       *  4. compute new element centroid by average of face centroids (or all edge centroids?)
       *
       *  alternative (avoids computing face centroids twice): 
       *
       *  0. cache set of nodes involved in moving/snapping to geometry
       *  1. get parent elements touching the boundary - just the set of elements with a node in the set from (0)
       *  2. for each face in NodeRegistry DB, if face has a node in set (0), move its centroid to the edge mid node average,
       *       but only if its node is not in set (0)
       *  4. loop over set (1), compute new element centroid by average of face centroids (or all edge centroids?)
       *  
       */
      static SubDimCellData empty_SubDimCellData;

      int debug = 1;  // or 2, 3,...

      if (m_doRemove)
        {
          throw std::runtime_error("Refiner::smoothGeometry: to use this method, you must call setRemoveOldElements(false)\n"
                                   "  or if using stk_adapt_exe, use option --remove_original_elements=0");
        }

      stk::mesh::Selector on_geometry_selector = getNodeWasSnappedSelector(mesh_geometry);

      std::string oldPartName = m_breakPattern[0]->getOldElementsPartName()+toString(m_eMesh.element_rank());
      mesh::Part *oldPart = m_eMesh.getFEM_meta_data()->get_part(oldPartName);
      //std::cout << "oldPartName= " << oldPartName << std::endl;
      if (!oldPart)
        {
          std::cout << "oldPartName= " << oldPartName << std::endl;
          throw std::runtime_error("oldpart is null in smoothGeometry");
        }
      stk::mesh::Selector oldPartSelector(*oldPart);

      int spatialDim = m_eMesh.getSpatialDim();
      int topoDim = -1;

      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.getBulkData()->buckets( m_eMesh.element_rank() );

      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
        {
          stk::mesh::Bucket & bucket = **k ;

          // only do "old" elements
          //if (!oldPartSelector(bucket))
          //  continue;

          const unsigned num_elements_in_bucket = bucket.size();
          for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
            {
              stk::mesh::Entity& element = bucket[iElement];

              const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);
                
              CellTopology cell_topo(cell_topo_data);

              if (m_eMesh.isParentElementLeaf(element, false))
                {
                  if (debug > 1) 
                    std::cout << "tmp debug smoothGeometry: parent elem id= " << element.identifier() << std::endl;
                  if (debug > 1)
                    {
                      std::cout << "tmp cell_topo = " << cell_topo.getName() << std::endl;
                      m_eMesh.printParentChildInfo(element, false);
                    }
                }
              else
                {
                  if (debug > 2) std::cout << "tmp debug smoothGeometry: child elem id= " << element.identifier() << std::endl;
                  continue;
                }


              // avoid ghosts
              if (m_eMesh.isGhostElement(element))
                {
                  continue;
                }


              // avoid elements that don't touch the boundary
              bool element_touches_snapped_geometry = false;
              stk::mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);
              for (unsigned inode=0; inode < elem_nodes.size(); inode++)
                {
                  stk::mesh::Entity *node = elem_nodes[inode].entity();
                  if (on_geometry_selector(*node))
                    {
                      element_touches_snapped_geometry = true;
                      break;
                    }
                }
              if (!element_touches_snapped_geometry)
                continue;


              // get sub-dim entities and look for centroid nodes
              //const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);

              // loop from top rank downwards (m_ranks has ranks sorted in top-down order)
              for (unsigned irank = 0; irank < m_ranks.size(); irank++)
                {
                  if (m_ranks[irank] != m_eMesh.element_rank())
                    continue;

                  vector<NeededEntityType> needed_entity_ranks;
                  m_breakPattern[irank]->fillNeededEntities(needed_entity_ranks);
                  unsigned elementType = m_breakPattern[irank]->getFromTypeKey();

                  topoDim = m_breakPattern[irank]->getTopoDim(cell_topo);

                  if (cell_topo.getKey() == elementType)
                    {

                      unsigned numFaces = 0;
                      unsigned numEdges = 0;

                      for (unsigned ineed_ent=0; ineed_ent < needed_entity_ranks.size(); ineed_ent++)
                        {
                          stk::mesh::EntityRank needed_entity_rank = needed_entity_ranks[ineed_ent].first;

                          if (needed_entity_rank == m_eMesh.edge_rank())
                            {
                              numEdges = cell_topo_data->edge_count;
                            }
                          else if (needed_entity_rank == m_eMesh.face_rank())
                            {
                              numFaces = cell_topo_data->side_count;
                            }
                        }

                      if (debug > 3)  std::cout << "tmp debug Refiner:: irank = " << irank << " ranks[irank] = " << m_ranks[irank]
                                       << " elementType= " << elementType
                                       << " cell_topo= " << cell_topo.getName()
                                       << " numFaces= " << numFaces 
                                       << " numEdges= " << numEdges
                                       << " topoDim= " << topoDim
                                       << std::endl;

                      if (topoDim == 3 && numFaces == 0)
                        continue;
                      if (topoDim == 2 && numEdges == 0)
                        continue;

                      // skip shells for now?
                      if (topoDim == 2 && spatialDim == 3)
                        continue;

                      // do this for either 2d or 3d (operate on edges to average their values to the face centroid - if in 2D, then we're done,
                      //   else in 3d have to average face centroids to element centroid as well)
                      {
                        EntityRank face_rank = (topoDim == 3? m_eMesh.face_rank() : m_eMesh.element_rank());
                        if (topoDim == 2) numFaces = 1;
                          
                        for (unsigned iSubDimOrd = 0; iSubDimOrd < numFaces; iSubDimOrd++)
                          {
                            // get the face (or edge in 2d)
                            SubDimCell_SDSEntityType subDimEntity;
                            m_nodeRegistry->getSubDimEntity(subDimEntity, element, face_rank, iSubDimOrd);

                            SubDimCellData* nodeId_elementOwnderId_ptr = m_nodeRegistry->getFromMapPtr(subDimEntity);
                            SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
                            bool is_empty = nodeId_elementOwnderId_ptr == 0;
                            if (is_empty)
                              {
                                if (1) std::cout << "tmp debug error 1 Refiner:: irank = " << irank << " ranks[irank] = " << m_ranks[irank]
                                                 << " elementType= " << elementType
                                                 << " cell_topo= " << cell_topo.getName()
                                                 << " numFaces= " << numFaces 
                                                 << " numEdges= " << numEdges
                                                 << " topoDim= " << topoDim
                                                 << std::endl;

                                throw std::logic_error("error1");
                              }

                            stk::mesh::Entity *node_at_centroid = 0;
                            double centroid[3] = {0,0,0};

                            NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
                            node_at_centroid = nodeIds_onSE[0]; // FIXME for quadratic elements

                            if (0 && on_geometry_selector(*node_at_centroid))
                              {
                                if (debug) std::cout << "tmp debug face centroid node is on geometry so we won't be moving it" << std::endl;
                                continue;
                              }
                            //stk::mesh::EntityRank erank = stk::mesh::entity_rank(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>());
                            //stk::mesh::EntityId owning_elementId = stk::mesh::entity_id(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>());
                                  
                            //if (owning_elementId == element.identifier())  // possible performance enhancement, but at the expense of visiting elements twice
                            {

                              const CellTopologyData *  face_topo_data = (topoDim == 2 ? cell_topo_data : cell_topo_data->side[iSubDimOrd].topology);
                              if (debug > 3) std::cout << "tmp debug face_topo_data->edge_count = " << face_topo_data->edge_count << std::endl;
                              double *edge_centroids[4];
                              for (unsigned ie=0; ie < face_topo_data->edge_count; ie++)
                                {
                                  SubDimCell_SDSEntityType edgeSubDim;
                                  if (topoDim == 3)
                                    {
                                      edgeSubDim.insert( elem_nodes[ cell_topo_data->side[iSubDimOrd].node[ face_topo_data->edge[ie].node[0] ] ].entity() );
                                      edgeSubDim.insert( elem_nodes[ cell_topo_data->side[iSubDimOrd].node[ face_topo_data->edge[ie].node[1] ] ].entity() );
                                    }
                                  else
                                    {
                                      edgeSubDim.insert( elem_nodes[  face_topo_data->edge[ie].node[0] ].entity() );
                                      edgeSubDim.insert( elem_nodes[  face_topo_data->edge[ie].node[1] ].entity() );
                                    }

                                  SubDimCellData* edge_data_ptr = m_nodeRegistry->getFromMapPtr(edgeSubDim);
                                  SubDimCellData& edge_data = (edge_data_ptr ? *edge_data_ptr : empty_SubDimCellData);
                                  bool is_empty_edge = edge_data_ptr == 0;
                                  if (is_empty_edge)
                                    {
                                      throw std::logic_error("couldn't find edge");
                                    }
                                  NodeIdsOnSubDimEntityType& nodeIds_onSE_edge = edge_data.get<SDC_DATA_GLOBAL_NODE_IDS>();

                                  stk::mesh::Entity *edge_node = nodeIds_onSE_edge[0]; // FIXME for quadratic elements
                                  double *fdata = stk::mesh::field_data( *m_eMesh.getCoordinatesField() , *edge_node );
                                  edge_centroids[ie] = fdata;
                                  if (debug > 3) 
                                    {
                                      std::cout << "tmp debug found edge node: " ; 
                                      m_eMesh.printEntity(std::cout, *edge_node);
                                    }
                                  for (int idim=0; idim < m_eMesh.getSpatialDim(); idim++)
                                    {
                                      centroid[idim] += fdata[idim];
                                    }
                                }

                              for (int idim=0; idim < m_eMesh.getSpatialDim(); idim++)
                                {
                                  centroid[idim] /= ((double)face_topo_data->edge_count);
                                }
                              double *fdata = stk::mesh::field_data( *m_eMesh.getCoordinatesField() , *node_at_centroid );

                              double dist_moved = 0.0;

                              // this can be commented out, which gives the original behavior (simple average all 4 edge nodes)
                              if (1)
                                move_quad_centroid_node_method1(spatialDim, fdata, centroid, edge_centroids);
                              
                              if (debug)
                                std::cout << "tmp debug moving face node with coords: ";
                              for (int idim=0; idim < m_eMesh.getSpatialDim(); idim++)
                                {
                                  if (debug)
                                    std::cout << fdata[idim] << " ";
                                  dist_moved += (fdata[idim] - centroid[idim])*(fdata[idim] - centroid[idim]);
                                  fdata[idim] = centroid[idim];
                                }
                              dist_moved = std::sqrt(dist_moved);
                              if (debug)
                                std::cout << " : a distance of " << dist_moved << std::endl;
                              
                            }
                          }

                        if (topoDim == 3)
                          {
                          
                            double element_centroid[3] = {0,0,0};
                            for (unsigned iSubDimOrd = 0; iSubDimOrd < numFaces; iSubDimOrd++)
                              {
                                SubDimCell_SDSEntityType subDimEntity;

                                m_nodeRegistry->getSubDimEntity(subDimEntity, element, m_eMesh.face_rank(), iSubDimOrd);

                                SubDimCellData* nodeId_elementOwnderId_ptr = m_nodeRegistry->getFromMapPtr(subDimEntity);
                                SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
                                bool is_empty = nodeId_elementOwnderId_ptr == 0;
                                if (is_empty)
                                  {
                                    //continue;
                                    throw std::logic_error("error2");
                                  }

                                NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
                                stk::mesh::Entity* face_node_at_centroid = nodeIds_onSE[0]; // FIXME for quadratic elements

                                double *fdata = stk::mesh::field_data( *m_eMesh.getCoordinatesField() , *face_node_at_centroid );
                                for (int idim=0; idim < m_eMesh.getSpatialDim(); idim++)
                                  {
                                    element_centroid[idim] += fdata[idim];
                                  }
                              }

                            // get the element centroid node

                            SubDimCell_SDSEntityType subDimEntity_element;
                            m_nodeRegistry->getSubDimEntity(subDimEntity_element, element, m_eMesh.element_rank(), 0);

                            SubDimCellData* nodeId_elementOwnderId_ptr = m_nodeRegistry->getFromMapPtr(subDimEntity_element);
                            SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
                            bool is_empty = nodeId_elementOwnderId_ptr == 0;
                            if (is_empty)
                              {
                                if (1) std::cout << "tmp debug error3 Refiner:: irank = " << irank << " ranks[irank] = " << m_ranks[irank]
                                                 << " elementType= " << elementType
                                                 << " cell_topo= " << cell_topo.getName()
                                                 << std::endl;

                                throw std::logic_error("error3");
                              }

                            NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();

                            stk::mesh::Entity *element_node_at_centroid = nodeIds_onSE[0]; // FIXME for quadratic elements

                            if (on_geometry_selector(*element_node_at_centroid))
                              {
                                std::cout << "tmp debug element centroid node is on geometry ... this is an error" << std::endl;
                                throw std::logic_error("tmp debug element centroid node is on geometry, error...");
                              }


                            for (int idim=0; idim < m_eMesh.getSpatialDim(); idim++)
                              {
                                element_centroid[idim] /= ((double)numFaces);
                              }
                            double *fdata = stk::mesh::field_data( *m_eMesh.getCoordinatesField() , *element_node_at_centroid );
                            double dist_moved = 0.0;
                            if (debug)
                              std::cout << "tmp debug moving element centroid node with coords: ";
                            for (int idim=0; idim < m_eMesh.getSpatialDim(); idim++)
                              {
                                if (debug)
                                  std::cout << fdata[idim] << " ";
                                dist_moved += (fdata[idim] - element_centroid[idim])*(fdata[idim] - element_centroid[idim]);
                                fdata[idim] = element_centroid[idim];
                              }
                            dist_moved = std::sqrt(dist_moved);
                            if (debug)
                              std::cout << " : a distance of " << dist_moved << std::endl;
                          }

                        }
                    }
                }
            }
        }
    }
#endif

    unsigned Refiner::
    doForAllElements(stk::mesh::EntityRank rank, NodeRegistry::ElementFunctionPrototype function,
                     vector< ColorerSetType >& elementColors, unsigned elementType,
                     vector<NeededEntityType>& needed_entity_ranks,
                     bool only_count, bool doAllElements)
    //bool only_count=false, bool doAllElements=true)
    {
      EXCEPTWATCH;
      unsigned num_elem = 0;

      int progress_meter_num_total = 0;
      if (m_doProgress)
        {
          m_doProgress = false;
          progress_meter_num_total = doForAllElements(rank, function, elementColors, elementType, needed_entity_ranks, true, doAllElements);
          m_doProgress = true;
          ProgressMeterData pd(ProgressMeterData::INIT, 0.0, "NodeRegistry passes");
          notifyObservers(&pd);
        }
      int progress_meter_when_to_post = progress_meter_num_total / m_progress_meter_frequency;
      if (0 == progress_meter_when_to_post)
        progress_meter_when_to_post = 1;
      double d_progress_meter_num_total = progress_meter_num_total;

      for (unsigned icolor = 0; icolor < elementColors.size(); icolor++)
        {
          if (elementColors[icolor].size() == 0)
            {
              std::cout << "tmp doForAllElements elementColors size = 0!!!" << std::endl;
              continue;
            }

          //stk::mesh::Entity* first_element_p = *(elementColors[icolor].begin());
          //const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(*first_element_p);

          // do in threaded mode FIXME
          for (ColorerSetType::iterator iele = elementColors[icolor].begin();
               iele !=  elementColors[icolor].end();
               iele++)
            {
              const stk::mesh::Entity * element_p =  *iele;
              const stk::mesh::Entity& element = * element_p;

              // FIXME
              // skip elements that are already a parent (if there's no family tree yet, it's not a parent, so avoid throwing an error if isParentElement)
              const bool check_for_family_tree = false;  
              bool isParent = m_eMesh.isParentElement(element, check_for_family_tree);
              
              if (isParent)
                continue;

              bool elementIsGhost = m_eMesh.isGhostElement(element);
              if (!elementIsGhost)
                ++num_elem;

              if (!only_count && (doAllElements || elementIsGhost))
                {
                  //m_nodeRegistry->doForAllSubEntities(function, element, needed_entity_ranks);
                  applyNodeRegistryFunctionForSubEntities(function, element, needed_entity_ranks);
                }

              if (m_doProgress && (num_elem % progress_meter_when_to_post == 0) )
                {
                  double progress_meter_percent = 100.0*((double)num_elem)/d_progress_meter_num_total;
                  ProgressMeterData pd(ProgressMeterData::RUNNING, progress_meter_percent, "NodeRegistry passes");
                  notifyObservers(&pd);
                  if (0) std::cout << "progress_meter_percent = " << progress_meter_percent << std::endl;
                }

            } // elements in this color
        } // icolor

      if (m_doProgress)
        {
          ProgressMeterData pd(ProgressMeterData::FINI, 0.0, "NodeRegistry passes");
          notifyObservers(&pd);
        }

      return num_elem;
    }

    void Refiner::
    createElementsAndNodesAndConnectLocal(stk::mesh::EntityRank rank, UniformRefinerPatternBase *breakPattern,
                                          vector< ColorerSetType >& elementColors,   vector<NeededEntityType>& needed_entity_ranks,
                                          vector<stk::mesh::Entity *>& new_elements_pool)
    {
      EXCEPTWATCH;
      static NewSubEntityNodesType s_new_sub_entity_nodes(stk::percept::EntityRankEnd);

      NewSubEntityNodesType& new_sub_entity_nodes = s_new_sub_entity_nodes;

      vector<stk::mesh::Entity *>::iterator element_pool_it = new_elements_pool.begin();

      int jele = 0;
      int numPrints = 20;

      // create new elements and connect them up

      for (unsigned icolor = 0; icolor < elementColors.size(); icolor++)
        {
          jele += elementColors[icolor].size();
        }

      int nele = jele;
      jele = 0;
      int printEvery = nele/numPrints;
      if (printEvery == 0) printEvery = 1;
      if (0)
        {
          std::cout << "Refiner::createElementsAndNodesAndConnectLocal: rank= " << rank
                    << " elementColors.size() = " << elementColors.size() << " num elements = " << nele
                    << " printEvery= " << printEvery
                    << std::endl;
        }
      if (m_doProgress)
        {
          ProgressMeterData pd(ProgressMeterData::INIT, 0.0, "createElementsAndNodesAndConnectLocal");
          notifyObservers(&pd);
        }

      for (unsigned icolor = 0; icolor < elementColors.size(); icolor++)
        {
          TRACE_PRINT(  "Refiner:createElementsAndNodesAndConnectLocal color= " + percept::toString(icolor) + " [ " +
                        (percept::toString (((double)icolor)/((double)elementColors.size())*100 )).substr(0,4) + " %] ");

          if (elementColors[icolor].size() == 0)
            {
              std::cout << "tmp elementColors size = 0!!!" << std::endl;
              continue;
            }

          stk::mesh::Entity* first_element_p = *(elementColors[icolor].begin());

          const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(*first_element_p);
          CellTopology cell_topo(cell_topo_data);

          for (ColorerSetType::iterator iele = elementColors[icolor].begin();  iele !=  elementColors[icolor].end();  iele++)
            {

              stk::mesh::Entity* element_p = *iele;
              if (!element_p)
                {
                  throw std::runtime_error("Refiner::createElementsAndNodesAndConnectLocal");
                }

              if (0 && (jele % printEvery == 0))
                {
                  std::cout << "Refiner::createElementsAndNodesAndConnectLocal: element # = " << jele << " ["
                            << (((double)jele)/((double)nele)*100.0) << " %]" << std::endl;
                }
              if (m_doProgress && (jele % printEvery == 0))
                {
                  ProgressMeterData pd(ProgressMeterData::RUNNING, 100.0*((double)jele)/((double)nele), "createElementsAndNodesAndConnectLocal RUN" );
                  notifyObservers(&pd);
                }

              stk::mesh::Entity& element = * element_p;

              if (m_proc_rank_field && rank == m_eMesh.element_rank())
                {
                  double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(m_proc_rank_field) , element );
                  fdata[0] = double(element.owner_rank());
                  //if (1 || element.owner_rank() == 3)
                  //  std::cout << "tmp element.owner_rank() = " << element.owner_rank() << std::endl;
                }
              // FIXME

              // skip elements that are already a parent (if there's no family tree yet, it's not a parent, so avoid throwing an error if isParentElement)
              const bool check_for_family_tree = false;  
              bool isParent = m_eMesh.isParentElement(element, check_for_family_tree);
              if (0)
                {
                  const unsigned FAMILY_TREE_RANK = m_eMesh.element_rank() + 1u;
                  stk::mesh::PairIterRelation element_to_family_tree_relations = element.relations(FAMILY_TREE_RANK);
                  if (element_to_family_tree_relations.size() == 1)
                    {
                      std::cout << "tmp isParent = " << isParent << " isChild = " << m_eMesh.isChildElement(element) << " element_to_family_tree_relations.size() = " << element_to_family_tree_relations.size() << std::endl;
                    }
                }
              
              if (isParent)
                continue;


              if (!m_eMesh.isGhostElement(element))
                {
                  //std::cout << "P["<< m_eMesh.getRank() << "] element.owner_rank() = " << element.owner_rank() << std::endl;
                  /**/                                                TRACE_CPU_TIME_AND_MEM_0(CONNECT_LOCAL_createNewNeededNodes);

                  if (createNewNeededNodeIds(cell_topo_data, element, needed_entity_ranks, new_sub_entity_nodes))
                    {
                      std::cout << "typeid= " << typeid(*breakPattern).name() << std::endl;
                      throw std::logic_error("needed_entity_ranks[ineed_ent].second");
                    }

                  /**/                                                TRACE_CPU_TIME_AND_MEM_1(CONNECT_LOCAL_createNewNeededNodes);

                  /**/                                                TRACE_CPU_TIME_AND_MEM_0(CONNECT_LOCAL_createNewElements);

                  breakPattern->createNewElements(m_eMesh, *m_nodeRegistry, element, new_sub_entity_nodes, element_pool_it, m_proc_rank_field);

                  /**/                                                TRACE_CPU_TIME_AND_MEM_1(CONNECT_LOCAL_createNewElements);
                }

              ++jele;
            }
        }

      if (m_doProgress)
        {
          ProgressMeterData pd(ProgressMeterData::FINI, 0.0, "createElementsAndNodesAndConnectLocal");
          notifyObservers(&pd);
        }

    }

    /// create a list of nodes from the new nodes that can be easily deciphered by the UniformRefinerPattern
    /// Returns the 3D array new_sub_entity_nodes[entity_rank][ordinal_of_sub_dim_entity][ordinal_of_node_on_sub_dim_entity]

    bool Refiner::
    createNewNeededNodeIds(const CellTopologyData * const cell_topo_data,
                           const stk::mesh::Entity& element, vector<NeededEntityType>& needed_entity_ranks, NewSubEntityNodesType& new_sub_entity_nodes)
    {
      EXCEPTWATCH;

      NodeRegistry& nodeRegistry = *m_nodeRegistry;

      const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);

      // CHECK - cache this
      for (unsigned ineed_ent=0; ineed_ent < needed_entity_ranks.size(); ineed_ent++)
        {
          unsigned numSubDimNeededEntities = 0;

          // special case of face in 3d or edge in 2d
          if (needed_entity_ranks[ineed_ent].first == element.entity_rank())
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
          else if (needed_entity_ranks[ineed_ent].first == m_eMesh.element_rank())
            {
              numSubDimNeededEntities = 1;
            }

          if (needed_entity_ranks[ineed_ent].first >= new_sub_entity_nodes.size())
            {
              throw std::logic_error("Refiner::createNewNeededNodeIds logic err #1");
            }
          new_sub_entity_nodes[needed_entity_ranks[ineed_ent].first].resize(numSubDimNeededEntities);

          if (0)
            {
              std::cout << "P[" << m_eMesh.getRank() << "]  needed_entity_ranks[ineed_ent]= " << needed_entity_ranks[ineed_ent].first
                        << " , " << needed_entity_ranks[ineed_ent].second << " numSubDimNeededEntities= " << numSubDimNeededEntities
                        << std::endl;
            }

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

              if (!nodeIds_onSE[0]) {

                if (nodeIds_onSE.m_entity_id_vector[0] == 0)
                  {
                    // for debugging cases that may have inconsistent edge marking schemes
#define DEBUG_ALLOW_0_ENTITY_ID_VECTOR 0
                    if (DEBUG_ALLOW_0_ENTITY_ID_VECTOR)
                      {
                        continue;
                      }
                    else
                      {
                        std::cout << "P[" << m_eMesh.getRank() << "] nodeId ## = 0 << "
                                  << " nodeIds_onSE.m_entity_id_vector[0] = " << nodeIds_onSE.m_entity_id_vector[0]
                                  << " element= " << element
                                  << " needed_entity_ranks= " << needed_entity_ranks[ineed_ent].first
                                  << " iSubDimOrd = " << iSubDimOrd
                                  <<  std::endl;
                        std::cout << " element= ";
                        m_eMesh.printEntity(std::cout, element, 0);

                        throw std::logic_error("Refiner::createNewNeededNodeIds logic err #5.0, nodeIds_onSE.m_entity_id_vector[i_new_node] == 0");
                      }
                  }

                stk::mesh::Entity * node1 = m_eMesh.getBulkData()->get_entity(stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE.m_entity_id_vector[0]);

                if (!node1)
                  {
                    std::cout << "P[" << m_eMesh.getRank() << "] nodeId ## = 0 << "
                              << " nodeIds_onSE.m_entity_id_vector[0] = " << nodeIds_onSE.m_entity_id_vector[0] << " node1= " << node1
                              << " element= " << element
                              << " needed_entity_ranks= " << needed_entity_ranks[ineed_ent].first
                              << " iSubDimOrd = " << iSubDimOrd
                              <<  std::endl;
                    throw std::logic_error("Refiner::createNewNeededNodeIds logic error #0");
                  }
              }

              unsigned num_new_nodes_needed = needed_entity_ranks[ineed_ent].second;
              if (0)
                {
                  const CellTopologyData * const cell_topo_data_0 = stk::percept::PerceptMesh::get_cell_topology(element);
                  CellTopology cell_topo_0(cell_topo_data_0);

                  std::cout << "tmp 43 cell_topo= " << cell_topo_0.getName() << " ineed_ent= " << ineed_ent << " needed_entity_ranks[ineed_ent].first/second = "
                            << needed_entity_ranks[ineed_ent].first << " "
                            << needed_entity_ranks[ineed_ent].second
                            << std::endl;
                }

              if (num_new_nodes_needed < 1)
                {
                  //std::cout << "needed_entity_ranks[ineed_ent].second = " << num_new_nodes_needed << std::endl;
                  //throw std::logic_error("needed_entity_ranks[ineed_ent].second");
                  return true;
                }

              if (iSubDimOrd >= new_sub_entity_nodes[needed_entity_ranks[ineed_ent].first].size())
                {
                  throw std::logic_error("Refiner::createNewNeededNodeIds logic err #2");
                }
              //std::cout << "tmp elementid, iSubDimOrd, num_new_nodes_needed = " << element.identifier() << " " << iSubDimOrd << " " << num_new_nodes_needed << std::endl;
              new_sub_entity_nodes[needed_entity_ranks[ineed_ent].first][iSubDimOrd].resize(num_new_nodes_needed);
              if (num_new_nodes_needed > nodeIds_onSE.size())
                {
                  std::cout << "Refiner::createNewNeededNodeIds logic err #3:  num_new_nodes_needed= " << num_new_nodes_needed
                            << " nodeIds_onSE.size() = " << nodeIds_onSE.size() << std::endl;
                  throw std::logic_error("Refiner::createNewNeededNodeIds logic err #3");
                }
              for (unsigned i_new_node = 0; i_new_node < num_new_nodes_needed; i_new_node++)
                {
                  if (!nodeIds_onSE[i_new_node])
                    {
                      if (nodeIds_onSE.m_entity_id_vector[i_new_node] == 0)
                        {
                          if (DEBUG_ALLOW_0_ENTITY_ID_VECTOR)
                            {
                              continue;
                            }
                          else
                            {
                              std::cout << "P[" << m_eMesh.getRank() << "] nodeId ## = 0 << "
                                        << " nodeIds_onSE.m_entity_id_vector[0] = " << nodeIds_onSE.m_entity_id_vector[0]
                                        << " element= " << element
                                        << " needed_entity_ranks= " << needed_entity_ranks[ineed_ent].first
                                        << " iSubDimOrd = " << iSubDimOrd
                                        <<  std::endl;
                              std::cout << " element= ";
                              m_eMesh.printEntity(std::cout, element, 0);

                            }
                          throw std::logic_error("Refiner::createNewNeededNodeIds logic err #5.1, nodeIds_onSE.m_entity_id_vector[i_new_node] == 0");
                        }
                      stk::mesh::Entity * node1 = m_eMesh.getBulkData()->get_entity(stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE.m_entity_id_vector[i_new_node]);

                      if (!node1)
                        {
                          throw std::logic_error("Refiner::createNewNeededNodeIds logic err #4");
                        }
                      nodeIds_onSE[i_new_node] = node1;
                    }
                  new_sub_entity_nodes[needed_entity_ranks[ineed_ent].first][iSubDimOrd][i_new_node] = nodeIds_onSE[i_new_node]->identifier();

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


    /** @deprecated */
    void Refiner::
    fixElementSides()
    {
    }

    void Refiner::
    fixElementSides1()
    {
      EXCEPTWATCH;
      if (getIgnoreSideSets()) return;

      if (m_eMesh.getSpatialDim() == 3)
        {
          fixElementSides1(m_eMesh.face_rank());
        }
      // FIXME
      else if (m_eMesh.getSpatialDim() == 2)
        {
          fixElementSides1(m_eMesh.edge_rank());
        }
    }



#if PERCEPT_USE_FAMILY_TREE == 0
    static const SameRankRelationValue * getChildVectorPtr(  SameRankRelation& repo , stk::mesh::Entity *parent)
    {
      SameRankRelation::const_iterator i = repo.find( parent );
      if (i != repo.end())
        return &i->second;
      else
        return 0;
    }
#endif

    /** Sets orientations and associativity of elements to sub-dimensional faces/edges after refinement.
     */
#define EXTRA_PRINT_UR_FES 0

#if PERCEPT_USE_FAMILY_TREE == 0
    void Refiner::
    fixElementSides1(stk::mesh::EntityRank side_rank)
    {
      EXCEPTWATCH;

      bool notFound = true;
      for (unsigned ibp = 0; ibp < m_breakPattern.size(); ibp++)
        {
          // only check the side elements
          if (m_breakPattern[ibp]->getPrimaryEntityRank() == side_rank)
            {
              notFound = false;
            }
        }
      if (notFound)
        {
          std::cout << "Refiner::fixElementSides1: missing sub-dim break pattern - logic error\n"
            " ---- for this refinement pattern to be able to handle sidesets and edgesets you must provide the sub-dim break pattern\n"
            " ---- or you must set the setIgnoreSideSets() flag " << std::endl;
          throw std::logic_error("Refiner::fixElementSides1: missing sub-dim break pattern - logic error");
          return;
        }

      SameRankRelation& parent_child = m_eMesh.adapt_parent_to_child_relations();

      //std::cout << "tmp parent_child.size() = " << parent_child.size() << std::endl;

      SameRankRelation::iterator pc_it;
      for (pc_it = parent_child.begin(); pc_it != parent_child.end(); pc_it++)
        {
          const SameRankRelationKey& parent = pc_it->first;
          SameRankRelationValue& child_vector = pc_it->second;

          if (0 == &parent)
            {
              throw std::logic_error("Refiner::fixElementSides1 parent is null");
            }

          if (0 == parent)
            {
              throw std::logic_error("Refiner::fixElementSides1 parent is null");
            }

          const CellTopologyData *parent_topo_data = stk::percept::PerceptMesh::get_cell_topology(*parent);
          if (0 == parent_topo_data)
            {
              throw std::logic_error("Refiner::fixElementSides1 parent_topo_data is null");
            }

          shards::CellTopology parent_topo(stk::percept::PerceptMesh::get_cell_topology(*parent));
          //unsigned parent_nsides = (unsigned)parent_topo.getSideCount();

          for (unsigned i_child = 0; i_child < child_vector.size(); i_child++)
            {
              stk::mesh::Entity *child = child_vector[i_child];
              //mesh::PairIterRelation child_sides = child->relations(side_rank);
              if (!child)
                {
                  std::cout << "fixElementSides1: child == null, i_child= " << i_child << " nchild= " << child_vector.size() << std::endl;
                  throw std::runtime_error("fixElementSides1: child == null");
                }

              shards::CellTopology child_topo(stk::percept::PerceptMesh::get_cell_topology(*child));
              unsigned child_nsides = (unsigned)child_topo.getSideCount();

              // if parent has any side relations, check if any of the sides' children match the parent's children's faces
              mesh::PairIterRelation parent_sides = parent->relations(side_rank);
              mesh::PairIterRelation side_to_parent = parent->relations(m_eMesh.element_rank());

              //std::cout << "tmp here 1 child_nsides= " << child_nsides
              //          << " parent_sides.size()=" << parent_sides.size() <<  " side_to_parent.size() = " << side_to_parent.size() << std::endl;

              for (unsigned i_parent_side = 0; i_parent_side < parent_sides.size(); i_parent_side++)
                {
                  stk::mesh::Entity *parent_side = parent_sides[i_parent_side].entity();
                  //unsigned local_parent_side_id = parent_sides[i_parent_side].identifier();

                  if (!parent_side)
                    {
                      throw std::logic_error("parent_side is null");
                    }
                  SameRankRelation& repo = m_eMesh.adapt_parent_to_child_relations();
                  //SameRankRelationValue& parent_side_children = m_eMesh.adapt_parent_to_child_relations()[parent_side];
                  const SameRankRelationValue* parent_side_children_ptr = getChildVectorPtr(repo, parent_side);
                  if (!parent_side_children_ptr)
                    continue;

                  //const SameRankRelationValue& parent_side_children = getChildVector(repo, parent_side);
                  const SameRankRelationValue& parent_side_children = *parent_side_children_ptr;

                  //std::cout << "tmp here 2 parent_side_children.size() = " << parent_side_children.size()
                  //          << std::endl;

                  for (unsigned i_parent_side_child = 0; i_parent_side_child < parent_side_children.size(); i_parent_side_child++)
                    {
                      stk::mesh::Entity *parent_side_child = parent_side_children[i_parent_side_child];

                      //std::cout << "tmp here 3 parent_side_child = " << *parent_side_child
                      //      << std::endl;

                      int permIndex = -1;
                      int permPolarity = 1;

                      // use of i_parent_side here implies that the children's sides match up with the parents, this could be untrue -
                      //  then will require a search through all child faces
                      // NOTE: have to search over child faces due to different topology cases - if parent & child have same topology,
                      //   we can save a few ops here TODO FIXME
                      unsigned k_child_side = 0;

#if 0
                      // FIXME - why is this #if'd out?
                      boolean sameTopology = false; // FIXME - get this from the break pattern
                      if (sameTopology)
                        {
                          PerceptMesh::element_side_permutation(*child, *parent_side_child, k_child_side, permIndex, permPolarity);
                        }
#endif

                      if (permIndex < 0)
                        {
                          // try search
                          for (unsigned j_child_side = 0; j_child_side < child_nsides; j_child_side++)
                            {
                              PerceptMesh::element_side_permutation(*child, *parent_side_child, j_child_side, permIndex, permPolarity);
                              if (0)
                                std::cout << "tmp j_child_side = " << j_child_side << " permIndex= " << permIndex
                                          << " child= " << *child
                                          << " parent_side_child= " << *parent_side_child
                                          <<  std::endl;

                              if (permIndex >= 0)
                                {
                                  k_child_side = j_child_side;
                                  break;
                                }
                            }
                        }

                      if (permIndex >= 0)
                        {
                          if (0)
                            std::cout << "tmp decl rel permIndex= " << permIndex
                                      << " child= " << *child
                                      << " parent_side_child= " << *parent_side_child
                                      <<  std::endl;
                          m_eMesh.getBulkData()->declare_relation(*child, *parent_side_child, k_child_side);
                          PerceptMesh::element_side_permutation(*child, *parent_side_child, k_child_side, permIndex, permPolarity);
                        }
                      else
                        {
                          // error condition?
                          //throw std::runtime_error("fixElementSides1: couldn't find a matching face");
                        }
                    }
                }
            }
        }
    }

#elif PERCEPT_USE_FAMILY_TREE == 1


    void Refiner::
    fixElementSides1(stk::mesh::EntityRank side_rank)
    {
      EXCEPTWATCH;

      bool notFound = true;
      for (unsigned ibp = 0; ibp < m_breakPattern.size(); ibp++)
        {
          // only check the side elements
          if (m_breakPattern[ibp]->getPrimaryEntityRank() == side_rank)
            {
              notFound = false;
            }
        }
      if (notFound)
        {
          std::cout << "Refiner::fixElementSides1: missing sub-dim break pattern - logic error\n"
            " ---- for this refinement pattern to be able to handle sidesets and edgesets you must provide the sub-dim break pattern\n"
            " ---- or you must set the setIgnoreSideSets() flag " << std::endl;
          throw std::logic_error("Refiner::fixElementSides1: missing sub-dim break pattern - logic error");
          return;
        }


      // FIXME
      const unsigned FAMILY_TREE_RANK = m_eMesh.element_rank() + 1u;
      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.getBulkData()->buckets( FAMILY_TREE_RANK );

      //std::cout << "tmp parent_child.size() = " << parent_child.size() << std::endl;

      // loop over all the available parent/child relations by looping over the super-nannies (entities of rank  = 1+element_rank() )
      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          const unsigned num_elements_in_bucket = bucket.size();
          for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
            {
              stk::mesh::Entity& family_tree = bucket[iElement];
              //for (unsigned rank = 1u; rank <= eMesh.element_rank(); rank++)
              //for (unsigned rank = m_eMesh.element_rank(); rank <= m_eMesh.element_rank(); rank++)
              unsigned element_rank = m_eMesh.element_rank();
                {
                  // only look at element rank
                  mesh::PairIterRelation family_tree_relations = family_tree.relations(element_rank);
                  if (family_tree_relations.size())
                    {
                      // get the parent from the family_tree
                      stk::mesh::Entity* parent = family_tree_relations[0].entity();


                      if (0)
                        {
                          std::cout << "tmp family_tree_relations, family_tree id= " << family_tree.identifier() << std::endl;
                          for (unsigned ipc = 0; ipc < family_tree_relations.size(); ipc++)
                            {
                              std::cout << "tmp ipc = " << ipc << " entity_rank = " << family_tree_relations[ipc].entity()->entity_rank()
                                        << " entity= " << *family_tree_relations[ipc].entity()
                                        << std::endl;
                            }
                        }

                      if (0 == parent)
                        {
                          throw std::logic_error("Refiner::fixElementSides1 parent is null");
                        }

                      const CellTopologyData *parent_topo_data = stk::percept::PerceptMesh::get_cell_topology(*parent);
                      if (0 == parent_topo_data)
                        {
                          throw std::logic_error("Refiner::fixElementSides1 parent_topo_data is null");
                        }

                      // skip non-leaf nodes of the tree
                      if (!m_eMesh.isParentElementLeaf(*parent))
                        continue;

                      shards::CellTopology parent_topo(stk::percept::PerceptMesh::get_cell_topology(*parent));
                      //unsigned parent_nsides = (unsigned)parent_topo.getSideCount();

                      // loop over each child associated with parent
                      for (unsigned i_child = 1; i_child < family_tree_relations.size(); i_child++)
                        {
                          stk::mesh::Entity *child = family_tree_relations[i_child].entity();

                          if (!child)
                            {
                              //std::cout << "fixElementSides1: child == null, i_child= " << i_child << " nchild= " << child_vector.size() << std::endl;
                              throw std::runtime_error("fixElementSides1: child == null");
                            }

                          shards::CellTopology child_topo(stk::percept::PerceptMesh::get_cell_topology(*child));
                          unsigned child_nsides = (unsigned)child_topo.getSideCount();

                          // if parent has any side relations, check if any of the sides' children match the parent's children's faces
                          mesh::PairIterRelation parent_sides = parent->relations(side_rank);

                          //??? mesh::PairIterRelation side_to_parent = parent->relations(m_eMesh.element_rank());

                          //std::cout << "tmp here 1 child_nsides= " << child_nsides
                          //          << " parent_sides.size()=" << parent_sides.size() <<  " side_to_parent.size() = " << side_to_parent.size() << std::endl;

                          for (unsigned i_parent_side = 0; i_parent_side < parent_sides.size(); i_parent_side++)
                            {
                              stk::mesh::Entity *parent_side = parent_sides[i_parent_side].entity();
                              //unsigned local_parent_side_id = parent_sides[i_parent_side].identifier();

                              if (!parent_side)
                                {
                                  throw std::logic_error("parent_side is null");
                                }
                              //SameRankRelation& repo = m_eMesh.adapt_parent_to_child_relations();
                              //SameRankRelationValue& parent_side_children = m_eMesh.adapt_parent_to_child_relations()[parent_side];
                              //const SameRankRelationValue* parent_side_children_ptr = getChildVectorPtr(repo, parent_side);

                              mesh::PairIterRelation parent_side_to_family_tree_relations = parent_side->relations(FAMILY_TREE_RANK);

                              if (! parent_side_to_family_tree_relations.size())
                                {
                                  //std::cout << "tmp found parent_side_to_family_tree_relations.size() == 0" << std::endl;
                                  continue;
                                }

                              //std::cout << "tmp found parent_side_to_family_tree_relations.size() != 0, = " << parent_side_to_family_tree_relations.size() << std::endl;
                              stk::mesh::Entity *sn_check = parent_side_to_family_tree_relations[0].entity();

                              // check
                              if (1)
                                {
                                  stk::mesh::Entity *ps_check = sn_check->relations(parent_side->entity_rank())[0].entity();
                                  if (ps_check != parent_side)
                                    {
                                      std::cout << "error: parent_side_to_family_tree_relations broken" << std::endl;
                                      throw std::logic_error( "error: parent_side_to_family_tree_relations broken");
                                    }
                                }

                              mesh::PairIterRelation family_tree_parent_side_relations = sn_check->relations(parent_side->entity_rank());

                              if (! family_tree_parent_side_relations.size())
                                {
                                  //std::cout << "tmp found family_tree_parent_side_relations.size() == 0" << std::endl;
                                  continue;
                                }

                              //std::cout << "tmp found family_tree_parent_side_relations.size() != 0, = " << family_tree_parent_side_relations.size() <<  std::endl;

                              for (unsigned i_parent_side_child = 1; i_parent_side_child < family_tree_parent_side_relations.size(); i_parent_side_child++)
                                {
                                  // the child of the parent's side
                                  stk::mesh::Entity *parent_side_child = family_tree_parent_side_relations[i_parent_side_child].entity();

                                  //std::cout << "tmp here 3 parent_side_child = " << *parent_side_child
                                  //      << std::endl;

                                  int permIndex = -1;
                                  int permPolarity = 1;

                                  // use of i_parent_side here implies that the children's sides match up with the parents, this could be untrue -
                                  //  then will require a search through all child faces
                                  // NOTE: have to search over child faces due to different topology cases - if parent & child have same topology,
                                  //   we can save a few ops here TODO FIXME
                                  unsigned k_child_side = 0;

#if 0
                                  // FIXME - why is this #if'd out?
                                  boolean sameTopology = false; // FIXME - get this from the break pattern
                                  if (sameTopology)
                                    {
                                      PerceptMesh::element_side_permutation(*child, *parent_side_child, k_child_side, permIndex, permPolarity);
                                    }
#endif

                                  if (permIndex < 0)
                                    {
                                      // try search
                                      for (unsigned j_child_side = 0; j_child_side < child_nsides; j_child_side++)
                                        {
                                          PerceptMesh::element_side_permutation(*child, *parent_side_child, j_child_side, permIndex, permPolarity);
                                          if (0)
                                            std::cout << "tmp j_child_side = " << j_child_side << " permIndex= " << permIndex
                                                      << " child= " << *child
                                                      << " parent_side_child= " << *parent_side_child
                                                      <<  std::endl;

                                          if (permIndex >= 0)
                                            {
                                              k_child_side = j_child_side;
                                              break;
                                            }
                                        }
                                    }

                                  if (permIndex >= 0)
                                    {
                                      if (0)
                                        std::cout << "tmp decl rel permIndex= " << permIndex
                                                  << " child= " << *child
                                                  << " parent_side_child= " << *parent_side_child
                                                  <<  std::endl;
                                      m_eMesh.getBulkData()->declare_relation(*child, *parent_side_child, k_child_side);
                                      PerceptMesh::element_side_permutation(*child, *parent_side_child, k_child_side, permIndex, permPolarity);
                                    }
                                  else
                                    {
                                      // error condition?
                                      //throw std::runtime_error("fixElementSides1: couldn't find a matching face");
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
#endif

    void Refiner::
    fixElementSides(stk::mesh::EntityRank side_rank)
    {

    }

#undef EXTRA_PRINT_UR_FES

    /** Sets orientations and associativity of elements to sub-dimensional faces/edges after refinement.
     */
    void Refiner::
    checkFixElementSides(stk::mesh::EntityRank side_rank, stk::mesh::EntityRank elem_rank)
    {
    }

    void Refiner::removeFamilyTrees()
    {
      EXCEPTWATCH;

      elements_to_be_destroyed_type elements_to_be_destroyed;

      const unsigned FAMILY_TREE_RANK = m_eMesh.element_rank() + 1u;
      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.getBulkData()->buckets( FAMILY_TREE_RANK );

      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          //if (removePartSelector(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_elements_in_bucket = bucket.size();

              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  stk::mesh::Entity& element = bucket[iElement];
                  stk::mesh::Entity* element_p = &element;

#if UNIFORM_REF_REMOVE_OLD_STD_VECTOR
                  elements_to_be_destroyed.push_back(element_p);
#else
                  elements_to_be_destroyed.insert(element_p);
#endif
                }
            }
        }
      //std::cout << "tmp P[" << m_eMesh.getRank() << "] removing family_trees, size() = "  << elements_to_be_destroyed.size() << std::endl;
      removeOldElements(elements_to_be_destroyed);
    }

    void Refiner::
    removeOldElements(stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern)
    {
      EXCEPTWATCH;

      const mesh::Part *oldPart = m_eMesh.getPart(breakPattern->getOldElementsPartName()+toString(rank));

      if (1 && oldPart)
        {
          const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(*oldPart);
          std::string ct_name = (cell_topo_data ? cell_topo_data->name : "");
          //std::cout << "tmp removeOldElements::name= " << oldPart->name() << " for rank= " << rank << " topology= " << ct_name << std::endl;
        }

      if (!oldPart)
        {
          std::cout << "name= " << breakPattern->getOldElementsPartName()+toString(rank) << std::endl;
          throw std::runtime_error("oldPart is null");
        }

      mesh::Selector removePartSelector (*oldPart);

      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.getBulkData()->buckets( rank );

      elements_to_be_destroyed_type elements_to_be_destroyed;

#if UNIFORM_REF_REMOVE_OLD_STD_VECTOR
      unsigned nel = 0u;
      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (removePartSelector(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_elements_in_bucket = bucket.size();

              nel += num_elements_in_bucket;
            }
        }
      elements_to_be_destroyed.reserve(nel);
#endif

      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (removePartSelector(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_elements_in_bucket = bucket.size();

              if (0)
                {
                  std::string str;
                  stk::mesh::PartVector pv;
                  bucket.supersets(pv);
                  for (unsigned ip = 0; ip < pv.size(); ip++)
                    {
                      str += " "+pv[ip]->name();
                    }
                  std::cout << "P[" << m_eMesh.getRank() << "] removing elements in bucket of parts: " << str << std::endl;
                }

              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  stk::mesh::Entity& element = bucket[iElement];
                  stk::mesh::Entity* element_p = &element;

                  if (!m_eMesh.isGhostElement(element))
                    {
#if UNIFORM_REF_REMOVE_OLD_STD_VECTOR
                      elements_to_be_destroyed.push_back(element_p);
#else
                      elements_to_be_destroyed.insert(element_p);
#endif
                      //std::cout << "tmp removing elem = " << *element_p << " ";
                      //m_eMesh.printEntity(std::cout, *element_p);
                    }
                }
            }
        }
      removeOldElements(elements_to_be_destroyed);

    }

    void Refiner::removeOldElements(elements_to_be_destroyed_type& elements_to_be_destroyed)
    {
      elements_to_be_destroyed_type elements_to_be_destroyed_pass2;

      for (elements_to_be_destroyed_type::iterator itbd = elements_to_be_destroyed.begin(); itbd != elements_to_be_destroyed.end();  ++itbd)
        {
          stk::mesh::Entity *element_p = *itbd;

          if (0)
            {
              std::cout << "tmp removeOldElements removing element_p = " << element_p << std::endl;
              if (element_p) std::cout << "tmp removeOldElements removing id= " << element_p->identifier() << std::endl;
            }

          if ( ! m_eMesh.getBulkData()->destroy_entity( element_p ) )
            {
#if UNIFORM_REF_REMOVE_OLD_STD_VECTOR
              elements_to_be_destroyed_pass2.push_back(element_p);
#else
              elements_to_be_destroyed_pass2.insert(element_p);
#endif
              //throw std::logic_error("Refiner::removeOldElements couldn't remove element");

            }
        }

      //std::cout << "tmp Refiner::removeOldElements pass2 size = " << elements_to_be_destroyed_pass2.size() << std::endl;
      for (elements_to_be_destroyed_type::iterator itbd = elements_to_be_destroyed_pass2.begin();
           itbd != elements_to_be_destroyed_pass2.end();  ++itbd)
        {
          stk::mesh::Entity *element_p = *itbd;
          if ( ! m_eMesh.getBulkData()->destroy_entity( element_p ) )
            {
              CellTopology cell_topo(stk::percept::PerceptMesh::get_cell_topology(*element_p));
              std::cout << "tmp Refiner::removeOldElements couldn't remove element in pass2,...\n tmp destroy_entity returned false: cell= " << cell_topo.getName() << std::endl;
              const mesh::PairIterRelation elem_relations = element_p->relations(element_p->entity_rank()+1);
              std::cout << "tmp elem_relations.size() = " << elem_relations.size() << std::endl;

              throw std::logic_error("Refiner::removeOldElements couldn't remove element, destroy_entity returned false.");
            }
        }
    }

    /// fix names of surfaces (changing for example surface_hex8_quad4 to surface_tet4_tri3)
    void Refiner::
    fixSurfaceAndEdgeSetNames(stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern)
    {
      EXCEPTWATCH;
      stk::mesh::PartVector toParts = breakPattern->getToParts();

      //std::cout << "toParts.size()= " << toParts.size() << " typeid= " << typeid(*breakPattern).name()  << std::endl;

      for (unsigned i_part = 0; i_part < toParts.size(); i_part++)
        {
          //const std::string & partName = toParts[i_part]->name();
          std::string * toPartName_p = const_cast<std::string *> (&toParts[i_part]->name());

          std::string toPartName = toParts[i_part]->name();
          if ( toPartName.find("surface_", 0) == std::string::npos)
            {
              if (0) std::cout << "tmp fixSurfaceAndEdgeSetNames:: skipping toPartName= " << toPartName << " typeid= " << typeid(*breakPattern).name()  << std::endl;
              continue;
            }

          std::string newToPartName = toPartName;

          StringStringMap::iterator map_it;
          StringStringMap str_map =  breakPattern->fixSurfaceAndEdgeSetNamesMap();
          if (0) std::cout << "tmp fixSurfaceAndEdgeSetNamesMap:: str_map.size()= " << str_map.size()
            //<< " " << breakPattern->getFromTopoPartName() << "__" << breakPattern->getToTopoPartName()
                           << " typeid= " << typeid(*breakPattern).name()
                           << std::endl;

          for (map_it = str_map.begin(); map_it != str_map.end(); map_it++)
            {
              std::string from_str = map_it->first;
              std::string to_str = map_it->second;
              Util::replace(newToPartName, from_str, to_str);
              if (0)
                std::cout << "tmp fixSurfaceAndEdgeSetNamesMap: old= " << toPartName << " new= " << newToPartName << std::endl;
            }

          *toPartName_p = newToPartName;

          if (0)
            std::cout << "tmp fixSurfaceAndEdgeSetNamesMap:: P[" << m_eMesh.getRank() << "] new part name= " << toParts[i_part]->name()
                      << " old part name = " << toPartName
                      << std::endl;
        }
    }

    // FIXME this is a hack to rename parts
    /// Renames as follows:
    ///   originalPartName -> originalPartName_uo_1000    The original part holds the elements to be converted, and is renamed to be the "old" part
    ///   originalPartName_urpconv -> originalPartName    The new part has the same name as the original part with urpconv appended, which
    ///                                                      is then changed back to the original part name
    ///
    /// So, after the renaming, the original part name holds the new elements, and the original elements are
    ///   in the part with the original name appended with _uo_1000.  These parts are ignored on subsequent input.
    ///
    void Refiner::
    renameNewParts(stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern)
    {
      EXCEPTWATCH;
      stk::mesh::PartVector toParts = breakPattern->getToParts();
      stk::mesh::PartVector fromParts = breakPattern->getFromParts();

      if (0)
        {
          for (unsigned i_part = 0; i_part < toParts.size(); i_part++)
            {
              std::cout << "tmp toParts[i_part]->name() = " << toParts[i_part]->name()
                        << " fromParts[i_part]->name() = " << fromParts[i_part]->name()  << std::endl;
            }

        }

      for (unsigned i_part = 0; i_part < toParts.size(); i_part++)
        {
          if (0) std::cout << "tmp before: fromPartName= " << fromParts[i_part]->name() << " toPartName= " << toParts[i_part]->name() << std::endl;

          std::string * toPartName_p = const_cast<std::string *> (&toParts[i_part]->name());

          std::string toPartName = toParts[i_part]->name();
          if (toParts[i_part]->name() == fromParts[i_part]->name())
            {
              continue;
            }
          std::string fromPartName = toPartName;
          int len = fromPartName.length();
          int clen = breakPattern->getAppendConvertString().length();
          fromPartName.erase(len - clen, clen);

          mesh::Part *fromPart = m_eMesh.getNonConstPart(fromPartName);
          std::string * fromPartName_p = const_cast<std::string *> (&fromPart->name());
          *toPartName_p = fromPartName;
          *fromPartName_p = fromPartName + breakPattern->getAppendOriginalString();

          if (0) std::cout << "tmp  after: fromPartName= " << fromParts[i_part]->name() << " toPartName= " << toParts[i_part]->name() << std::endl;

          if (0)
            std::cout << "tmp P[" << m_eMesh.getRank() << "] fromPartName: " << fromPartName << " part= " << toParts[i_part]->name()
                      << " old part name = " << fromPart->name()
                      << std::endl;
        }
    }


    std::vector< RefinementInfoByType >&
    Refiner::
    getRefinementInfoByType()
    {
      return m_refinementInfoByType;
    }

    void
    Refiner::
    setQueryPassOnly(bool doQueryOnly)
    {
      m_doQueryOnly = doQueryOnly;
    }

    // ====================================================================================================
    // ====================================================================================================
    // ====================================================================================================

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
                // skip elements that are already a parent (if there's no family tree yet, it's not a parent, so avoid throwing an error is isParentElement)
                const bool check_for_family_tree = false;  
                bool isParent = m_eMesh.isParentElement(element, check_for_family_tree);
              
                if (isParent)
                  continue;

                const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);

                if (elem_nodes.size() && m_eMesh.isChildWithoutNieces(element, false) )
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

    void Refiner::
    filterUnrefSet(ElementUnrefineCollection& elements_to_unref)
    {
      const unsigned FAMILY_TREE_RANK = m_eMesh.element_rank() + 1u;
      ElementUnrefineCollection elements_to_unref_copy;

      for (ElementUnrefineCollection::iterator u_iter = elements_to_unref.begin();
           u_iter != elements_to_unref.end(); ++u_iter)
        {
          //stk::mesh::Entity * element_ptr = *u_iter;
          stk::mesh::Entity& element = **u_iter;
          
          const bool check_for_family_tree = false;  
          bool isParent = m_eMesh.isParentElement(element, check_for_family_tree);
              
          if (isParent)
            continue;

          const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);

          if (elem_nodes.size() && m_eMesh.isChildWithoutNieces(element, false) )
            {
              //bool elementIsGhost = m_eMesh.isGhostElement(element);

              //std::vector<stk::mesh::Entity *> siblings;
              stk::mesh::PairIterRelation child_to_family_tree_relations = element.relations(FAMILY_TREE_RANK);

              // look for level 0 only - these are children with no children
              unsigned child_ft_level_0 = m_eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, element);

              stk::mesh::Entity *family_tree = child_to_family_tree_relations[child_ft_level_0].entity();
              stk::mesh::PairIterRelation family_tree_relations = family_tree->relations(m_eMesh.element_rank());
              if (family_tree_relations.size() == 0)
                {
                  throw std::logic_error("Refiner::filterUnrefSet family_tree_relations.size() == 0");
                }

              for (unsigned ichild=1; ichild < family_tree_relations.size(); ichild++)
                {
                  stk::mesh::Entity *child = family_tree_relations[ichild].entity();
                  if (m_eMesh.isParentElement(*child))
                    {
                      throw std::logic_error("Refiner::filterUnrefSet isParentElement not expected");
                    }
              
                  elements_to_unref_copy.insert(child);
                }
            }
        }
      //std::cout << "tmp filterUnrefSet::elements_to_unref.size = " << elements_to_unref.size() 
      //          <<  " filtered size= " << elements_to_unref_copy.size() << std::endl;
      elements_to_unref = elements_to_unref_copy;
    }

    void Refiner::
    getKeptNodes(NodeSetType& kept_nodes, ElementUnrefineCollection& elements_to_unref)
    {
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

                //if (m_eMesh.isLeafElement(element) && !m_eMesh.isGhostElement(element))
                {
                  bool in_unref_set = elements_to_unref.find( &element ) != elements_to_unref.end();
                  if (!in_unref_set)
                    {
                      const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);

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
      // remove elements asked to be unrefined if they have a kept node
      // FIXME

      // mark deleted nodes that aren't marked 'kept'
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
      const unsigned FAMILY_TREE_RANK = m_eMesh.element_rank() + 1u;

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

      if (0 && DEBUG_UNREF)
        {
          NodeSetType d1;
          std::cout << "tmp dump for node14: " << std::endl;
          m_nodeRegistry->dumpDB();

          stk::mesh::Entity *node14= m_eMesh.getBulkData()->get_entity(stk::mesh::fem::FEMMetaData::NODE_RANK, 14u);
          std::cout << "node14= " ; 
          m_eMesh.printEntity(std::cout, *node14);
          d1.insert(node14);
          m_nodeRegistry->cleanDeletedNodes(d1, true);
          //exit(123);
        }

      // filter unref set
      filterUnrefSet(elements_to_unref);

      // get kept nodes
      getKeptNodes(kept_nodes, elements_to_unref);

      // get deleted nodes
      getDeletedNodes(deleted_nodes, kept_nodes, elements_to_unref);

      // remove deleted nodes and their associated sub-dim entities
      m_nodeRegistry->cleanDeletedNodes(deleted_nodes);

      // remove elements to be unrefined
      ElementUnrefineCollection copied_children_to_be_removed = elements_to_unref;

#if DEBUG_UNREF
      std::cout << "tmp copied_children_to_be_removed.size() [= num elements to be urefined that are children and !ghosts]= " << copied_children_to_be_removed.size() << std::endl;
#endif

      unsigned nchild_removed = 0;

      typedef std::set<stk::mesh::Entity *> SetOfEntities;
      SetOfEntities family_trees_to_be_removed;
      SetOfEntities children_to_be_removed;
      SetOfEntities children_to_be_removed_with_ghosts;

      ElementUnrefineCollection parent_elements;
      ElementUnrefineCollection parent_elements_copy;

      // remove elements marked for unrefine (make sure they are children and not ghosts)
      for (ElementUnrefineCollection::iterator u_iter = elements_to_unref.begin();
           u_iter != elements_to_unref.end(); ++u_iter)
        {
          stk::mesh::Entity * element_p = *u_iter;
          bool isGhostElement = m_eMesh.isGhostElement(*element_p);
          //if (isGhostElement)
          //  continue;

          if (!m_eMesh.isChildElement(*element_p))
            continue;

#if DEBUG_UNREF
          //std::cout << "tmp element to be removed id= " << element_p->identifier() << " " << std::endl;
#endif
          if (copied_children_to_be_removed.find(element_p) != copied_children_to_be_removed.end())
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
              parent_elements_copy.insert(parent);

              family_trees_to_be_removed.insert(family_tree);

              for (unsigned ichild=0; ichild < siblings.size(); ichild++)
                {
                  stk::mesh::Entity *child = siblings[ichild];

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

      m_nodeRegistry->clear_element_owner_data(children_to_be_removed_with_ghosts);

      for(SetOfEntities::iterator family_tree_it = family_trees_to_be_removed.begin();
          family_tree_it != family_trees_to_be_removed.end(); ++family_tree_it)
        {
          stk::mesh::Entity *family_tree = *family_tree_it;
          if ( ! m_eMesh.getBulkData()->destroy_entity( family_tree ) )
            {
              throw std::logic_error("Refiner::unrefineTheseElements couldn't remove element, destroy_entity returned false for family_tree.");
            }
        }

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

      //m_nodeRegistry->cleanDeletedNodes(deleted_nodes);

      //
      // FIXME refactor to a generic function operating on a collection of elements; incorporate with the doBreak() calls above
      //
      // re-mesh
      if (1)
      {
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

                if (!m_eMesh.isGhostElement(element))
                  {
                    ++num_elem_not_ghost;
                  }
              }

            unsigned num_elem_needed = num_elem_not_ghost * m_breakPattern[irank]->getNumNewElemPerElem();
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
                stk::mesh::Entity *parent = *p_iter;

                const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(*parent);
                CellTopology cell_topo(cell_topo_data);
                unsigned elementType = cell_topo.getKey();
                unsigned bpElementType = m_breakPattern[irank]->getFromTypeKey();
                if (elementType == bpElementType)
                  {
                    stk::mesh::Entity& element = *parent;

                    if (!m_eMesh.isGhostElement(element))
                      {
#if DEBUG_UNREF
                        //std::cout << "P["<< m_eMesh.getRank() << "] element.owner_rank() = " << element.owner_rank() << std::endl;
                        std::cout << "tmp Parent to be remeshed = ";
                        m_eMesh.printEntity(std::cout, element);
#endif
                        if (createNewNeededNodeIds(cell_topo_data, element, needed_entity_ranks, new_sub_entity_nodes))
                          {
                            //std::cout << "typeid= " << typeid(*breakPattern).name() << std::endl;
                            throw std::logic_error("unrefineTheseElements:: createNewNeededNodeIds failed");
                          }

                        vector<stk::mesh::Entity *>::iterator element_pool_it_b4 = element_pool_it;
                        m_breakPattern[irank]->createNewElements(m_eMesh, *m_nodeRegistry, element, new_sub_entity_nodes, element_pool_it, m_proc_rank_field);
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

      m_eMesh.getBulkData()->modification_end();

#if CHECK_DEBUG
      check_db("after unrefineTheseElements");
#endif
    }

    //    ========================================================================================================================
    //    ========================================================================================================================
    //    ========================================================================================================================


    void Refiner::check_db(std::string msg)
    {
      std::cout << "tmp check_db msg= " << msg << std::endl;
      check_db_ownership_consistency();
      //check_db_hanging_nodes();
    }

    void Refiner::check_db_ownership_consistency()
    {
      SubDimCellToDataMap& cell_2_data_map = m_nodeRegistry->getMap();

      for (SubDimCellToDataMap::iterator cell_iter = cell_2_data_map.begin(); cell_iter != cell_2_data_map.end(); ++cell_iter)
        {
          SubDimCellData& nodeId_elementOwnderId = (*cell_iter).second;
          stk::mesh::EntityId owning_elementId = stk::mesh::entity_id(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>());
          NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
          unsigned owning_elementRank = stk::mesh::entity_rank(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>());

          if (nodeIds_onSE.size())
            {

              if (!owning_elementId)
                throw std::logic_error("check_db_ownership_consistency:: error #1");

              stk::mesh::Entity * owning_element = m_eMesh.getBulkData()->get_entity(owning_elementRank, owning_elementId);

              if (!owning_element)
                throw std::logic_error("check_db_ownership_consistency:: error #2");

              if (!m_eMesh.isGhostElement(*owning_element))
                {
                
                  for (unsigned inode = 0; inode < nodeIds_onSE.size(); inode++)
                    {
                      stk::mesh::Entity *node = nodeIds_onSE[inode];
                      if (!node)
                        throw std::logic_error("check_db_ownership_consistency:: error #3");

                      stk::mesh::Entity * node1 = m_eMesh.getBulkData()->get_entity(stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE.m_entity_id_vector[inode]);
                      if (!node1)
                        throw std::logic_error("check_db_ownership_consistency:: error #3a");

                      stk::mesh::Entity * node2 = m_eMesh.getBulkData()->get_entity(stk::mesh::fem::FEMMetaData::NODE_RANK, node->identifier() );
                      if (!node2)
                        throw std::logic_error("check_db_ownership_consistency:: error #3b");
                      if (node != node2)
                        throw std::logic_error("check_db_ownership_consistency:: error #3c");
              
                    }
                }
            }
        }
    }

    void Refiner::check_db_hanging_nodes()
    {
      std::set<stk::mesh::Entity *> node_set;

      // check for hanging nodes - ensure all parents have their sub-entities in the DB
      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.getBulkData()->buckets( m_eMesh.element_rank() );

      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;

          const unsigned num_entity_in_bucket = bucket.size();
          for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
            {
              stk::mesh::Entity& element = bucket[ientity];
              if (m_eMesh.hasFamilyTree(element) && m_eMesh.isParentElement(element, false))
                {
                  for (unsigned irank=0; irank < m_ranks.size(); irank++)
                    {
                      vector<NeededEntityType> needed_entity_ranks;
                      m_breakPattern[irank]->fillNeededEntities(needed_entity_ranks);

                      const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);
                      CellTopology cell_topo(cell_topo_data);
                      unsigned elementType = cell_topo.getKey();
                      unsigned bpElementType = m_breakPattern[irank]->getFromTypeKey();
                      if (elementType == bpElementType)
                        {
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
                              else if (needed_entity_rank == m_eMesh.element_rank())
                                {
                                  numSubDimNeededEntities = 1;
                                }

                              for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                                {
                                  static SubDimCellData empty_SubDimCellData;
                                  SubDimCell_SDSEntityType subDimEntity;
                                  m_nodeRegistry->getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);

                                  SubDimCellData* nodeId_elementOwnderId_ptr = m_nodeRegistry->getFromMapPtr(subDimEntity);
                                  SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
                                  bool is_empty = nodeId_elementOwnderId_ptr == 0;
                                  if (!is_empty)
                                    {
                                      NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
                                      if (nodeIds_onSE.size() == 0)
                                        {
                                          if (1) std::cout << "error check_db_hanging_nodes  irank = " << irank << " ranks[irank] = " << m_ranks[irank]
                                                           << " elementType= " << elementType
                                                           << " cell_topo= " << cell_topo.getName()
                                                           << std::endl;

                                          throw std::logic_error("check_db_hanging_nodes:: error #1");
                                        }
                                    }
                              
                                  /*if (nodeIds_onSE.size() != 0)
                                    {
                                    //std::cout << "tmp" << std::endl;
                                  
                                    }
                                  */
                          
                                }
                            }
                        }
                    }
                }
            }
        }
    }


  } // namespace adapt
} // namespace stk
