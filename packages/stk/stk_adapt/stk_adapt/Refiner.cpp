#include <exception>
#include <fstream>
#include <set>
#include <typeinfo>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <stk_adapt/Refiner.hpp>


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
      m_doProgress(true && (0 == eMesh.getRank()) )
    {
      bp.setSubPatterns(m_breakPattern, eMesh);
      m_nodeRegistry = new NodeRegistry (m_eMesh);
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
    addOldElementsToPart(stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern, unsigned *elementType)
    {
      EXCEPTWATCH;
      //m_eMesh.getBulkData()->modification_begin();
      std::string oldPartName = breakPattern->getOldElementsPartName()+toString(rank);
      mesh::Part *oldPart = m_eMesh.getFEM_meta_data()->get_part(oldPartName);
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
          //std::cout << "tmp element = " << *elems[ielem] << std::endl;
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

      std::vector<stk::mesh::EntityRank> ranks;

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

          m_nodeRegistry->initialize();                           /**/  TRACE_PRINT("Refiner: beginRegistration (top-level rank)... ");

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
                unsigned num_elem_not_ghost_0_incr = doForAllElements(ranks[irank], &NodeRegistry::registerNeedNewNode, elementColors, elementType, needed_entity_ranks);

                num_elem_not_ghost_0 += num_elem_not_ghost_0_incr;
              }
            }

          m_nodeRegistry->endRegistration();                    /**/  TRACE_PRINT("Refiner: endRegistration (top-level rank)... ");
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

                num_elem = doForAllElements(ranks[irank], &NodeRegistry::checkForRemote, elementColors, elementType, needed_entity_ranks, false, false);
              }
            }
          m_nodeRegistry->endCheckForRemote();                /**/   TRACE_PRINT("Refiner: endCheckForRemote (top-level rank)... ");

          if (0)
            {
              std::cout << "num_elem= " << num_elem << std::endl;
              MPI_Barrier( MPI_COMM_WORLD );
              std::cout << "P["<< m_eMesh.getRank() 
                        <<"] ========================================================================================================================" << std::endl;
              m_nodeRegistry->checkDB();
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

                num_elem = doForAllElements(ranks[irank], &NodeRegistry::getFromRemote, elementColors, elementType, needed_entity_ranks, false, false);
              }
            }

          m_nodeRegistry->endGetFromRemote();                    /**/  TRACE_PRINT("Refiner: endGetFromRemote (top-level rank)... ");

          //stk::diag::printTimersTable(std::cout, perceptTimer(), stk::diag::METRICS_ALL, false);

          if (0)
            {
              std::cout << "num_elem= " << num_elem << std::endl;
              MPI_Barrier( MPI_COMM_WORLD );
              std::cout << "P["<< m_eMesh.getRank() 
                        <<"] ========================================================================================================================" << std::endl;
              m_nodeRegistry->checkDB();
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
            /**/                                                TRACE_PRINT("Refiner: registerNeedNewNode count_only(true) ranks[irank]==ranks[0]... ");
            unsigned num_elem_not_ghost = doForAllElements(ranks[irank], &NodeRegistry::registerNeedNewNode, elementColors, elementType, needed_entity_ranks, count_only);
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

          /**/                                                TRACE_PRINT("Refiner: addToExistingParts [etc.]... ");
#if !STK_ADAPT_URP_LOCAL_NODE_COMPS
          if (0)
          {
            EXCEPTWATCH;
            //if (ranks[irank] == ranks[0])
            // only need to do this once: the map is fully built and we loop over the map's faces/edges, which are fixed after the getFromRemote step
            if (irank == 0) 
              {
                if (0)
                  {
                    shards::CellTopology cell_topo(m_breakPattern[irank]->getFromTopology());

                    if (1) std::cout << "tmp Refiner:: calling addToExistingPartsNew() irank = " << irank << " ranks[irank] = " << ranks[irank] 
                                     << " cell_topo= " << cell_topo.getName()
                                     << std::endl;
                  }

                m_nodeRegistry->addToExistingPartsNew();
                //std::cout << "tmp makeCentroid... " << std::endl;
                m_nodeRegistry->makeCentroid(m_eMesh.getCoordinatesField());
                //std::cout << "tmp makeCentroid...done " << std::endl;
                //std::cout << "tmp interpolateFields... " << std::endl;
                //FIXME m_nodeRegistry->interpolateFields();
                //std::cout << "tmp interpolateFields...done " << std::endl;
              }
          }
#endif
          /**/                                                TRACE_PRINT("Refiner: addToExistingParts [etc.] ...done ");

          // this is for testing removing old elements as early as possible for memory reasons
          // FIXME - remove old elements on the fly?
          if (0 && m_doRemove)
            {
              EXCEPTWATCH;

              /**/                                                TRACE_PRINT( "Refiner: remove old elements...start " );

              removeOldElements(ranks[irank], m_breakPattern[irank]);
              renameNewParts(ranks[irank], m_breakPattern[irank]);
              fixSurfaceAndEdgeSetNames(ranks[irank], m_breakPattern[irank]);

              /**/                                                TRACE_PRINT( "Refiner: remove old elements...done " );
            }

          if (TRACE_STAGE_PRINT && !m_eMesh.getRank()) {
            Util::trace_cpu_time_and_mem_print(CONNECT_LOCAL, "CONNECT_LOCAL");
            Util::trace_cpu_time_and_mem_print(CONNECT_LOCAL_createNewNeededNodes, "CONNECT_LOCAL_createNewNeededNodes");
            Util::trace_cpu_time_and_mem_print(CONNECT_LOCAL_createNewElements, "CONNECT_LOCAL_createNewElements");
            Util::trace_cpu_time_and_mem_print(CONNECT_LOCAL_URP_createOrGetNode, "CONNECT_LOCAL_URP_createOrGetNode");
            Util::trace_cpu_time_and_mem_print(CONNECT_LOCAL_URP_declare_relation, "CONNECT_LOCAL_URP_declare_relation");
          }

          /**/                                                TRACE_PRINT("Refiner: modification_end...start... ");
          //FIXME FIXME FIXME 
          //bulkData.modification_end();
          /**/                                                TRACE_PRINT("Refiner: modification_end...done ");

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
        }
#endif

      /**/                                                TRACE_PRINT( "Refiner:doBreak ... done");

      //std::cout << "tmp m_nodeRegistry.m_gee_cnt= " << m_nodeRegistry->m_gee_cnt << std::endl;
      //std::cout << "tmp m_nodeRegistry.m_gen_cnt= " << m_nodeRegistry->m_gen_cnt << std::endl;

    } // doBreak


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

              if (!nodeIds_onSE[0]) {
                
                if (nodeIds_onSE.m_entity_id_vector[0] == 0)
                  {
                    continue;
                    //throw std::logic_error("Refiner::createNewNeededNodeIds logic err #5.0, nodeIds_onSE.m_entity_id_vector[i_new_node] == 0");
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
                          continue;
                          //throw std::logic_error("Refiner::createNewNeededNodeIds logic err #5, nodeIds_onSE.m_entity_id_vector[i_new_node] == 0");
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
              unsigned rank = m_eMesh.element_rank();
                {
                  // only look at element rank
                  mesh::PairIterRelation family_tree_relations = family_tree.relations(rank);
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
                      //std::cout << "tmp removing elem = " << *element_p << std::endl;
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
      
    void 
    Refiner::
    unrefineTheseElements(ElementUnrefineCollection& elements_to_unref)
    {
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

    }


  } // namespace adapt
} // namespace stk
