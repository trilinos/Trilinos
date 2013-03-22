#include <exception>
#include <fstream>
#include <set>
#include <typeinfo>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <stk_adapt/Refiner.hpp>

#include <stk_percept/MeshUtil.hpp>

#include <stk_adapt/NodeRegistryDef.hpp>

#if defined( STK_PERCEPT_HAS_GEOMETRY )
#include <stk_percept/mesh/geometry/kernel/GeometryKernelOpenNURBS.hpp>
#include <stk_percept/mesh/geometry/kernel/MeshGeometry.hpp>
#include <stk_percept/mesh/geometry/kernel/GeometryFactory.hpp>
#endif

#if defined ( STK_PERCEPT_HAS_MESQUITE )
#define StackTraceTmp StackTrace
#undef StackTrace
#include <stk_percept/mesh/mod/mesquite-interface/PerceptMesquiteMesh.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PerceptMesquiteMeshDomain.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PMMLaplaceSmoother.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PMMLaplaceSmoother1.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PMMShapeImprover.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PMMParallelShapeImprover.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PMMShapeSizeOrientImprover.hpp>
#define StackTrace StackTraceTmp

#endif

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
      m_doProgress(true && (0 == eMesh.get_rank()) ),
      m_alwaysInitNodeRegistry(true),
      m_doSmoothGeometry(true)
    {
      bp.setSubPatterns(m_breakPattern, eMesh);
      m_nodeRegistry = new NodeRegistry (m_eMesh);
      m_nodeRegistry->initialize();
      m_nodeRegistry->init_comm_all();
      m_allocated = false;
    }
    
    //NLM new constructor for Python users
    //takes in a Pattern refine_type to determine which UniformRefinerPatternBase will be used
    
    Refiner::Refiner(percept::PerceptMesh& eMesh, Pattern refine_type, stk::mesh::FieldBase *proc_rank_field) :
      m_eMesh(eMesh), m_breakPattern(),
      m_nodeRegistry(0),
      m_proc_rank_field(proc_rank_field), m_doRemove(true), m_ranks(), m_ignoreSideSets(false),
      m_geomFile(""), m_geomSnap(false),
      m_doQueryOnly(false),
      m_progress_meter_frequency(20),
      m_doProgress(true && (0 == eMesh.get_rank()) ),
      m_alwaysInitNodeRegistry(true),
      m_doSmoothGeometry(true)
    {
      m_nodeRegistry = new NodeRegistry (m_eMesh);
      m_nodeRegistry->initialize();
      m_nodeRegistry->init_comm_all();
      m_allocated = true;
      UniformRefinerPatternBase* bp;
      switch(refine_type)
        {
        case LINE2_LINE2_2: 
          bp = (UniformRefinerPatternBase*) new Line2_Line2_2(eMesh); 
          break;
        case BEAM2_BEAM2_2: 
          bp = (UniformRefinerPatternBase*) new Beam2_Beam2_2(eMesh); 
          break;
        case SHELLLINE2_SHELLLINE2_2:
          bp = (UniformRefinerPatternBase*) new ShellLine2_ShellLine2_2(eMesh);
          break;
        case SHELLLINE3_SHELLLINE3_2:
          bp = (UniformRefinerPatternBase*) new ShellLine3_ShellLine3_2(eMesh);
          break;
        case QUAD4_QUAD4_4_OLD:
          bp = (UniformRefinerPatternBase*) new Quad4_Quad4_4_Old(eMesh);
          break;
        case QUAD4_QUAD4_4:
          bp = (UniformRefinerPatternBase*) new Quad4_Quad4_4(eMesh);
          break;
        case QUAD4_QUAD4_4_SIERRA:
          bp = (UniformRefinerPatternBase*) new Quad4_Quad4_4_Sierra(eMesh);
          break;
        case TRI3_TRI3_4:
          bp = (UniformRefinerPatternBase*) new Tri3_Tri3_4(eMesh);
          break;
        case SHELLTRI3_SHELLTRI3_4:
          bp = (UniformRefinerPatternBase*) new ShellTri3_ShellTri3_4(eMesh);
          break;
        case SHELLTRI6_SHELLTRI6_4:
          bp = (UniformRefinerPatternBase*) new ShellTri6_ShellTri6_4(eMesh);
          break;
        case SHELLQUAD4_SHELLQUAD4_4:
          bp = (UniformRefinerPatternBase*) new ShellQuad4_ShellQuad4_4(eMesh);
          break;
        case SHELLQUAD8_SHELLQUAD8_4:
          bp = (UniformRefinerPatternBase*) new ShellQuad8_ShellQuad8_4(eMesh);
          break;
        case TET4_TET4_8:
          bp = (UniformRefinerPatternBase*) new Tet4_Tet4_8(eMesh);
          break;
        case HEX8_HEX8_8:
          bp = (UniformRefinerPatternBase*) new Hex8_Hex8_8(eMesh);
          break;
        case WEDGE6_WEDGE6_8:
          bp = (UniformRefinerPatternBase*) new Wedge6_Wedge6_8(eMesh);
          break;
        case PYRAMID5_PYRAMID5_10:
          bp = (UniformRefinerPatternBase*) new Pyramid5_Pyramid5_10(eMesh);
          break;
        case LINE3_LINE3_2:
          bp = (UniformRefinerPatternBase*) new Line3_Line3_2(eMesh);
          break;
        case BEAM3_BEAM3_2:
          bp = (UniformRefinerPatternBase*) new Beam3_Beam3_2(eMesh);
          break;
        case TRI6_TRI6_4:
          bp = (UniformRefinerPatternBase*) new Tri6_Tri6_4(eMesh);
          break;
        case QUAD9_QUAD9_4:
          bp = (UniformRefinerPatternBase*) new Quad9_Quad9_4(eMesh);
          break;
        case QUAD8_QUAD8_4:
          bp = (UniformRefinerPatternBase*) new Quad8_Quad8_4(eMesh);
          break;
        case HEX27_HEX27_8:
          bp = (UniformRefinerPatternBase*) new Hex27_Hex27_8(eMesh);
          break;
        case HEX20_HEX20_8:
          bp = (UniformRefinerPatternBase*) new Hex20_Hex20_8(eMesh);
          break;
        case TET10_TET10_8:
          bp = (UniformRefinerPatternBase*) new Tet10_Tet10_8(eMesh);
          break;
        case WEDGE15_WEDGE15_8:
          bp = (UniformRefinerPatternBase*) new Wedge15_Wedge15_8(eMesh);
          break;
        case WEDGE18_WEDGE18_8:
          bp = (UniformRefinerPatternBase*) new Wedge18_Wedge18_8(eMesh);
          break;
        case PYRAMID13_PYRAMID13_10:
          bp = (UniformRefinerPatternBase*) new Pyramid13_Pyramid13_10(eMesh);
          break;
        case QUAD4_QUAD9_1:
          bp = (UniformRefinerPatternBase*) new Quad4_Quad9_1(eMesh);
          break;
        case QUAD4_QUAD8_1:
          bp = (UniformRefinerPatternBase*) new Quad4_Quad8_1(eMesh);
          break;
        case BEAM2_BEAM3_1:
          bp = (UniformRefinerPatternBase*) new Beam2_Beam3_1(eMesh);
          break;
        case SHELLQUAD4_SHELLQUAD8_1:
          bp = (UniformRefinerPatternBase*) new ShellQuad4_ShellQuad8_1(eMesh);
          break;
        case TRI3_TRI6_1:
          bp = (UniformRefinerPatternBase*) new Tri3_Tri6_1(eMesh);
          break;
        case TET4_TET10_1:
          bp = (UniformRefinerPatternBase*) new Tet4_Tet10_1(eMesh);
          break;
        case HEX8_HEX27_1:
          bp = (UniformRefinerPatternBase*) new Hex8_Hex27_1(eMesh);
          break;
        case HEX8_HEX20_1:
          bp = (UniformRefinerPatternBase*) new Hex8_Hex20_1(eMesh);
          break;
        case WEDGE6_WEDGE15_1:
          bp = (UniformRefinerPatternBase*) new Wedge6_Wedge15_1(eMesh);
          break;
        case WEDGE6_WEDGE18_1:
          bp = (UniformRefinerPatternBase*) new Wedge6_Wedge18_1(eMesh);
          break;
        case PYRAMID5_PYRAMID13_1:
          bp = (UniformRefinerPatternBase*) new Pyramid5_Pyramid13_1(eMesh);
          break;
        case QUAD4_TRI3_2:
          bp = (UniformRefinerPatternBase*) new Quad4_Tri3_2(eMesh);
          break;
        case QUAD4_TRI3_4:
          bp = (UniformRefinerPatternBase*) new Quad4_Tri3_4(eMesh);
          break;
        case QUAD4_TRI3_6:
          bp = (UniformRefinerPatternBase*) new Quad4_Tri3_6(eMesh);
          break;
        case HEX8_TET4_24:
          bp = (UniformRefinerPatternBase*) new Hex8_Tet4_24(eMesh);
          break;
        case HEX8_TET4_6_12:
          bp = (UniformRefinerPatternBase*) new Hex8_Tet4_6_12(eMesh);
          break;
        default:
          throw std::runtime_error("Refiner::Refiner Unrecognized refinement pattern");
          break;
        }

      bp->setSubPatterns(m_breakPattern, eMesh);
      
      m_proc_rank_field = proc_rank_field;
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
      mesh::Part *oldPart = m_eMesh.get_fem_meta_data()->get_part(oldPartName);
      //std::cout << "tmp removeFromOldPart:: oldPartName= " << oldPartName << std::endl;
      if (!oldPart)
        {
          std::cout << "oldPartName= " << oldPartName << std::endl;
          throw std::runtime_error("oldpart is null");
        }

      mesh::PartVector remove_parts(1, oldPart);
      mesh::PartVector add_parts;
      mesh::Selector on_locally_owned_part =  ( m_eMesh.get_fem_meta_data()->locally_owned_part() );

      std::vector<stk::mesh::Entity*> elems;
      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( rank );

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
          m_eMesh.get_bulk_data()->change_entity_parts( *elems[ielem], add_parts, remove_parts );
        }
    }

    void Refiner::
    addOldElementsToPart(stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern, unsigned *elementType)
    {
      EXCEPTWATCH;
      //m_eMesh.get_bulk_data()->modification_begin();
      std::string oldPartName = breakPattern->getOldElementsPartName()+toString(rank);
      mesh::Part *oldPart = m_eMesh.get_fem_meta_data()->get_part(oldPartName);

#define DEBUG_REMOVE_OLD_PARTS 0

      if (DEBUG_REMOVE_OLD_PARTS) std::cout << "tmp addOldElementsToPart:: oldPartName= " << oldPartName << std::endl;

      if (!oldPart)
        {
          std::cout << "oldPartName= " << oldPartName << std::endl;
          throw std::runtime_error("oldpart is null");
        }

      mesh::PartVector add_parts(1, oldPart);
      mesh::PartVector remove_parts;
      mesh::Selector on_locally_owned_part =  ( m_eMesh.get_fem_meta_data()->locally_owned_part() );

      // The list of Parts that this break pattern will refine.  Only remove elements belonging to these parts.
      mesh::Selector fromPartsSelector = mesh::selectUnion( breakPattern->getFromParts() );

      std::vector<stk::mesh::Entity*> elems;
      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( rank );

      unsigned nele=0;
      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (on_locally_owned_part(**k) && fromPartsSelector(**k) )
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
                      if (DEBUG_REMOVE_OLD_PARTS) std::cout << "tmp adding to oldParts = " << element << std::endl;
                    }
                }
            }
        }


      for (unsigned ielem=0; ielem < elems.size(); ielem++)
        {
          if (DEBUG_REMOVE_OLD_PARTS) std::cout << "tmp addOldElementsToPart element = " << *elems[ielem] << std::endl;
          m_eMesh.get_bulk_data()->change_entity_parts( *elems[ielem], add_parts, remove_parts );
        }

      //m_eMesh.get_bulk_data()->modification_end();
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

      m_eMesh.get_bulk_data()->modification_begin();
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
      m_eMesh.get_bulk_data()->modification_end();

    }

    void Refiner::
    refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity& element, vector<NeededEntityType>& needed_entity_ranks)
    {
      const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);

      shards::CellTopology cell_topo(cell_topo_data);
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

              (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, true);

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

      CommDataType buffer_entry;

      stk::mesh::BulkData& bulkData = *m_eMesh.get_bulk_data();
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

          unsigned num_orig_nodes=0;
          {
            std::vector<unsigned> count1 ;

            stk::mesh::count_entities(stk::mesh::Selector(m_eMesh.get_fem_meta_data()->universal_part()) , *m_eMesh.get_bulk_data(), count1 );
            if (count1.size() < 3)
              {
                throw std::logic_error("logic error in Refiner m_refinementInfoByType");
              }
            num_orig_nodes = count1[0];
          }

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

              mesh::Selector selector(m_eMesh.get_fem_meta_data()->locally_owned_part());
              if (fromPartsAll.size())
                {
                  selector = mesh::Selector();
                  for (unsigned ipart = 0; ipart < fromPartsAll.size(); ipart++)
                    {
                      mesh::Part *part = fromPartsAll[ipart];
                      const CellTopologyData * part_cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(*part);
                      if (part_cell_topo_data)
                        {
                          shards::CellTopology part_cell_topo(part_cell_topo_data);
                          if (part_cell_topo.getKey() == elementType)
                            {
                              selector = selector | *part;
                            }
                        }
                    }
                  selector = selector & (mesh::Selector(m_eMesh.get_fem_meta_data()->locally_owned_part()));
                }
              std::vector<unsigned> count ;
              stk::mesh::count_entities( selector, *m_eMesh.get_bulk_data(), count );
              if (count.size() < 3)
                {
                  throw std::logic_error("logic error in Refiner m_refinementInfoByType");
                }
              unsigned n_ele = count[ ranks[irank] ];

              m_refinementInfoByType[irank].m_rank = ranks[irank];
              m_refinementInfoByType[irank].m_numOrigElems = n_ele;

              m_refinementInfoByType[irank].m_numNewElems = n_ele * m_breakPattern[irank]->getNumNewElemPerElem();
              m_refinementInfoByType[irank].m_topology = cell_topo;
              m_refinementInfoByType[irank].m_numOrigNodes = num_orig_nodes;
              m_refinementInfoByType[irank].m_numNewNodes = 0; // can't predict this
            }

          // sum info from all procs
          {
            stk::ParallelMachine pm = m_eMesh.get_bulk_data()->parallel();

            for (unsigned irank = 0; irank < ranks.size(); irank++)
              {
                stk::all_reduce( pm, stk::ReduceSum<1>( &m_refinementInfoByType[irank].m_numOrigElems ) );
                stk::all_reduce( pm, stk::ReduceSum<1>( &m_refinementInfoByType[irank].m_numNewElems ) );
                stk::all_reduce( pm, stk::ReduceSum<1>( &m_refinementInfoByType[irank].m_numOrigNodes ) );
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

                unsigned num_elem_not_ghost_0_incr = doForAllElements(irank, "Register New Nodes",
                                                                      ranks[irank], &NodeRegistry::registerNeedNewNode, elementColors, elementType, needed_entity_ranks,
                                                                      count_only, doAllElements);

                num_elem_not_ghost_0 += num_elem_not_ghost_0_incr;

                if (0) std::cout << "tmp irank= " << irank << " ranks[irank]= " << ranks[irank] << " nodeRegistry size= " << m_nodeRegistry->getMap().size() << std::endl;
                
              }
            }

          m_nodeRegistry->endRegistration();                    /**/  TRACE_PRINT("Refiner: endRegistration (top-level rank)... ");
        }
        m_nodeRegistry->dumpDB("after registration");

#define CHECK_DEBUG 0
        if (CHECK_DEBUG)
          {
            MPI_Barrier( MPI_COMM_WORLD );
            std::cout << "P["<< m_eMesh.get_rank()
                      <<"] ========================================================================================================================" << std::endl;
            m_nodeRegistry->checkDB("after registerNeedNewNode");
            check_db("after registerNeedNewNode");
            MPI_Barrier( MPI_COMM_WORLD );
            std::cout << "P["<< m_eMesh.get_rank()
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
                num_elem = doForAllElements(irank, "Check For Comm", ranks[irank], &NodeRegistry::checkForRemote, elementColors, elementType, needed_entity_ranks, count_only, doAllElements);
              }
            }
          m_nodeRegistry->endCheckForRemote();                /**/   TRACE_PRINT("Refiner: endCheckForRemote (top-level rank)... ");

          if (1 && CHECK_DEBUG)
            {
              std::cout << "num_elem= " << num_elem << std::endl;
              MPI_Barrier( MPI_COMM_WORLD );
              std::cout << "P["<< m_eMesh.get_rank()
                        <<"] ========================================================================================================================" << std::endl;
              m_nodeRegistry->checkDB("after checkForRemote");
              check_db("after checkForRemote");
              MPI_Barrier( MPI_COMM_WORLD );
              std::cout << "P["<< m_eMesh.get_rank()
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
                num_elem = doForAllElements(irank, "Get From Remote", ranks[irank], &NodeRegistry::getFromRemote, elementColors, elementType, needed_entity_ranks,  count_only, doAllElements);
              }
            }

          m_nodeRegistry->endGetFromRemote();                    /**/  TRACE_PRINT("Refiner: endGetFromRemote (top-level rank)... ");
          m_nodeRegistry->dumpDB("after endGetFromRemote");

          //stk::diag::printTimersTable(std::cout, perceptTimer(), stk::diag::METRICS_ALL, false);

          if (CHECK_DEBUG)
            {
              std::cout << "num_elem= " << num_elem << std::endl;
              MPI_Barrier( MPI_COMM_WORLD );
              std::cout << "P["<< m_eMesh.get_rank()
                        <<"] ========================================================================================================================" << std::endl;
              m_nodeRegistry->checkDB("end getFromRemote");
              check_db("end getFromRemote");
              MPI_Barrier( MPI_COMM_WORLD );
              
              std::cout << "P["<< m_eMesh.get_rank()
                        <<"] ========================================================================================================================" << std::endl;
            }
        }  // get from remote
      } // start top-level ranks

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // for each element type, in top-down rank order, do the rest of the refinement operations
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      bulkData.modification_begin();

      if (CHECK_DEBUG)
        {
          m_nodeRegistry->checkDB("after mod begin() after end getFromRemote");
        }

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
            unsigned num_elem_not_ghost = doForAllElements(irank, "Register New Nodes", ranks[irank], &NodeRegistry::registerNeedNewNode, elementColors, elementType, needed_entity_ranks, count_only, doAllElements);
            /**/                                                TRACE_PRINT("Refiner: registerNeedNewNode count_only(true) ranks[irank]==ranks[0]... done ");

            unsigned num_elem_needed = num_elem_not_ghost * m_breakPattern[irank]->getNumNewElemPerElem();

            // FIXME TMP
            if (CHECK_DEBUG)
              { 
                unsigned nele_col=0;
                for (unsigned icolor = 0; icolor < elementColors.size(); icolor++)
                  {
                    nele_col += elementColors[icolor].size();
                  }
                std::cout << "tmp Refiner::doBreak: irank= " << irank << " ranks[irank]= " << ranks[irank] << " bp= [" 
                          << m_breakPattern[irank]->getFromTopoPartName() << " ==> "
                          << m_breakPattern[irank]->getToTopoPartName() << "] num_elem_needed= " << num_elem_needed 
                          << " num_elem_not_ghost= " << num_elem_not_ghost
                          << " nelementColors= " << nele_col << " elementColors.size()= " << elementColors.size()
                          << std::endl;
              }

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
          if (CHECK_DEBUG)
            {
              m_nodeRegistry->checkDB("after endLocalMeshMods");
            }

          /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          ///  Global element ops: here's where we e.g. connect the new elements by declaring new relations
          /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          /**/                                                TRACE_PRINT("Refiner: createElementsAndNodesAndConnectLocal... ");
          /**/                                                TRACE_CPU_TIME_AND_MEM_0(CONNECT_LOCAL);

          createElementsAndNodesAndConnectLocal(irank, ranks[irank], m_breakPattern[irank], elementColors, needed_entity_ranks, new_elements);

          /**/                                                TRACE_CPU_TIME_AND_MEM_1(CONNECT_LOCAL);
          /**/                                                TRACE_PRINT("Refiner: createElementsAndNodesAndConnectLocal...done ");

          /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          ///  Global node loop operations:  this is where we perform ops like adding new nodes to the right parts, interpolating fields, etc.
          /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

          if (TRACE_STAGE_PRINT && !m_eMesh.get_rank()) {
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
          m_nodeRegistry->makeCentroid(m_eMesh.get_coordinates_field());
          //std::cout << "tmp makeCentroid...done " << std::endl;
          //std::cout << "tmp interpolateFields... " << std::endl;
          m_nodeRegistry->interpolateFields();
          //std::cout << "tmp interpolateFields...done " << std::endl;
        }
      //std::cout << "tmp dump_elements 1" << std::endl;
      // m_eMesh.dump_elements();
#endif
      /**/                                                TRACE_PRINT("Refiner: addToExistingParts [etc.] ...done ");

      /***********************/                           TRACE_PRINT("Refiner: fixElementSides1 ");
      fixElementSides1();
      m_eMesh.adapt_parent_to_child_relations().clear();
      /***********************/                           TRACE_PRINT("Refiner: fixElementSides1...done ");

      //std::cout << "tmp dump_elements 2" << std::endl;
      //m_eMesh.dump_elements();

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

          for (unsigned irank = 0; irank < ranks.size(); irank++)
            {

#if PERCEPT_USE_FAMILY_TREE
              if (irank == 0)
                removeFamilyTrees();
#endif
              removeOldElements(irank, ranks[irank], m_breakPattern[irank]);
              renameNewParts(ranks[irank], m_breakPattern[irank]);
              fixSurfaceAndEdgeSetNames(ranks[irank], m_breakPattern[irank]);
            }
        }

      // remove any elements that are empty (these can exist when doing local refinement)
      removeEmptyElements();

      // remove nodes not referred to by elements
      removeDanglingNodes();

      // remove pseudo elements
      //m_nodeRegistry->removePseudoEntities();

      set_active_part();

      /**/                                                TRACE_PRINT("Refiner: modification_end...start... ");
      bulkData.modification_end();
      // force a flush of all pending deletes, etc
      bulkData.modification_begin();
      bulkData.modification_end();
      /**/                                                TRACE_PRINT("Refiner: modification_end...done ");

      // remove pseudo elements
      if (0)
        {
          bulkData.modification_begin();
          // FIXME  - remove only those that are not shared?
          //removeNodesWithOnlyPseudoNodeRelations();
          m_nodeRegistry->removePseudoEntities();
          bulkData.modification_end();
        }

      //std::cout << "tmp dump_elements 3" << std::endl;
      //m_eMesh.dump_elements();

      snapAndSmooth(m_geomSnap, m_geomFile);

      /**/                                                TRACE_PRINT( "Refiner:doBreak ... done");

      //std::cout << "tmp m_nodeRegistry.m_gee_cnt= " << m_nodeRegistry->m_gee_cnt << std::endl;
      //std::cout << "tmp m_nodeRegistry.m_gen_cnt= " << m_nodeRegistry->m_gen_cnt << std::endl;
      RefinementInfoByType::countCurrentNodes(m_eMesh, getRefinementInfoByType());

      getNodeRegistry().init_entity_repo();
      //getNodeRegistry().clear_dangling_elements();

      m_nodeRegistry->dumpDB("after doBreak");
#if CHECK_DEBUG
      //check_db("after doBreak");
      m_nodeRegistry->checkDB("after doBreak");
#endif

    } // doBreak
    
    // FIXME - temp until we figure out what to do with parent/child, persistence, etc.
    // FIXME - just deletes elements, not family trees for now

    /// Delete all elements that aren't child elements
    void Refiner::deleteParentElements()
    {
      //check_sidesets_2(" deleteParentElements:: start");
      //check_sidesets(" deleteParentElements:: start");
      //check_sidesets_1(" deleteParentElements:: start");

      std::vector<stk::mesh::EntityRank> ranks_to_be_deleted;
      ranks_to_be_deleted.push_back(m_eMesh.element_rank());
      ranks_to_be_deleted.push_back(m_eMesh.side_rank());

      //std::cout << "tmp srk ranks_to_be_deleted= " << ranks_to_be_deleted << std::endl;

      elements_to_be_destroyed_type parents;
      for (unsigned irank=0; irank < ranks_to_be_deleted.size(); irank++)
        {

          const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( ranks_to_be_deleted[irank] );
          int npar=0;
          int nchild=0;
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
                  if (!m_eMesh.isParentElement(element, false))
                  //if (!m_eMesh.hasFamilyTree(element) || m_eMesh.isChildElement(element, true))
                  //if (!m_eMesh.hasFamilyTree(element) || m_eMesh.isChildElementLeaf(element, true))
                    {
                      // it has no family tree, so it's a leaf, or it has no children
                      ++nchild;
                    }
                  else
                    {
                      ++npar;
#if UNIFORM_REF_REMOVE_OLD_STD_VECTOR
                      parents.push_back(&element);
#else
                      parents.insert(&element);
#endif
                    }
                }
            }
          //std::cout << "tmp removeElements(parents) irank, size= " << ranks_to_be_deleted[irank] << " " << npar << " nchild= " << nchild << std::endl;

        }

      m_eMesh.get_bulk_data()->modification_begin();
      //SidePartMap side_part_map;
      //get_side_part_relations(false, side_part_map);
#if PERCEPT_USE_FAMILY_TREE
      removeFamilyTrees();
#endif
      //std::cout << "tmp removeElements(parents) size= " << parents.size() << std::endl;
      removeElements(parents);
      //fix_side_sets_3(false, side_part_map);
      fix_side_sets_2();

      m_eMesh.get_bulk_data()->modification_end();

      //check_sidesets_2(" deleteParentElements:: end");
      //check_sidesets(" deleteParentElements:: end");
      //check_sidesets_1(" deleteParentElements:: end");

    }

    void Refiner::removeEmptyElements()
    {

      elements_to_be_destroyed_type list;

      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.element_rank() );

      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
        {
          stk::mesh::Bucket & bucket = **k ;

          const unsigned num_elements_in_bucket = bucket.size();
          for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
            {
              stk::mesh::Entity& element = bucket[iElement];
              if (0 == element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK).size())
                {
#if UNIFORM_REF_REMOVE_OLD_STD_VECTOR
                  list.push_back(&element);
#else
                  list.insert(&element);
#endif
                }
            }
        }

      //m_eMesh.get_bulk_data()->modification_begin();
      //std::cout << "tmp removeElements(parents) " << std::endl;
      removeElements(list);
      //m_eMesh.get_bulk_data()->modification_end();

    }

    void Refiner::removeDanglingNodes()
    {
      SetOfEntities node_list;
      SetOfEntities pseudos;

      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.node_rank() );

      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
        {
          stk::mesh::Bucket & bucket = **k ;

          const unsigned num_nodes_in_bucket = bucket.size();
          for (unsigned iElement = 0; iElement < num_nodes_in_bucket; iElement++)
            {
              stk::mesh::Entity& node = bucket[iElement];
              if (0 == node.relations().size())
                {
                  node_list.insert(&node);
                }
              else if (1 == node.relations().size() && node.relations()[0].entity()->entity_rank() == m_eMesh.element_rank() + PSEUDO_ELEMENT_RANK_SHIFT)
              {
                pseudos.insert( node.relations()[0].entity() );
                node_list.insert(&node);
              }
            }
        }

      //if (1 && node_list.size()) std::cout << "P[" << m_eMesh.get_rank() << "] tmp number of dangling nodes = " << node_list.size() << " and pseudos= " << pseudos.size() << std::endl;
      //!srk
      getNodeRegistry().clear_dangling_nodes(&node_list);

      if (1)
        {
          for (SetOfEntities::iterator itbd = pseudos.begin(); itbd != pseudos.end();  ++itbd)
            {
              stk::mesh::Entity *pseudo_p = *itbd;

              if ( ! m_eMesh.get_bulk_data()->destroy_entity( pseudo_p ) )
                {
                  throw std::logic_error("Refiner::removeDanglingNodes couldn't remove pseudo");

                }
            }
        }

      for (SetOfEntities::iterator itbd = node_list.begin(); itbd != node_list.end();  ++itbd)
        {
          stk::mesh::Entity *node_p = *itbd;

          if ( ! m_eMesh.get_bulk_data()->destroy_entity( node_p ) )
            {
              throw std::logic_error("Refiner::removeDanglingNodes couldn't remove node");

            }
        }

      // check for any null entities
      //std::cout << "check for any null entities..." << std::endl;
      const vector<stk::mesh::Bucket*> & elem_buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.element_rank() );

      for ( vector<stk::mesh::Bucket*>::const_iterator k = elem_buckets.begin() ; k != elem_buckets.end() ; ++k ) 
        {
          stk::mesh::Bucket & bucket = **k ;

          const unsigned num_nodes_in_bucket = bucket.size();
          for (unsigned iElement = 0; iElement < num_nodes_in_bucket; iElement++)
            {
              stk::mesh::Entity& elem = bucket[iElement];
              const stk::mesh::PairIterRelation& rels = elem.relations(m_eMesh.node_rank());
              for (unsigned j=0; j < rels.size(); j++)
                {
                  if (rels[j].entity() == 0) throw std::runtime_error("bad node in an element");
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

    void Refiner::snapAndSmooth(bool geomSnap, std::string geomFile)
    {
      std::cout << " geomFile= " << geomFile << " geomSnap= " << geomSnap << std::endl;
#if defined( STK_PERCEPT_HAS_GEOMETRY )
      if (geomFile == "") return;
      std::cout << " 2 geomFile= " << geomFile << " geomSnap= " << geomSnap << std::endl;

      //SMOOTHING_OPTIONS option = SNAP_PLUS_SMOOTH;
      SMOOTHING_OPTIONS option = USE_LINE_SEARCH_WITH_MULTIPLE_STATES;

      GeometryKernelOpenNURBS gk;
      // set to 0.0 for no checks, > 0.0 for a fixed check delta, < 0.0 (e.g. -0.5) to check against local edge length average times this |value|
      double doCheckMovement = 0.0; 
      //double doCheckMovement = -1.0; 

      // anything exceeding a value > 0.0 will be printed
      double doCheckCPUTime = 0.0;  
      //double doCheckCPUTime = 0.1;

      MeshGeometry mesh_geometry(&gk, doCheckMovement, doCheckCPUTime);
      GeometryFactory factory(&gk, &mesh_geometry);
      factory.read_file(geomFile, &m_eMesh);

      switch(option) {
      case SNAP_PLUS_SMOOTH:
        {
          mesh_geometry.snap_points_to_geometry(&m_eMesh);
          if (doCheckMovement != 0.0) 
            mesh_geometry.print_node_movement_summary();

#if defined (STK_PERCEPT_HAS_MESQUITE)
          if (m_doSmoothGeometry)
            {
              smoothGeometry(mesh_geometry,option);
              mesh_geometry.snap_points_to_geometry(&m_eMesh);
            }
#endif
        }
        break;
      case USE_LINE_SEARCH_WITH_MULTIPLE_STATES:
        {
          //VERIFY_OP_ON(m_eMesh.get_coordinates_field()->number_of_states(), ==, 3, "Must use PerceptMesh::set_num_coordinate_field_states(3) to use new smoothing.");
          stk::mesh::FieldBase *nm1_field = m_eMesh.get_field("coordinates_NM1");

#ifdef STK_PERCEPT_HAS_MESQUITE
          if (m_doSmoothGeometry && nm1_field)
            {
              // make a copy of current non-snapped state (dst,src)
              m_eMesh.copy_field(m_eMesh.get_field("coordinates_NM1"), m_eMesh.get_coordinates_field() );
            }
#endif

          // do the snap
          if (geomSnap)
            mesh_geometry.snap_points_to_geometry(&m_eMesh);
          if (doCheckMovement != 0.0) 
            mesh_geometry.print_node_movement_summary();

#ifdef STK_PERCEPT_HAS_MESQUITE
          if (m_doSmoothGeometry && nm1_field)
            {
              // make a copy of current snapped state
              m_eMesh.copy_field(m_eMesh.get_field("coordinates_N"), m_eMesh.get_coordinates_field() );

              // reset current state to non-snapped state
              m_eMesh.copy_field(m_eMesh.get_coordinates_field(), m_eMesh.get_field("coordinates_NM1") );

              smoothGeometry(mesh_geometry,option);
              //mesh_geometry.snap_points_to_geometry(&m_eMesh);
            }
#endif
              
        }
        break;
      }

#endif
    }

#if defined( STK_PERCEPT_HAS_MESQUITE ) && defined(STK_PERCEPT_HAS_GEOMETRY)
    void Refiner::smoothGeometry(MeshGeometry& mesh_geometry, SMOOTHING_OPTIONS option)
    {
      bool do_mesquite_smoothing = true;
      if (do_mesquite_smoothing)
        {
          int  msq_debug             = 2; // 1,2,3 for more debug info
          bool always_smooth         = true;
          bool do_laplace            = false;
          bool do_jacobi             = true;

          PerceptMesquiteMeshDomain pmd(&m_eMesh, &mesh_geometry);

          switch(option) {
          case SNAP_PLUS_SMOOTH:
            {
#define ALWAYS_PMM_PARALLEL 1
#if ALWAYS_PMM_PARALLEL
              PerceptMesquiteMesh pmm0(&m_eMesh, &pmd);
              PerceptMesquiteMesh::PMMParallelMesh pmm(&pmm0);
              Mesquite::MsqError err;
              //pmm.helper.set_communication_model(Mesquite::ParallelHelperImpl::Blocking, err);
#else
              PerceptMesquiteMesh pmm(&m_eMesh, &pmd);
#endif
              if (do_laplace)
                {
                  PMMLaplaceSmoother1 ls;
                  if (do_jacobi) ls.get_smoother().do_jacobi_optimization();
                  ls.run(pmm, pmd, always_smooth, msq_debug);
                }
              else
                {
                  //PMMShapeImprover(int innerIter=100, double gradNorm = 1.e-8, int parallelIterations=20) : 
                  PMMShapeImprover si(100, 1.e-8, 20);
                  //PMMShapeImprover si(5, 1.e-8, 20);
                  //const double max_vertex_movement_term_crit=10;
                  //PMMShapeSizeOrientImprover si(10);
                  si.run(pmm, pmd, always_smooth, msq_debug);
                }
            }
            break;
          case USE_LINE_SEARCH_WITH_MULTIPLE_STATES:
            {
#if 1
              // geometry used for classification of fixed/non-fixed nodes
              PerceptMesquiteMesh pmm(&m_eMesh, &pmd);

              //PMMShapeImprover(int innerIter=100, double gradNorm = 1.e-8, int parallelIterations=20) : 

              percept::PMMParallelShapeImprover pmmpsi(1001, 1.e-4, 1);
              //pmmpsi.run(pmm, &pmd, always_smooth, msq_debug);
              pmmpsi.run(pmm, &pmd, always_smooth, msq_debug);


#endif
            }
            break;
          }
          return;
        }
    }
#endif

    unsigned Refiner::
    doForAllElements(unsigned irank, std::string function_info,
                     stk::mesh::EntityRank rank, NodeRegistry::ElementFunctionPrototype function,
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
          progress_meter_num_total = doForAllElements(irank, function_info, rank, function, elementColors, elementType, needed_entity_ranks,  true, doAllElements);
          m_doProgress = true;
          std::ostringstream oss; oss << function_info <<" [" << 100.0*((double)irank)/((double)m_ranks.size()) << " %]";
          ProgressMeterData pd(ProgressMeterData::INIT, 0.0, oss.str());
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
                  refineMethodApply(function, element, needed_entity_ranks);
                }

              if (m_doProgress && (num_elem % progress_meter_when_to_post == 0) )
                {
                  double progress_meter_percent = 100.0*((double)num_elem)/std::max(d_progress_meter_num_total,1.0);
                  std::ostringstream oss; oss << function_info << " [" << 100.0*((double)irank)/((double)m_ranks.size()) << " %]";
                  ProgressMeterData pd(ProgressMeterData::RUNNING, progress_meter_percent, oss.str());
                  notifyObservers(&pd);
                  if (0) std::cout << "progress_meter_percent = " << progress_meter_percent << std::endl;
                }

            } // elements in this color
        } // icolor

      if (m_doProgress)
        {
          std::ostringstream oss; oss << function_info << " [" << 100.0*((double)irank)/((double)m_ranks.size()) << " %]";
          ProgressMeterData pd(ProgressMeterData::FINI, 0.0, oss.str());
          notifyObservers(&pd);
        }

      return num_elem;
    }

    void Refiner::
    createElementsAndNodesAndConnectLocal(unsigned irank, stk::mesh::EntityRank rank, UniformRefinerPatternBase *breakPattern,
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
          std::ostringstream oss; oss << "Create Elements pass [" << 100.0*((double)irank)/((double)m_ranks.size()) << " %]";
          ProgressMeterData pd(ProgressMeterData::INIT, 0.0, oss.str());
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
          shards::CellTopology cell_topo(cell_topo_data);

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
                  std::ostringstream oss; oss << "Create Elements pass [" << 100.0*((double)irank)/((double)m_ranks.size()) << " %]";
                  ProgressMeterData pd(ProgressMeterData::RUNNING, 100.0*((double)jele)/((double)std::max(nele,1)), oss.str());
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
                  //std::cout << "P["<< m_eMesh.get_rank() << "] element.owner_rank() = " << element.owner_rank() << std::endl;
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
          std::ostringstream oss; oss << "Create Elements pass [" << 100.0*((double)irank)/((double)m_ranks.size()) << " %]";
          ProgressMeterData pd(ProgressMeterData::FINI, 0.0, oss.str());
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
              std::cout << "P[" << m_eMesh.get_rank() << "]  needed_entity_ranks[ineed_ent]= " << needed_entity_ranks[ineed_ent].first
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
                        std::cout << "P[" << m_eMesh.get_rank() << "] nodeId ## = 0 << "
                                  << " nodeIds_onSE.m_entity_id_vector[0] = " << nodeIds_onSE.m_entity_id_vector[0]
                                  << " element= " << element
                                  << " needed_entity_ranks= " << needed_entity_ranks[ineed_ent].first
                                  << " iSubDimOrd = " << iSubDimOrd
                                  <<  std::endl;
                        std::cout << " element= ";
                        m_eMesh.print_entity(std::cout, element, 0);

                        throw std::logic_error("Refiner::createNewNeededNodeIds logic err #5.0, nodeIds_onSE.m_entity_id_vector[i_new_node] == 0");
                      }
                  }

                stk::mesh::Entity * node1 = m_eMesh.get_bulk_data()->get_entity(stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE.m_entity_id_vector[0]);
                nodeIds_onSE[0] = node1;

                if (!node1)
                  {
                    if (!m_nodeRegistry->getUseCustomGhosting())
                    {
                      static stk::mesh::PartVector empty_parts;
                      node1 = & m_eMesh.get_bulk_data()->declare_entity(stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE.m_entity_id_vector[0], empty_parts);
                    }

                    if (!node1)
                    {
                      std::cout << "P[" << m_eMesh.get_rank() << "] nodeId ## = 0 << "
                              << " nodeIds_onSE.m_entity_id_vector[0] = " << nodeIds_onSE.m_entity_id_vector[0] << " node1= " << node1
                              << " element= " << element
                              << " needed_entity_ranks= " << needed_entity_ranks[ineed_ent].first
                              << " iSubDimOrd = " << iSubDimOrd
                              <<  std::endl;
                      throw std::logic_error("Refiner::createNewNeededNodeIds logic error #0");
                    }
                  }
                
              }

              unsigned num_new_nodes_needed = needed_entity_ranks[ineed_ent].second;
              if (0)
                {
                  const CellTopologyData * const cell_topo_data_0 = stk::percept::PerceptMesh::get_cell_topology(element);
                  shards::CellTopology cell_topo_0(cell_topo_data_0);

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
                              std::cout << "P[" << m_eMesh.get_rank() << "] nodeId ## = 0 << "
                                        << " nodeIds_onSE.m_entity_id_vector[0] = " << nodeIds_onSE.m_entity_id_vector[0]
                                        << " element= " << element
                                        << " needed_entity_ranks= " << needed_entity_ranks[ineed_ent].first
                                        << " iSubDimOrd = " << iSubDimOrd
                                        <<  std::endl;
                              std::cout << " element= ";
                              m_eMesh.print_entity(std::cout, element, 0);

                            }
                          throw std::logic_error("Refiner::createNewNeededNodeIds logic err #5.1, nodeIds_onSE.m_entity_id_vector[i_new_node] == 0");
                        }
                      stk::mesh::Entity * node1 = m_eMesh.get_bulk_data()->get_entity(stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE.m_entity_id_vector[i_new_node]);

                      if (!node1)
                        {
                          if (!m_nodeRegistry->getUseCustomGhosting())
                          {
                            static stk::mesh::PartVector empty_parts;
                            node1 = & m_eMesh.get_bulk_data()->declare_entity(stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE.m_entity_id_vector[i_new_node], empty_parts);
                          }
                          if (!node1)
                          { 
                            throw std::logic_error("Refiner::createNewNeededNodeIds logic err #4");
                          }
                        }
                      nodeIds_onSE[i_new_node] = node1;
                      VERIFY_OP_ON(node1->identifier(), ==, nodeIds_onSE.m_entity_id_vector[i_new_node], "Refiner::createNewNeededNodeIds logic err #4.1");
                    }
                  new_sub_entity_nodes[needed_entity_ranks[ineed_ent].first][iSubDimOrd][i_new_node] = nodeIds_onSE[i_new_node]->identifier();

#ifndef NDEBUG
                  stk::mesh::Entity * node2 = m_eMesh.get_bulk_data()->get_entity(stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE[i_new_node]->identifier() );
                  if (!node2)
                    {
                      std::cout << "P[" << m_eMesh.get_rank() << "] element is ghost = " << m_eMesh.isGhostElement(element)
                                << " needed_entity_ranks= " << needed_entity_ranks[ineed_ent].first << " iSubDimOrd= " << iSubDimOrd
                                << " i_new_node= " << i_new_node 
                                << " id= " << nodeIds_onSE[i_new_node]->identifier() 
                                << " entity_vec_id= " << nodeIds_onSE.m_entity_id_vector[i_new_node] 
                                << std::endl;
                        
                      VERIFY_OP_ON(node2, !=, 0, "Refiner::createNewNeededNodeIds logic err #6 - node2 is null");
                    }
#endif

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

      if (m_eMesh.get_spatial_dim() == 3)
        {
          fixElementSides1(m_eMesh.face_rank());
        }
      // FIXME
      else if (m_eMesh.get_spatial_dim() == 2)
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
                          m_eMesh.get_bulk_data()->declare_relation(*child, *parent_side_child, k_child_side);
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
    check_sidesets(std::string msg)
    {
      EXCEPTWATCH;

      std::cout << "tmp check_sidesets start..." << std::endl;

      stk::mesh::EntityRank node_rank = m_eMesh.node_rank();
      stk::mesh::EntityRank side_rank = m_eMesh.side_rank();
      stk::mesh::EntityRank element_rank = m_eMesh.element_rank();

      const vector<stk::mesh::Bucket*> & side_buckets = m_eMesh.get_bulk_data()->buckets( side_rank );
      for ( vector<stk::mesh::Bucket*>::const_iterator it_side_bucket = side_buckets.begin() ; it_side_bucket != side_buckets.end() ; ++it_side_bucket )
        {
          stk::mesh::Bucket & side_bucket = **it_side_bucket ;
          const unsigned num_elements_in_side_bucket = side_bucket.size();
          for (unsigned i_side = 0; i_side < num_elements_in_side_bucket; i_side++)
            {
              stk::mesh::Entity& side = side_bucket[i_side];
              if (side.relations(node_rank).size() == 0)
                continue;
              if (!m_eMesh.hasFamilyTree(side) || m_eMesh.isChildElement(side, true))
                {
                  bool found = false;
                  const vector<stk::mesh::Bucket*> & element_buckets = m_eMesh.get_bulk_data()->buckets( element_rank );
                  for ( vector<stk::mesh::Bucket*>::const_iterator it_element_bucket = element_buckets.begin() ; 
                        it_element_bucket != element_buckets.end() ; ++it_element_bucket )
                    {
                      stk::mesh::Bucket & element_bucket = **it_element_bucket ;
                      const unsigned num_elements_in_element_bucket = element_bucket.size();
                      for (unsigned i_element = 0; i_element < num_elements_in_element_bucket; i_element++)
                        {
                          stk::mesh::Entity& element = element_bucket[i_element];
                          if (element.relations(node_rank).size() == 0)
                            continue;

                          shards::CellTopology element_topo(stk::percept::PerceptMesh::get_cell_topology(element));
                          unsigned element_nsides = (unsigned)element_topo.getSideCount();

                          for (unsigned j_element_side = 0; j_element_side < element_nsides; j_element_side++)
                            {
                              int permIndex = -1;
                              int permPolarity = 1;
                              PerceptMesh::element_side_permutation(element, side, j_element_side, permIndex, permPolarity);
                              if (permIndex >= 0)
                                {
                                  found = true;
                                  break;
                                }
                            }
                        }
                    } // elements...
                  if (!found)
                    {
                      std::cout << "ERROR: side = " << side << std::cout;
                      throw std::logic_error("check_sidesets error");
                    }
                }
            }
        }
      std::cout << "tmp check_sidesets ...end" << std::endl;
    }

    // check for two sides sharing all nodes
    void Refiner::
    check_sidesets_1(std::string msg)
    {
      EXCEPTWATCH;

      std::cout <<  "P["<< m_eMesh.get_rank() << "] tmp check_sidesets_1 start... " << msg << std::endl;

      stk::mesh::EntityRank node_rank = m_eMesh.node_rank();
      stk::mesh::EntityRank side_rank = m_eMesh.side_rank();

      typedef std::set<stk::mesh::Entity *> SetOfEntities;
      SetOfEntities side_set;

      const vector<stk::mesh::Bucket*> & side_buckets = m_eMesh.get_bulk_data()->buckets( side_rank );
      for ( vector<stk::mesh::Bucket*>::const_iterator it_side_bucket = side_buckets.begin() ; it_side_bucket != side_buckets.end() ; ++it_side_bucket )
        {
          stk::mesh::Bucket & side_bucket = **it_side_bucket ;
          const unsigned num_elements_in_side_bucket = side_bucket.size();
          for (unsigned i_side = 0; i_side < num_elements_in_side_bucket; i_side++)
            {
              stk::mesh::Entity& side = side_bucket[i_side];

              if (m_eMesh.isGhostElement(side))
                continue;

              if (side.relations(node_rank).size() == 0)
                continue;

              if (m_eMesh.check_entity_duplicate(side))
                {
                  throw std::logic_error("found two sides with same nodes");
                }
            }
        }
      std::cout <<  "P["<< m_eMesh.get_rank() << "] tmp check_sidesets_1 ... end " << msg << std::endl;
    }

    // fast check if side elems have relations to elements and vice versa
    void Refiner::
    check_sidesets_2(std::string msg)
    {
      EXCEPTWATCH;


      stk::mesh::EntityRank node_rank = m_eMesh.node_rank();
      stk::mesh::EntityRank side_rank = m_eMesh.side_rank();
      stk::mesh::EntityRank element_rank = m_eMesh.element_rank();

      std::cout << "tmp check_sidesets_2 start... " << msg << " side_rank= " << side_rank << " element_rank= " << element_rank << std::endl;

      typedef std::set<stk::mesh::Entity *> SetOfEntities;
      SetOfEntities side_set;

      const vector<stk::mesh::Bucket*> & side_buckets = m_eMesh.get_bulk_data()->buckets( side_rank );
      for ( vector<stk::mesh::Bucket*>::const_iterator it_side_bucket = side_buckets.begin() ; it_side_bucket != side_buckets.end() ; ++it_side_bucket )
        {
          stk::mesh::Bucket & side_bucket = **it_side_bucket ;
          const unsigned num_elements_in_side_bucket = side_bucket.size();
          for (unsigned i_side = 0; i_side < num_elements_in_side_bucket; i_side++)
            {
              stk::mesh::Entity& side = side_bucket[i_side];
              if (m_eMesh.isGhostElement(side))
                continue;

              if (side.relations(node_rank).size() == 0)
                {
                  continue;
                }

              if (side.relations(element_rank).size() > 1)
                {
                  throw std::logic_error("check_sidesets_2: too many side relations");
                }

              if (side.relations(element_rank).size() < 1)
                {
                  throw std::logic_error("check_sidesets_2: too few side relations");
                }

              //if (!m_eMesh.isLeafElement(side))
              //  continue;

              bool found = false;

              stk::mesh::PairIterRelation side_nodes = side.relations(node_rank);

              for (unsigned isnode=0; isnode < side_nodes.size(); isnode++)
                {
                  stk::mesh::PairIterRelation node_elements = side_nodes[isnode].entity()->relations(element_rank);
                  for (unsigned ienode=0; ienode < node_elements.size(); ienode++)
                    {
                      stk::mesh::Entity& element = *node_elements[ienode].entity();
              
                      if (element.relations(node_rank).size() == 0)
                        continue;
                      if (m_eMesh.isGhostElement(element))
                        continue;

                      // FIXME
                      //if (m_eMesh.isLeafElement(element))
                      {
                        shards::CellTopology element_topo(stk::percept::PerceptMesh::get_cell_topology(element));
                        unsigned element_nsides = (unsigned)element_topo.getSideCount();

                        for (unsigned j_element_side = 0; j_element_side < element_nsides; j_element_side++)
                          {
                            int permIndex = -1;
                            int permPolarity = 1;
                            PerceptMesh::element_side_permutation(element, side, j_element_side, permIndex, permPolarity);
                            if (permIndex >= 0)
                              {
                                found = true;
                                //std::cout << "found side element needing fixing, id= " << side.identifier() 
                                //          <<  " ele id= " << element.identifier() <<  std::endl;
                                //m_eMesh.get_bulk_data()->declare_relation(element, side, j_element_side);
                              }
                            if (found) break;
                          }
                        if (found) break;
                      }
                      if (found) break;
                    }
                  if (found) break;
                }

              if (!found)
                {
                  std::cout << "ERROR: side = " << side << std::cout;
                  throw std::logic_error("check_sidesets_2 error");
                }
            }
        }
      std::cout << "tmp check_sidesets_2 ...end" << std::endl;
    }

    // find remaining missing relations
    void Refiner::
    fix_side_sets_1()
    {
      EXCEPTWATCH;

      //std::cout << "tmp fix_side_sets_1 start..." << std::endl;

      stk::mesh::EntityRank node_rank = m_eMesh.node_rank();
      stk::mesh::EntityRank side_rank = m_eMesh.side_rank();
      stk::mesh::EntityRank element_rank = m_eMesh.element_rank();

      SetOfEntities side_set_with_empty_relations;

      const vector<stk::mesh::Bucket*> & side_buckets = m_eMesh.get_bulk_data()->buckets( side_rank );
      for ( vector<stk::mesh::Bucket*>::const_iterator it_side_bucket = side_buckets.begin() ; it_side_bucket != side_buckets.end() ; ++it_side_bucket )
        {
          stk::mesh::Bucket & side_bucket = **it_side_bucket ;
          const unsigned num_elements_in_side_bucket = side_bucket.size();
          for (unsigned i_side = 0; i_side < num_elements_in_side_bucket; i_side++)
            {
              stk::mesh::Entity& side = side_bucket[i_side];

              if (side.relations(node_rank).size() == 0)
                continue;

              //if (!m_eMesh.hasFamilyTree(side) || m_eMesh.isChildElement(side, true))
              {
                if (side.relations(element_rank).size() > 1)
                  {
                    throw std::logic_error("fix_side_sets_1: too many side relations");
                  }

                // already found
                if (side.relations(element_rank).size())
                  continue;

                side_set_with_empty_relations.insert(&side);
              }
            }
        }

      for (SetOfEntities::iterator side_iter = side_set_with_empty_relations.begin(); side_iter != side_set_with_empty_relations.end(); ++side_iter)
        {
          stk::mesh::Entity& side = **side_iter;

          if (side.relations(node_rank).size() == 0)
            continue;

          if (side.relations(element_rank).size() > 1)
            {
              throw std::logic_error("fix_side_sets_1: too many side relations");
            }

          // already found
          if (side.relations(element_rank).size())
            continue;
          
          bool found = false;
          const vector<stk::mesh::Bucket*> & element_buckets = m_eMesh.get_bulk_data()->buckets( element_rank );
          for ( vector<stk::mesh::Bucket*>::const_iterator it_element_bucket = element_buckets.begin() ; 
                it_element_bucket != element_buckets.end() ; ++it_element_bucket )
            {
              stk::mesh::Bucket & element_bucket = **it_element_bucket ;
              const unsigned num_elements_in_element_bucket = element_bucket.size();
              for (unsigned i_element = 0; i_element < num_elements_in_element_bucket; i_element++)
                {
                  stk::mesh::Entity& element = element_bucket[i_element];
                  if (element.relations(node_rank).size() == 0)
                    continue;

                  if (!m_eMesh.hasFamilyTree(element) || m_eMesh.isChildElement(element, true))
                    {
                      shards::CellTopology element_topo(stk::percept::PerceptMesh::get_cell_topology(element));
                      unsigned element_nsides = (unsigned)element_topo.getSideCount();

                      for (unsigned j_element_side = 0; j_element_side < element_nsides; j_element_side++)
                        {
                          int permIndex = -1;
                          int permPolarity = 1;
                          PerceptMesh::element_side_permutation(element, side, j_element_side, permIndex, permPolarity);
                          if (permIndex >= 0)
                            {
                              found = true;
                              if (1)
                                {
                                  std::cout << "found side element needing fixing, id= " << side.identifier() <<  std::endl;
                                  std::cout << "found side element needing fixing, ele id= " << element.identifier() <<  std::endl;
                                  //exit(123);
                                  throw std::logic_error("fix_side_sets_2 error 1");
                                }
                              

                              mesh::PairIterRelation rels = element.relations(side_rank);

                              //std::cout << "found 1 side element needing fixing, id= " << side.identifier() <<  std::endl;

                              bool found_existing_rel = false;
                              for (unsigned irels=0; irels < rels.size(); irels++)
                                {
                                  if (rels[irels].entity() == &side)
                                    {
                                      //std::cout << "found 2 side element needing fixing, id= " << side.identifier() <<  std::endl;
                                      found_existing_rel = true;
                                      break;
                                    }
                                }
                              if (found_existing_rel)
                                {
                                  throw std::logic_error("found_existing_rel where none expected");
                                }
                              else
                                {
                                  //std::cout << "found 3 side element needing fixing, id= " << side.identifier() <<  std::endl;
                                  m_eMesh.get_bulk_data()->declare_relation(element, side, j_element_side);
                                  //std::cout << "found 4 side element needing fixing, id= " << side.identifier() <<  std::endl;
                                }

                              break;
                            }
                          if (found) break;
                        }
                      if (found) break;
                    }
                  if (found) break;
                }
              if (found) break;
            } // element buckets...

          if (!found)
            {
              std::cout << "ERROR: side = " << side << std::cout;
              throw std::logic_error("fix_side_sets_1 error");
            }
        }
      //std::cout << "tmp fix_side_sets_1 ...end" << std::endl;
    }

    // determine side part to elem part relations
    void Refiner::
    get_side_part_relations(bool checkParentChild, SidePartMap& side_part_map)
    {
      EXCEPTWATCH;
      std::cout << "get_side_part_relations start...\n";

      stk::mesh::EntityRank node_rank = m_eMesh.node_rank();
      stk::mesh::EntityRank edge_rank = m_eMesh.edge_rank();
      stk::mesh::EntityRank side_rank = m_eMesh.side_rank();
      stk::mesh::EntityRank element_rank = m_eMesh.element_rank();

      int spatialDim = m_eMesh.get_spatial_dim();

      unsigned side_rank_iter_begin = side_rank;
      unsigned side_rank_iter_end = side_rank;
      if (spatialDim == 3)
        {
          side_rank_iter_begin = edge_rank;
        }

      // get super-relations (side_part.name() --> elem_part.name())
      for (unsigned side_rank_iter = side_rank_iter_begin; side_rank_iter <= side_rank_iter_end; side_rank_iter++)
        {
          const vector<stk::mesh::Bucket*> & side_buckets = m_eMesh.get_bulk_data()->buckets( side_rank_iter );
          for ( vector<stk::mesh::Bucket*>::const_iterator it_side_bucket = side_buckets.begin() ; it_side_bucket != side_buckets.end() ; ++it_side_bucket )
            {
              stk::mesh::Bucket & side_bucket = **it_side_bucket ;
              stk::mesh::PartVector side_parts;
              side_bucket.supersets(side_parts);
              if (0)
                for (unsigned isp=0; isp < side_parts.size(); isp++)
                  {
                    std::cout << "side_part= " << side_parts[isp]->name() << std::endl;
                  }

              const unsigned num_elements_in_side_bucket = side_bucket.size();
              for (unsigned i_side = 0; i_side < num_elements_in_side_bucket; i_side++)
                {
                  stk::mesh::Entity& side = side_bucket[i_side];

                  if (m_eMesh.isGhostElement(side))
                    continue;

                  if (side.relations(node_rank).size() == 0)
                    continue;

                  if (side.relations(element_rank).size() > 1)
                    {
                      std::cout << "get_side_part_relations: too many side relations" << std::endl;
                      throw std::logic_error("get_side_part_relations: too many side relations");
                    }

                  bool isLeafElement = !checkParentChild || m_eMesh.isLeafElement(side);
                  if (isLeafElement)
                    {
                      stk::mesh::PairIterRelation side_to_elem_rels = side.relations(element_rank);
                      for (unsigned irel = 0; irel < side_to_elem_rels.size(); irel++)
                        {
                          stk::mesh::Entity& elem = *side_to_elem_rels[irel].entity();
                          stk::mesh::PartVector elem_parts;
                          elem.bucket().supersets(elem_parts);
                          for (unsigned isp=0; isp < side_parts.size(); isp++)
                            {
                              if ( stk::mesh::is_auto_declared_part(*side_parts[isp]) )
                                continue;
                              const CellTopologyData *const topology = stk::percept::PerceptMesh::get_cell_topology(*side_parts[isp]);
                              if (!topology)
                                continue;

                              for (unsigned iep=0; iep < elem_parts.size(); iep++)
                                {
                                  if ( stk::mesh::is_auto_declared_part(*elem_parts[iep]) )
                                    continue;

                                  //if (elem_parts[iep]->name().find(UniformRefinerPatternBase::getOldElementsPartName()) != std::string::npos)
                                  //  continue;

                                  const STK_Adapt_Auto_Part *auto_part = elem_parts[iep]->attribute<STK_Adapt_Auto_Part>();
                                  if (elem_parts[iep]->name().find(UniformRefinerPatternBase::getOldElementsPartName()) != std::string::npos)
                                    {
                                      if (!auto_part) throw std::runtime_error("Refiner::get_side_part_relations: bad old part attribute for auto");
                                    }

                                  if (auto_part) 
                                    {
                                      continue;
                                    }

                                  if (elem_parts[iep] != side_parts[isp])
                                    {
                                      SidePartMap::iterator found = side_part_map.find(side_parts[isp]->name());
                                      if (found != side_part_map.end())
                                        {
                                          if (found->second != elem_parts[iep]->name())
                                            {
                                              std::cout << "side_part = " << side_parts[isp]->name() 
                                                        << " elem_part = " << elem_parts[iep]->name()
                                                        << " found_elem_part = " << found->second
                                                        << std::endl;
                                              throw std::runtime_error("Refiner::get_side_part_relations: too many side_part to elem_part relations");
                                            }
                                        }
                                      side_part_map[side_parts[isp]->name()] = elem_parts[iep]->name();
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

      if (1)
        {
          SidePartMap::iterator iter;
          for (iter = side_part_map.begin(); iter != side_part_map.end(); iter++)
            {
              std::cout << "Refiner::get_side_part_relations: side_part = " << iter->first << " elem_part= " << iter->second << std::endl;
            }
          //exit(1);
        }
      std::cout << "get_side_part_relations ... done\n";
    }

    // fast reconnector - checks for which super-connections are valid (from the original mesh, look
    //   at which side-part touches which element-part, save the info), and reconnects found side/elem matches.
    void Refiner::
    fix_side_sets_3(bool checkParentChild, SidePartMap& side_part_map)
    {
      EXCEPTWATCH;

      std::cout << "tmp fix_side_sets_3 start... checkParentChild= " <<  checkParentChild << std::endl;

      stk::mesh::EntityRank node_rank = m_eMesh.node_rank();
      stk::mesh::EntityRank edge_rank = m_eMesh.edge_rank();
      stk::mesh::EntityRank side_rank = m_eMesh.side_rank();
      stk::mesh::EntityRank element_rank = m_eMesh.element_rank();

      int spatialDim = m_eMesh.get_spatial_dim();

      //const unsigned FAMILY_TREE_RANK = m_eMesh.element_rank() + 1u;
      //const vector<stk::mesh::Bucket*> & family_tree_buckets = m_eMesh.get_bulk_data()->buckets( FAMILY_TREE_RANK );
      //bool have_family_tree = family_tree_buckets.size() > 0;

      //std::cout << "tmp fix_side_sets_3 side_rank= " << side_rank << " element_rank= " << element_rank << std::endl;

      // loop over all sides that are leaves (not parent or have no family tree), 
      //   loop over their nodes and their associated elements,
      //     connect element and side if they share a face

      unsigned side_rank_iter_begin = side_rank;
      unsigned side_rank_iter_end = side_rank;
      if (spatialDim == 3)
        {
          side_rank_iter_begin = edge_rank;
        }

      for (unsigned side_rank_iter = side_rank_iter_begin; side_rank_iter <= side_rank_iter_end; side_rank_iter++)
        {
          SetOfEntities side_set;
          int connections=0;

          const vector<stk::mesh::Bucket*> & side_buckets = m_eMesh.get_bulk_data()->buckets( side_rank_iter );
          for ( vector<stk::mesh::Bucket*>::const_iterator it_side_bucket = side_buckets.begin() ; it_side_bucket != side_buckets.end() ; ++it_side_bucket )
            {
              stk::mesh::Bucket & side_bucket = **it_side_bucket ;
              const unsigned num_elements_in_side_bucket = side_bucket.size();
              for (unsigned i_side = 0; i_side < num_elements_in_side_bucket; i_side++)
                {
                  stk::mesh::Entity& side = side_bucket[i_side];

                  if (m_eMesh.isGhostElement(side))
                    continue;

                  if (side.relations(node_rank).size() == 0)
                    continue;

                  if (side.relations(element_rank).size() > 1)
                    {
                      std::cout << "fix_side_sets_3: too many side relations" << std::endl;
                      throw std::logic_error("fix_side_sets_3: too many side relations");
                    }

                  bool isLeafElement = !checkParentChild || m_eMesh.isLeafElement(side);
                  if (!isLeafElement)
                    continue;

                  side_set.insert(&side);
                }
            }
      
          for (SetOfEntities::iterator it_side=side_set.begin(); it_side != side_set.end(); ++it_side)
            {
              stk::mesh::Entity& side = **it_side;

              bool found = false;

              stk::mesh::PairIterRelation side_nodes = side.relations(node_rank);

              for (unsigned isnode=0; isnode < side_nodes.size(); isnode++)
                {
                  stk::mesh::PairIterRelation node_elements = side_nodes[isnode].entity()->relations(element_rank);
                  for (unsigned ienode=0; ienode < node_elements.size(); ienode++)
                    {
                      stk::mesh::Entity& element = *node_elements[ienode].entity();

                      if (element.relations(node_rank).size() == 0)
                        continue;
                      if (m_eMesh.isGhostElement(element))
                        continue;

                      // FIXME
                      bool isLeafElement = !checkParentChild || m_eMesh.isLeafElement(element);
                      if (isLeafElement)
                        {
                          if (connectSides(&element, &side, &side_part_map))
                            {
                              ++connections;
                              found = true;
                            }
                        }
                      if (found) break;
                    }
                  if (found) break;
                }

              if (!found)
                {
                  std::cout << "ERROR: side = " << side << std::endl;
                  throw std::logic_error("fix_side_sets_3 error 2");
                }
            }
          std::cout << "Refiner::fix_side_sets_3 number of connections for side_rank= " << side_rank_iter << " = " << connections << std::endl;
          
        }

      std::cout << "tmp fix_side_sets_3 ...end" << std::endl;
    }

    // fast check if side elems have relations to elements and vice versa
    void Refiner::
    fix_side_sets_2()
    {
      EXCEPTWATCH;

      //std::cout << "tmp fix_side_sets_2 start... " << std::endl;

      stk::mesh::EntityRank node_rank = m_eMesh.node_rank();
      stk::mesh::EntityRank side_rank = m_eMesh.side_rank();
      stk::mesh::EntityRank element_rank = m_eMesh.element_rank();

      int spatialDim = m_eMesh.get_spatial_dim();

      //std::cout << "tmp fix_side_sets_2 side_rank= " << side_rank << " element_rank= " << element_rank << std::endl;

      // loop over all sides that are leaves (not parent or have no family tree), 
      //   loop over their nodes and their associated elements,
      //     connect element and side if they share a face

      unsigned side_rank_iter_begin = side_rank;
      unsigned side_rank_iter_end = side_rank;
      if (spatialDim == 3)
        {
          side_rank_iter_begin = m_eMesh.edge_rank();
        }
      for (unsigned side_rank_iter = side_rank_iter_begin; side_rank_iter <= side_rank_iter_end; side_rank_iter++)
        {
          SetOfEntities side_set;

          const vector<stk::mesh::Bucket*> & side_buckets = m_eMesh.get_bulk_data()->buckets( side_rank_iter );
          for ( vector<stk::mesh::Bucket*>::const_iterator it_side_bucket = side_buckets.begin() ; it_side_bucket != side_buckets.end() ; ++it_side_bucket )
            {
              stk::mesh::Bucket & side_bucket = **it_side_bucket ;
              const unsigned num_elements_in_side_bucket = side_bucket.size();
              for (unsigned i_side = 0; i_side < num_elements_in_side_bucket; i_side++)
                {
                  stk::mesh::Entity& side = side_bucket[i_side];

                  if (m_eMesh.isGhostElement(side))
                    continue;

                  if (side.relations(node_rank).size() == 0)
                    continue;

                  if (side.relations(element_rank).size() > 1)
                    {
                      std::cout << "fix_side_sets_2: too many side relations" << std::endl;
                      throw std::logic_error("fix_side_sets_2: too many side relations");
                    }
                  // if we already have a connection, skip this side
                  if (side.relations(element_rank).size() == 1)
                    continue;

                  if (!m_eMesh.isLeafElement(side))
                    continue;

                  side_set.insert(&side);
                }
            }
      
          for (SetOfEntities::iterator it_side=side_set.begin(); it_side != side_set.end(); ++it_side)
            {
              stk::mesh::Entity& side = **it_side;

              bool found = false;

              stk::mesh::PairIterRelation side_nodes = side.relations(node_rank);

              for (unsigned isnode=0; isnode < side_nodes.size(); isnode++)
                {
                  stk::mesh::PairIterRelation node_elements = side_nodes[isnode].entity()->relations(element_rank);
                  for (unsigned ienode=0; ienode < node_elements.size(); ienode++)
                    {
                      stk::mesh::Entity& element = *node_elements[ienode].entity();
              
                      if (element.relations(node_rank).size() == 0)
                        continue;
                      if (m_eMesh.isGhostElement(element))
                        continue;

                      // FIXME
                      if (m_eMesh.isLeafElement(element))
                        {
                          if (connectSides(&element, &side))
                            found = true;
                        }
                      if (found) break;
                    }
                  if (found) break;
                }

              if (!found)
                {
                  std::cout << "ERROR: side = " << side << std::endl;
                  throw std::logic_error("fix_side_sets_2 error 2");
                }
            }
        }

      //std::cout << "tmp fix_side_sets_2 ...end" << std::endl;
    }

    void Refiner::
    fixElementSides1(stk::mesh::EntityRank side_rank)
    {
      EXCEPTWATCH;

      // FIXME - cleanup, remove all old code

      //std::cout << "fixElementSides1 " << std::endl;
      // we don't check parent child here so we can pick up the original elements
      // FIXME - should we have option of checking parent, child, parent|child, etc?
      //SidePartMap side_part_map;
      //get_side_part_relations(false, side_part_map);
      //fix_side_sets_3(true, side_part_map);
      fix_side_sets_2();

    }

    // if the element (element) has a side that matches  the given side (side_elem), connect them but first delete old connections
    bool Refiner::connectSides(stk::mesh::Entity *element, stk::mesh::Entity *side_elem, SidePartMap *side_part_map)
    {
      EXCEPTWATCH;
      shards::CellTopology element_topo(stk::percept::PerceptMesh::get_cell_topology(*element));
      unsigned element_nsides = (unsigned)element_topo.getSideCount();

      // check validity of connection 
      if (side_part_map)
        {
          bool valid = false;
          stk::mesh::PartVector side_parts, elem_parts;
          element->bucket().supersets(elem_parts);
          side_elem->bucket().supersets(side_parts);
          for (unsigned isp = 0; isp < side_parts.size(); isp++)
            {
              if ( stk::mesh::is_auto_declared_part(*side_parts[isp]) )
                continue;
              const CellTopologyData *const topology = stk::percept::PerceptMesh::get_cell_topology(*side_parts[isp]);
              if (!topology)
                continue;

              SidePartMap::iterator found = side_part_map->find(side_parts[isp]->name());
              if (found == side_part_map->end())
                {
                  //std::cout << "side_part = " << side_parts[isp]->name() << std::endl;
                  //throw std::runtime_error("Refiner::connectSides: couldn't find side map part");
                  continue;
                }
              std::string& elem_part_name = found->second;
              for (unsigned iep = 0; iep < elem_parts.size(); iep++)
                {
                  if ( stk::mesh::is_auto_declared_part(*elem_parts[iep]) )
                    continue;
                  if (elem_parts[iep]->name() == elem_part_name)
                    //if (elem_parts[iep] == found->second)
                    {
                      valid = true;
                      break;
                    }
                }
            }
          if (!valid) return false;
        }

      // special case for shells
      int topoDim = UniformRefinerPatternBase::getTopoDim(element_topo);

      bool isShell = false;
      if (topoDim < (int)element->entity_rank())
        {
          isShell = true;
        }
      int spatialDim = m_eMesh.get_spatial_dim();
      if (spatialDim == 3 && isShell && side_elem->entity_rank() == m_eMesh.edge_rank())
        {
          element_nsides = (unsigned) element_topo.getEdgeCount();
        }

      int permIndex = -1;
      int permPolarity = 1;

      unsigned k_element_side = 0;

      // try search
      for (unsigned j_element_side = 0; j_element_side < element_nsides; j_element_side++)
        {
          PerceptMesh::element_side_permutation(*element, *side_elem, j_element_side, permIndex, permPolarity);
          if (permIndex >= 0)
            {
              k_element_side = j_element_side;
              break;
            }
        }

      if (permIndex >= 0)
        {
          mesh::PairIterRelation rels = side_elem->relations(m_eMesh.element_rank());

          if (rels.size() > 1)
            {
              throw std::logic_error("rels.size() > 1");
            }

          if (rels.size())
            {
              stk::mesh::Entity *to_rel = rels[0].entity();
              stk::mesh::RelationIdentifier to_id = rels[0].identifier();
              bool del = m_eMesh.get_bulk_data()->destroy_relation( *to_rel, *side_elem, to_id);
              if (!del)
                throw std::logic_error("connectSides:: destroy_relation failed");
            }

          // special case for shells
          if (isShell)
            {
              // FIXME for 2D
              if (side_elem->entity_rank() == m_eMesh.face_rank())
                {
                  stk::mesh::PairIterRelation elem_sides = element->relations(side_elem->entity_rank());
                  unsigned elem_sides_size= elem_sides.size();
                  //std::cout << "tmp srk found shell, elem_sides_size= " << elem_sides_size << std::endl;
                  if (elem_sides_size == 1)
                    {
                      stk::mesh::RelationIdentifier rel_id = elem_sides[0].identifier();
                      if (rel_id > 1) 
                        throw std::logic_error("connectSides:: logic 1");
                      k_element_side = (rel_id == 0 ? 1 : 0);
                      //std::cout << "tmp srk k_element_side= " << k_element_side << " rel_id= " << rel_id << std::endl;
                    }
                }
            }

          m_eMesh.get_bulk_data()->declare_relation(*element, *side_elem, k_element_side);
          return true;
        }
      else
        {
          // error condition?
          //throw std::runtime_error("fixElementSides2: couldn't find a matching face");
          return false;
        }
    }


    void Refiner::
    fixElementSides2()
    {
      EXCEPTWATCH;
      //std::cout << "fixElementSides2 start... " << std::endl;

      stk::mesh::EntityRank side_rank = m_eMesh.side_rank();

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
          std::cout << "Refiner::fixElementSides2: missing sub-dim break pattern - logic error\n"
            " ---- for this refinement pattern to be able to handle sidesets and edgesets you must provide the sub-dim break pattern\n"
            " ---- or you must set the setIgnoreSideSets() flag " << std::endl;
          throw std::logic_error("Refiner::fixElementSides2: missing sub-dim break pattern - logic error");
          return;
        }

      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.element_rank() );

      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          const unsigned num_elements_in_bucket = bucket.size();
          for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
            {
              stk::mesh::Entity& element = bucket[iElement];
              if (m_eMesh.hasFamilyTree(element) && m_eMesh.isParentElementLeaf(element) && !m_eMesh.isGhostElement(element))
                {
                  stk::mesh::Entity* parent = &element;

                  fixSides(parent);
                }
            }
        }

      //std::cout << "fixElementSides2 ...end " << std::endl;

       fix_side_sets_1();

    }

    void Refiner::fixSides(stk::mesh::Entity *parent)
    {

      std::vector<stk::mesh::Entity*> children;
      VERIFY_OP_ON(m_eMesh.hasFamilyTree(*parent), == , true, "error 101");
      m_eMesh.getChildren(*parent, children);

      // if parent has any side relations, check if any of the sides' children match the parent's children's faces
      mesh::PairIterRelation parent_sides = parent->relations(m_eMesh.side_rank());

      for (unsigned i_parent_side = 0; i_parent_side < parent_sides.size(); i_parent_side++)
        {
          stk::mesh::Entity *parent_side = parent_sides[i_parent_side].entity();

          //VERIFY_OP_ON(parent_side->entity_rank(), ==, side_rank, "side ranks mismatch");

          // parent_side has no children
          if (m_eMesh.isLeafElement(*parent_side))
            {
              // find a child of parent to hook to the parent_side, enforce it must find one, replace old ones
              for (unsigned i_child = 0; i_child < children.size(); i_child++)
                {
                  stk::mesh::Entity *child = children[i_child];

                  connectSides(child, parent_side);
                }          
            }
          else
            {
              std::vector<stk::mesh::Entity*> parent_side_children;
              m_eMesh.getChildren(*parent_side, parent_side_children);

              for (unsigned i_parent_side_child = 0; i_parent_side_child < parent_side_children.size(); i_parent_side_child++)
                {
                  // the child of the parent's side
                  stk::mesh::Entity *parent_side_child = parent_side_children[i_parent_side_child];

                  // loop over each child associated with parent
                  for (unsigned i_child = 0; i_child < children.size(); i_child++)
                    {
                      stk::mesh::Entity *child = children[i_child];

                      connectSides(child, parent_side_child);
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
      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( FAMILY_TREE_RANK );

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
      //std::cout << "tmp P[" << m_eMesh.get_rank() << "] removing family_trees, size() = "  << elements_to_be_destroyed.size() << std::endl;
      removeElements(elements_to_be_destroyed);
    }

    void Refiner::
    removeOldElements(unsigned irank, stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern)
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

      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( rank );

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
                  std::cout << "P[" << m_eMesh.get_rank() << "] removing elements in bucket of parts: " << str << std::endl;
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
                      //m_eMesh.print_entity(std::cout, *element_p);
                    }
                }
            }
        }
      removeElements(elements_to_be_destroyed, irank);

    }

    void Refiner::removeElements(elements_to_be_destroyed_type& elements_to_be_destroyed, unsigned irank)
    {
      elements_to_be_destroyed_type elements_to_be_destroyed_pass2;

      if (m_doProgress)
        {
          std::ostringstream oss; oss << "Delete Original Elements pass [" << 100.0*((double)irank)/((double)m_ranks.size()) << " %]";
          ProgressMeterData pd(ProgressMeterData::INIT, 0.0, oss.str());
          notifyObservers(&pd);
        }
      unsigned num_elem = 0;
      int progress_meter_num_total = elements_to_be_destroyed.size();
      int progress_meter_when_to_post = progress_meter_num_total / m_progress_meter_frequency;
      if (0 == progress_meter_when_to_post)
        progress_meter_when_to_post = 1;
      double d_progress_meter_num_total = progress_meter_num_total;

      for (elements_to_be_destroyed_type::iterator itbd = elements_to_be_destroyed.begin(); itbd != elements_to_be_destroyed.end();  ++itbd)
        {
          stk::mesh::Entity *element_p = *itbd;

          if (m_doProgress && (num_elem % progress_meter_when_to_post == 0) )
            {
              double progress_meter_percent = 100.0*((double)num_elem)/std::max(d_progress_meter_num_total,1.0);
              std::ostringstream oss; oss << "Delete Original Elements pass [" << 100.0*((double)irank)/((double)m_ranks.size()) << " %]";
              ProgressMeterData pd(ProgressMeterData::RUNNING, progress_meter_percent, oss.str());
              notifyObservers(&pd);
            }

          ++num_elem;

          if (0)
            {
              std::cout << "tmp removeElements removing element_p = " << element_p << std::endl;
              if (element_p) std::cout << "tmp removeElements removing id= " << element_p->identifier() << std::endl;
            }

          if ( ! m_eMesh.get_bulk_data()->destroy_entity( element_p ) )
            {
#if UNIFORM_REF_REMOVE_OLD_STD_VECTOR
              elements_to_be_destroyed_pass2.push_back(element_p);
#else
              elements_to_be_destroyed_pass2.insert(element_p);
#endif
              //throw std::logic_error("Refiner::removeElements couldn't remove element");

            }
        }

      if (m_doProgress)
        {
          std::ostringstream oss; oss << "Delete Original Elements pass [" << 100.0*((double)irank)/((double)m_ranks.size()) << " %]";
          ProgressMeterData pd(ProgressMeterData::FINI, 0.0, oss.str());
          notifyObservers(&pd);
        }

      //std::cout << "tmp Refiner::removeElements pass2 size = " << elements_to_be_destroyed_pass2.size() << std::endl;
      for (elements_to_be_destroyed_type::iterator itbd = elements_to_be_destroyed_pass2.begin();
           itbd != elements_to_be_destroyed_pass2.end();  ++itbd)
        {
          stk::mesh::Entity *element_p = *itbd;
          if ( ! m_eMesh.get_bulk_data()->destroy_entity( element_p ) )
            {
              shards::CellTopology cell_topo(stk::percept::PerceptMesh::get_cell_topology(*element_p));
              std::cout << "tmp Refiner::removeElements couldn't remove element in pass2,...\n tmp destroy_entity returned false: cell= " << cell_topo.getName() << std::endl;
              const mesh::PairIterRelation elem_relations = element_p->relations(element_p->entity_rank()+1);
              std::cout << "tmp elem_relations.size() = " << elem_relations.size() << std::endl;

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
    void Refiner::
    renameNewParts(stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern)
    {
      EXCEPTWATCH;
      
      stk::mesh::PartVector toParts = breakPattern->getToParts();
      stk::mesh::PartVector fromParts = breakPattern->getFromParts();
      
      if (DEBUG_RENAME_NEW_PARTS)
        {
          for (unsigned i_part = 0; i_part < toParts.size(); i_part++)
            {
              std::cout << "tmp toParts[i_part]->name() = " << toParts[i_part]->name() << std::endl;
            }
          for (unsigned i_part = 0; i_part < fromParts.size(); i_part++)
            {
              std::cout << " fromParts[i_part]->name() = " << fromParts[i_part]->name()  << std::endl;
            }
        }
      
      if (fromParts.size() == toParts.size())
        {
          for (unsigned i_part = 0; i_part < toParts.size(); i_part++)
            {
              if (DEBUG_RENAME_NEW_PARTS) std::cout << "tmp before: fromPartName= " << fromParts[i_part]->name() 
                                                    << " toPartName= " << toParts[i_part]->name() << std::endl;

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
              //mesh::Part *fromPart = m_eMesh.get_non_const_part(fromPartName);
              mesh::Part *fromPart = fromParts[i_part];
              VERIFY_OP_ON(fromPart, !=, 0, std::string("Refiner::renameNewParts null fromPart found, fromPart= ")+fromPartName);
              std::string * fromPartName_p = const_cast<std::string *> (&fromPart->name());
              *toPartName_p = fromPartName;
              *fromPartName_p = fromPartName + breakPattern->getAppendOriginalString();

              if (DEBUG_RENAME_NEW_PARTS) {
                std::cout << "tmp  after: fromPartName= " << fromParts[i_part]->name() << " toPartName= " << toParts[i_part]->name() << std::endl;
                std::cout << "tmp P[" << m_eMesh.get_rank() << "] fromPartName: " << fromPartName << " part= " << toParts[i_part]->name()
                          << " old part name = " << fromPart->name()
                          << std::endl;
              }
            }
        }
      else
        {
          for (unsigned i_part = 0; i_part < fromParts.size(); i_part++)
            {
              if (DEBUG_RENAME_NEW_PARTS) std::cout << "tmp before: fromPartName= " << fromParts[i_part]->name() << std::endl;
              std::string * fromPartName_p = const_cast<std::string *> (&fromParts[i_part]->name());
              *fromPartName_p = fromParts[i_part]->name() + breakPattern->getAppendOriginalString();
              if (DEBUG_RENAME_NEW_PARTS) {
                std::cout << "tmp  after: fromPartName= " << fromParts[i_part]->name() << std::endl;
              }
            }

          for (unsigned i_part = 0; i_part < toParts.size(); i_part++)
            {
              if (DEBUG_RENAME_NEW_PARTS) std::cout << "tmp before: toPartName= " << toParts[i_part]->name() << std::endl;

              std::string * toPartName_p = const_cast<std::string *> (&toParts[i_part]->name());
              std::string toPartName = toParts[i_part]->name();
              if (toPartName.find(breakPattern->getAppendConvertString()) == std::string::npos)
                continue;
              int len = toPartName.length();
              int clen = breakPattern->getAppendConvertString().length();
              toPartName.erase(len - clen, clen);
              *toPartName_p = toPartName;

              if (DEBUG_RENAME_NEW_PARTS) {
                std::cout << "tmp  after: toPartName= " << toParts[i_part]->name() << std::endl;
              }
            }
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


    void Refiner::set_active_part()
    {
      // deal with parts
      stk::mesh::Part* child_elements_part = m_eMesh.get_non_const_part("refine_active_elements_part");
      stk::mesh::Part* parent_elements_part = m_eMesh.get_non_const_part("refine_inactive_elements_part");
      if (child_elements_part && parent_elements_part)
        {
          //m_eMesh.get_bulk_data()->modification_begin();
          std::vector<stk::mesh::Part*> child_parts(1, child_elements_part);
          std::vector<stk::mesh::Part*> parent_parts(1, parent_elements_part);
          //mesh::Selector in_child_part(*child_elements_part);
          mesh::Selector on_locally_owned_part =  ( m_eMesh.get_fem_meta_data()->locally_owned_part() );

          std::vector<stk::mesh::Entity *> child_entities;
          std::vector<stk::mesh::Entity *> parent_entities;

          const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.element_rank() );

          for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;

              if (on_locally_owned_part(bucket))
                {

                  const unsigned num_entity_in_bucket = bucket.size();
                  for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                    {
                      stk::mesh::Entity& element = bucket[ientity];
                      if (m_eMesh.hasFamilyTree(element) && m_eMesh.isParentElement(element, true))
                        {
                          //if (in_child_part(element))
                            {
                              parent_entities.push_back(&element);
                            }
                        }
                      else
                        {
                          //if (!in_child_part(element))
                            {
                              child_entities.push_back(&element);
                            }
                        }
                    }

                }
            }

          //std::cout << "tmp Refiner::set_active_part: child_entities= " << child_entities.size() << " parent_entities= " << parent_entities.size() << std::endl;
          for (unsigned iv=0; iv < child_entities.size(); iv++)
            {
              m_eMesh.get_bulk_data()->change_entity_parts( *child_entities[iv],   child_parts, parent_parts );
            }
          for (unsigned iv=0; iv < parent_entities.size(); iv++)
            {
              m_eMesh.get_bulk_data()->change_entity_parts( *parent_entities[iv],  parent_parts, child_parts );
            }

          /* for future
          if (!m_doIOSaveInactiveElements) 
            {
              m_eMesh.set_io_omitted_parts(parent_parts);
            }
          */
          //m_eMesh.get_bulk_data()->modification_end();
        }

    }

    //    ========================================================================================================================
    //    ========================================================================================================================
    //    ========================================================================================================================


    void Refiner::check_db(std::string msg)
    {
      std::cout << "P[" << m_eMesh.get_rank() << "] tmp check_db msg= " << msg << std::endl;
      check_db_entities_exist(msg);
      check_db_ownership_consistency(msg);
      std::cout << "P[" << m_eMesh.get_rank() << "] tmp check_db done msg= " << msg << std::endl;
      //check_db_hanging_nodes();
    }


    void Refiner::check_db_entities_exist(std::string msg)
    {
      std::vector<stk::mesh::EntityRank> ranks_to_check;
      ranks_to_check.push_back(m_eMesh.node_rank());
      ranks_to_check.push_back(m_eMesh.element_rank());
      for (unsigned irank=0; irank < ranks_to_check.size(); irank++)
        {

          const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( ranks_to_check[irank] );

          for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;

              const unsigned num_entity_in_bucket = bucket.size();
              for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                {
                  stk::mesh::Entity& element = bucket[ientity];
                  stk::mesh::Entity *element_1 = m_eMesh.get_bulk_data()->get_entity(ranks_to_check[irank], element.identifier());
                  if (&element != element_1 || element.identifier() != element_1->identifier())
                    {
                      std::cout << "msg= " << msg << " error element, element_1, ids= " 
                                << &element << " " << element_1 << " " << element.identifier() << " " << element_1->identifier() << std::endl;
                      throw std::logic_error("check_db_entities_exist:: error #1");
                    }

                }
            }
        }
    }

    void Refiner::check_db_ownership_consistency(std::string msg)
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
                throw std::logic_error("check_db_ownership_consistency:: error #1 msg= "+msg);

              stk::mesh::Entity * owning_element = m_eMesh.get_bulk_data()->get_entity(owning_elementRank, owning_elementId);

              if (!owning_element)
                {
                  std::cout << "P[" << m_eMesh.get_rank() << "] "
                            << " error check_db_ownership_consistency: msg= " << msg << " owning_elementId= " << owning_elementId << std::endl;
                  throw std::logic_error("check_db_ownership_consistency:: error #2, msg= "+msg);
                }

              if (stk::mesh::EntityLogDeleted == owning_element->log_query() )
                {
                  std::cout << "P[" << m_eMesh.get_rank() << "] "
                            << " error check_db_ownership_consistency: EntityLogDeleted, msg= " << msg << " owning_elementId= " << owning_elementId << std::endl;
                  throw std::logic_error("check_db_ownership_consistency:: error #2.0, msg= "+msg);
                }


              if (owning_element->identifier() != owning_elementId)
                {
                  std::cout << "msg= " << msg << " check_db_ownership_consistency error element, element_1, ids= " 
                            << owning_elementId << " " << owning_element->identifier() << std::endl;
                  throw std::logic_error("check_db_ownership_consistency:: error #1.1, msg= "+msg);
                }

              if (!m_eMesh.isGhostElement(*owning_element))
                {
                
                  for (unsigned inode = 0; inode < nodeIds_onSE.size(); inode++)
                    {
                      stk::mesh::Entity *node = nodeIds_onSE[inode];
                      if (!node)
                        throw std::logic_error("check_db_ownership_consistency:: error #3, msg= "+msg);

                      stk::mesh::Entity * node1 = m_eMesh.get_bulk_data()->get_entity(stk::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE.m_entity_id_vector[inode]);
                      if (!node1)
                        throw std::logic_error("check_db_ownership_consistency:: error #3a, msg= "+msg);

                      stk::mesh::Entity * node2 = m_eMesh.get_bulk_data()->get_entity(stk::mesh::fem::FEMMetaData::NODE_RANK, node->identifier() );
                      if (!node2)
                        throw std::logic_error("check_db_ownership_consistency:: error #3b, msg= "+msg);
                      if (node != node2)
                        throw std::logic_error("check_db_ownership_consistency:: error #3c, msg= "+msg);
              
                    }
                }
            }
        }
    }

    void Refiner::check_db_hanging_nodes()
    {
      std::set<stk::mesh::Entity *> node_set;

      // check for hanging nodes - ensure all parents have their sub-entities in the DB
      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.element_rank() );

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
                      shards::CellTopology cell_topo(cell_topo_data);
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
