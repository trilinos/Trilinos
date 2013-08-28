
/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <utility>
#include <stdint.h>
#include <map>
#include <list>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/RunEnvironment.hpp>
#include <stk_percept/ProgressMeter.hpp>
#include <stk_percept/Histograms.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/util/memory_util.hpp>

#include <stk_mesh/base/MemoryUsage.hpp>

#include <stk_adapt/RefinerUtil.hpp>

#include <stk_adapt/UniformRefiner.hpp>
#include <stk_adapt/UniformRefinerPattern.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/static_assert.hpp>
#include <Ioss_Utils.h>
#include <Ioss_SerializeIO.h>

#include <stk_adapt/SerializeNodeRegistry.hpp>

#include <stk_percept/mesh/mod/smoother/SpacingFieldUtil.hpp>

#if defined( STK_PERCEPT_HAS_GEOMETRY )
#include <stk_percept/mesh/geometry/kernel/GeometryKernelOpenNURBS.hpp>
#include <stk_percept/mesh/geometry/stk_geom/LocalCubicSplineFit.hpp>
#endif

#define ALLOW_MEM_TEST 1
#define DEBUG_ADAPT_MAIN 0

#include "AdaptMain.hpp"
namespace stk {

  namespace adapt {

    static int s_spatialDim=0;

    typedef size_t MemorySizeType;

    static double MegaByte(MemorySizeType x) { return  ((double)x/1024.0/1024.0); }

    struct MemoryMultipliers
    {
      MemorySizeType num_hex8;
      MemorySizeType num_tet4;
      MemorySizeType num_nodes;

      typedef MemorySizeType MemMultType;
      MemMultType mult_hex8;
      MemMultType mult_tet4;
      MemMultType mult_nodes;

      MemoryMultipliers(MemMultType mult_hex8=1490, MemMultType mult_tet4=702, MemMultType mult_nodes=0):
        //MemoryMultipliers(MemMultType mult_hex8=381, MemMultType mult_tet4=538, MemMultType mult_nodes=1017):
        num_hex8(0ul),
        num_tet4(0ul),
        num_nodes(0ul),
        mult_hex8(mult_hex8),
        mult_tet4(mult_tet4),
        mult_nodes(mult_nodes)
      {
      }

      void read_simple(std::string file_name)
      {
        std::ifstream file(file_name.c_str());
        if (file.good())
          file >> mult_hex8 >> mult_tet4 >> mult_nodes;
        //std::string line1;
        //file >> line1;
        //std::cout << "mult_hex8= " << mult_hex8 << " mult_tet4= " << mult_tet4 << " mult_nodes=" << mult_nodes << std::endl;
      }

      MemorySizeType estimate_memory()
      {
        return mult_nodes*num_nodes + mult_hex8*num_hex8 + mult_tet4*num_tet4;
      }

      MemorySizeType estimate_memory(std::vector<RefinementInfoByType>& refInfo, bool use_new=true)
      {
        num_hex8=0ul;
        num_tet4=0ul;
        num_nodes=0ul;

        for (unsigned i = 0; i < refInfo.size(); i++)
          {
            num_nodes= refInfo[0].m_numNewNodes;
            //std::cout << "irank, rank, m_numNewNodes, m_numNewElems= " << i << " " << refInfo[i].m_rank << " " << refInfo[i].m_numNewNodes
            //<< " " << refInfo[i].m_numNewElemsLast
            //<< std::endl;

            //             if (refInfo[i].m_rank == 0)
            //               {
            //                 num_nodes += refInfo[i].m_numNewNodes;
            //               }
            //             else
              {
                switch(refInfo[i].m_topology.getKey())
                  {
                  case shards::Hexahedron<8>::key:
                    if (use_new)
                      num_hex8 += refInfo[i].m_numNewElemsLast;
                    else
                      num_hex8 += refInfo[i].m_numOrigElemsLast;

                    break;
                  case shards::Tetrahedron<4>::key:
                    if (use_new)
                      num_tet4 += refInfo[i].m_numNewElemsLast;
                    else
                      num_tet4 += refInfo[i].m_numOrigElemsLast;
                    break;
                  default:
                    break;
                  }
              }
          }

        return estimate_memory();
      }

      static void process_estimate(MemorySizeType tot_mem, PerceptMesh& eMesh, std::vector<RefinementInfoByType>& refInfo, std::string memory_multipliers_file, std::string input_file, bool use_new=true)
      {
        //const stk::ParallelMachine& comm = eMesh.get_bulk_data()->parallel();

        // this is a data gather pass
        if (tot_mem)
          {
            /*
              mesh::Selector sel_locally_owned(eMesh.get_fem_meta_data()->locally_owned_part());
              mesh::Selector sel_globally_shared(eMesh.get_fem_meta_data()->globally_shared_part());
              mesh::Selector sel_universal(eMesh.get_fem_meta_data()->universal_part());

              std::vector<unsigned> count ;
              stk::mesh::count_entities( sel_universal, *eMesh.get_bulk_data(), count );

              unsigned nnodes = count[0];

              stk::ParallelMachine pm = eMesh.get_bulk_data()->parallel();
              stk::all_reduce( pm, stk::ReduceSum<1>( &nnodes ) );
            */
            MemoryMultipliers memMults;
            // FIXME, here's where we would read in some values for memMults from memory_multipliers_file
            if (memory_multipliers_file.size())
              memMults.read_simple(memory_multipliers_file);
            RefinementInfoByType::countCurrentNodes(eMesh, refInfo);
            MemorySizeType estMem = memMults.estimate_memory(refInfo, use_new);
            //std::cout << "tmp srk tot_mem = " << MegaByte(tot_mem) << " estMem= " << MegaByte(estMem) << std::endl;
            if (eMesh.get_rank() == 0)
              {
                if (1)
                std::cout << "MemEst: num_nodes= " << memMults.num_nodes << " num_tet4= " << memMults.num_tet4 << " num_hex8= " << memMults.num_hex8 << " actualMem[MB]= " << MegaByte(tot_mem)
                          << " estMem[MB]= " << MegaByte(estMem)
                          << " mult_hex8= " << memMults.mult_hex8 << " mult_tet4= " << memMults.mult_tet4 << " mult_nodes=" << memMults.mult_nodes << std::endl;

                std::cout << "(*MemEstMM: {nn,ntet,nhex,mem,estMem} " << input_file << " *) ,{" << memMults.num_nodes << ", " << memMults.num_tet4 << ", " << memMults.num_hex8 << ", " << MegaByte(tot_mem)
                          << ", " << MegaByte(estMem) << "}" << std::endl;
              }

          }
        else
          {
            // this is an estimate multipliers pass (computes memory using current multipliers)
            MemoryMultipliers memMults;
            // FIXME, here's where we would read in some values for memMults from memory_multipliers_file
            if (memory_multipliers_file.size())
              memMults.read_simple(memory_multipliers_file);
            RefinementInfoByType::countCurrentNodes(eMesh, refInfo);
            MemorySizeType estMem = memMults.estimate_memory(refInfo);
            //std::cout << "tmp srk tot_mem = " << MegaByte(tot_mem) << " estMem= " << MegaByte(estMem) << std::endl;
            if (eMesh.get_rank() == 0)
              {
                std::cout << "MemEst: num_nodes= " << memMults.num_nodes << " num_tet4= " << memMults.num_tet4 << " num_hex8= " << memMults.num_hex8
                  //<< " memory[MB]= " << MegaByte(tot_mem)
                          << " estimatedMem[MB]= " << MegaByte(estMem)
                          << " mult_hex8= " << memMults.mult_hex8 << " mult_tet4= " << memMults.mult_tet4 << " mult_nodes=" << memMults.mult_nodes << std::endl;

                if (0)
                  std::cout << "(*MemEstMM: " << input_file << " *) ,{" << memMults.num_nodes << ", " << memMults.num_tet4 << "," << memMults.num_hex8 << "," << MegaByte(tot_mem)
                            << ", " << MegaByte(estMem) << "}" << std::endl;
              }

          }
      }
    };


    stk::mesh::Part& create_clone_part(PerceptMesh& eMesh, const std::string& clone_name, const stk::mesh::EntityRank rank_to_clone, bool make_part_io_part=true )
    {
      stk::mesh::Part& clone = eMesh.get_fem_meta_data()->declare_part(clone_name, rank_to_clone);
      if (make_part_io_part && clone.attribute<Ioss::GroupingEntity>() == NULL) {
        stk::io::put_io_part_attribute(clone);
      }
      return clone;
    }

    stk::mesh::Part& clone_part_bulk(PerceptMesh& eMesh, const stk::mesh::Part& part, const std::string& clone_name, const stk::mesh::EntityRank rank_to_clone, bool make_part_io_part=true )
    {
      stk::mesh::Part * clone_p = eMesh.get_fem_meta_data()->get_part(clone_name);
      if (!clone_p)
        throw std::runtime_error("AdaptMain::clone_part: no part named: "+clone_name);
      stk::mesh::Part& clone = *clone_p;

      stk::mesh::Selector this_part(part);
      std::vector<stk::mesh::Entity> entities;
      stk::mesh::PartVector add_parts(1,&clone), remove_parts;
      const std::vector<stk::mesh::Bucket*> & entity_buckets = eMesh.get_bulk_data()->buckets( rank_to_clone );
      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = entity_buckets.begin() ; k != entity_buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          if (this_part(bucket))
            {
              const unsigned num_entities_in_bucket = bucket.size();
              for (unsigned i_entity = 0; i_entity < num_entities_in_bucket; i_entity++)
                {
                  stk::mesh::Entity entity = bucket[i_entity];
                  entities.push_back(entity);
                }
            }
        }

      for (unsigned ii=0; ii < entities.size(); ii++)
        {
          eMesh.get_bulk_data()->change_entity_parts( entities[ii], add_parts, remove_parts );
        }
      return clone;
    }

    /// returns true if closed found
    bool get_sorted_curve_node_entities(PerceptMesh& eMesh, stk::mesh::Part& part, std::vector<stk::mesh::Entity>& sorted_entities)
    {
      bool debug_print = false;
#define APRINTLN(a) do { if (debug_print) std::cout << #a << " = " << a << std::endl; } while(0)
#define APRINTLN2(a,b) do { if (debug_print) std::cout << #a << " = " << a << " " << #b << " = " << b << std::endl; } while(0)

      bool closed = false;
      if (debug_print) std::cout << "AdaptMain::get_sorted_curve_node_entities: processing part = " << part.name() << std::endl;
      VERIFY_OP_ON(part.primary_entity_rank(), ==, eMesh.edge_rank(), "bad part");

      typedef std::set<stk::mesh::Entity, stk::mesh::EntityLess> SetOfEntities;
      SetOfEntities node_set(*eMesh.get_bulk_data());

      //typedef std::list<stk::mesh::Entity, stk::mesh::EntityLess> ListOfEntities;
      //ListOfEntities node_list(*eMesh.get_bulk_data());
      typedef std::list<stk::mesh::Entity> ListOfEntities;
      ListOfEntities node_list;

      stk::mesh::Selector this_part(part);
      stk::mesh::Entity edge_first = stk::mesh::Entity();
      const std::vector<stk::mesh::Bucket*> & edge_buckets = eMesh.get_bulk_data()->buckets( eMesh.edge_rank() );
      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = edge_buckets.begin() ; k != edge_buckets.end() ; ++k )
        {
          if (this_part(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              edge_first = bucket[0];
              break;
            }
        }

      int node_set_nnodes=0;
      if (debug_print)
        {
          const std::vector<stk::mesh::Bucket*> & node_buckets = eMesh.get_bulk_data()->buckets( eMesh.node_rank() );
          for ( std::vector<stk::mesh::Bucket*>::const_iterator k = node_buckets.begin() ; k != node_buckets.end() ; ++k )
            {
              if (this_part(**k))
                {
                  stk::mesh::Bucket & bucket = **k ;
                  const unsigned num_nodes_in_bucket = bucket.size();
                  node_set_nnodes += num_nodes_in_bucket;
                  for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
                    {
                      stk::mesh::Entity node = bucket[iNode];
                      node_set.insert(node);
                    }
                }
            }
        }

      const MyPairIterRelation edge_nodes(*eMesh.get_bulk_data(), edge_first, eMesh.node_rank() );
      VERIFY_OP_ON(edge_nodes.size(), ==, 2, "bad edge");

      for (unsigned idir=0; idir < 2; idir++)
        {
          stk::mesh::Entity current_edge = edge_first;
          stk::mesh::Entity last_node = stk::mesh::Entity();
          bool continue_proc=true;
          while(continue_proc)
            {
              const MyPairIterRelation current_edge_nodes(*eMesh.get_bulk_data(), current_edge, eMesh.node_rank() );
              unsigned JDIR=0;
              if (current_edge == edge_first)
                JDIR = idir;
              else
                {
                  for (unsigned jdir=0; jdir < 2; jdir++)
                    {
                      if (current_edge_nodes[jdir].entity() != last_node)
                        {
                          JDIR = jdir;
                          break;
                        }
                    }
                }
              last_node = current_edge_nodes[JDIR].entity();

              bool is_not_in_list = std::find(node_list.begin(), node_list.end(), last_node) == node_list.end();
              VERIFY_OP_ON(is_not_in_list, ==, true, "bad node list");

              bool is_node_in_part = this_part(eMesh.bucket(last_node));
              VERIFY_OP_ON(is_node_in_part, ==, true, "bad node not in part");
              if (idir==0)
                node_list.push_back(last_node);
              else
                node_list.push_front(last_node);

              is_not_in_list = std::find(node_list.begin(), node_list.end(), last_node) == node_list.end();
              VERIFY_OP_ON(is_not_in_list, ==, false, "bad node list 2");

              const MyPairIterRelation last_node_edges(*eMesh.get_bulk_data(), last_node, eMesh.edge_rank() );
              bool found_new_edge = false;
              for (unsigned kdir=0; kdir < 2; kdir++)
                {
                  if (last_node_edges[kdir].entity() != current_edge && this_part(eMesh.bucket(last_node_edges[kdir].entity() )))
                    {
                      current_edge = last_node_edges[kdir].entity();
                      // check for closed
                      if (idir == 0 && current_edge == edge_first)
                        {
                          if (debug_print) std::cout << "AdaptMain::get_sorted_curve_node_entities: found closed condition..." << std::endl;
                          VERIFY_OP_ON(last_node, ==, edge_nodes[1].entity(), "integrity check failed");
                          node_list.push_front(edge_nodes[1].entity());
                          ++idir; //force loop exit
                        }
                      else
                        {
                          found_new_edge = true;
                        }
                      break;
                    }
                }
              if (!found_new_edge)
                continue_proc = false;

              //VERIFY_OP_ON(found, ==, true, "hmmm");
            }
        }
      APRINTLN(node_list.size());
      APRINTLN2(node_set_nnodes,node_set.size());

      if (node_list.front() == node_list.back())
        {
          // closed case
          closed = true;
          sorted_entities.resize(0);
          sorted_entities.assign(node_list.begin(), node_list.end());
        }
      else
        {
          sorted_entities.resize(0);
          sorted_entities.assign(node_list.begin(), node_list.end());

          if (debug_print) VERIFY_OP_ON(node_list.size(), ==, node_set.size(), "node set/list error for part= "+part.name());
        }
      return closed;
    }

    static void fit_geometry_create_parts_meta(PerceptMesh& eMesh, stk::mesh::PartVector& surface_parts, stk::mesh::PartVector& topo_parts)
    {
      //bool debug_print = false;

      const stk::mesh::PartVector & parts = eMesh.get_fem_meta_data()->get_parts();
      surface_parts.resize(0);
      topo_parts.resize(0);

      int n_topo = 10000;
      unsigned nparts = parts.size();
      for (unsigned ipart=0; ipart < nparts; ipart++)
        {
          stk::mesh::Part& part = *parts[ipart];
          if (stk::mesh::is_auto_declared_part(part) )
            //|| (part.name().find("oldElem") != std::string::npos))
            continue;

          const STK_Adapt_Auto_Part *side_auto_part = part.attribute<STK_Adapt_Auto_Part>();
          if (side_auto_part)
            continue;

          if (part.primary_entity_rank() != eMesh.edge_rank())
            continue;
          if (part.subsets().size() == 0)  // skip parts like surface_quad4_edge2_4
            continue;

          std::string clone_name = "tbc_curve_block_"+toString(n_topo++);
          stk::mesh::Part& clone = create_clone_part(eMesh, clone_name, eMesh.node_rank());
          topo_parts.push_back(&clone);
          surface_parts.push_back(&part);
        }
    }


    static void fit_geometry(PerceptMesh& eMesh, std::string filename, stk::mesh::PartVector& surface_parts, stk::mesh::PartVector& topo_parts)
    {
#if defined(STK_PERCEPT_HAS_GEOMETRY)
      using namespace stk::geom;

      bool debug_print = false;
      SplineFit::s_debug_print = debug_print;
      ON::Begin();

      FILE* file = ON::OpenFile( filename.c_str(), "wb" );
      if (!file)
        throw std::runtime_error("couldn't open file: "+filename+" in fit_geometry");

      // allow additional parts to be added
      ONX_Model model;

      // set revision history information
      model.m_properties.m_RevisionHistory.NewRevision();

      // set application information
      model.m_properties.m_Application.m_application_name = "Percept";
      model.m_properties.m_Application.m_application_URL = "http://www.sandia.gov";
      model.m_properties.m_Application.m_application_details = "Fit cubic splines.";

      // some notes
      model.m_properties.m_Notes.m_notes = "notes";
      model.m_properties.m_Notes.m_bVisible = true;

      // file settings (units, tolerances, views, ...)
      model.m_settings.m_ModelUnitsAndTolerances.m_unit_system = ON::inches; //ON::meters; //ON::inches;
      model.m_settings.m_ModelUnitsAndTolerances.m_absolute_tolerance = 0.001;
      model.m_settings.m_ModelUnitsAndTolerances.m_angle_tolerance = ON_PI/180.0; // radians
      model.m_settings.m_ModelUnitsAndTolerances.m_relative_tolerance = 0.01; // 1%

      // layer table
      {
        // OPTIONAL - define some layers
        ON_Layer layer[1];

        layer[0].SetLayerName("Default");
        layer[0].SetVisible(true);
        layer[0].SetLocked(false);
        layer[0].SetColor( ON_Color(0,0,0) );

        model.m_layer_table.Append(layer[0]);
      }

      // object table
      {
        //const stk::mesh::PartVector & parts = eMesh.get_fem_meta_data()->get_parts();
        //stk::mesh::PartVector surface_parts;
        //stk::mesh::PartVector topo_parts;

        eMesh.get_bulk_data()->modification_begin();
        unsigned nparts = topo_parts.size();
        for (unsigned ipart=0; ipart < nparts; ipart++)
          {
            clone_part_bulk(eMesh, *surface_parts[ipart], topo_parts[ipart]->name(), eMesh.node_rank());
          }
        eMesh.get_bulk_data()->modification_end();

        nparts = topo_parts.size();
        for (unsigned ipart=0; ipart < nparts; ipart++)
          {
            stk::mesh::Part& part = *topo_parts[ipart];

            std::vector<stk::mesh::Entity> sorted_entities;
            bool isClosed = get_sorted_curve_node_entities(eMesh, *surface_parts[ipart], sorted_entities);

            LocalCubicSplineFit cf;
            if (isClosed) {
              cf.setIsPeriodic(true);  // default is to assume a smooth seam at the closed/repeated node
            }
            int n = sorted_entities.size();
            Vectors2D Q(n);
            for (int i=0; i < n; i++)
              {
                double *cdata = eMesh.field_data(eMesh.get_coordinates_field(), sorted_entities[i]);
                Q[i] = Vector2D(cdata[0], cdata[1]);
              }

            DPRINTLN(Q);
            ON_Curve *curve = cf.fit(Q);
            if (debug_print)
              {
                std::cout << "Part = " << part.name() << std::endl;
                cf.print();
              }

            if ( curve->IsValid() )
              {
                ONX_Model_Object& mo = model.m_object_table.AppendNew();
                mo.m_object = curve;
                mo.m_bDeleteObject = true;
                mo.m_attributes.m_layer_index = 0;
                mo.m_attributes.m_name = part.name().c_str();
              }
            else
              delete curve;
          }
      }

      // start section comments
      const char* sStartSectionComment = __FILE__ "PerceptMesh::fit_geometry" __DATE__;

      ON_BinaryFile archive( ON::write3dm, file );

      // Set uuid's, indices, etc.
      model.Polish();
      // writes model to archive
      int version = 5; // File can be read by Rhino 5
      ON_TextLog error_log;
      bool ok = model.Write(archive, version, sStartSectionComment, &error_log );
      VERIFY_OP_ON(ok,==,true,"failed write of 3dm file");

      ON::End();

#endif
    }


    static MemorySizeType memory_dump(int dump_level, const stk::ParallelMachine& comm, stk::mesh::BulkData& bulkData, NodeRegistry* node_reg, std::string msg)
    {
      MemorySizeType returned_total_memory;
      MemorySizeType rss_current_0 = 0;
      MemorySizeType rss_high_water_mark_0 = 0;
      get_memory_usage(rss_current_0, rss_high_water_mark_0);

      const MemorySizeType MB = 1024*1024;

      stk::all_reduce( comm, stk::ReduceSum<1>( &rss_current_0 ) );
      stk::all_reduce( comm, stk::ReduceSum<1>( &rss_high_water_mark_0 ) );

      stk::mesh::MemoryUsage mem_usage;
      stk::mesh::compute_memory_usage(bulkData, mem_usage);

      MemorySizeType node_reg_mem = 0;
      if (node_reg)
        node_reg_mem = node_reg->get_memory_usage();
      MemorySizeType node_reg_mem_sum = node_reg_mem;
      MemorySizeType node_reg_mem_max = node_reg_mem;
      stk::all_reduce( comm, stk::ReduceSum<1>( &node_reg_mem_sum ) );
      stk::all_reduce( comm, stk::ReduceMax<1>( &node_reg_mem_max ) );

      MemorySizeType mem_total_bytes_sum = mem_usage.total_bytes;
      MemorySizeType mem_total_bytes_max = mem_usage.total_bytes;
      stk::all_reduce( comm, stk::ReduceSum<1>( &mem_total_bytes_sum ) );
      stk::all_reduce( comm, stk::ReduceMax<1>( &mem_total_bytes_max ) );
      if (bulkData.parallel_rank() == 0)
        {
          if (dump_level > 1)

            {
              std::cout << "P[" << bulkData.parallel_rank() << "] AdaptMain::memory_dump stk_mesh counted memory usage at stage [" << msg << "] "
                " parallel sum, max memory [MB]= " << ((double)mem_total_bytes_sum)/MB << " , " << ((double)mem_total_bytes_max)/MB << std::endl;
              if (dump_level > 2)
                stk::mesh::print_memory_usage(mem_usage, std::cout);

              std::cout << "P[" << bulkData.parallel_rank() << "] AdaptMain::memory_dump rss total (sum all proc) at stage [" << msg << "] = "
                        << " current high-water-mark [MB]= " << ((double)rss_current_0)/MB
                        << " , " << ((double)rss_high_water_mark_0)/MB
                        << std::endl;
            }

          {
            std::cout << "AdaptMain::memory_dump summary for " << msg << " : stk_mesh [sum], NodeRegistry [sum], rss[sum] [MB] "
                      << ((double)mem_total_bytes_sum)/MB << " , "
                      << ((double)node_reg_mem_sum)/MB << " , "
                      << ((double)rss_current_0)/MB
                      << std::endl;
          }

        }
      if (rss_current_0)
        returned_total_memory = rss_current_0;
      else
        returned_total_memory = mem_total_bytes_sum+node_reg_mem_sum;

      return returned_total_memory;
    }

    //extern void test_memory(int, int);
    void test_memory(percept::PerceptMesh& eMesh, MemorySizeType n_elements, MemorySizeType n_nodes)
    {
      vector<stk::mesh::Entity> new_elements;
      vector<stk::mesh::Entity> new_nodes;

      eMesh.get_bulk_data()->modification_begin();

      std::cout << "creating " << n_elements << " elements..." <<std::endl;
      eMesh.createEntities( stk::mesh::MetaData::ELEMENT_RANK, n_elements, new_elements);
      std::cout << "... done creating " << n_elements << " elements" << std::endl;

      std::cout << "creating " << n_nodes << " nodes..." <<std::endl;
      eMesh.createEntities( stk::mesh::MetaData::NODE_RANK, n_nodes, new_nodes);
      std::cout << "... done creating " << n_nodes << " nodes" << std::endl;

      MemorySizeType num_prints = std::min(static_cast<MemorySizeType>(100ul), n_elements);
      MemorySizeType print_mod = n_elements/num_prints;
      MemorySizeType i_node = 0;
      int n_node_per_element = 8;
      for (MemorySizeType i_element = 0; i_element < n_elements; i_element++)
        {
          if (!i_element || (i_element % print_mod == 0))
            {
              std::cout << "declare_relation for i_element = " << i_element << " [" << n_elements << "] = " << ((double)i_element)/((double)n_elements)*100 << "%"
                        << std::endl;
            }
          stk::mesh::Entity element = new_elements[i_element];

          for (int j_node = 0; j_node < n_node_per_element; j_node++)
            {
              stk::mesh::Entity node = new_nodes[i_node];

              eMesh.get_bulk_data()->declare_relation(element, node, j_node);

              i_node++;
              if (i_node >= n_nodes-1)
                i_node = 0;
            }
        }

      std::cout << " doing modification_end ... " << std::endl;
      eMesh.get_bulk_data()->modification_end();
      std::cout << " done modification_end ... " << std::endl;


    }

    static void checkInput(std::string option, std::string value, std::string allowed_values, RunEnvironment& run_environment)
    {
      //if (value.length() == 0) return;
      std::vector<std::string> vals = Util::split(allowed_values, ", ");
      for (unsigned i = 0; i < vals.size(); i++)
        {
          if (vals[i] == value)
            return;
        }
      std::ostringstream oss;
      oss << "\nCommand line syntax error: bad option for " << option << " (= " << value << ") \n allowed values = " << allowed_values;
      //throw std::runtime_error(oss.str());
      {
        std::cout << oss.str() << std::endl;
        run_environment.printHelp();
        std::cout << oss.str() << std::endl;
        exit(1);
      }

    }

    static void print_simple_usage(int argc, char **argv)
    {
      std::cout << "AdaptMain::print_simple_usage number of arguments = " << argc << std::endl;
      for (int i = 0; i < argc; i++)
        {
          std::cout << "AdaptMain::print_simple_usage arg[" << i << "]= " << argv[i] << std::endl;
        }
      std::cout << "usage: exe_name [convert|enrich|refine] input_file_name [output_file_name] [number_refines]" << std::endl;
    }

    int adapt_main(int argc, char **argv) ;
    int adapt_main_full_options(int argc, char **argv) ;

    // FIXME
    static int check_for_simple_options(int argc, char **argv)
    {
      int simple = 0;
      for (int i = 1; i < argc; i++)
        {
          if (std::string(argv[i]) == "refine" || std::string(argv[i]) == "enrich" || std::string(argv[i]) == "convert")
            return i;
        }
      return simple;
    }

    static bool debug = false;
    int adapt_main_simple_options(int argc_in, char **argv_in)
    {

      // format: exe_name [convert|enrich|refine] input_file_name
      int simple_options_index = check_for_simple_options(argc_in, argv_in);

      //       if (!simple_options_index)
      //         {
      //           print_simple_usage();
      //           return 1;
      //         }

      int argc = argc_in;
      if (argc != 2 + simple_options_index && argc != 3 + simple_options_index && argc != 4 + simple_options_index )
        {
          print_simple_usage(argc_in, argv_in);
          return 1;
        }

      std::vector<std::string> argv(argc);
      for (int i = 0; i < argc; i++)
        {
          argv[i] = (const char *)argv_in[i];
        }
      if (debug)
        std::cout << "argc = " << argc << " argv= \n" << argv << std::endl;

      std::string exe_name        = argv[0];
      std::string option          = argv[0+simple_options_index];
      if (option != "refine" && option != "enrich" && option != "convert")
        {
          print_simple_usage(argc_in, argv_in);
          return 1;
        }
      std::string input_mesh = argv[1+simple_options_index];
      std::string number_refines = (argc == 3+simple_options_index? argv[2+simple_options_index] : "1");
      int nref=0;

      bool isInt = false;
      try {
        nref = boost::lexical_cast<int>(number_refines);
        (void)nref;
        isInt = true;
      }
      catch( ... )
        {
        }

      std::string output_mesh = input_mesh;
      std::string extension = input_mesh.substr(input_mesh.length()-2,input_mesh.length());
      if (debug) std::cout << " extension= " << extension << std::endl;
      std::string new_ext = "";
      new_ext += "_";
      if (option == "refine")
        new_ext += "refined_"+number_refines+extension;
      else
        new_ext += option+"ed_"+extension;

      if (!isInt && (argc == 3+simple_options_index))
        {
          output_mesh = number_refines;
          number_refines = (argc == 4+simple_options_index? argv[3+simple_options_index] : "1");
          //std::cout << "tmp output_mesh= " << output_mesh << std::endl;
          //std::cout << "tmp number_refines= " << number_refines << std::endl;
        }
      else
        {
          Util::replace(output_mesh, extension, new_ext);
        }
      if (debug) std::cout << " output_mesh= " << output_mesh << std::endl;

      std::vector<std::string> argv_new;
      for (int i = 0; i < simple_options_index; i++)
        argv_new.push_back(argv[i]);
      argv_new.push_back("--input_mesh="+input_mesh);
      argv_new.push_back("--output_mesh="+output_mesh);
      if (option == "refine")
        argv_new.push_back("--refine=DEFAULT");
      else if (option == "enrich")
        argv_new.push_back("--enrich=DEFAULT");
      else
        argv_new.push_back("--convert=Hex8_Tet4_24");
      argv_new.push_back("--load_balance=1");
      argv_new.push_back("--remove_original_elements=1");
      argv_new.push_back("--number_refines="+number_refines);

      if ( debug) std::cout << "new argv = \n" << argv_new << std::endl;
      int argc_new = argv_new.size();
      char **argv_new_cstr = new char*[argc_new];
      for (int i = 0; i < argc_new; i++)
        {
          argv_new_cstr[i] = (char *)argv_new[i].c_str();
        }
      int ret_val = adapt_main_full_options(argc_new, argv_new_cstr);
      delete[] argv_new_cstr;
      return ret_val;
    }

    static void dump_args(int argc, char **argv)
    {
      std::cout << "argc = " << argc << std::endl;
      for (int i = 0; i < argc; i++)
        {
          std::cout << "argv[" << i << "]= " << argv[i] << std::endl;
        }
    }

    static void check_args(int argc, char **argv)
    {
      std::string errors="";
      for (int i = 1; i < argc; i++)
        {
          std::string av(argv[i]);
          if (av == "--help" || av == "--Help" || av == "-h" || av == "--h" || av == "-H" || av == "--H") continue;
          size_t equal_pos = av.find("=",0);
          if (equal_pos == std::string::npos)
            {
              errors += "ERROR in options: no '=' found in option: "+av+"\n";
            }
          if (av.length() < equal_pos+2)
            {
              errors += "ERROR in options: no value given after '=', found in option: "+av+"\n";
            }
        }
      if (errors.length())
        {
          std::cout << "ERRORS found in options: debug dump of options= " << std::endl;
          dump_args(argc, argv);
          std::cout << errors << std::endl;
          throw std::runtime_error(errors);
        }
    }

    int adapt_main(int argc, char **argv)
    {
      if (debug)
        dump_args(argc, argv);
      // allow positional arguments, etc.
      if (check_for_simple_options(argc, argv))
        return adapt_main_simple_options(argc, argv);
      else
        {
          check_args(argc, argv);
          return adapt_main_full_options(argc, argv);
        }
    }

    int adapt_main_full_options(int argc, char **argv)
    {
      EXCEPTWATCH;
      bool debug_re = false;

      RunEnvironment run_environment(&argc, &argv, debug_re);
      unsigned p_rank = stk::parallel_machine_rank(run_environment.m_comm);
      unsigned p_size = stk::parallel_machine_size(run_environment.m_comm);

      bool found_Help = false;
      for (int i = 0; i < argc; ++i) {
        const std::string s(argv[i]);
        if ( s == "-H" || s == "-Help" || s == "--Help" || s == "--H") {
            found_Help = true;
            //std::cout << "Found Help:: Usage: " << (*argv)[0] << " [options...]" << std::endl;
          }
      }

      std::string options_description_desc = "stk_adapt options";

      // NOTE: Options --directory --output-log --runtest are handled/defined in RunEnvironment
      std::string input_mesh="";
      std::string input_geometry="";
      std::string fit_geometry_file="";
      std::string output_mesh="";
      std::string block_name_inc = "";
      std::string block_name_exc = "";

      // for Salinas
#if defined(STK_BUILT_IN_SIERRA)
      std::string rbar_blocks= "";
#endif
      // for Salinas and other codes
      //std::string ignore_blocks = "";
      // just use block_name_inc to exclude....

      std::string convert="";
      std::string refine="";
      //std::string refine="";
      std::string enrich="";
      bool doRefineMesh = true;
      int load_balance = 1;
      std::string convert_Hex8_Tet4_24 = "Hex8_Tet4_24";
      int print_info=0;
      int serialized_io_group_size = 0;
      int remove_original_elements = 1;
      int verify_meshes = 0;
      int number_refines = 1;
      int proc_rank_field = 0;
      int query_only = 0;
      int progress_meter = 0;
      int smooth_geometry = 0;
      int smooth_use_reference_mesh = 1;
      int fix_all_block_boundaries = 0;
      std::string ioss_write_options = "";
      std::string ioss_read_options = "";
      int snap_geometry = 0;
      std::string internal_test = "";
      int respect_spacing = 1;
      int smooth_surfaces = 0;
      //double min_spacing_factor = 0.25; // range [0,0.5]
      int remove_geometry_blocks = 0;
      int dump_geometry_file = 0;
      int sync_io_regions = 0;
      int delete_parents = 1;
      int print_memory_usage = 0;
      // a list of comma-separated names like Entity, Relation, Field, etc.
      std::string memory_multipliers_file="";
      int estimate_memory_usage=0;
      int streaming_size=0;
      int streaming_rank=0;
      int streaming_pass_start= -2;  // FIXME - change to not start from -1 below
      int streaming_pass_end= -2;
      //std::string streaming_instruction="";
      int streaming_W = 0;
      int streaming_iW = 0;
      std::string compute_hmesh = "";
      int print_hmesh_surface_normal = 0;
      int save_internal_fields = 0;

      double hmesh_factor = 0.0;
      double hmesh_min_max_ave_factor[3] = {0,0,0};
      //std::string histograms_root="histograms_root";
      std::string histograms_root="cout";
      //std::string histogram_options = "{mesh: [edge_length, quality_edge, quality_vol_edge_ratio, volume] }";
      std::string histogram_options = "";

      //  Hex8_Tet4_24 (default), Quad4_Quad4_4, Qu
      std::string block_name_desc =
        "block name(s) to convert: there are several options\n"
        "  (1) empty string or option not specified: convert all blocks in the input mesh file\n"
        "  (2) file:my_filename.my_ext (e.g. file:filelist.dat) which will read input block names\n"
        "            from the given file\n"
        "  (3) [+]block_name_1,[+]block_name_2, etc ,block_name_n to include only these blocks, plus sign is optional\n"
        "  (4) a single input block name (e.g. block_3) to be converted \n"
        "  (5) -block_3,-block_5 to exclude blocks from those included (all blocks or include-only blocks), minus sign is mandatory\n"
        "  (6) block_1..block_10 include the range of blocks #1 to #10 \n"
        "  (7) any combination of [+] and - options and range (..) option can be specified \n"
        "Note: wherever you specify block_# this can be replaced with just the #, e.g. \"1,2,4,5\"";

      std::string convert_options = UniformRefinerPatternBase::s_convert_options;
      std::string refine_options  = UniformRefinerPatternBase::s_refine_options;
      std::string enrich_options  = UniformRefinerPatternBase::s_enrich_options;

      std::string def1= Util::split(convert_options, ", ")[0] ;
      if (0) std::cout << def1 << "tmp split = " << Util::split(convert_options, ", ") << std::endl;
      int test_memory_elements = 0;
      int test_memory_nodes = 0;

      //convert_options = "DEFAULT or one of "+convert_options;
      //refine_options = "DEFAULT or one of "+refine_options;
      //enrich_options = "DEFAULT or one of "+enrich_options;

      // : if not specified, use input mesh name appended with _{converted,refined_#refines,enriched}");

      std::string block_name_desc_inc = "which blocks to include, specified as: "+block_name_desc;
      std::string block_name_desc_exc = "which blocks to exclude, specified as: "+block_name_desc;

      int help = 0;

      run_environment.clp.setOption("help"                     , &help                     , "print this usage message");

      // files
      run_environment.clp.setOption("input_mesh"               , &input_mesh               , "input mesh name");
      run_environment.clp.setOption("output_mesh"              , &output_mesh              , "output mesh name");
      run_environment.clp.setOption("load_balance"             , &load_balance             , "load balance (decomp/slice/spread) input mesh file");

      // operation type
      run_environment.clp.setOption("refine"                   , &refine                   , refine_options.c_str());
      run_environment.clp.setOption("number_refines"           , &number_refines           , "number of refinement passes");
      run_environment.clp.setOption("convert"                  , &convert                  , convert_options.c_str());
      run_environment.clp.setOption("enrich"                   , &enrich                   , enrich_options.c_str());

      // spacing
#if !defined(__IBMCPP__)
      run_environment.clp.setOption("respect_spacing"          , &respect_spacing          , "respect the initial mesh spacing during refinement");
#endif

      // query/control
      run_environment.clp.setOption("verify_meshes"            , &verify_meshes            , "verify positive volumes for input and output meshes");
      run_environment.clp.setOption("query_only"               , &query_only               , "query only, no refinement done");
      run_environment.clp.setOption("progress_meter"           , &progress_meter           , "progress meter on or off");
      run_environment.clp.setOption("print_info"               , &print_info               , ">= 0  (higher values print more info)");
      run_environment.clp.setOption("print_memory_usage"       , &print_memory_usage       , "print memory usage");
      run_environment.clp.setOption("estimate_memory_usage"    , &estimate_memory_usage    ,
                                    " use internal or memory_multipliers_file values to estimate memory needed.\n"
                                    "   if query_only=1, use multipliers from memory_multipliers_file to estimate memory to be used in refinements, if memory_multipliers_file is set.\n"
                                    "   if query_only=1, and no memory_multipliers_file is set, use internal values for memory_multipliers.\n"
                                    "   If query_only=0, print actual memory data and estimates.");

      // geometry
      run_environment.clp.setOption("input_geometry"           , &input_geometry           , "input geometry name");
      run_environment.clp.setOption("smooth_geometry"          , &smooth_geometry          , "smooth geometry - moves nodes after geometry projection to try to avoid bad meshes");
      run_environment.clp.setOption("smooth_surfaces"          , &smooth_surfaces          , "allow nodes to move on surfaces when smoothing");
      run_environment.clp.setOption("dump_geometry_file"       , &dump_geometry_file       , "debug print geometry (OpenNURBS 3dm) file contents");
      run_environment.clp.setOption("fit_geometry_file"        , &fit_geometry_file        , "for 2D meshes, create an OpenNURBS 3dm file fitting cubic splines to all boundaries");

      // smoothing
      run_environment.clp.setOption("smooth_use_reference_mesh", &smooth_use_reference_mesh, "for most cases, set to 1 (default) - can be used for smoothing with no reference mesh");
      run_environment.clp.setOption("fix_all_block_boundaries" , &fix_all_block_boundaries , "when smoothing without geometry, fix all inner and outer block boundaries");

      // mesh query
      run_environment.clp.setOption("compute_hmesh"            , &compute_hmesh            , "compute mesh parameter using method eigens|edges");
      run_environment.clp.setOption("print_hmesh_surface_normal"  , &print_hmesh_surface_normal            , "prints a table of normal mesh spacing at each surface");

      // subsetting
      run_environment.clp.setOption("blocks"                   , &block_name_inc           , block_name_desc_inc.c_str());

      // histograms
#if STK_ADAPT_HAVE_YAML_CPP
      run_environment.clp.setOption("histogram"  , &histogram_options  ,
                                    "\n  either a single filename, which reads further commands in YAML format (yaml.org) from that file, or\n"
                                    "  a string of the form \"{ fields: [field_1,...,field_n], file_root: my_histograms, \n"
                                    "     mesh: [edge_length, quality_edge, quality_vol_edge_ratio, volume], time: 0.1, step: 2 }\" \n"
                                    "  where field_i are field names to get stats for, file_root is a root filename, \n"
                                    "  and mesh: gives options for mesh quality histograms.\n"
                                    "  time: or step: options allow specifying which timestep in the database should be used.\n"
                                    "  If read from a file, file format is like this: \n"
                                    "    fields:\n"
                                    "      - pressure\n"
                                    "      - velocity\n"
                                    "      - temperature\n"
                                    "    file_root: my_histograms\n"
                                    "    mesh:\n"
                                    "      - edge_length\n"
                                    "      - quality_edge\n"
                                    "      - quality_vol_edge_ratio\n"
                                    "      - volume");
#endif
      run_environment.clp.setOption("histogram_file_root"      , &histograms_root          , " if cout, use screen, else use this as the root name of histogram files.");

      // ioss options
      run_environment.clp.setOption("ioss_read_options"        , &ioss_read_options        ,
                                    "options to IOSS/Exodus for e.g. large files | auto-decomp | auto-join\n"
                                    "to use, set the string to a combination of \n"
                                    "{\"large\", \"auto-decomp:yes\",  \"auto-decomp:no\", \n"
                                    "   \"auto-join:yes\", \"auto-join:no\" }\n"
                                    "   e.g. \"large,auto-decomp:yes\" \n"
                                    " Note: set options for read and/or write (ioss_write_options)");
      run_environment.clp.setOption("ioss_write_options"       , &ioss_write_options       , "see ioss_read_options");

      // debugging/advanced
      run_environment.clp.setOption("proc_rank_field"          , &proc_rank_field          , " add an element field to show processor rank");
      run_environment.clp.setOption("remove_original_elements" , &remove_original_elements , " remove original (converted) elements (default=true)");
      run_environment.clp.setOption("remove_geometry_blocks"   , &remove_geometry_blocks   , "remove geometry blocks from output Exodus file after refinement/geometry projection");

      // internal
      run_environment.clp.setOption("delete_parents"           , &delete_parents           , "DEBUG: delete parents from a nested, multi-refine mesh - used for debugging");
      run_environment.clp.setOption("snap_geometry"            , &snap_geometry            , "project nodes to geometry - used for internal testing only");
      run_environment.clp.setOption("internal_test"            , &internal_test            , "run the specified internal test");
      run_environment.clp.setOption("sync_io_regions"          , &sync_io_regions          , "synchronize input/output region's Exodus id's");
      run_environment.clp.setOption("save_internal_fields"     , &save_internal_fields     , "save internally created fields to the output file");

#if defined(STK_BUILT_IN_SIERRA)
      // Salinas
      run_environment.clp.setOption("rbar_blocks"              , &rbar_blocks              , "list of blocks to treat in special Salinas fashion for RBARs - see block_name description for format.");
#endif
      //run_environment.clp.setOption("exclude"                  , &block_name_exc           , block_name_desc_exc.c_str());

      run_environment.clp.setOption("memory_multipliers_file"  , &memory_multipliers_file  ,
                                    "[experimental]  filename with 3 space-separated entries, with estimate for bytes-per-hex8 tet4 and nodes, e.g. 300 280 200\n"
                                    "  If not set, use internal estimates for memory multipliers.");

#if ALLOW_MEM_TEST
      run_environment.clp.setOption("test_memory_elements"     , &test_memory_elements     , " give a number of elements");
      run_environment.clp.setOption("test_memory_nodes"        , &test_memory_nodes        , " give a number of nodes");
#endif
      run_environment.clp.setOption("serialized_io_group_size" , &serialized_io_group_size , "[experimental] set to non-zero to use this many i/o groups to minimize disk contention");

      run_environment.clp.setOption("streaming_size"           , &streaming_size      ,
                                    "INTERNAL use only by python script streaming refinement interface:\n"
                                    "  run in streaming mode - this number specifies how many virtual procs the mesh is split into\n"
                                    "    i.e. we expect to see files like file.e.N.iN where N = streaming_size iN=0..N");
      run_environment.clp.setOption("streaming_rank"           , &streaming_rank     ,
                                    "INTERNAL use only by python script streaming refinement interface:\n"
                                    "  run in streaming mode - this number specifies which virtual proc this is.");
      run_environment.clp.setOption("streaming_pass_start"     , &streaming_pass_start           ,
                                    "INTERNAL use only by python script streaming refinement interface:");
      run_environment.clp.setOption("streaming_pass_end"       , &streaming_pass_end           ,
                                    "INTERNAL use only by python script streaming refinement interface:");
      run_environment.clp.setOption("streaming_W"              , &streaming_W           ,
                                    "INTERNAL use only by python script streaming refinement interface:");
      run_environment.clp.setOption("streaming_iW"             , &streaming_iW          ,
                                    "INTERNAL use only by python script streaming refinement interface:");

      // Detailed doc
      std::string docString = s_docString;

#ifdef STK_BUILT_IN_SIERRA
      std::string docStringSalinas =
        "Salinas Special Command --rbar_blocks\n"
        "\n"
        "Percept can refine a mesh that contains Salinas RBAR elements (beams connecting nodes of one surface\n"
        "to another, e.g. to model a joint or spring).  The new mesh will contain RBARs connecting new nodes\n"
        "on one surface to new nodes on the other.  Specify a list of block names that contain the RBAR\n"
        "elements (see the --blocks command for format).\n";

      docString = docString + docStringSalinas;
#endif

      if (!found_Help) docString = "";
      run_environment.clp.setDocString(docString.c_str());

      int err_clp = run_environment.processCommandLine();
      if (err_clp) return err_clp;

      if (found_Help) {
        run_environment.printHelp();
#if defined( STK_HAS_MPI )
          MPI_Barrier( run_environment.m_comm );
          MPI_Finalize();
#endif
          std::exit(0);
      }

      std::string histogram_basic_options = "{file_root: "+histograms_root + ", mesh: [edge_length, quality_edge, quality_vol_edge_ratio, volume] }";
      int result = 0;
      unsigned failed_proc_rank = 0u;

      double t0   = 0.0;
      double t1   = 0.0;
      double cpu0 = 0.0;
      double cpu1 = 0.0;

      if (serialized_io_group_size)
      {
        std::cout << "Info: found non-zero serialized_io_group_size on command-line= "
                              << serialized_io_group_size << std::endl;
        if (serialized_io_group_size < 0 || serialized_io_group_size > (int)p_size || (int)p_size % serialized_io_group_size != 0)
          {
            if (p_rank==0)
              std::cout << "Error: Job requested serialized_io_group_size of " << serialized_io_group_size
                   << " which is incompatible with MPI size= " << p_size
                   << "... shutting down." << std::endl;
            throw std::runtime_error("bad value for serialized_io_group_size");
          }
        Ioss::SerializeIO::setGroupFactor(serialized_io_group_size);
      }

      if (convert.length()+enrich.length()+refine.length() == 0)
        {
          std::cout << "\nCommand line syntax error: you must give a value for one (and only one) of the options: refine, enrich, or convert.\n" << std::endl;
          run_environment.printHelp();
          std::cout << "\nCommand line syntax error: you must give a value for one (and only one) of the options: refine, enrich, or convert.\n" << std::endl;
          exit(1);
        }

      if (convert.length())
        checkInput("convert", convert, convert_options, run_environment);

      if (enrich.length())
        checkInput("enrich", enrich, enrich_options, run_environment);

      if (refine.length())
        checkInput("refine", refine, refine_options, run_environment);

      if (print_info)
        {
          doRefineMesh = false;
        }

      if (help
          || input_mesh.length() == 0
          || output_mesh.length() == 0
          || ((convert.length() == 0 && refine.length()==0 && enrich.length()==0) && number_refines)
          //||  not (convert == "Hex8_Tet4_24" || convert == "Quad4_Quad4_4" || convert == "Quad4_Tri3_6")
          )
        {
          run_environment.printHelp();
          return 1;
        }

#if defined( STK_HAS_MPI )
      MPI_Barrier( MPI_COMM_WORLD );
#endif

#if defined( STK_PERCEPT_HAS_GEOMETRY )
      if (dump_geometry_file)
        {
          GeometryKernelOpenNURBS gko;
          gko.debug_dump_file(input_geometry);
        }
#endif

      // ========================================================================================================================================================================================================
      //    START
      // ========================================================================================================================================================================================================
      std::string input_mesh_save = input_mesh;
      std::string output_mesh_save = output_mesh;

      // streaming
      int m_M = 1;
      int m_W = 1;
      int m_iW = 0;
      int m_M_0 = 0;
      int m_M_1 = 0;
#if STK_ADAPT_HAVE_YAML_CPP
      if (streaming_size)
        {
          m_M = streaming_size;
          m_W = streaming_W ? streaming_W : 1;
          m_iW = streaming_W ? streaming_iW : 0;
          SerializeNodeRegistry::getStreamingPiece(m_M, m_W, m_iW, m_M_0, m_M_1);
        }
      //std::cout << "tmp srk AdaptMain: " << PERCEPT_OUT(streaming_size) << PERCEPT_OUT(streaming_W) << PERCEPT_OUT(streaming_iW) << PERCEPT_OUT(m_M) << PERCEPT_OUT(m_W) << PERCEPT_OUT(m_iW) << PERCEPT_OUT(m_M_0) << PERCEPT_OUT(m_M_1) << std::endl;
#endif

      // FIXME - starting from -1 pass is bogus
#if STK_ADAPT_HAVE_YAML_CPP
      int remove_original_elements_save = remove_original_elements;
      int delete_parents_save = delete_parents;

      if ((streaming_pass_start == -2 && streaming_pass_end != -2) ||
          (streaming_pass_start != -2 && streaming_pass_end == -2))
        {
          throw std::runtime_error("must specify both streaming_pass_start and streaming_pass_end");
        }

      if (streaming_pass_start != -2 && (streaming_pass_start < -1  || streaming_pass_start > SerializeNodeRegistry::MaxPass))
        {
          throw std::runtime_error("streaming_pass_start bad value");
        }
      if (streaming_pass_end != -2 && (streaming_pass_end < -1  || streaming_pass_end > SerializeNodeRegistry::MaxPass))
        {
          throw std::runtime_error("streaming_pass_end bad value");
        }
      if (streaming_pass_end != -2 && (streaming_pass_end < streaming_pass_start))
        {
          throw std::runtime_error("streaming_pass_start > streaming_pass_end");
        }

      if (streaming_pass_end != -2)
        {
          std::cout << "\n\nWARNING: running passes from command line: streaming_pass_start,end= [" << streaming_pass_start << ", " << streaming_pass_end << "]\n\n";
        }

      // allow for driving this from a script
      if (streaming_pass_start == -2 && streaming_pass_end == -2)
        {
          streaming_pass_start = streaming_size ? -1 : 0;
          streaming_pass_end = streaming_size ? SerializeNodeRegistry::MaxPass : 0;
        }
#else
      streaming_pass_start = 0;
      streaming_pass_end = 0;
#endif

      bool do_normal_pass = true;
      if (m_W > 1 && (m_iW == -1 || m_iW == m_W))
        {
          do_normal_pass = false;
        }

      if (streaming_size && s_spatialDim == 0)
        {
          PerceptMesh eMesh(0);
          std::string mesh_name = Ioss::Utils::decode_filename(input_mesh_save, 0, m_M);
          eMesh.open(mesh_name);
          if (smooth_geometry == 1) eMesh.add_coordinate_state_fields();
#if !defined(__IBMCPP__)
          if (respect_spacing == 1) {
            eMesh.set_respect_spacing(true);
            eMesh.add_spacing_fields();
          }
#endif
          if (smooth_surfaces == 1) eMesh.set_smooth_surfaces(true);
          s_spatialDim = eMesh.get_spatial_dim();
          VERIFY_OP_ON(s_spatialDim, >=, 2, "AdaptMain bad spatial_dim");
        }

      for (int i_pass=streaming_pass_start; i_pass <= streaming_pass_end; i_pass++)
        {

#if STK_ADAPT_HAVE_YAML_CPP
          if (streaming_size)
            {

              // init pass
              if (m_W > 1 && m_iW == -1)
                {
                  PerceptMesh eMesh(3);
                  eMesh.openEmpty();
                  SerializeNodeRegistry snr(eMesh, 0, input_mesh, output_mesh, m_M, m_M_0,  m_W, m_iW, m_M_0, m_M_1);
                  snr.pass_init(i_pass);
                  VERIFY_OP_ON(streaming_pass_start, ==, streaming_pass_end, "must have one pass at a time with W > 1");
                  continue;
                }

              remove_original_elements = 0;
              delete_parents = 0;
              if (i_pass == 1)
                {
                  remove_original_elements = remove_original_elements_save;
                  delete_parents = delete_parents_save;
                }
            }
          //std::cout << "tmp srk i_pass= " << i_pass << " delete_parents= " << delete_parents << " remove_original_elements= " << remove_original_elements << std::endl;
#endif

          if (do_normal_pass)
            {
              if (streaming_size)
                {
                  // special cases

                  // read all files, get global parts and node info
                  if (i_pass == -1)
                    {
#if STK_ADAPT_HAVE_YAML_CPP
                      PerceptMesh eMesh(s_spatialDim);
                      eMesh.openEmpty();
                      SerializeNodeRegistry snr(eMesh, 0, input_mesh, output_mesh, m_M, 0, m_W, m_iW, m_M_0, m_M_1);
                      snr.pass(i_pass);
#else
                      throw std::runtime_error("must have YAML for streaming refine");
#endif
                      continue;
                    }

                  if (i_pass == 2)
                    {
#if STK_ADAPT_HAVE_YAML_CPP
                      if (m_W == 1)
                        {
                          //  no exodus files i/o
                          PerceptMesh eMesh(s_spatialDim);
                          eMesh.openEmpty();
                          SerializeNodeRegistry snr(eMesh, 0, input_mesh, output_mesh, m_M, 0, m_W, m_iW, m_M_0, m_M_1);
                          snr.pass(i_pass);
                        }

#else
                      throw std::runtime_error("must have YAML for streaming refine");
#endif
                      continue;
                    }
                }

              for (int m_iM = m_M_0; m_iM <= m_M_1; m_iM++)
                {
                  if (streaming_size)
                    {
                      input_mesh = Ioss::Utils::decode_filename(input_mesh_save, m_iM, m_M);
                      output_mesh = Ioss::Utils::decode_filename(output_mesh_save, m_iM, m_M);

                      if (i_pass == 1)
                        {
                          output_mesh = output_mesh_save+"-pass1";
                          output_mesh = Ioss::Utils::decode_filename(output_mesh, m_iM, m_M);
                        }
                    }

                  try {

                    if (load_balance)
                      {
                        if (streaming_size) throw std::runtime_error("can't load balance and stream");
                        RunEnvironment::doLoadBalance(run_environment.m_comm, input_mesh);
                      }

                    percept::PerceptMesh eMesh(0);  // FIXME
                    //percept::PerceptMesh eMesh;  // FIXME

                    if (streaming_size)
                      {
                        // special cases
                        if (i_pass == 3)
                          {
#if STK_ADAPT_HAVE_YAML_CPP
                            // fini
                            //if ((m_W > 1 && m_iW == m_W) || (m_W == 1 && m_iM == m_M_1))
                            if (m_iM == m_M_1)
                              {
                                NodeRegistry *some_nr = 0;
                                PerceptMesh eMeshEmpty(0);
                                eMeshEmpty.openEmpty();
                                SerializeNodeRegistry snr(eMeshEmpty, some_nr, input_mesh_save, output_mesh_save, m_M, 0, m_W, m_iW, m_M_0, m_M_1);
                                if (remove_geometry_blocks) snr.set_geometry_file(input_geometry);
                                snr.pass3_new();
                              }
#else
                            throw std::runtime_error("must have YAML for streaming refine");
#endif
                            continue;
                          }
                      }

                    // ==============  START Normal Pass  ==================================================================

                    if (do_normal_pass)
                      {
                        if (ioss_read_options.length() || ioss_write_options.length())
                          {
                            if (!eMesh.get_rank())
                              {
                                std::cout << "INFO: ioss_read_options= " << ioss_read_options << " ioss_write_options= " << ioss_write_options << std::endl;
                                stk::percept::pout() << "INFO: ioss_read_options= " << ioss_read_options << " ioss_write_options= " << ioss_write_options << std::endl;
                              }
                          }

                        if (ioss_read_options.length())  eMesh.set_ioss_read_options(ioss_read_options);
                        if (ioss_write_options.length()) eMesh.set_ioss_write_options(ioss_write_options);

                        eMesh.open(input_mesh);
                        eMesh.set_save_internal_fields(save_internal_fields);
                        if (smooth_geometry == 1) eMesh.add_coordinate_state_fields();
#if !defined(__IBMCPP__)
                        if (respect_spacing == 1) {
                          eMesh.set_respect_spacing(true);
                          eMesh.add_spacing_fields();
                        }
#endif
                        if (smooth_surfaces == 1) eMesh.set_smooth_surfaces(true);
                        if (!sync_io_regions) eMesh.set_sync_io_regions(false);
                        if (!s_spatialDim) s_spatialDim = eMesh.get_spatial_dim();

                        Util::setRank(eMesh.get_rank());

                        Teuchos::RCP<UniformRefinerPatternBase> pattern;

                        if (doRefineMesh)
                          {
                            // FIXME move this next block of code to a method on UniformRefiner
                            BlockNamesType block_names(stk::percept::EntityRankEnd+1u);

#if defined(STK_BUILT_IN_SIERRA)
                            if (rbar_blocks.length())
                              {
                                BlockNamesType rbar_names(stk::percept::EntityRankEnd+1u);
                                std::string block_name_inc_orig = block_name_inc;
                                if (rbar_blocks.length())
                                  {
                                    rbar_names = RefinerUtil::getBlockNames(rbar_blocks, eMesh.get_rank(), eMesh);
                                    std::cout << "rbar_blocks= " << rbar_blocks << " rbar_names= " << rbar_names << std::endl;
                                  }
                                for (unsigned ii=0; ii < rbar_names[eMesh.element_rank()].size(); ii++)
                                  {
                                    std::string srb = rbar_names[eMesh.element_rank()][ii];
                                    Util::replace(srb, "+", "-");
                                    block_name_inc = block_name_inc+(block_name_inc.length()?",":"")+srb;
                                  }
                                if (!eMesh.get_rank())
                                  std::cout << "rbar: original block_name option = " << block_name_inc_orig << " new = " << block_name_inc << std::endl;
                              }
#endif

                            if (block_name_inc.length())
                              {
                                block_names = RefinerUtil::getBlockNames(block_name_inc, eMesh.get_rank(), eMesh);
                                if (1)
                                  {
                                    eMesh.commit();
                                    block_names = RefinerUtil::correctBlockNamesForPartPartConsistency(eMesh, block_names);

                                    eMesh.close();
                                    eMesh.open(input_mesh);
                                    if (smooth_geometry == 1) eMesh.add_coordinate_state_fields();
#if !defined(__IBMCPP__)
                                    if (respect_spacing == 1) {
                                      eMesh.set_respect_spacing(true);
                                      eMesh.add_spacing_fields();
                                    }
#endif
                                    if (smooth_surfaces == 1) eMesh.set_smooth_surfaces(true);

                                  }
                                if (!eMesh.get_rank()) std::cout << "block_names after processing: " << block_names << std::endl;
                              }

                            pattern = UniformRefinerPatternBase::createPattern(refine, enrich, convert, eMesh, block_names);

                            if (0)
                              {
                                run_environment.printHelp();
                                exit(1);
                              }
                          }

                        if (0)
                          {
                            std::cout << "tmp convert = " << convert << std::endl;
                            std::cout << "tmp refine = " << refine << std::endl;
                            std::cout << "tmp enrich = " << enrich << std::endl;
                          }

                        int scalarDimension = 0; // a scalar

                        stk::mesh::FieldBase* proc_rank_field_ptr = 0;
                        if (proc_rank_field)
                          {
                            proc_rank_field_ptr = eMesh.add_field("proc_rank", stk::mesh::MetaData::ELEMENT_RANK, scalarDimension);
                          }

                        if (fix_all_block_boundaries)
                          {
                            bool make_part_io_part=true;
                            eMesh.add_part("inner_skin_part", make_part_io_part);
                          }

#if STK_ADAPT_HAVE_YAML_CPP
                        // FIXME - this is this needed? see above
                        // add global parts not found in this file
                        if (streaming_size)
                          {
                            NodeRegistry *some_nr = 0;
                            SerializeNodeRegistry snr(eMesh, some_nr, input_mesh, output_mesh, m_M, m_iM,  m_W, m_iW, m_M_0, m_M_1);
                            snr.declareGlobalParts();
                          }
#endif
                        stk::mesh::PartVector fg_surface_parts;
                        stk::mesh::PartVector fg_topo_parts;
                        if (fit_geometry_file != "")
                          {
                            fit_geometry_create_parts_meta(eMesh, fg_surface_parts, fg_topo_parts);
                          }

                        eMesh.commit();

                        //eMesh.delete_side_sets();

                        if (fix_all_block_boundaries)
                          {
                            eMesh.get_skin_part("inner_skin_part", true);
                          }

                        if (print_info)
                          {
                            eMesh.print_info("PerceptMesh info:", print_info);
                          }

                        if (fit_geometry_file != "")
                          {
                            fit_geometry(eMesh, fit_geometry_file, fg_surface_parts, fg_topo_parts);
                          }

                        // print message about rbars being treated and beams being refined
//                         if (!eMesh.get_rank())
//                           std::cout << "P[" << eMesh.get_rank() << "] Adding rbar elements as requested by user for block[" << ipart << "]= " << part.name()
//                                     << "\n  NOTE:  This block is automatically ignored during refinement."
//                                     << std::endl;

                        if (verify_meshes)
                          {
                            bool print_table=true;
                            double badJac=1.e-10;
                            int dump_all_elements = verify_meshes - 1;
                            if (!eMesh.get_rank()) std::cout << "Verify input mesh..." << std::endl;
                            if (eMesh.check_mesh_volumes(print_table, badJac, dump_all_elements))
                              {
                                throw std::runtime_error("ERROR: verify_meshes shows a bad input mesh");
                              }
                          }

#if STK_ADAPT_HAVE_YAML_CPP
                        if (histogram_options.size() != 0)
                          {
                            Histograms<double> histograms;
                            HistogramsParser<double> hparser(histogram_options);
                            hparser.create(histograms);

                            double hopt_time = histograms.m_database_time;
                            int hopt_step = histograms.m_database_step;
                            int current_step = eMesh.get_current_database_step();
                            if (hopt_time >= 0.0)
                              {
                                eMesh.read_database_at_time(hopt_time);
                              }
                            else if (hopt_step >= 0)
                              {
                                eMesh.read_database_at_step(hopt_step);
                              }
                            else
                              {
                                // get last step
                                int step = eMesh.get_database_time_step_count();
                                eMesh.read_database_at_step(step);
                                //std::cout << "step= " << step << " current_step= " << current_step << std::endl;
                                //eMesh.read_database_at_step(step?step:1);
                              }

                            eMesh.mesh_field_stats(&histograms);
                            histograms.compute_uniform_bins(10);

                            if (!p_rank) {
                              std::cout << "Before refine, user-requested histograms= " << std::endl;
                              histograms.print(true);
                              //histograms.print(stk::percept::pout());
                            }
                            // reset
                            eMesh.read_database_at_step(current_step);
                          }
#endif

                        if (compute_hmesh.size() != 0)
                          {
                            double hmesh=0.0;
                            Histograms<double> histograms(histograms_root);
                            //Histograms<double> histograms;
#if STK_ADAPT_HAVE_YAML_CPP
                            HistogramsParser<double> hparser(histogram_basic_options);
                            hparser.create(histograms);
#endif
                            if (compute_hmesh == "eigens")
                              {
                                hmesh = eMesh.hmesh_stretch_eigens(hmesh_min_max_ave_factor, &histograms["stretch_eigens"], &histograms["quality_edge_eigens"]);
                                histograms["stretch_eigens"].set_titles("Stretch Eigens Histogram");
                                histograms["quality_edge_eigens"].set_titles("Stretch Eigens Max/Min Quality Histogram");
                              }
                            else if (compute_hmesh == "edges")
                              {
                                hmesh = eMesh.hmesh_edge_lengths(hmesh_min_max_ave_factor, &histograms["edge_length"], &histograms["quality_edge"]);
                                histograms["edge_length"].set_titles("Edge Length Histogram");
                                histograms["quality_edge_eigens"].set_titles("Edge Max/Min Quality Histogram");
                              }
                            else
                              {
                                throw std::runtime_error("unknown option for compute_hmesh: "+compute_hmesh);
                              }
                            histograms.compute_uniform_bins(10);

                            if (!p_rank) {
                              std::cout << "Before refine, Mesh size (h-parameter) = " << hmesh
                                        << " (min = " << hmesh_min_max_ave_factor[0]
                                        << " max = " << hmesh_min_max_ave_factor[1]
                                        << " ave = " << hmesh_min_max_ave_factor[2]
                                        << ") "
                                        << std::endl;
                              stk::percept::pout() << "Before refine, Mesh size (h-parameter) = " << hmesh
                                                   << " (min = " << hmesh_min_max_ave_factor[0]
                                                   << " max = " << hmesh_min_max_ave_factor[1]
                                                   << " ave = " << hmesh_min_max_ave_factor[2]
                                                   << ")\n " ;
                              histograms.print(true);
                              //histograms.print(stk::percept::pout());
                            }
                            hmesh_factor = hmesh;
                          }
                        if (print_hmesh_surface_normal)
                          {
                            std::string msg="before refine";
                            eMesh.print_hmesh_surface_normal(msg, std::cout);
                            eMesh.print_hmesh_surface_normal(msg, stk::percept::pout());
                          }
#if !defined(__IBMCPP__)
                        if (respect_spacing)
                          {
                            SpacingFieldUtil sfu(eMesh);
                            sfu.compute_spacing_field();
                          }
#endif

                        if (print_memory_usage)
                          memory_dump(print_memory_usage, run_environment.m_comm, *eMesh.get_bulk_data(), 0, "after file open");

                        if (test_memory_nodes && test_memory_elements)
                          {
                            std::cout << "test_memory_elements and nodes are nonzero, will not refine exodus files." << std::endl;

                            test_memory(eMesh, test_memory_elements, test_memory_nodes);

                            if (print_memory_usage)
                              memory_dump(print_memory_usage, run_environment.m_comm, *eMesh.get_bulk_data(), 0, "after test memory");

                            if (estimate_memory_usage && !query_only)
                              {
                                MemorySizeType tot_mem = memory_dump(false, run_environment.m_comm, *eMesh.get_bulk_data(), 0, "after test memory");

                                //std::cout << "MemEst: num_nodes= " << test_memory_nodes << " num_tet4=0 hum_hex8= " << test_memory_elements << " memory= " << MegaByte(tot_mem) << std::endl;
                                //MemoryMultipliers::process_estimate(tot_mem, eMesh, breaker.getRefinementInfoByType(), memory_multipliers_file);
                                MemoryMultipliers memMults;
                                if (memory_multipliers_file.size())
                                  memMults.read_simple(memory_multipliers_file);
                                memMults.num_hex8=test_memory_elements;
                                memMults.num_nodes=test_memory_nodes;
                                MemorySizeType estMem = memMults.estimate_memory();
                                //                 std::cout << "MemEst: num_nodes= " << memMults.num_nodes << " num_tet4= " << memMults.num_tet4 << " num_hex8= " << memMults.num_hex8 << " memory= " << MegaByte(tot_mem)
                                //                           << " estMem= " << MegaByte(estMem) << std::endl;
                                std::cout << "MemEst: num_nodes= " << memMults.num_nodes << " num_tet4= " << memMults.num_tet4 << " num_hex8= " << memMults.num_hex8 << " memory[MB]= " << MegaByte(tot_mem)
                                          << " estMem[MB]= " << MegaByte(estMem)
                                          << " mult_hex8= " << memMults.mult_hex8 << " mult_tet4= " << memMults.mult_tet4 << " mult_nodes=" << memMults.mult_nodes << std::endl;
                                std::cout << "(*MemEstMM: " << input_mesh << " *) ,{" << memMults.num_nodes << ", " << memMults.num_tet4 << "," << memMults.num_hex8 << "," << MegaByte(tot_mem)
                                          << ", " << MegaByte(estMem) << "}" << std::endl;

                              }
                            if (estimate_memory_usage && query_only)
                              {
                                //MemoryMultipliers::process_estimate(0, eMesh, breaker.getRefinementInfoByType(), memory_multipliers_file, input_mesh);
                              }

                            return 0;
                          }

                        // FIXME
                        if (0)
                          {
                            eMesh.save_as("outtmp.e");
                            exit(1);
                          }

                        // FIXME
                        if (0)
                          {
                            eMesh.print_info("before convert", 2);
                            exit(1);
                          }

                        if (doRefineMesh)
                          {
                            t0 =  stk::wall_time();
                            cpu0 = stk::cpu_time();

                            UniformRefiner breaker(eMesh, *pattern, proc_rank_field_ptr);

                            ProgressMeter pm(breaker);
                            //pm.setActive(true);

                            //std::cout << "P[" << p_rank << ", " << p_size << "] input_geometry = " << input_geometry << std::endl;

                            if (input_geometry != "")
                              {
                                breaker.setGeometryFile(input_geometry);
                                breaker.setSmoothGeometry(smooth_geometry == 1);
                                //breaker.setRemoveGeometryBlocks(remove_geometry_blocks == 1);
                              }
                            breaker.setRemoveOldElements(remove_original_elements);
                            breaker.setFixAllBlockBoundaries(fix_all_block_boundaries);
                            breaker.setQueryPassOnly(query_only == 1);
                            breaker.setDoProgressMeter(progress_meter == 1 && 0 == p_rank);
                            //breaker.setIgnoreSideSets(true);
#if defined(STK_BUILT_IN_SIERRA)
                            if (rbar_blocks.length())
                              {
                                BlockNamesType rbar_names(stk::percept::EntityRankEnd+1u);
                                if (rbar_blocks.length())
                                  {
                                    rbar_names = RefinerUtil::getBlockNames(rbar_blocks, eMesh.get_rank(), eMesh);
                                    std::cout << "rbar_blocks= " << rbar_blocks << " rbar_names= " << rbar_names << std::endl;
                                  }
                                breaker.set_rbar_special_treatment(rbar_names);
                              }
#endif

                            for (int iBreak = 0; iBreak < number_refines; iBreak++)
                              {
                                if (!eMesh.get_rank())
                                  {
                                    std::cout << "Refinement pass # " << (iBreak+1) << " start..." << std::endl;
                                  }
                                //breaker.setPassNumber(iBreak);
                                breaker.doBreak();
                                //RefinementInfoByType::countCurrentNodes(eMesh, breaker.getRefinementInfoByType());
                                if (!eMesh.get_rank())
                                  {
                                    std::cout << std::endl;
                                    int ib = iBreak;
                                    if (!query_only) ib = 0;
                                    bool printAllTopologies = false;
                                    RefinementInfoByType::printTable(std::cout, breaker.getRefinementInfoByType(), ib , printAllTopologies);
                                    RefinementInfoByType::printTable(stk::percept::pout(), breaker.getRefinementInfoByType(), ib , printAllTopologies);
                                    std::cout << std::endl;
                                  }
                                if (print_memory_usage)
                                  {
                                    memory_dump(print_memory_usage, run_environment.m_comm, *eMesh.get_bulk_data(), &breaker.getNodeRegistry(),
                                                std::string("after refine pass: ")+toString(iBreak));
                                  }

                                if (estimate_memory_usage)
                                  {
                                    if (query_only)
                                      {
                                        if (number_refines == 1)
                                          {
                                            MemorySizeType tot_mem = memory_dump(false, run_environment.m_comm, *eMesh.get_bulk_data(), &breaker.getNodeRegistry(),
                                                                                 std::string("after mesh read"));
                                            if (!p_rank)
                                              std::cout << "P[" << p_rank << "] tmp srk mem after mesh read= " << MegaByte(tot_mem) << std::endl;
                                            bool use_new = false;
                                            MemoryMultipliers::process_estimate(tot_mem, eMesh, breaker.getRefinementInfoByType(), memory_multipliers_file, input_mesh, use_new);
                                          }

                                        RefinementInfoByType::estimateNew(breaker.getRefinementInfoByType(), iBreak);
                                        MemoryMultipliers::process_estimate(0, eMesh, breaker.getRefinementInfoByType(), memory_multipliers_file, input_mesh);
                                      }
                                    else
                                      {
                                        MemorySizeType tot_mem = memory_dump(false, run_environment.m_comm, *eMesh.get_bulk_data(), &breaker.getNodeRegistry(),
                                                                             std::string("after refine pass: ")+toString(iBreak));
                                        if (!p_rank)
                                          std::cout << "P[" << p_rank << "] tmp srk tot_mem= " << MegaByte(tot_mem) << std::endl;
                                        MemoryMultipliers::process_estimate(tot_mem, eMesh, breaker.getRefinementInfoByType(), memory_multipliers_file, input_mesh);
                                      }
                                  }

                              } // iBreak

                            if (number_refines == 0 && (smooth_geometry == 1 || snap_geometry == 1))
                              {
#if defined(STK_PERCEPT_HAS_GEOMETRY)
                                breaker.setSmoothGeometry((smooth_geometry == 1));
                                breaker.snapAndSmooth((snap_geometry == 1), input_geometry, (smooth_use_reference_mesh == 1) );
#elif defined(__IBMCPP__)
      				std::ostringstream oss;
      				oss << "\nERROR: Geometry and/or smoothing is not currently supported on this platform. Try running with geometry turned off.";
      				throw std::runtime_error(oss.str());
#endif
                              }

                            if (streaming_size)
                              {
#if STK_ADAPT_HAVE_YAML_CPP
                                {
                                  SerializeNodeRegistry snr(eMesh, &breaker.getNodeRegistry(), input_mesh, output_mesh, m_M, m_iM,  m_W, m_iW, m_M_0, m_M_1);
                                  snr.pass(i_pass);
                                }

                                if (i_pass == 0)
                                  {
                                    NodeRegistry *some_nr = 0;
                                    SerializeNodeRegistry snr(eMesh, some_nr, input_mesh, output_mesh, m_M, m_iM,  m_W, m_iW, m_M_0, m_M_1);
                                    snr.getGlobalPartMap();
                                    snr.passM1_mergeGlobalParts();
                                    snr.setGlobalPartMap();
                                  }
#else
                                throw std::runtime_error("must have YAML for streaming refine");
#endif
                              }

                            if (delete_parents)
                              breaker.deleteParentElements();

                            t1 =  stk::wall_time();
                            cpu1 = stk::cpu_time();

                            if (verify_meshes)
                              {
                                bool print_table=true;
                                double badJac=1.e-10;
                                int dump_all_elements = verify_meshes - 1;
                                if (!eMesh.get_rank()) std::cout << "Verify output mesh..." << std::endl;
                                if (eMesh.check_mesh_volumes(print_table, badJac, dump_all_elements))
                                  {
                                    throw std::runtime_error("ERROR: verify_meshes shows a bad output mesh");
                                  }
                              }
#if STK_ADAPT_HAVE_YAML_CPP
                            if (histogram_options.size() != 0)
                              {
                                Histograms<double> histograms;
                                HistogramsParser<double> hparser(histogram_options);
                                hparser.create(histograms);

                                eMesh.mesh_field_stats(&histograms);
                                histograms.compute_uniform_bins(10);

                                if (!p_rank) {
                                  std::cout << "After refine, user-requested histograms= " << std::endl;
                                  histograms.print(true);
                                  //histograms.print(stk::percept::pout());
                                }
                              }
#endif

                            if (compute_hmesh.size() != 0)
                              {
                                double hmesh=0.0;
                                double min_max_ave[3];
                                Histograms<double> histograms(histograms_root);
                                //Histograms<double> histograms;
#if STK_ADAPT_HAVE_YAML_CPP
                                HistogramsParser<double> hparser(histogram_basic_options);
                                hparser.create(histograms);
#endif

                                if (compute_hmesh == "eigens")
                                  {
                                    hmesh = eMesh.hmesh_stretch_eigens(min_max_ave, &histograms["stretch_eigens"], &histograms["quality_edge_eigens"]);
                                    histograms["stretch_eigens"].set_titles("Stretch Eigens Histogram");
                                    histograms["quality_edge_eigens"].set_titles("Stretch Eigens Max/Min Quality Histogram");
                                  }
                                else if (compute_hmesh == "edges")
                                  {
                                    hmesh = eMesh.hmesh_edge_lengths(min_max_ave, &histograms["edge_length"], &histograms["quality_edge"]);
                                    histograms["edge_length"].set_titles("Edge Length Histogram");
                                    histograms["quality_edge_eigens"].set_titles("Edge Max/Min Quality Histogram");
                                  }
                                else
                                  {
                                    throw std::runtime_error("unknown option for compute_hmesh: "+compute_hmesh);
                                  }
                                histograms.compute_uniform_bins(10);

                                hmesh_factor /= hmesh;
                                hmesh_min_max_ave_factor[0] /= min_max_ave[0];
                                hmesh_min_max_ave_factor[1] /= min_max_ave[1];
                                hmesh_min_max_ave_factor[2] /= min_max_ave[2];
                                if (!p_rank) {
                                  std::cout << "After refine, Mesh size (h-parameter) = " << hmesh << " oldH/newH factor= " << hmesh_factor
                                            << "\n (new min = " << min_max_ave[0]
                                            << " max = " << min_max_ave[1]
                                            << " ave = " << min_max_ave[2]
                                            << ") "
                                            << "\n (old/new min = " << hmesh_min_max_ave_factor[0]
                                            << " max = " << hmesh_min_max_ave_factor[1]
                                            << " ave = " << hmesh_min_max_ave_factor[2]
                                            << ") "
                                            << std::endl;
                                  stk::percept::pout() << "After refine, Mesh size (h-parameter) = " << hmesh << " oldH/newH factor= " << hmesh_factor
                                            << "\n (new min = " << min_max_ave[0]
                                            << " max = " << min_max_ave[1]
                                            << " ave = " << min_max_ave[2]
                                            << ") "
                                            << "\n (old/new min = " << hmesh_min_max_ave_factor[0]
                                            << " max = " << hmesh_min_max_ave_factor[1]
                                            << " ave = " << hmesh_min_max_ave_factor[2]
                                                       << ")\n";
                                  histograms.print(true);
                                  //histograms.print(stk::percept::pout());
                                }
                              }
                            if (print_hmesh_surface_normal)
                              {
                                std::string msg="after refine";
                                eMesh.print_hmesh_surface_normal(msg, std::cout);
                                eMesh.print_hmesh_surface_normal(msg, stk::percept::pout());
                              }

                            if (DEBUG_ADAPT_MAIN || 0 == p_rank) {
                              stk::percept::pout() << "P[" << p_rank << "] AdaptMain::  saving mesh... \n";
                              std::cout << "P[" << p_rank << "]  AdaptMain:: saving mesh... " << std::endl;
                            }
                            if (streaming_size) eMesh.setStreamingSize(m_M);
                            if (remove_geometry_blocks) eMesh.remove_geometry_blocks_on_output(input_geometry);
                            if (0) eMesh.dump_vtk(output_mesh+".vtk",false);
                            eMesh.save_as(output_mesh);
                            if (DEBUG_ADAPT_MAIN || 0 == p_rank) {
                              stk::percept::pout() << "P[" << p_rank << "] AdaptMain:: ... mesh saved\n";
                              std::cout << "P[" << p_rank << "]  AdaptMain:: mesh saved" << std::endl;
                            }

                            if (print_memory_usage)
                              memory_dump(print_memory_usage, run_environment.m_comm, *eMesh.get_bulk_data(), &breaker.getNodeRegistry(), "after final save mesh");

                          } // doRefineMesh
                      } // do_normal_pass

#if 0
                    for (int itime=0; itime < 10; itime++)
                      {
                        std::cout << "tmp timer[" << itime << "]= " << s_timers[itime] << " " << s_timers[itime]/s_timers[3]*100 << " %" << std::endl;
                      }
#endif

                  }
                  catch ( const std::exception * X ) {
                    std::cout << "AdaptMain::  unexpected exception POINTER: " << X->what() << std::endl;
                    failed_proc_rank = p_rank+1u;
                  }
                  catch ( const std::exception & X ) {
                    std::cout << "AdaptMain:: unexpected exception: " << X.what() << std::endl;
                    stk::percept::pout() << "AdaptMain:: unexpected exception: " << X.what() << "\n";
                    failed_proc_rank = p_rank+1u;
                  }
                  catch( ... ) {
                    std::cout << "AdaptMain::  ... exception" << std::endl;
                    failed_proc_rank = p_rank+1u;
                  }

                  //stk::all_reduce( run_environment.m_comm, stk::ReduceSum<1>( &failed_proc_rank ) );
                  if (failed_proc_rank)
                    {
                      std::cout << "P[" << p_rank << "] AdaptMain::exception found on processor " << (failed_proc_rank-1) << std::endl;
                      stk::percept::pout() << "P[" << p_rank << "]  exception found on processor " << (failed_proc_rank-1) << "\n";
                      exit(1);
                    }

                  if (DEBUG_ADAPT_MAIN)
                    {
                      stk::percept::pout() << "P[" << p_rank << ", " << p_size << "]  wall clock time on processor [" << p_rank << ", " << p_size << "]= " << (t1-t0) << " (sec) "
                                           << " cpu time= " << (cpu1 - cpu0) << " (sec)\n";
                      std::cout << "P[" << p_rank << ", " << p_size << "]  wall clock time on processor [" << p_rank << ", " << p_size << "]= " << (t1-t0) << " (sec) "
                                << " cpu time= " << (cpu1 - cpu0) << " (sec) " << std::endl;
                    }

                  double cpuMax = (cpu1-cpu0);
                  double wallMax = (t1-t0);
                  double cpuSum = (cpu1-cpu0);

                  stk::all_reduce( run_environment.m_comm, stk::ReduceSum<1>( &cpuSum ) );
                  stk::all_reduce( run_environment.m_comm, stk::ReduceMax<1>( &cpuMax ) );
                  stk::all_reduce( run_environment.m_comm, stk::ReduceMax<1>( &wallMax ) );

                  if (0 == p_rank)
                    {
                      stk::percept::pout() << "P[" << p_rank << ", " << p_size << "]  max wall clock time = " << wallMax << " (sec)\n";
                      stk::percept::pout() << "P[" << p_rank << ", " << p_size << "]  max cpu  clock time = " << cpuMax << " (sec)\n";
                      stk::percept::pout() << "P[" << p_rank << ", " << p_size << "]  sum cpu  clock time = " << cpuSum << " (sec)\n";
                      std::cout << "P[" << p_rank << ", " << p_size << "]  max wall clock time = " << wallMax << " (sec)" << std::endl;
                      std::cout << "P[" << p_rank << ", " << p_size << "]  max cpu  clock time = " << cpuMax << " (sec)" << std::endl;
                      std::cout << "P[" << p_rank << ", " << p_size << "]  sum cpu  clock time = " << cpuSum << " (sec)" << std::endl;
                    }
                } // m_iM

            } // do_normal_pass

#if STK_ADAPT_HAVE_YAML_CPP
          // fini pass
          if (m_W > 1 && m_iW == m_W)
            {
              VERIFY_OP_ON(s_spatialDim, >=, 2, "AdaptMain bad spatial_dim");
              PerceptMesh eMesh(s_spatialDim);
              eMesh.openEmpty();
              SerializeNodeRegistry snr(eMesh, 0, input_mesh, output_mesh, m_M, 0,  m_W, m_iW, m_M_0, m_M_1);
              snr.pass_final(i_pass);
            }
#endif

        } // i_pass

      return result;
    }

  }
}

//#include "pyencore.h"
//#if !PY_PERCEPT
int main(int argc, char **argv) {

  int res=0;
  res = stk::adapt::adapt_main(argc, argv);
  stk::adapt::ParallelMachineFinalize pm(true);
  return res;
}
//#endif
