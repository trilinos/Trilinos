// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/structured/StructuredGridRefiner.hpp>
#include <percept/PerceptUtils.hpp>
#include <adapt/UniformRefinerPattern.hpp>
#include <adapt/UniformRefiner.hpp>
#include <sys/time.h>
#include <stk_util/parallel/ParallelReduce.hpp>

namespace percept {

  void get_new_sizes(const SGridSizes& input_sizes,SGridSizes& output_sizes  ) {
    // new sizes
    const unsigned m_index_base =0;
    UInt sizes[3], new_sizes[3];
    for (UInt ii=0; ii < 3; ++ii)
      {
        sizes[ii] = input_sizes.node_max[ii] - input_sizes.node_min[ii] + 1;
        new_sizes[ii] = 2 * sizes[ii] - 1;
        output_sizes.node_min[ii] = 2*(input_sizes.node_min[ii] - m_index_base) + m_index_base;
        output_sizes.node_max[ii] = output_sizes.node_min[ii] + new_sizes[ii] - 1;
        output_sizes.cell_min[ii] = 2*(input_sizes.cell_min[ii] - m_index_base) + m_index_base;
        output_sizes.cell_max[ii] = output_sizes.cell_min[ii] + new_sizes[ii] - 2;
        output_sizes.node_size[ii] = new_sizes[ii];
        output_sizes.cell_size[ii] = new_sizes[ii] - 1;
        output_sizes.node_size_global[ii] = 2 * input_sizes.node_size_global[ii] - 1;
      }
  }

  unsigned StructuredGridRefiner::do_refine_structured()
  {
    unsigned num_refined_cells = 0;

#if HAVE_CGNS
    m_output->m_sblocks.resize(0);
    const int my_rank = stk::parallel_machine_rank(m_output->m_comm);

    double runTime = 0;
    for (unsigned iblock=0; iblock < m_input->m_sblocks.size(); ++iblock)
      {
        std::shared_ptr<StructuredBlock> sgi = m_input->m_sblocks[iblock];

        std::array<unsigned, 9> nijk{{0,0,0,0,0,0,0,0,0}};
        std::array<unsigned, 3> node_size_global{{sgi->m_sizes.node_size_global[0],
                sgi->m_sizes.node_size_global[1], sgi->m_sizes.node_size_global[2]}};

        std::shared_ptr<StructuredBlock> sgiNew ( new StructuredBlock(sgi->m_comm, iblock, nijk, node_size_global, sgi->m_name+"_refined", sgi->m_base, sgi->m_zone, m_output.get()) );
        m_output->m_sblocks.push_back(sgiNew);

        get_new_sizes(sgi->m_sizes,sgiNew->m_sizes);

        num_refined_cells += (sgiNew->m_sizes).cell_size[0]*
                             (sgiNew->m_sizes).cell_size[1]*
                             (sgiNew->m_sizes).cell_size[2];

        bool debug = m_debug > 0;
        if (debug)
          std::cout << "StructuredGridRefiner::do_refine_structured: iblock= " << iblock << " name= " << sgi->m_name
                    << " new sizes= " << sgiNew->m_sizes.node_size[0] << " " << sgiNew->m_sizes.node_size[1] << " " << sgiNew->m_sizes.node_size[2] << std::endl;

        const int A0 = sgi->m_access_ordering[0], A1 = sgi->m_access_ordering[1], A2 = sgi->m_access_ordering[2];
        if (debug)
          std::cout << "A0= " << A0 << " A1= " << A1 << " A2= " << A2
                    << " sizes= "
                    << sgiNew->m_sizes.node_size[A0] << " "
                    << sgiNew->m_sizes.node_size[A1] << " "
                    << sgiNew->m_sizes.node_size[A2]
                    << std::endl;

        Kokkos::resize(sgiNew->m_sgrid_coords,
                                     sgiNew->m_sizes.node_size[A0],
                                     sgiNew->m_sizes.node_size[A1],
                                     sgiNew->m_sizes.node_size[A2], 3);

        StructuredGridRefinerImpl<unsigned, uint64_t, MTSGridField::Array4D> refiner(sgi,sgiNew,m_debug);

        if (debug && my_rank == 0)
          std::cout << "StructuredGridRefiner: start block " << sgi->m_name << " refine..." << std::endl;

        struct timeval begin, end;
        gettimeofday(&begin,NULL);
        {
          stk::diag::Timer     my_timer(sgi->m_name, rootTimerStructured());
          stk::diag::Timer     my_refine_timer("refine", my_timer);
          stk::diag::TimeBlock my_timeblock(my_refine_timer);

          refiner.refine();
        }
        gettimeofday(&end,NULL);
        runTime += 1.0*(end.tv_sec-begin.tv_sec) +
          1.0e-6*(end.tv_usec-begin.tv_usec);

        // deal with zone connectivity
        {
          int nconn = sgi->m_zoneConnectivity.size();
          for (int i = 0; i < nconn; i++) {
            Ioss::ZoneConnectivity& zc = sgi->m_zoneConnectivity[i];

            std::string connectname = zc.m_connectionName;
            std::string donorname = zc.m_donorName;

            std::array<cgsize_t, 3> range_beg = zc.m_ownerRangeBeg;
            std::array<cgsize_t, 3> range_end = zc.m_ownerRangeEnd;
            std::array<cgsize_t, 3> donor_beg = zc.m_donorRangeBeg;
            std::array<cgsize_t, 3> donor_end = zc.m_donorRangeEnd;
            for (unsigned j=0; j < 3; ++j)
              {
                    range_beg[j] = 2*(range_beg[j] - m_index_base) + m_index_base;

                    range_end[j] = 2*(range_end[j] - m_index_base) + m_index_base;

                    donor_beg[j] = 2*(donor_beg[j] - m_index_base) + m_index_base;

                    donor_end[j] = 2*(donor_end[j] - m_index_base) + m_index_base;
              }
            std::array<int, 3> transform = zc.m_transform;
            int owner_zone = zc.m_ownerZone;
            (void) owner_zone;
            int donor_zone = zc.m_donorZone;

            sgiNew->m_zoneConnectivity.emplace_back(connectname, sgi->m_zone, donorname, donor_zone, transform,
                                                    range_beg, range_end, donor_beg, donor_end);
          }
        }

        // deal with boundary conditions
        {
          int nbc = sgi->m_boundaryConditions.size();
          for (int i = 0; i < nbc; i++) {
            Ioss::BoundaryCondition& bc = sgi->m_boundaryConditions[i];

            std::string bcname = bc.m_bcName;
            //std::string donorname = zc.m_donorName;
            // FIXME
            CG_BCType_t bctype = CGNS_ENUMV( BCTypeNull );
            (void)bctype;

            Ioss::IJK_t range_beg = bc.m_rangeBeg;
            Ioss::IJK_t range_end = bc.m_rangeEnd;
            for (unsigned j=0; j < 3; ++j)
              {
                    range_beg[j] = 2*(range_beg[j] - m_index_base) + m_index_base;

                    range_end[j] = 2*(range_end[j] - m_index_base) + m_index_base;
              }

            sgiNew->m_boundaryConditions.emplace_back(bcname, range_beg, range_end);
          }
        }
      }

    stk::all_reduce(m_output->m_comm, stk::ReduceSum<1>(&num_refined_cells));

    if (my_rank == 0) {
        std::cout << "Total number of cells after refinement = " << num_refined_cells << std::endl;
        std::cout << "Runtime (seconds) over " <<  m_input->m_sblocks.size() << " blocks was " << runTime << std::endl;
    }
#endif

    return num_refined_cells;
  }

  unsigned StructuredGridRefiner::do_refine()
  {
    return do_refine_structured();
  }

  void StructuredGridRefiner::print(std::ostream& out, int level)
  {
    out << "initial grid: " << std::endl;
    m_input->print(out, level);
    out << "refined grid: " << std::endl;
    m_output->print(out,level);
  }

  // Create a structured block with twice the cells in each direction
  //   - for now, should get called in the same order as the create_structured_block
  //       was called to get the node_ and cell_offsets correct.

  // This is just prototype code... FIXME
  void StructuredGridRefiner::double_structured_block(Ioss::Region *new_region,
                                                      Ioss::StructuredBlock *oblock, size_t& num_node, size_t& num_cell)
  {
#if 0
    // FIXME
    //typedef size_t cgsize_t;

    cgsize_t base;
    int zone;

    cgsize_t size[9]{0,0,0,0,0,0,0,0,0};
    std::string zone_name;

    base = oblock->get_property("base").get_int();
    zone = oblock->get_property("zone").get_int();


    size_t node_offset = oblock->get_node_offset();
    (void) node_offset;
    size_t cell_offset = oblock->get_cell_offset();
    (void) cell_offset;
    size_t nnode = oblock->get_property("node_count").get_int();
    (void) nnode;
    size_t ncell = oblock->get_property("cell_count").get_int();
    (void) ncell;
    int ni = 1+ oblock->get_property("ni").get_int();
    int nj = 1+ oblock->get_property("nj").get_int();
    int nk = 1+ oblock->get_property("nk").get_int();

    Iocgns::DatabaseIO *dbio = reinterpret_cast< Iocgns::DatabaseIO *>(oblock->get_database());
    VERIFY_OP_ON(dbio, !=, 0, "must pass in a CGNS type DatabaseIO");
    auto zmap = dbio->getZoneNameMap();
    for (auto& zz : zmap)
      {
        if (zz.second == zone)
          {
            zone_name = zz.first;
            break;
          }
      }
    VERIFY_OP_ON(zone_name.size(), !=, 0, "bad zone_name");

    size[0] = 2*ni-1;
    size[1] = 2*nj-1;
    size[2] = 2*nk-1;
    size[3] = size[0]-1;
    size[4] = size[1]-1;
    size[5] = size[2]-1;

    cgsize_t index_dim = oblock->get_property("component_degree").get_int();
    assert(index_dim == 3);  //only for 3D

    // An Ioss::StructuredBlock corresponds to a CG_Structured zone...
    Ioss::StructuredBlock *block =
      new Ioss::StructuredBlock(dbio, zone_name.c_str(), index_dim, size[3], size[4], size[5]);

    new_region->add(block);

    // get_region()->remove(oblock);

    block->property_add(Ioss::Property("base", base));
    block->property_add(Ioss::Property("zone", zone));
    dbio->get_region()->add(block);

    block->set_node_offset(num_node);
    block->set_cell_offset(num_cell);
    num_node += block->get_property("node_count").get_int();
    num_cell += block->get_property("cell_count").get_int();

    // Handle zone-grid-connectivity...
    int nconn = oblock->m_zoneConnectivity.size();
    //cg_n1to1(cgnsFilePtr, base, zone, &nconn);
    for (int i = 0; i < nconn; i++) {
      Ioss::ZoneConnectivity& zc = oblock->m_zoneConnectivity[i];
      //(connectname, zone, donorname, donor_zone, transform,
      //range_beg, range_end, donor_beg, donor_end);

      std::string connectname = zc.m_connectionName;
      std::string donorname = zc.m_donorName;

      std::array<cgsize_t, 3> range_beg = zc.m_ownerRangeBeg;
      std::array<cgsize_t, 3> range_end = zc.m_ownerRangeEnd;
      std::array<cgsize_t, 3> donor_beg = zc.m_donorRangeBeg;
      std::array<cgsize_t, 3> donor_end = zc.m_donorRangeEnd;
      for (unsigned i=0; i < 3; ++i)
        {
          range_beg[i] = 2*range_beg[i] - 1;
          range_end[i] = 2*range_end[i] - 1;
          donor_beg[i] = 2*donor_beg[i] - 1;
          donor_end[i] = 2*donor_end[i] - 1;
        }
      std::array<int, 3> transform = zc.m_transform;
      int owner_zone = zc.m_ownerZone;
      (void) owner_zone;
      int donor_zone = zc.m_donorZone;

      block->m_zoneConnectivity.emplace_back(connectname, zone, donorname, donor_zone, transform,
                                             range_beg, range_end, donor_beg, donor_end);
    }
#endif
  }



}
