// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/Percept.hpp>
#if !STK_PERCEPT_LITE

#include <percept/structured/StructuredBlock.hpp>
#include <percept/structured/BlockStructuredGrid.hpp>
#include <percept/PerceptUtils.hpp>


namespace percept {


  StructuredBlock::StructuredBlock(stk::ParallelMachine comm, unsigned iblock, Ioss::StructuredBlock *sb, BlockStructuredGrid *bsg) :
    m_comm(comm), m_iblock(iblock), m_sblock(sb), m_parent(bsg), m_name(sb->name()),
    m_sgrid_coords_impl(new Array4D),
    m_sgrid_coords(*m_sgrid_coords_impl.get()),
    m_loop_ordering{{0,1,2}}, m_access_ordering{{0,1,2}}
    {
      m_fields["coordinates"] = m_sgrid_coords_impl;

      // Ioss::StructuredBlock stores the "number of intervals" in each direction - add one to get number of nodes
      unsigned nijk[3] = {0,0,0};
      nijk[0] = 1 + sb->get_property("ni").get_int();
      nijk[1] = 1 + sb->get_property("nj").get_int();
      nijk[2] = 1 + sb->get_property("nk").get_int();
      m_base = sb->get_property("base").get_int();
      m_zone = sb->get_property("zone").get_int();

      m_sizes.node_size_global[0] = 1 + sb->get_property("ni_global").get_int();
      m_sizes.node_size_global[1] = 1 + sb->get_property("nj_global").get_int();
      m_sizes.node_size_global[2] = 1 + sb->get_property("nk_global").get_int();

      for (unsigned ii=0; ii < 3; ++ii)
        {
          m_sizes.node_min[ii] = 0;
          m_sizes.node_max[ii] = nijk[ii] - 1;
          m_sizes.cell_min[ii] = 0;
          m_sizes.cell_max[ii] = nijk[ii] - 2;
          m_sizes.node_size[ii] = m_sizes.node_max[ii] - m_sizes.node_min[ii] + 1;
          m_sizes.cell_size[ii] = m_sizes.cell_max[ii] - m_sizes.cell_min[ii] + 1;
        }
      m_zoneConnectivity = m_sblock->m_zoneConnectivity;
      m_boundaryConditions = m_sblock->m_boundaryConditions;
    }

  StructuredBlock::StructuredBlock(stk::ParallelMachine comm, unsigned iblock, const std::array<unsigned,9>& nijk, const std::array<unsigned,3>& node_size_global, const std::string& name, int base, int zone, BlockStructuredGrid *bsg) :
    m_comm(comm), m_iblock(iblock), m_sblock(0), m_parent(bsg), m_name(name), m_base(base), m_zone(zone),
    m_sgrid_coords_impl(new Array4D),
    m_sgrid_coords(*m_sgrid_coords_impl.get()),
    m_loop_ordering{{0,1,2}}, m_access_ordering{{0,1,2}}
    {
      //VERIFY_OP_ON(stk::parallel_machine_size(m_comm), ==, 1, "only for serial");
      m_fields["coordinates"] = m_sgrid_coords_impl;

      // Ioss::StructuredBlock stores the "number of intervals" in each direction - add one to get number of nodes

      for (unsigned ii=0; ii < 3; ++ii)
        {
          m_sizes.node_min[ii] = 0;
          m_sizes.node_max[ii] = nijk[ii] - 1;
          m_sizes.cell_min[ii] = 0;
          m_sizes.cell_max[ii] = nijk[ii] - 2;
          m_sizes.node_size[ii] = m_sizes.node_max[ii] - m_sizes.node_min[ii] + 1;
          m_sizes.node_size_global[ii] = node_size_global[ii];
          m_sizes.cell_size[ii] = m_sizes.cell_max[ii] - m_sizes.cell_min[ii] + 1;
        }
    }

  bool StructuredBlock::is_empty()
  {
    return (m_sizes.node_size[0] <= 1 || m_sizes.node_size[1] <= 1 || m_sizes.node_size[2] <= 1);
  }

  std::shared_ptr<StructuredBlock::Array4D> StructuredBlock::register_field(const std::string& field, unsigned fourth_dim)
  {
    const unsigned A0 = m_access_ordering[0], A1 = m_access_ordering[1], A2 = m_access_ordering[2];
    unsigned Asizes[3] = {m_sizes.node_size[A0], m_sizes.node_size[A1], m_sizes.node_size[A2]};
    auto iter = m_fields.find(field);
    std::shared_ptr<Array4D> nfield;
    if (iter == m_fields.end())
      {
        nfield.reset(new Array4D);
        m_fields[field] = nfield;
      }
    else
      {
        nfield = iter->second;
      }
    Kokkos::resize(*nfield.get(), Asizes[0], Asizes[1], Asizes[2], fourth_dim);
    return nfield;
  }

  // static
  std::shared_ptr<StructuredBlock>
  StructuredBlock::fixture_1(stk::ParallelMachine comm, std::array<unsigned,3> sizes, unsigned iblock, int base, int zone, BlockStructuredGrid *bsg, std::array<double,3> dim_width,std::array<double,3> dim_offset)
  {
    std::array<unsigned,9> nijk{{sizes[0],sizes[1],sizes[2], 0,0,0, 0,0,0}};
    std::array<unsigned, 3> node_size_global{{sizes[0],sizes[1],sizes[2]}};

    std::shared_ptr<StructuredBlock> sgi(new StructuredBlock(comm, iblock, nijk, node_size_global, "test", base, zone, bsg));

    const unsigned L0 = sgi->m_loop_ordering[0], L1 = sgi->m_loop_ordering[1], L2 = sgi->m_loop_ordering[2];
    const unsigned A0 = sgi->m_access_ordering[0], A1 = sgi->m_access_ordering[1], A2 = sgi->m_access_ordering[2];

    for (unsigned ii=0; ii < 3; ++ii)
      {
        sgi->m_sizes.node_min[ii] = 0;
        sgi->m_sizes.node_max[ii] = sizes[ii]-1;
        sgi->m_sizes.cell_min[ii] = 0;
        sgi->m_sizes.cell_max[ii] = sizes[ii]-2;
        sgi->m_sizes.node_size[ii] = sizes[ii];
        sgi->m_sizes.cell_size[ii] = sizes[ii]-1;
      }

    unsigned Asizes[3] = {sgi->m_sizes.node_size[A0], sgi->m_sizes.node_size[A1], sgi->m_sizes.node_size[A2]};

    Kokkos::resize(sgi->m_sgrid_coords, Asizes[0], Asizes[1], Asizes[2], 3);

    unsigned indx[3]{0,0,0};

    // FIXME indexing from 0

    Array4D::HostMirror interimNodes;
    interimNodes = Kokkos::create_mirror_view(*sgi->m_sgrid_coords_impl);

    for (indx[L2] = 0; indx[L2] < sizes[L2]; ++indx[L2])
      {
        for (indx[L1] = 0; indx[L1] < sizes[L1]; ++indx[L1])
          {
            for (indx[L0] = 0; indx[L0] < sizes[L0]; ++indx[L0])
              {
                //auto ilo = sgi->local_offset(indx[A0], indx[A1], indx[A2], &Asizes[0]);
                double xyz[3] = {  ( double(indx[A0])/double(sizes[A0]-1) ), double(indx[A1])/double(sizes[A1]-1), double(indx[A2])/double(sizes[A2]-1)};
                for (unsigned ic = 0; ic < 3; ++ic)
                  {
                    interimNodes(indx[A0], indx[A1], indx[A2], ic) = dim_width[ic]*xyz[ic] + dim_offset[ic];
                  }
              }
          }
      }
    Kokkos::deep_copy( *sgi->m_sgrid_coords_impl, interimNodes);

    // bc's
    int isize[3] = {(int)sizes[0], (int)sizes[1], (int)sizes[2]};
    //madbrew: this was inconsistent with imported cgns meshes (which are 1 based indexing) so I switched from 0 based to be consistent
    sgi->m_boundaryConditions.emplace_back("BCInflow",  Ioss::IJK_t({{1,                   1,          1}}), Ioss::IJK_t({{1,          isize[1], isize[2]}}) );
    sgi->m_boundaryConditions.emplace_back("BCOutflow", Ioss::IJK_t({{isize[0],          1,          1}}), Ioss::IJK_t({{isize[0], isize[1], isize[2]}}) );
    sgi->m_boundaryConditions.emplace_back("BCWall",    Ioss::IJK_t({{1,                   1,          1}}), Ioss::IJK_t({{isize[0],          1, isize[2]}}) );
    sgi->m_boundaryConditions.emplace_back("BCWall",    Ioss::IJK_t({{1,          isize[1],          1}}), Ioss::IJK_t({{isize[0], isize[1], isize[2]}}) );
    sgi->m_boundaryConditions.emplace_back("BCWall",    Ioss::IJK_t({{1,          1,                   1}}), Ioss::IJK_t({{isize[0], isize[1],          1}}) );
    sgi->m_boundaryConditions.emplace_back("BCWall",    Ioss::IJK_t({{1,          1,          isize[2]}}), Ioss::IJK_t({{isize[0], isize[1], isize[2]}}) );

#if HAVE_CGNS
    sgi->m_bc_types.push_back(CGNS_ENUMV( BCInflow ));
    sgi->m_bc_types.push_back(CGNS_ENUMV( BCOutflow ));
    sgi->m_bc_types.push_back(CGNS_ENUMV( BCWall ));
    sgi->m_bc_types.push_back(CGNS_ENUMV( BCWall ));
    sgi->m_bc_types.push_back(CGNS_ENUMV( BCWall ));
    sgi->m_bc_types.push_back(CGNS_ENUMV( BCWall ));
#endif

    return sgi;
  }

  uint64_t StructuredBlock::local_offset(uint64_t i, uint64_t j, uint64_t k, unsigned *sizes)
  {
    // assumes 1-based i, j, k - returns 0-based offset
    //(k - 1) * (sizes[0] + 1) * (sizes[1] + 1) + (j - 1) * (sizes[0] + 1) + i - 1;

    // assume 0-based ijk, returns 0-based offset
    return k * sizes[0] * sizes[1] + j * sizes[0] + i;
  }


  void StructuredBlock::multi_dim_indices_from_local_offset(uint64_t local_offset, std::array<unsigned,3>& indx)
  {
    const int L0 = m_loop_ordering[0], L1 = m_loop_ordering[1], L2 = m_loop_ordering[2];
    const unsigned sizes[3] = {
        1+ m_sizes.node_max[L0] - m_sizes.node_min[L0],
        1+ m_sizes.node_max[L1] - m_sizes.node_min[L1],
        1+ m_sizes.node_max[L2] - m_sizes.node_min[L2]
        };

    indx[L2] = m_sizes.node_min[L2] + (local_offset / (sizes[L0]*sizes[L1] ));
    indx[L1] = m_sizes.node_min[L1] + ((local_offset / sizes[L0]) % sizes[L1] );
    indx[L0] = m_sizes.node_min[L0] + (local_offset % sizes[L0]);
  }

  void StructuredBlock::read_cgns()
  {
    stk::diag::Timer     my_timer(m_name, rootTimerStructured());
    stk::diag::Timer     my_read_timer("read_cgns", my_timer);
    stk::diag::TimeBlock my_timeblock(my_read_timer);

    const unsigned L0 = m_loop_ordering[0], L1 = m_loop_ordering[1], L2 = m_loop_ordering[2];
    const unsigned A0 = m_access_ordering[0], A1 = m_access_ordering[1], A2 = m_access_ordering[2];

    unsigned sizes[3] = {m_sizes.node_size[0], m_sizes.node_size[1], m_sizes.node_size[2]};
    unsigned Asizes[3] = {m_sizes.node_size[A0], m_sizes.node_size[A1], m_sizes.node_size[A2]};

    Kokkos::resize(m_sgrid_coords, Asizes[0], Asizes[1], Asizes[2], 3);

    auto &block = *m_sblock;
    std::vector<double> coord_tmp[3];
    block.get_field_data("mesh_model_coordinates_x", coord_tmp[0]);
    block.get_field_data("mesh_model_coordinates_y", coord_tmp[1]);
    block.get_field_data("mesh_model_coordinates_z", coord_tmp[2]);

    if (is_empty())
      return;

    unsigned indx[3]{0,0,0};

    Array4D::HostMirror interimNodes;
    interimNodes = Kokkos::create_mirror_view(*m_sgrid_coords_impl);

    for (indx[L2] = 0; indx[L2] < sizes[L2]; ++indx[L2])
      {
        for (indx[L1] = 0; indx[L1] < sizes[L1]; ++indx[L1])
          {
            for (indx[L0] = 0; indx[L0] < sizes[L0]; ++indx[L0])
              {
                auto ilo = local_offset(indx[A0], indx[A1], indx[A2], &Asizes[0]);
                for (unsigned ic = 0; ic < 3; ++ic)
                  {
                    interimNodes(indx[A0], indx[A1], indx[A2], ic) = coord_tmp[ic][ilo];
                  }
              }
          }
      }

    {
      stk::diag::Timer     my_copy_timer("deep_copy_to_device", my_timer);
      stk::diag::TimeBlock my_timeblock2(my_copy_timer);

      Kokkos::deep_copy(*m_sgrid_coords_impl,interimNodes);
    }
  }

  void StructuredBlock::print(std::ostream& out, int level)
  {
    Array4D::HostMirror interimNodes;
    interimNodes = Kokkos::create_mirror_view(*m_sgrid_coords_impl);
    Kokkos::deep_copy(interimNodes,*m_sgrid_coords_impl);

    out << "StructuredBlock: iblock= " << m_iblock << " name= " << m_name << " base= " << m_base << " zone= " << m_zone << std::endl;
    out << "node_min         = " << m_sizes.node_min[0] << " " << m_sizes.node_min[1] << " " << m_sizes.node_min[2] << std::endl;
    out << "node_max         = " << m_sizes.node_max[0] << " " << m_sizes.node_max[1] << " " << m_sizes.node_max[2] << std::endl;
    out << "node_size        = " << m_sizes.node_size[0] << " " << m_sizes.node_size[1] << " " << m_sizes.node_size[2] << std::endl;
    out << "node_size_global = " << m_sizes.node_size_global[0] << " " << m_sizes.node_size_global[1] << " " << m_sizes.node_size_global[2] << std::endl;
    out << "loop_ordering    = " << m_loop_ordering[0] << " " << m_loop_ordering[1] << " " << m_loop_ordering[2] << std::endl;
    out << "access_ordering  = " << m_access_ordering[0] << " " << m_access_ordering[1] << " " << m_access_ordering[2] << std::endl;
    if (level > 0)
      {
        out << "ZoneConnectivity: num= " << m_zoneConnectivity.size() << std::endl;
        for (unsigned i=0; i < m_zoneConnectivity.size(); ++i)
          {
            out << "zc[" << i << "]= " << m_zoneConnectivity[i] << std::endl;
          }
        out << "BoundaryCondition: num= " << m_boundaryConditions.size() << std::endl;
        for (unsigned i=0; i < m_boundaryConditions.size(); ++i)
          {
            out << "bc[" << i << "]= " << m_boundaryConditions[i] << std::endl;
          }
      }
    if (level > 2 && !is_empty())
      {
        const unsigned L0 = m_loop_ordering[0], L1 = m_loop_ordering[1], L2 = m_loop_ordering[2];
        const unsigned A0 = m_access_ordering[0], A1 = m_access_ordering[1], A2 = m_access_ordering[2];

        unsigned sizes[3] = {m_sizes.node_size[0], m_sizes.node_size[1], m_sizes.node_size[2]};

        unsigned indx[3]{0,0,0};

        for (indx[L2] = 0; indx[L2] < sizes[L2]; ++indx[L2])
          {
            for (indx[L1] = 0; indx[L1] < sizes[L1]; ++indx[L1])
              {
                for (indx[L0] = 0; indx[L0] < sizes[L0]; ++indx[L0])
                  {
                    for (unsigned ic = 0; ic < 3; ++ic)
                      {
                    	out << "coord: " << indx[A0] << " " << indx[A1] << " " << indx[A2] << " "
                    	                            <<interimNodes(indx[A0], indx[A1], indx[A2], ic) << std::endl;
                      }
                  }
              }
          }
      }
  }

  void StructuredBlock::dump_vtk(const std::string& file_prefix)
  {
    Array4D::HostMirror interimNodes;
    interimNodes = Kokkos::create_mirror_view(*m_sgrid_coords_impl);
    Kokkos::deep_copy(interimNodes,*m_sgrid_coords_impl);

    const unsigned L0 = m_loop_ordering[0], L1 = m_loop_ordering[1], L2 = m_loop_ordering[2];
    const unsigned A0 = m_access_ordering[0], A1 = m_access_ordering[1], A2 = m_access_ordering[2];

    std::ostringstream fstr;
    fstr << file_prefix << "." << m_iblock;
    if (stk::parallel_machine_size(m_comm) == 1)
      fstr << ".vts";
    else
      fstr << "." << stk::parallel_machine_size(m_comm) << "." << stk::parallel_machine_rank(m_comm) << ".vts";

    std::ofstream out(fstr.str().c_str());
    //std::string
    char buf[1024];
    sprintf(buf,
            "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"
            "  <StructuredGrid WholeExtent=\"0 %u 0 %u 0 %u\">\n"
            "  <Piece Extent=\"0 %u 0 %u 0 %u\">\n"
            "    <PointData>\n"
            "    </PointData>\n",
            m_sizes.node_size_global[A0]-1, m_sizes.node_size_global[A1]-1, m_sizes.node_size_global[A2]-1,
            m_sizes.node_size[A0]-1, m_sizes.node_size[A1]-1, m_sizes.node_size[A2]-1
            );
    out << buf << std::endl;

    if (1)
      {
        out <<
          "    <Points>\n"
          "       <DataArray type=\"Float64\" Name=\"Points\"  NumberOfComponents=\"3\" format=\"ascii\">\n";

        {

          std::ios_base::fmtflags old_flags = out.flags();
          int old_precision = out.precision();
          out.precision(4);
          out << std::fixed;

          unsigned sizes[3] = {m_sizes.node_size[0], m_sizes.node_size[1], m_sizes.node_size[2]};
          unsigned indx[3]{0,0,0};

          for (indx[L2] = 0; indx[L2] < sizes[L2]; ++indx[L2])
            {
              for (indx[L1] = 0; indx[L1] < sizes[L1]; ++indx[L1])
                {
                  for (indx[L0] = 0; indx[L0] < sizes[L0]; ++indx[L0])
                    {
                      out << "         ";
                      for (unsigned ic = 0; ic < 3; ++ic)
                        {
                    	  double val = interimNodes(indx[A0], indx[A1], indx[A2], ic);

                    	  if (std::fabs(val)<std::numeric_limits<double>::epsilon()) val = 0.0;
                    	  out << std::setw(16) << val;
                        }
                      out << "\n";
                    }
                }
            }

          out.flags(old_flags);
          out.precision(old_precision);
        }

        out <<
          "       </DataArray>\n"
          "    </Points>\n";
      }

    if (!is_empty())
      {
        bool saveProcId = true;
        bool saveBlockId = true;
        if (saveProcId || saveBlockId)
          {
            out <<
              "    <CellData Scalars=\"" << (saveProcId ? " proc_id " : "") << (saveBlockId ? " block_id" : "") << "\">\n";
          }
        if (saveProcId)
          {
            out <<
              "       <DataArray Name=\"proc_id\" type=\"Int32\" format=\"ascii\">\n";

            unsigned sizes[3] = {m_sizes.node_size[0]-1, m_sizes.node_size[1]-1, m_sizes.node_size[2]-1 };
            unsigned indx[3]{0,0,0};
            size_t count=0;
            out << "         ";
            for (indx[L2] = 0; indx[L2] < sizes[L2]; ++indx[L2])
              {
                for (indx[L1] = 0; indx[L1] < sizes[L1]; ++indx[L1])
                  {
                    for (indx[L0] = 0; indx[L0] < sizes[L0]; ++indx[L0])
                      {
                        out << " " << stk::parallel_machine_rank(m_comm);
                        ++count;
                        if (count % 10 == 0)
                          out << "\n         ";
                      }
                  }
              }
            out << "\n"
              "       </DataArray>\n";
          }

        if (saveBlockId)
          {
            out <<
              "       <DataArray Name=\"block_id\" type=\"Float32\" format=\"ascii\">\n";

            size_t nblocks = m_parent->m_sblocks.size();
            unsigned sizes[3] = {m_sizes.node_size[0]-1, m_sizes.node_size[1]-1, m_sizes.node_size[2]-1 };
            unsigned indx[3]{0,0,0};
            size_t count=0;
            out << "         ";
            for (indx[L2] = 0; indx[L2] < sizes[L2]; ++indx[L2])
              {
                for (indx[L1] = 0; indx[L1] < sizes[L1]; ++indx[L1])
                  {
                    for (indx[L0] = 0; indx[L0] < sizes[L0]; ++indx[L0])
                      {
                        // workaround a bug in paraview - it appears to scale the mesh using all data, not just coords
                        out << " " << double(m_iblock) / double(nblocks);
                        ++count;
                        if (count % 10 == 0)
                          out << "\n         ";
                      }
                  }
              }
            out << "\n"
              "       </DataArray>\n";

          }
        if (saveProcId || saveBlockId)
          {
            out <<
              "    </CellData>\n";
          }
      }

    out <<
      "  </Piece>\n"
      "  </StructuredGrid>\n"
      "</VTKFile>\n" << std::endl;
  }//dump_vtk

}


#endif

