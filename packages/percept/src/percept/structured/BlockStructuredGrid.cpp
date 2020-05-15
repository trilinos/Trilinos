// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/structured/BlockStructuredGrid.hpp>
#include <percept/MeshType.hpp>

#include <stk_util/parallel/CommSparse.hpp>

#include <Ioss_Region.h>

namespace percept {

  using Array4D = typename MTSGridField::Array4D;

  static const int permute_to_unstructured_local[8] = {0,1,3,2, 4,5,7,6};
  const int *BlockStructuredGrid::permute_to_unstructured = &permute_to_unstructured_local[0];

  BlockStructuredGrid::BlockStructuredGrid(stk::ParallelMachine comm, Ioss::Region *region) :  m_comm(comm), m_region(region)
  {}

  void BlockStructuredGrid::print(std::ostream& out, int level)
  {
    for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
      {
        m_sblocks[iblock]->print(out, level);
      }
  }

  void BlockStructuredGrid::read_cgns()
  {
    auto& blocks = m_region->get_structured_blocks();
    unsigned iblock=0;
    m_sblocks.resize(0);

    for (auto &block : blocks) {
      std::shared_ptr<StructuredBlock> ptr ( new StructuredBlock(m_comm, iblock, block, this) );
      ptr->read_cgns();    
      m_sblocks.push_back(ptr);
      ++iblock;
    }
    std::cout << stk::parallel_machine_rank(m_comm) << " nblocks= " << m_sblocks.size() << std::endl;
  }


  void BlockStructuredGrid::register_field(const std::string& field, unsigned fourth_dim)
  {
    std::shared_ptr<MTSGridField> nf (new MTSGridField(field));
    m_fields[field] = nf;

    nf->m_block_fields.resize(m_sblocks.size());
    for (size_t ib=0; ib < m_sblocks.size(); ++ib)
      {
        nf->m_block_fields[ib] = m_sblocks[ib]->register_field(field, fourth_dim);
      }
  }

  void BlockStructuredGrid::create_pvd(const std::string& file_prefix, bool paraviewBroken)
  {
    if (stk::parallel_machine_rank(m_comm))
      return;

    std::ofstream out(file_prefix+".pvd");
    out << "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    out << "  <Collection>\n";

    if (paraviewBroken)
      {
        int p_size = stk::parallel_machine_size(m_comm);
        for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
          {
            for (int p_rank=0; p_rank < p_size; ++p_rank)
              {
                out << "      <DataSet part=\"" << iblock << "\" file=\"" << file_prefix << "." << iblock;
                if (p_size > 1 )
                  out << "." << p_size << "." << p_rank;
                out << ".vts\"/>" << std::endl;
              }
          }
      }
    else
      {
        std::string par_prefix = "";
        if (stk::parallel_machine_size(m_comm) > 1)
          par_prefix = "p";
        for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
          {
            out << "      <DataSet part=\"" << iblock << "\" file=\"" << file_prefix << "." << iblock << "." + par_prefix + "vts\"/>" << std::endl;
            m_sblocks[iblock]->dump_vtk(file_prefix);
          }
      }

    out << "  </Collection>\n";
    out << "</VTKFile>\n";
  }

  /// creates a file, one per block, with pointers to all the parallel parts 
  ///   of the block in individual prefix_block.proc.vts files

  void BlockStructuredGrid::create_pvts(const std::string& file_prefix)
  {
    if (stk::parallel_machine_size(m_comm) == 1)
      return;

    stk::CommSparse comm_sparse(m_comm);
    //unsigned proc_rank = comm_sparse.parallel_rank();
    unsigned proc_size = comm_sparse.parallel_size();

    // send grid sizes to proc 0
    for (int stage=0; stage < 2; ++stage)
      {
        for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
          {
            std::shared_ptr<StructuredBlock> sb = m_sblocks[iblock];

            for (int ii=0; ii < 3; ++ii)
              comm_sparse.send_buffer( 0 ).pack< unsigned > (sb->m_sizes.node_size[ii]);
          }

        if (stage == 0)
          {
            comm_sparse.allocate_buffers();
          }
        else
          {
            comm_sparse.communicate();
          }
      }

    if (stk::parallel_machine_rank(m_comm) == 0)
      {
        MDArrayUInt block_info(stk::parallel_machine_size(m_comm), m_sblocks.size(), 3);

        for(unsigned from_proc = 0; from_proc < proc_size; ++from_proc )
          {
            //if (from_proc == proc_rank)
            //  continue;
            stk::CommBuffer & recv_buffer = comm_sparse.recv_buffer( from_proc );

            for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
              {
                for (unsigned ii=0; ii < 3; ++ii)
                  recv_buffer.unpack< unsigned >( block_info(from_proc, iblock, ii) );
              }
            VERIFY_OP_ON(recv_buffer.remaining(), ==, false, "bad unpack");
          }

        for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
          {
            std::shared_ptr<StructuredBlock> sb = m_sblocks[iblock];

            std::ofstream out1(file_prefix+"."+toString(iblock)+".pvts");
            char buf0[1024];
            sprintf(buf0,
                    "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >\n"
                    "  <PStructuredGrid WholeExtent=\"0 %u 0 %u 0 %u\" GhostLevel=\"0\" >\n"
                    "      <PPoints>\n"
                    "         <PDataArray NumberOfComponents=\"3\" type=\"Float32\" />\n"
                    "      </PPoints>\n",
                    sb->m_sizes.node_size_global[0]-1, sb->m_sizes.node_size_global[1]-1, sb->m_sizes.node_size_global[2]-1);
            out1 << buf0;

            for(unsigned from_proc = 0; from_proc < proc_size; ++from_proc )
              {
                char buf1[1024];
                sprintf(buf1,
                        "      <PPiece Extent=\"0 %u 0 %u 0 %u\" Source=\"%s.%u.%u.vts\" />\n",
                        block_info(from_proc, iblock, 0) - 1,
                        block_info(from_proc, iblock, 1) - 1,
                        block_info(from_proc, iblock, 2) - 1,
                        file_prefix.c_str(), iblock, from_proc);
                out1 << buf1;
              }

            out1 <<  "  </PStructuredGrid>\n";
            out1 << "</VTKFile>​\n";
          }
      }
  }

  void BlockStructuredGrid::dump_vtk(const std::string& file_prefix)
  {
    if (stk::parallel_machine_rank(m_comm) == 0)
      std::cout << "BlockStructuredGrid: saving VTK format file= " << file_prefix << " ... " << std::endl;

    bool paraviewBroken = true;
    create_pvd(file_prefix, paraviewBroken);
    // Paraview is broken, doesn't seem to like pvts, so, we just create more .vts pieces
    if (!paraviewBroken)
      create_pvts(file_prefix);

    for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
      {
        m_sblocks[iblock]->dump_vtk(file_prefix);
      }
  }

  void BlockStructuredGrid::get_nodes(std::vector<StructuredCellIndex>& nodes, unsigned offset)
  {
    nodes.resize(0);
    for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
      {
        std::shared_ptr<StructuredBlock> sgrid = m_sblocks[iblock];
        const unsigned L0 = sgrid->m_loop_ordering[0], L1 = sgrid->m_loop_ordering[1], L2 = sgrid->m_loop_ordering[2];
        if (sgrid->is_empty())
          continue;

        unsigned sizes[3] = {sgrid->m_sizes.node_size[0], sgrid->m_sizes.node_size[1], sgrid->m_sizes.node_size[2]};

        unsigned indx[3]{0,0,0};

        for (indx[L2] = 0; indx[L2] < sizes[L2]-offset; ++indx[L2])
          {
            for (indx[L1] = 0; indx[L1] < sizes[L1]-offset; ++indx[L1])
              {
                for (indx[L0] = 0; indx[L0] < sizes[L0]-offset; ++indx[L0])
                  {
                    StructuredCellIndex node{{indx[0], indx[1], indx[2], iblock}};
                    nodes.push_back(node);
                  }
              }
          }
      }
  }

  void BlockStructuredGrid::get_element_nodes(const StructuredCellIndex& element, std::vector<StructuredCellIndex>& nodes)
  {
    nodes.resize(0);
    for (unsigned i2 = element[2]; i2 <= element[2]+1; ++i2)
      {
        for (unsigned i1 = element[1]; i1 <= element[1]+1; ++i1)
          {
            for (unsigned i0 = element[0]; i0 <= element[0]+1; ++i0)
              {
                nodes.push_back(StructuredCellIndex{{i0,i1,i2,element[3]}});
              }
          }
      }
  }

void BlockStructuredGrid::get_elements(
        std::vector<StructuredCellIndex>& elements) { //writes all the elements of a structured block at certain index inside m-sblocks to an std::vector
    get_nodes(elements, 1);
}

void BlockStructuredGrid::get_nodes_of_sb(
        std::vector<StructuredCellIndex>& nodes, unsigned iBlock,
        unsigned offset) { //writes all the ijkb components of the nodes or elements of a structured block at certain index inside m_sblocks to an std::vector of SCell indices
    if (iBlock > m_sblocks.size()) {
        std::cout << "Block index out of range, exiting get_nodes_of_sb"
                << std::endl;
        return;
    }

    std::shared_ptr<StructuredBlock> sgrid = m_sblocks[iBlock];
    const unsigned L0 = sgrid->m_loop_ordering[0], L1 =
            sgrid->m_loop_ordering[1], L2 = sgrid->m_loop_ordering[2];

    if (sgrid->is_empty()) {
        std::cout << "Block is empty, exiting get_nodes_of_sb" << std::endl;
        return;
    }

    unsigned sizes[3] = { sgrid->m_sizes.node_size[0] - offset,
            sgrid->m_sizes.node_size[1] - offset, sgrid->m_sizes.node_size[2]
                    - offset };

    unsigned total_entities = sizes[0] * sizes[1] * sizes[2];

    nodes.resize(total_entities);

    unsigned indx[3] { 0, 0, 0 };

    unsigned iEnt = 0;

    for (indx[L2] = 0; indx[L2] < sizes[L2]; ++indx[L2]) {
        for (indx[L1] = 0; indx[L1] < sizes[L1]; ++indx[L1]) {
            for (indx[L0] = 0; indx[L0] < sizes[L0]; ++indx[L0]) {
                StructuredCellIndex node { { indx[0], indx[1], indx[2], iBlock } };
                nodes[iEnt][0] = node[0];
                nodes[iEnt][1] = node[1];
                nodes[iEnt][2] = node[2];
                nodes[iEnt][3] = node[3];
                iEnt++;
            }
        }
    }
} //get_nodes_sb



void BlockStructuredGrid::get_elements_of_sb(
        std::vector<StructuredCellIndex>&elements, unsigned iBlock) { //writes all the ijkb components of the elements of a structured block at certain index inside m_sblocks to an std::vector
    get_nodes_of_sb(elements, iBlock, 1);
}

  struct SB_copy_field {
    Array4D m_field_dest;
    const Array4D m_field_src;

    SB_copy_field(Array4D field_dest, const Array4D field_src) : m_field_dest(field_dest), m_field_src(field_src) {}

    void run()
    {
        unsigned sz = m_field_dest.size();
        VERIFY_OP_ON(sz, ==, (unsigned)m_field_src.size(), "bad size");
        Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,(unsigned)sz),*this);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(unsigned index) const
    { //madbrew: faster but doesn't leverage kokkos to full advantage; it ignores layout
      m_field_dest.data()[index] = m_field_src.data()[index];
    }
  };

  void BlockStructuredGrid::copy_field(typename StructuredGrid::MTField* field_dest, typename StructuredGrid::MTField* field_src)
  {
    for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
      {
        std::shared_ptr<StructuredBlock> sgrid = m_sblocks[iblock];
        if (sgrid->is_empty())
          continue;
        Array4D dest = (*field_dest->m_block_fields[iblock]);
        Array4D src = (*field_src->m_block_fields[iblock]);
        SB_copy_field cf(dest, src);
        cf.run();
      }
  }


  /// axpbypgz calculates: z = alpha*x + beta*y + gamma*z
  struct SB_nodal_field_axpbypgz {
    const double m_alpha, m_beta, m_gamma;
    const Array4D m_field_x;
    const Array4D m_field_y;
    Array4D m_field_z;

    SB_nodal_field_axpbypgz(double alp, const Array4D field_x,
                            double bet, const Array4D field_y,
                            double gam,  Array4D field_z) :
      m_alpha(alp), m_beta(bet), m_gamma(gam),
      m_field_x(field_x), m_field_y(field_y), m_field_z(field_z) {}

    void run()
    {
      unsigned sz = m_field_x.size();
      VERIFY_OP_ON(sz, ==, (unsigned)m_field_y.size(), "bad size");
      VERIFY_OP_ON(sz, ==, (unsigned)m_field_z.size(), "bad size");
      Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,sz),*this);
    }
    KOKKOS_INLINE_FUNCTION
    void operator()(const unsigned& index) const
    {
      m_field_z.data()[index] = m_alpha*m_field_x.data()[index] + m_beta*m_field_y.data()[index] + m_gamma*m_field_z.data()[index];
    }
  };

  void BlockStructuredGrid::nodal_field_axpbypgz(double alpha, typename StructuredGrid::MTField* field_x,
                                                 double beta, typename StructuredGrid::MTField* field_y,
                                                 double gamma, typename StructuredGrid::MTField* field_z)
  {
    for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
      {
        std::shared_ptr<StructuredBlock> sgrid = m_sblocks[iblock];
        if (sgrid->is_empty())
          continue;
        Array4D fx = (*field_x->m_block_fields[iblock]);
        Array4D fy = (*field_y->m_block_fields[iblock]);
        Array4D fz = (*field_z->m_block_fields[iblock]);
        SB_nodal_field_axpbypgz na(alpha, fx,
                                   beta,  fy,
                                   gamma, fz);
        na.run();
      }
  }

  /// axpby calculates: y = alpha*x + beta*y
  struct SB_nodal_field_axpby {
    const double m_alpha, m_beta;
    const Array4D m_field_x;
    Array4D m_field_y;

    SB_nodal_field_axpby(double alp, const Array4D field_x,
                         double bet, Array4D field_y) :
      m_alpha(alp), m_beta(bet),
      m_field_x(field_x), m_field_y(field_y) {}

    void run()
    {
      unsigned sz = m_field_x.size();
      VERIFY_OP_ON(sz, ==, (unsigned)m_field_y.size(), "bad size");
      Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,sz),*this);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const unsigned& index) const
    {
      m_field_y.data()[index] = m_alpha*m_field_x.data()[index] + m_beta*m_field_y.data()[index];
    }
  };

  void BlockStructuredGrid::nodal_field_axpby(double alpha, typename StructuredGrid::MTField* field_x, double beta, typename StructuredGrid::MTField* field_y)
  {
    for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
      {
        std::shared_ptr<StructuredBlock> sgrid = m_sblocks[iblock];
        if (sgrid->is_empty())
          continue;
        Array4D fx = (*field_x->m_block_fields[iblock]);
        Array4D fy = (*field_y->m_block_fields[iblock]);
        SB_nodal_field_axpby na(alpha, fx, beta, fy);
        na.run();
      }
  }


  struct SB_nodal_field_dot {
    DeviceSafeSGridBoundaryNotInterfaceSelector boundSelector;
    DetermineIfLowestRankedOwner rank_determinator;

    unsigned spatialDim;

    Double m_sum;
    const Array4D m_field_x;
    const Array4D m_field_y;
    SGridSizes m_sizes;
    Kokkos::Array<unsigned int, 3> loop_orderings;

    SB_nodal_field_dot(const Array4D field_x,
                       const Array4D field_y,
                       std::shared_ptr<StructuredBlock> sgrid) :
                           boundSelector(sgrid), rank_determinator(sgrid),
                           m_sum(0.0),
                           m_field_x(field_x), m_field_y(field_y)
    {
        m_sizes = sgrid->m_sizes;
        for(unsigned ijk=0;ijk<3;ijk++)
            loop_orderings[ijk]=sgrid->m_loop_ordering[ijk];

        spatialDim = 3;
    }

    void run()
    {
        Double tot=0;
        unsigned sz = m_sizes.node_size[0]*m_sizes.node_size[1]*m_sizes.node_size[2];
        Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0,sz),*this,tot);
        m_sum+=tot;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const unsigned& iNode, Double& loc_sum) const
    {
        Kokkos::Array<unsigned, 3> ijk;
        sgrid_multi_dim_indices_from_index_node(m_sizes,loop_orderings,iNode,ijk);
        if(    !(ijk[0] == 0 || ijk[0] == m_sizes.node_size[0]-1 ||
                ijk[1] == 0 || ijk[1] == m_sizes.node_size[1]-1 ||
                ijk[2] == 0 || ijk[2] == m_sizes.node_size[2]-1 ) //if the node is purely interior
                || boundSelector(ijk) //or if it is on a boundary and not on an interface
                || rank_determinator(ijk) ) //or if it is on an interface and this block is the lowest ranking block sharing that node
            for(unsigned iDim=0;iDim<spatialDim;iDim++)
                loc_sum += m_field_x(ijk[0],ijk[1],ijk[2],iDim)*m_field_y(ijk[0],ijk[1],ijk[2],iDim); //then dot it
    }
  };

  struct SB_parallel_count_nodes
  {
      DeviceSafeSGridBoundaryNotInterfaceSelector boundSelector;
      DetermineIfLowestRankedOwner rank_determinator;

      unsigned spatialDim;

      unsigned num_nodes;
      const Array4D m_field;
      SGridSizes m_sizes;
      Kokkos::Array<unsigned int, 3> loop_orderings;

      SB_parallel_count_nodes(const Array4D field,
                         std::shared_ptr<StructuredBlock> sgrid) :
                             boundSelector(sgrid), rank_determinator(sgrid),
                             num_nodes(0.0),
                             m_field(field)
      {
          m_sizes = sgrid->m_sizes;
          for(unsigned ijk=0;ijk<3;ijk++)
              loop_orderings[ijk]=sgrid->m_loop_ordering[ijk];

          spatialDim = 3;
      }

      void run()
      {
          unsigned tot=0;
          unsigned sz = m_sizes.node_size[0]*m_sizes.node_size[1]*m_sizes.node_size[2];//m_field_x.size();
          Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0,sz),*this,tot);
          num_nodes+=tot;
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const unsigned& iNode, unsigned& loc_sum) const
      {
          Kokkos::Array<unsigned, 3> ijk;
          sgrid_multi_dim_indices_from_index_node(m_sizes,loop_orderings,iNode,ijk);
          if(    !(ijk[0] == 0 || ijk[0] == m_sizes.node_size[0]-1 ||
                  ijk[1] == 0 || ijk[1] == m_sizes.node_size[1]-1 ||
                  ijk[2] == 0 || ijk[2] == m_sizes.node_size[2]-1 ) //if the node is purely interior
                  || boundSelector(ijk) //or if it is on a boundary and not on an interface
                  || rank_determinator(ijk) ) //or if it is on an interface and this block is the lowest ranking block sharing that node
              loc_sum++;
      }
  };

  long double BlockStructuredGrid::nodal_field_dot(typename StructuredGrid::MTField* field_x, typename StructuredGrid::MTField* field_y)
  {
    Double sum=0.0;
    for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
      {
        std::shared_ptr<StructuredBlock> sgrid = m_sblocks[iblock];
        if (sgrid->is_empty())
          continue;
        Array4D fx = (*field_x->m_block_fields[iblock]);
        Array4D fy = (*field_y->m_block_fields[iblock]);

        unsigned spatDim = fx.extent(3);

        SB_nodal_field_dot na(fx, fy, sgrid);
        na.spatialDim = spatDim;

        na.run();
        sum += na.m_sum;
      }
    return sum;
  }

  /// set field to constant value
  struct SB_nodal_field_set_value {
    const Double m_value;
    Array4D m_field_x;
    SGridSizes m_block_sizes;
   Kokkos::Array<unsigned int, 3> m_loop_orderings;
    SB_nodal_field_set_value(Double val, Array4D field_x,SGridSizes block_sizes_in,
            std::array<unsigned int, 3> loop_orderings_in) :
      m_value(val), m_field_x(field_x) {m_block_sizes=block_sizes_in;
                              m_loop_orderings[0]=loop_orderings_in[0];
                              m_loop_orderings[1]=loop_orderings_in[1];
                              m_loop_orderings[2]=loop_orderings_in[2];}

    void run()
    {
      unsigned sz =  m_field_x.size();
      Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,(unsigned)sz),*this);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const unsigned& index) const
    {
      m_field_x.data()[index] = m_value;
    }
  };


  void BlockStructuredGrid::nodal_field_set_value(typename StructuredGrid::MTField* field_x, double value )
  {
    for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
      {
        std::shared_ptr<StructuredBlock> sgrid = m_sblocks[iblock];
        if (sgrid->is_empty())
          continue;
        Array4D fx = (*field_x->m_block_fields[iblock]);
        SB_nodal_field_set_value na(value, fx,sgrid->m_sizes,sgrid->m_loop_ordering);
        na.run();
      }
  }


  void BlockStructuredGrid::comm_fields(std::vector<const typename StructuredGrid::MTField*>& fields)
  {
    if (stk::parallel_machine_rank(m_comm) == 1) // only if parallel run
      std::cout << "comm_fields not yet impl" << std::endl;
  }


  void BlockStructuredGrid::sum_fields(std::vector<const typename StructuredGrid::MTField*>& fields)
  {
    if (stk::parallel_machine_rank(m_comm) == 1) // only if parallel run
      std::cout << "sum_fields not yet impl" << std::endl;
  }

  size_t BlockStructuredGrid::parallel_count_nodes()
  {
      size_t sum = 0;
      for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
        {
          std::shared_ptr<StructuredBlock> sgrid = m_sblocks[iblock];
          if (sgrid->is_empty())
            continue;
          Array4D f = *(m_sblocks[iblock]->m_sgrid_coords_impl.get());

          SB_parallel_count_nodes counter(f,  sgrid);

          counter.run();
          sum += counter.num_nodes;
        }
      return sum;
  }

  std::shared_ptr<BlockStructuredGrid>   BlockStructuredGrid::
  fixture_1(stk::ParallelMachine comm, std::array<unsigned,3> sizes, std::array<double,3> dim_widths,std::array<double,3> dim_offsets)
  {
    std::shared_ptr<BlockStructuredGrid> bsg(new BlockStructuredGrid(comm,0));
    std::shared_ptr<StructuredBlock> sbi = StructuredBlock::fixture_1(comm, sizes, 0, 0, 0, bsg.get(), dim_widths ,dim_offsets);
    bsg->m_sblocks.push_back(sbi);
    // FIXME - make this automatic, always have coordinates registered
    bsg->register_field("coordinates",3);
    return bsg;
  }

}//percept

