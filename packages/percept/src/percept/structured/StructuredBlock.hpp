// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_StructuredBlock_hpp
#define percept_StructuredBlock_hpp

#include <percept/Percept.hpp>
#if !STK_PERCEPT_LITE

#include <Ioss_StructuredBlock.h>

#if HAVE_CGNS
#include <cgnslib.h>
#endif

#include <percept/function/MDArray.hpp>
#include <percept/MeshType.hpp>

#include <array>
#include <unordered_map>

#include <stdio.h>

namespace percept {

  class BlockStructuredGrid;

  class StructuredBlock {
  public:
    //typedef MDArray Array4D;
    using Array4D = StructuredGrid::MTField::Array4D;
    stk::ParallelMachine m_comm;

  protected:

    unsigned m_iblock;
    Ioss::StructuredBlock *m_sblock;

  public:

    BlockStructuredGrid *m_parent;

    std::string m_name;
    int m_base;
    int m_zone;

    std::shared_ptr<Array4D> m_sgrid_coords_impl;
    Array4D& m_sgrid_coords;
    std::map<std::string, std::shared_ptr<Array4D> > m_fields;

    SGridSizes m_sizes;
    std::array<unsigned,3> m_loop_ordering; // loop ordering
    std::array<unsigned,3> m_access_ordering; // access ordering

    // boundary condition info - either (ilo,jlo,klo),(ihi,jhi,khi) or specify each point (i,j,k)
    std::vector<Ioss::BoundaryCondition> m_boundaryConditions;
    //enum FaceDir { FACE_I_MIN, FACE_I_MAX, FACE_J_MIN, FACE_J_MAX, FACE_K_MIN, FACE_K_MAX };
#if HAVE_CGNS
    std::vector<CG_BCType_t> m_bc_types;
#endif

    std::vector<Ioss::ZoneConnectivity> m_zoneConnectivity;

  public:

    // create a StructuredBlock from an Ioss CGNS-based Ioss::StructuredBlock
    StructuredBlock(stk::ParallelMachine comm, unsigned iblock, Ioss::StructuredBlock *sb, BlockStructuredGrid *bsg);

    // create a StructuredBlock with the given sizes @param njik (this mirrors the data in CGNS - nijk[0..2] give the node sizes in i,j,k, other entries ignored
    StructuredBlock(stk::ParallelMachine comm, unsigned iblock, const std::array<unsigned,9>& nijk, const std::array<unsigned,3>& node_size_global, const std::string& name, int base, int zone, BlockStructuredGrid *bsg);

    // can be "empty" if the block is not on this processor
    bool is_empty();

    // register a nodal field of name @param fieldName, with @param fourth_dim components
    std::shared_ptr<Array4D> register_field(const std::string& fieldName, unsigned fourth_dim);

    // create a 1x1x1 cube of given sizes
    static std::shared_ptr<StructuredBlock>
    fixture_1(stk::ParallelMachine comm, std::array<unsigned,3> sizes=std::array<unsigned,3>{{2u,2u,2u}}, unsigned iblock = 0, int base = 0, int zone = 0, BlockStructuredGrid *bsg = 0,
            std::array<double,3> dim_widths=std::array<double,3>{{1.0,1.0,1.0}},std::array<double,3> dim_offsets=std::array<double,3>{{0.0,0.0,0.0}});

    // return the ordinal of the given i,j,k in the list of all nodes
    uint64_t local_offset(uint64_t i, uint64_t j, uint64_t k, unsigned *sizes);

    // return the i,j,k given the local offset
    void multi_dim_indices_from_local_offset(uint64_t local_offset, std::array<unsigned,3>& indx);

    // extract data from the CGNS file if constructed using an Ioss::StructuredBlock
    void read_cgns();

    void print(std::ostream& out = std::cout, int level=0);

    // create a VTK file representation of this grid
    void dump_vtk(const std::string& file_prefix);

    unsigned get_block_id() const
    {return m_iblock;}

  };



}


#endif
#endif
