// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_StructuredGridRefiner_hpp
#define percept_StructuredGridRefiner_hpp

#include <percept/structured/BlockStructuredGrid.hpp>
#if HAVE_CGNS
#  if defined(STK_BUILT_IN_SIERRA)
#    include <cgns/Iocgns_DatabaseIO.h>
#  else
#    include <Iocgns_DatabaseIO.h>
#  endif
#endif

namespace percept {

  class StructuredGridRefiner {

    // allowing for future 1-based indexing
    const unsigned m_index_base = 1;//0;

    unsigned do_refine_structured();

  public:
    std::shared_ptr<BlockStructuredGrid> m_input, m_output;
    int m_debug;

    StructuredGridRefiner(std::shared_ptr<BlockStructuredGrid> input, int debug=0) :
      m_debug(debug)
    {
      m_input = input;
      m_output.reset(new BlockStructuredGrid(input->m_comm, 0));
    }

    unsigned do_refine();

    void print(std::ostream& out = std::cout, int level=0);

    // tmp - shows next steps for pulling in sideset info, etc.
    void double_structured_block(Ioss::Region *new_region,
                                 Ioss::StructuredBlock *oblock, size_t& num_node, size_t& num_cell);
  };


}

#include <percept/structured/StructuredGridRefinerDef.hpp>


#endif
