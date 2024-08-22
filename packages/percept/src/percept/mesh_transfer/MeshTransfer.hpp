// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_MeshTransfer_hpp
#define percept_MeshTransfer_hpp

#include <stdexcept>
#include <sstream>
#include <vector>
#include <iostream>

#include <stk_mesh/base/Field.hpp>

#include "Teuchos_CommandLineProcessor.hpp"

namespace stk
{
  namespace io
  {
    class StkMeshIoBroker;
  }
}

  namespace percept
  {

    class MeshTransfer
    {
    public:
      MeshTransfer(stk::ParallelMachine comm_in) : 
    comm(comm_in),
    src_mesh(),
	dst_mesh(),
	target_mesh(),
	dst_entity(),
	coarse_search_expansion_factor(1.5),
	field_name(),
	thvec_name(),
	rzvec_name(),
	dst_field_name(),
	xrot(0),
	yrot(0),
	zrot(0),
	xtrans(0),
	ytrans(0),
	ztrans(0),
	clp(false)
	{}

      void run(int argc,  char** argv);

    private:
      void process_options();

      stk::ParallelMachine comm;
      std::string src_mesh, dst_mesh, target_mesh;
      std::string dst_entity;
      double coarse_search_expansion_factor;
      std::string field_name, thvec_name, rzvec_name, dst_field_name;

      double xrot, yrot, zrot;
      double xtrans, ytrans, ztrans;

      Teuchos::CommandLineProcessor clp;
    };
      
  }//namespace percept

#endif
