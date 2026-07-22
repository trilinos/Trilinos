// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_Verifier_hpp
#define percept_Verifier_hpp

#include <stdexcept>
#include <sstream>
#include <vector>
#include <iostream>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/BroadcastArg.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/Comm.hpp>


#include <percept/TopologyVerifier.hpp>
#include <percept/GeometryVerifier.hpp>
#include <percept/PerceptMesh.hpp>

#include <percept/RunEnvironment.hpp>


  namespace percept
  {

    class Verifier
    {
    public:
      Verifier() : printTable_opt(1), fullPrint_opt(0) {}
      void verify(int argc,  char** argv);

      void process_options(RunEnvironment& re);
      // command line
      std::string mesh_opt;
      int printTable_opt;
      int fullPrint_opt;
    };

    class MeshVerifier : public Verifier
    {
    };

    class CodeVerifier : public Verifier
    {
    };

    class MeshGeometryVerifier : public MeshVerifier
    {

    };

    class MeshTopologyVerifier : public MeshVerifier
    {
      
    };

      
  }//namespace percept

#endif
