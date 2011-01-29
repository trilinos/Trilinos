#ifndef stk_percept_Verifier_hpp
#define stk_percept_Verifier_hpp

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
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/fem/FieldTraits.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/Comm.hpp>

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/Stencils.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

#include <stk_percept/TopologyVerifier.hpp>
#include <stk_percept/GeometryVerifier.hpp>
#include <stk_percept/PerceptMesh.hpp>

#include <stk_percept/RunEnvironment.hpp>


namespace stk
{
  namespace percept
  {

    class Verifier
    {
    public:
      Verifier() {}
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
}//namespace stk

#endif
