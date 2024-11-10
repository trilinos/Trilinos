/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <KrinoMeshAdapt.hpp>
#include <KrinoMeshAdaptAlgorithmParameters.hpp>
#include <Akri_Refinement.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/diag/Timer.hpp>
#include <stk_util/environment/Env.hpp>
#include "Akri_TransitionElementEdgeMarker.hpp"

namespace krino {

void refine_mesh_with_params(stk::mesh::BulkData & mesh, const MeshAdaptAlgorithmParameters &algParams, const stk::ParallelMachine comm)
{
  stk::diag::Timer refinementTimer("Refinement", sierra::Diag::sierraTimer());
  stk::mesh::Part * activePart = nullptr;
  Refinement refinement(mesh.mesh_meta_data(), activePart, algParams.force64Bit, algParams.assert32Bit, refinementTimer);

  std::vector<size_t> counts;
  stk::mesh::comm_mesh_counts(mesh, counts);

  sierra::Env::outputP0() << "Uniform refinement: before refine, mesh has  " << counts[0] << " nodes, " << counts[1]
           << " edges, " << counts[2] << " faces, " << counts[3] << " elements" << std::endl;

  refinement.do_uniform_refinement(algParams.numUniformRefinementLevels);
  refinement.delete_parent_elements();

  stk::mesh::comm_mesh_counts(mesh, counts);

  sierra::Env::outputP0() << "Uniform refinement: after refine, mesh has  " << counts[0] << " nodes, " << counts[1]
           << " edges, " << counts[2] << " faces, " << counts[3] << " elements" << std::endl;
}

}

