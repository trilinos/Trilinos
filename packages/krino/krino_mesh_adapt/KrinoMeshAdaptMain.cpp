/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <Ioss_PropertyManager.h>
#include <KrinoMeshAdapt.hpp>
#include <mpi.h>                           // for MPI_Comm_rank, MPI_Comm_size
#include <KrinoMeshAdaptInputData.hpp>
#include <KrinoMeshAdaptParser.hpp>
#include <iostream>                        // for operator<<, basic_ostream:...
#include <memory>                          // for unique_ptr
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_finalize
#include <stk_io/FillMesh.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_util/environment/Env.hpp>
#include <stk_util/environment/EnvData.hpp>
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/ProgramOptions.hpp>
#include <stk_util/registry/ProductRegistry.hpp>  // for ProductRegistry

void krino_mesh_adapt(const krino::MeshAdaptInputData &inputData, const stk::ParallelMachine comm)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(comm).create();
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();
  meta.enable_late_fields();

  stk::io::fill_mesh_with_auto_decomp(inputData.meshIn, *bulk);

  krino::refine_mesh_with_params(*bulk, inputData.algorithmParams, bulk->parallel());

  Ioss::PropertyManager properties;
  if (inputData.algorithmParams.autoCompose)
    properties.add(Ioss::Property("COMPOSE_RESULTS", 1));

  stk::io::StkMeshIoBroker io_broker;
  io_broker.set_bulk_data(bulk);
  size_t resultFileIndex = io_broker.create_output_mesh(inputData.meshOut, stk::io::WRITE_RESULTS, properties);
  io_broker.write_output_mesh(resultFileIndex);
}

bool run_krino_mesh_adapt(int argc,  char **argv, stk::ParallelMachine comm)
{
    const double startWallTime = stk::wall_time();
    const double startCpuTime = stk::cpu_time();

    int proc = stk::parallel_machine_rank(comm);
    if(proc != 0)
        stk::EnvData::instance().m_outputP0 = &stk::EnvData::instance().m_outputNull;

    krino::MeshAdaptParser parser;
    bool didParseOk = parser.read_command_line(argc, argv, comm);

    bool successfullyMeshed{false};
    double elapsedWallTime{0.0};
    double elapsedCpuTime{0.0};

    if (didParseOk)
    {
      const krino::MeshAdaptInputData &inputData = parser.get_input_data();

      krino_mesh_adapt(inputData, comm);

      stk::parallel_machine_barrier(comm);

      elapsedWallTime = stk::wall_time() - startWallTime;
      elapsedCpuTime = stk::cpu_time() - startCpuTime;

      if (0 == stk::parallel_machine_rank(comm))
          std::cout << "Elapsed wall time: " + std::to_string(elapsedWallTime) + " seconds" << std::endl;

      if (0 == stk::parallel_machine_rank(comm))
          std::cout << "Elapsed cpu time: " + std::to_string(elapsedCpuTime) + " seconds" << std::endl;

      successfullyMeshed = true;
    }

//    const double startTime{static_cast<double>(::time(nullptr))};
//    common_utils::write_audit_file(
//                            "krino_mesh_adapt",
//                            "mesh refinement",
//                            startTime,
//                            elapsedWallTime,
//                            elapsedCpuTime,
//                            successfullyMeshed,
//                            [&](const auditdata &data) {OutputAuditLog(&data);},
//                            comm);
    return successfullyMeshed;
}

int main(int argc,  char **argv)
{
    stk::ParallelMachine comm{stk::parallel_machine_init(&argc, &argv)};

    bool successfullyMeshed{false};

    try
    {
        successfullyMeshed = run_krino_mesh_adapt(argc, argv, comm);
    }
    catch(std::exception &e)
    {
        std::cerr << "Proc " << stk::parallel_machine_rank(comm) << ": " << e.what() << std::endl;
        const int errorCode{-1};
        MPI_Abort(comm, errorCode);
    }

    stk::parallel_machine_finalize();

    return successfullyMeshed ? 0 : 1;
}
