#include "mpi.h"
#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/balanceMtoN.hpp>
#include <stk_balance/internal/Inputs.hpp>
#include <stk_balance/internal/M2NDecomposer.hpp>
#include <stk_balance/setup/M2NParser.hpp>

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

#include <stk_util/command_line/CommandLineParser.hpp>
#include <stk_util/command_line/CommandLineParserParallel.hpp>
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <stk_util/environment/Env.hpp>
#include <stk_util/environment/EnvData.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/environment/FileUtils.hpp>
#include <stk_util/util/string_utils.hpp>
#include <stk_util/util/human_bytes.hpp>

#include <string>
#include <iostream>
#include <fstream>
#include <memory>

namespace
{

void set_output_streams(MPI_Comm comm)
{
  if (stk::parallel_machine_rank(comm) != 0) {
    stk::EnvData::instance().m_outputP0 = &stk::EnvData::instance().m_outputNull;
  }
  Ioss::Utils::set_output_stream(sierra::Env::outputP0());
}

void rebalance_m_to_n(stk::balance::M2NParsedOptions &parsedOptions, MPI_Comm comm)
{
    stk::mesh::MetaData meta;
    stk::mesh::BulkData bulk(meta, comm);

    stk::mesh::Field<double> &field = meta.declare_field<stk::mesh::Field<double> >(stk::topology::ELEMENT_RANK, "TargetDecomp", 1);
    stk::mesh::put_field_on_mesh(field, meta.universal_part(), (double*)nullptr);

    stk::io::StkMeshIoBroker ioBroker;
    stk::io::fill_mesh_preexisting(ioBroker, parsedOptions.inFile, bulk);

    stk::balance::internal::rebalanceMtoN(ioBroker, field, parsedOptions);
}

}

int main(int argc, const char**argv)
{
    MPI_Init(&argc, const_cast<char***>(&argv));
    MPI_Comm comm = MPI_COMM_WORLD;

    stk::balance::M2NParser parser(comm);
    stk::balance::M2NParsedOptions parsedOptions;
    parser.parse_command_line_options(argc, argv, parsedOptions);

    set_output_streams(comm);
    rebalance_m_to_n(parsedOptions, comm);

    size_t hwmMax = 0, hwmMin = 0, hwmAvg = 0;
    stk::get_memory_high_water_mark_across_processors(comm, hwmMax, hwmMin, hwmAvg);
    sierra::Env::outputP0() << "Memory HWM across procs, max/min/avg: "
            << stk::human_bytes(hwmMax) << " / " << stk::human_bytes(hwmMin) << " / "<< stk::human_bytes(hwmAvg) << std::endl;
    MPI_Finalize();
    return 0;
}
