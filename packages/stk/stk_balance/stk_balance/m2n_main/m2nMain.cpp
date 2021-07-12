#include <stk_util/parallel/Parallel.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/m2n/m2nRebalance.hpp>
#include <stk_balance/setup/M2NParser.hpp>

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <stk_util/environment/Env.hpp>
#include <stk_util/environment/EnvData.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/util/human_bytes.hpp>

namespace
{

void set_output_streams(MPI_Comm comm)
{
  if (stk::parallel_machine_rank(comm) != 0) {
    stk::EnvData::instance().m_outputP0 = &stk::EnvData::instance().m_outputNull;
  }
  Ioss::Utils::set_output_stream(sierra::Env::outputP0());
}

void rebalance_m2n(stk::balance::M2NBalanceSettings &balanceSettings, MPI_Comm comm)
{
  stk::mesh::MetaData meta;
  stk::mesh::BulkData bulk(meta, comm);
  stk::io::StkMeshIoBroker ioBroker;
  stk::io::fill_mesh_preexisting(ioBroker, balanceSettings.get_input_filename(), bulk);

  stk::balance::m2n::m2nRebalance(ioBroker, balanceSettings);
}

}

int main(int argc, const char**argv)
{
    MPI_Init(&argc, const_cast<char***>(&argv));
    MPI_Comm comm = MPI_COMM_WORLD;

    stk::balance::M2NBalanceSettings balanceSettings;

    stk::balance::M2NParser parser(comm);
    parser.parse_command_line_options(argc, argv, balanceSettings);

    set_output_streams(comm);
    rebalance_m2n(balanceSettings, comm);

    size_t hwmMax = 0, hwmMin = 0, hwmAvg = 0;
    stk::get_memory_high_water_mark_across_processors(comm, hwmMax, hwmMin, hwmAvg);
    sierra::Env::outputP0() << "Memory HWM across procs, max/min/avg: "
                            << stk::human_bytes(hwmMax) << " / "
                            << stk::human_bytes(hwmMin) << " / "
                            << stk::human_bytes(hwmAvg) << std::endl;
    MPI_Finalize();
    return 0;
}
