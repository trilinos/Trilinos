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
#include <stk_util/environment/FileUtils.hpp>
#include <stk_util/util/string_utils.hpp>

#include <string>
#include <iostream>
#include <fstream>
#include <memory>

namespace
{

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

    rebalance_m_to_n(parsedOptions, comm);

    MPI_Finalize();
    return 0;
}
