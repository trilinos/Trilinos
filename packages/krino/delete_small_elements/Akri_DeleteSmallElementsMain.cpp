// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_MeshHelpers.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_util/command_line/CommandLineParser.hpp>
#include <stk_util/command_line/CommandLineParserParallel.hpp>
#include <stk_util/diag/WriterRegistry.hpp>
#include <stk_util/environment/Env.hpp>
#include <stk_util/environment/EnvData.hpp>

namespace krino {

struct DeleteSmallElementsInputData
{
    static constexpr const char * mDefaultString{"None"};
    std::string meshIn{mDefaultString};
    std::string meshOut{mDefaultString};
    bool blockNameSpecified{false};
    std::string blockName{mDefaultString};
    double minNodalVolume{0.0};
    bool minNodalVolumeSpecified{false};
    double minRelativeNodalVolume{-1.0};
    bool minRelativeNodalVolumeSpecified{false};
};

bool read_command_line( int argc, char *argv[], DeleteSmallElementsInputData & inputData, stk::ParallelMachine comm)
{
    DeleteSmallElementsInputData defaultValues;

    const stk::CommandLineOption inmeshOption{"inmesh", "i", "Filename of input genesis mesh."};
    const stk::CommandLineOption outmeshOption{"outmesh", "o", "Filename of output genesis mesh to write."};

    const stk::CommandLineOption minNodalVolumeOption{"min_nodal_volume", "m", "Remove all nodes and attached entities smaller than this size. Either --min_nodal_volume or --min_relative_nodal_volume can be specified, but not both."};
    const stk::CommandLineOption minRelativeNodalVolumeOption{"min_relative_nodal_volume", "r", "Remove all all nodes and attached entities smaller than this factor times the maximum element volume in the input mesh. Either --min_nodal_volume or --min_relative_nodal_volume can be specified, but not both"};
    const stk::CommandLineOption blockNameOption{"block", "b", "Confine consideration to the specified block."};

    stk::CommandLineParserParallel commandLine(comm);
    commandLine.disallow_unrecognized();

    commandLine.add_required<std::string>(inmeshOption);
    commandLine.add_required<std::string>(outmeshOption);
    commandLine.add_optional<std::string>(blockNameOption, inputData.blockName);
    commandLine.add_optional<double>(minNodalVolumeOption, inputData.minNodalVolume);
    commandLine.add_optional<double>(minRelativeNodalVolumeOption, inputData.minRelativeNodalVolume);


    stk::CommandLineParser::ParseState state = commandLine.parse(argc, const_cast<const char**>(argv));
    if(state == stk::CommandLineParser::ParseComplete)
    {
        inputData.meshIn = commandLine.get_option_value<std::string>(inmeshOption.name);
        inputData.meshOut = commandLine.get_option_value<std::string>(outmeshOption.name);

        inputData.blockNameSpecified = commandLine.is_option_parsed(blockNameOption.name);
        inputData.minNodalVolumeSpecified = commandLine.is_option_parsed(minNodalVolumeOption.name);
        inputData.minRelativeNodalVolumeSpecified = commandLine.is_option_parsed(minRelativeNodalVolumeOption.name);

        if(inputData.blockNameSpecified)
        {
            inputData.blockName = commandLine.get_option_value<std::string>(blockNameOption.name);
        }

        if(inputData.minNodalVolumeSpecified && inputData.minRelativeNodalVolumeSpecified)
        {
            sierra::Env::outputP0() << "ERROR: You cannot specify both --"+minNodalVolumeOption.name+" and --"+minRelativeNodalVolumeOption.name+"." << std::endl;
            state = stk::CommandLineParser::ParseError;
        }

        if(inputData.minNodalVolumeSpecified)
        {
            inputData.minNodalVolume = commandLine.get_option_value<double>(minNodalVolumeOption.name);
            if (inputData.minNodalVolume < 0.0)
            {
                sierra::Env::outputP0() << "ERROR: size specified via --"+minNodalVolumeOption.name+" must be >= 0." << std::endl;
                state = stk::CommandLineParser::ParseError;
            }
        }

        if(inputData.minRelativeNodalVolumeSpecified)
        {
            inputData.minRelativeNodalVolume = commandLine.get_option_value<double>(minRelativeNodalVolumeOption.name);
            if (inputData.minRelativeNodalVolume < 0.0)
            {
                sierra::Env::outputP0() << "ERROR: edge length ratio specified via --"+minRelativeNodalVolumeOption.name+" must be >= 0." << std::endl;
                state = stk::CommandLineParser::ParseError;
            }
        }
    }

    if(state != stk::CommandLineParser::ParseComplete)
        sierra::Env::outputP0() << commandLine.get_usage() << std::endl;

    return state == stk::CommandLineParser::ParseComplete;
}

static bool delete_small_elements(const DeleteSmallElementsInputData& inputData,
                const stk::ParallelMachine comm)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(comm).create();
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  stk::io::fill_mesh_with_auto_decomp(inputData.meshIn, *bulk);

  // delete infinitesimal elements
  double minEdgeLength, maxEdgeLength, minElementVolume, maxElementVolume;
  compute_element_quality(*bulk, minEdgeLength, maxEdgeLength, minElementVolume, maxElementVolume);
  sierra::Env::outputP0() << "Overall mesh size results: minEdgeLength=" << minEdgeLength << ", maxEdgeLength=" << maxEdgeLength << ", minElementVolume=" << minElementVolume << ", maxElementVolume=" << maxElementVolume << std::endl;

  stk::mesh::Selector blockSelector = meta.universal_part();
  if (inputData.blockNameSpecified)
  {
    stk::mesh::Part * block = meta.get_part(inputData.blockName);
    if (!block)
    {
      sierra::Env::outputP0() << "Did not find part with name " << inputData.blockName << "." << std::endl;
      return false;
    }
    blockSelector = *block;
  }

  const double minRetainedElementVolume = inputData.minNodalVolumeSpecified ? inputData.minNodalVolume : (inputData.minRelativeNodalVolume*maxElementVolume);
  if (minElementVolume < minRetainedElementVolume)
    delete_all_entities_using_nodes_with_nodal_volume_below_threshold(*bulk, blockSelector, minRetainedElementVolume);
  else
    sierra::Env::outputP0() << "All nodes already have nodal volume larger than " << minRetainedElementVolume << "." << std::endl;

  stk::io::StkMeshIoBroker io_broker;
  io_broker.set_bulk_data(bulk);
  size_t resultFileIndex = io_broker.create_output_mesh(inputData.meshOut, stk::io::WRITE_RESULTS);
  io_broker.write_output_mesh(resultFileIndex);

  return true;
}

static bool run_delete_small_elements(int argc,  char **argv, stk::ParallelMachine comm)
{
    int proc = stk::parallel_machine_rank(comm);
    if(proc != 0)
        stk::EnvData::instance().m_outputP0 = &stk::EnvData::instance().m_outputNull;

    DeleteSmallElementsInputData inputData;
    bool didParseOk = read_command_line(argc, argv, inputData, comm);

    bool successfullyRun{false};

    if (didParseOk)
    {
      successfullyRun = delete_small_elements(inputData, comm);

      stk::parallel_machine_barrier(comm);
    }

    return successfullyRun;
}

}

int main(int argc,  char **argv)
{
    stk::ParallelMachine comm{stk::parallel_machine_init(&argc, &argv)};

    bool successfullyRun{false};

    try
    {
        successfullyRun = krino::run_delete_small_elements(argc, argv, comm);
    }
    catch(std::exception &e)
    {
        std::cerr << "Proc " << stk::parallel_machine_rank(comm) << ": " << e.what() << std::endl;
        const int errorCode{-1};
        MPI_Abort(comm, errorCode);
    }

    stk::parallel_machine_finalize();

    return successfullyRun ? 0 : 1;
}
