#include <mpi.h>
#include <stk_balance/balance.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/Inputs.hpp>
#include <stk_balance/internal/LastStepFieldWriter.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

#include <stk_util/parallel/ParallelReduceBool.hpp>

#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>

#include <string>
#include <iostream>
#include <fstream>

bool create_path(const std::string &filename);
bool am_proc0();
void print_usage_msg(const std::string &programName);
void print_file_not_found_msg(const std::string &filename);
void print_running_msg(const stk::balance::Inputs& parsedInputs);
void print_directory_does_not_exist_msg(const std::string& directoryName);
void run_stk_rebalance(const stk::balance::Inputs& parsedInputs, MPI_Comm comm);

enum { OK=0, NOT_OK = 1 };

int handle_usage_info(const std::string& exe)
{
    print_usage_msg(exe);
    return NOT_OK;
}

int exodus_file_does_not_exist(const std::string& exe, const std::string& filename)
{
    print_usage_msg(exe);
    print_file_not_found_msg(filename);
    return NOT_OK;
}

int handle_output_directory(const std::string& outputDirectory)
{
    bool path_ok = true;
    if(am_proc0())
        path_ok = stk::balance::create_path(outputDirectory);

    if(!stk::is_true_on_all_procs(MPI_COMM_WORLD, path_ok))
    {
        print_directory_does_not_exist_msg(outputDirectory);
        return NOT_OK;
    }
    return OK;
}

bool are_inputs_valid(stk::balance::Inputs& parsedInputs)
{
    int returnCode = OK;

    if(stk::balance::should_write_usage_info(parsedInputs.get_exodus_filename()))
    {
        returnCode = handle_usage_info(parsedInputs.get_executable_name());
    }
    else
    {
        if(!stk::balance::does_file_exist(parsedInputs.get_exodus_filename()))
        {
            returnCode = exodus_file_does_not_exist(parsedInputs.get_executable_name(), parsedInputs.get_exodus_filename());
        }
        else
        {
            returnCode = handle_output_directory(parsedInputs.get_output_directory());
        }
    }
    return returnCode == OK;
}


int main(int argc, const char**argv)
{
    MPI_Init(&argc, const_cast<char***>(&argv));
    MPI_Comm comm = MPI_COMM_WORLD;

    stk::balance::Inputs parsedInputs(argc, argv);

    int returnCode = NOT_OK;

    if(are_inputs_valid(parsedInputs))
    {
        print_running_msg(parsedInputs);
        run_stk_rebalance(parsedInputs, comm);
        returnCode = OK;
    }

    MPI_Finalize();
    return returnCode;
}

bool am_proc0()
{
    return stk::parallel_machine_rank(MPI_COMM_WORLD)==0;
}

void print_usage_msg(const std::string &programName)
{
    if(am_proc0())
    {
        std::cerr << "\n";
        std::cerr << "\tUsage:\n";
        std::cerr << "\t\t" << programName << " exodus_file [output_directory]\n";
        std::cerr << "\n";
    }
}

void print_file_not_found_msg(const std::string &filename)
{
    if(am_proc0())
    {
        std::cerr << "\n";
        std::cerr << "\tError:\n";
        std::cerr << "\t\tFile " << filename << " does not exist. Exiting.\n";
        std::cerr << "\n";
    }
}

void print_running_msg(const stk::balance::Inputs& parsedInputs)
{
    if(am_proc0())
    {
        std::cerr << "\n";
        std::cerr << "\tRunning: " << parsedInputs.get_executable_name() << " " << parsedInputs.get_exodus_filename() << " " << parsedInputs.get_output_directory() << std::endl;
        std::cerr << "\n";
    }
}

void print_directory_does_not_exist_msg(const std::string& directoryName)
{
    if(am_proc0())
        std::cerr << "Could not create or find directory: " <<  directoryName<< ". Exiting." << std::endl;
}



void run_stk_rebalance(const stk::balance::Inputs& parsedInputs, MPI_Comm comm)
{
    stk::mesh::MetaData meta;
    stk::mesh::BulkData bulk(meta, comm);

    stk::balance::internal::LastStepFieldWriterAutoDecomp fieldWriter(bulk, parsedInputs.get_exodus_filename());

    stk::balance::GraphCreationSettings graphOptions;
    stk::balance::balanceStkMesh(graphOptions, bulk);

    std::string outputFilename = parsedInputs.get_output_directory() + "/" + parsedInputs.get_exodus_filename();
    fieldWriter.write_mesh(outputFilename);
}



