#include <gtest/gtest.h>
#include <string>
#include <mpi.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/IossBridge.hpp>
#include <Ioss_SubSystem.h>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Types.hpp>
#include <fieldNameTestUtils.hpp>
#include <restartTestUtils.hpp>

namespace {

TEST(StkMeshIoBrokerHowTo, writeTwoResultFiles)
{
    const std::string resultsFilename1 = "output1.results";
    const std::string resultsFilename2 = "output2.results";
    const std::string displacementFieldName = "displacement";
    const std::string velocityFieldName = "velocity";

    MPI_Comm communicator = MPI_COMM_WORLD;
    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    setupMeshAndFieldsForTest(stkMeshIoBroker, displacementFieldName, velocityFieldName);
    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();

    // Now that input mesh with field data have been set, start output

    stk::mesh::FieldBase *displacementField = stkMeshMetaData.get_field(displacementFieldName);

    // For file one, set up results and global variables
    size_t file1Handle = stkMeshIoBroker.create_output_mesh(resultsFilename1);
    stkMeshIoBroker.add_results_field(file1Handle, *displacementField);
    std::string globalVarNameFile1 = "eigenValue";
    stkMeshIoBroker.add_global(file1Handle, globalVarNameFile1, Ioss::Field::REAL);

    // For file two, set up results and global variables
    size_t file2Handle = stkMeshIoBroker.create_output_mesh(resultsFilename2);
    std::string nameOnOutputFile("deformations");
    stkMeshIoBroker.add_results_field_with_alternate_name(file2Handle, *displacementField, nameOnOutputFile);
    stk::mesh::FieldBase *velocityField = stkMeshMetaData.get_field(velocityFieldName);
    stkMeshIoBroker.add_results_field(file2Handle, *velocityField);
    std::string globalVarNameFile2 = "kineticEnergy";
    stkMeshIoBroker.add_global(file2Handle, globalVarNameFile2, Ioss::Field::REAL);

    // Write output
    double time = 0.0;
    stkMeshIoBroker.begin_output_step(time, file1Handle);
    stkMeshIoBroker.process_output_request(file1Handle);
    const double globalVarValue1 = 13.0;
    stkMeshIoBroker.write_global(file1Handle, globalVarNameFile1, globalVarValue1);
    stkMeshIoBroker.end_output_step(file1Handle);

    stkMeshIoBroker.begin_output_step(time, file2Handle);
    stkMeshIoBroker.process_output_request(file2Handle);
    const double globalVarValue2 = 14.0;
    stkMeshIoBroker.write_global(file2Handle, globalVarNameFile2, globalVarValue2);
    stkMeshIoBroker.end_output_step(file2Handle);

    // Code below is for testing that above code is working
    std::vector<std::string> goldNodalVarNamesInFile1;
    goldNodalVarNamesInFile1.push_back(displacementFieldName);
    checkFileForNodalVarNames(resultsFilename1, goldNodalVarNamesInFile1);
    checkFileForGlobal(resultsFilename1, globalVarNameFile1, globalVarValue1);

    std::vector<std::string> goldNodalVarNamesInFile2;
    goldNodalVarNamesInFile2.push_back(nameOnOutputFile);
    goldNodalVarNamesInFile2.push_back(velocityFieldName);
    checkFileForNodalVarNames(resultsFilename2, goldNodalVarNamesInFile2);
    checkFileForGlobal(resultsFilename2, globalVarNameFile2, globalVarValue2);

    unlink(resultsFilename1.c_str());
    unlink(resultsFilename2.c_str());
}

}
