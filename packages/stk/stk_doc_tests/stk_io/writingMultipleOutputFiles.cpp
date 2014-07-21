#include <gtest/gtest.h>
#include <string>
#include <stk_io/StkMeshIoBroker.hpp>
#include <Ioss_SubSystem.h>
#include <stk_mesh/base/MetaData.hpp>
#include <restartTestUtils.hpp>

namespace {

  TEST(StkMeshIoBrokerHowTo, writingMultipleOutputFiles)
  {
    // ============================================================
    //+ INITIALIZATION
    const std::string resultsFilename1 = "output1.results";
    const std::string resultsFilename2 = "output2.results";
    const std::string displacementFieldName = "displacement";
    const std::string velocityFieldName = "velocity";

    MPI_Comm communicator = MPI_COMM_WORLD;
    stk::io::StkMeshIoBroker stkIo(communicator);
    setupMeshAndFieldsForTest(stkIo, displacementFieldName, velocityFieldName);
    stk::mesh::MetaData &meta_data = stkIo.meta_data();

    //-BEGIN
    // ============================================================
    //+ EXAMPLE -- Two results output files
    stk::mesh::FieldBase *displacementField =
      meta_data.get_field(stk::topology::NODE_RANK, displacementFieldName);

    //+ For file one, set up results and global variables
    size_t file1Handle = stkIo.create_output_mesh(resultsFilename1,
						  stk::io::WRITE_RESULTS);
    stkIo.add_field(file1Handle, *displacementField);
    std::string globalVarNameFile1 = "eigenValue";
    stkIo.add_global(file1Handle, globalVarNameFile1, Ioss::Field::REAL);

    //+ For file two, set up results and global variables
    size_t file2Handle = stkIo.create_output_mesh(resultsFilename2,
						  stk::io::WRITE_RESULTS);
    std::string nameOnOutputFile("deformations");
    stkIo.add_field(file2Handle, *displacementField, nameOnOutputFile);
    stk::mesh::FieldBase *velocityField = meta_data.get_field(stk::topology::NODE_RANK, velocityFieldName);
    stkIo.add_field(file2Handle, *velocityField);
    std::string globalVarNameFile2 = "kineticEnergy";
    stkIo.add_global(file2Handle, globalVarNameFile2, Ioss::Field::REAL);

    //+ Write output
    double time = 0.0;
    stkIo.begin_output_step(file1Handle, time);
    stkIo.write_defined_output_fields(file1Handle);
    const double globalVarValue1 = 13.0;
    stkIo.write_global(file1Handle, globalVarNameFile1, globalVarValue1);
    stkIo.end_output_step(file1Handle);

    stkIo.begin_output_step(file2Handle, time);
    stkIo.write_defined_output_fields(file2Handle);
    const double globalVarValue2 = 14.0;
    stkIo.write_global(file2Handle, globalVarNameFile2, globalVarValue2);
    stkIo.end_output_step(file2Handle);

    //-END
    // ============================================================
    //+ Validation
    std::vector<std::string> goldNodalVarNamesInFile1;
    goldNodalVarNamesInFile1.push_back(displacementFieldName);
    checkFileForNodalVarNames(resultsFilename1, goldNodalVarNamesInFile1);
    checkFileForGlobal(resultsFilename1, globalVarNameFile1, globalVarValue1);

    std::vector<std::string> goldNodalVarNamesInFile2;
    goldNodalVarNamesInFile2.push_back(nameOnOutputFile);
    goldNodalVarNamesInFile2.push_back(velocityFieldName);
    checkFileForNodalVarNames(resultsFilename2, goldNodalVarNamesInFile2);
    checkFileForGlobal(resultsFilename2, globalVarNameFile2, globalVarValue2);

    // Cleanup
    unlink(resultsFilename1.c_str());
    unlink(resultsFilename2.c_str());
  }

}
