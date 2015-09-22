/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/


#include <fei_sstream.hpp>
#include <snl_fei_Utils.hpp>
#include <test_utils/test_FEI.hpp>

#include <test_utils/fei_test_utils.hpp>
#include <test_utils/LibraryFactory.hpp>

#include <fei_ParameterSet.hpp>
#include <fei_utils.hpp>
#include <test_utils/DataReader.hpp>
#include <test_utils/FEI_tester.hpp>     //FEI_tester tests the "old" FEI.
#include <test_utils/snl_fei_tester.hpp> //snl_fei_tester tests the "new" fei

#undef fei_file
#define fei_file "test_FEI.cpp"
#include <fei_ErrMacros.hpp>

test_FEI::test_FEI(MPI_Comm comm)
  : tester(comm),
    fileName_()
{
}

test_FEI::~test_FEI()
{
}

int test_FEI::runtests()
{
  if (fileName_.empty()) return(-1);

  CHK_ERR( test1() );
  return(0);
}

int test_FEI::test1()
{
  fei::SharedPtr<fei::ParameterSet> test_params = get_test_parameters();

  int errcode = 0;

  int numSolves = 0;
  errcode = test_params->getIntParamValue("NUM_SOLVES", numSolves);

  std::string solverName;
  errcode = test_params->getStringParamValue("SOLVER_LIBRARY", solverName);
  if (errcode != 0) {
    throw std::runtime_error(".input file doesn't contain 'SOLVER_LIBRARY'");
  }

  std::string whichFEI;
  errcode = test_params->getStringParamValue("WHICH_FEI", whichFEI);
  if (errcode != 0) {
    throw std::runtime_error(".input file doesn't contain 'WHICH_FEI'");
  }

  std::string inputFileName;
  errcode = test_params->getStringParamValue("INPUT_FILE",inputFileName);
  if (errcode != 0) {
    throw std::runtime_error(".input file doesn't contain 'INPUT_FILE'");
  }

  std::string paramFile;
  errcode = test_params->getStringParamValue("PARAM_FILE", paramFile);
  if (errcode != 0) {
    throw std::runtime_error(".input file doesn't contain 'PARAM_FILE'");
  }

  std::string solnFile;
  errcode = test_params->getStringParamValue("SOLN_FILE", solnFile);
  if (errcode != 0) {
    throw std::runtime_error(".input file doesn't contain 'SOLN_FILE'");
  }

  std::string checkFile;
  errcode = test_params->getStringParamValue("CHECK_FILE", checkFile);
  if (errcode != 0) {
    throw std::runtime_error(".input file doesn't contain 'CHECK_FILE'");
  }

  std::string readerType;
  errcode = test_params->getStringParamValue("READER_TYPE", readerType);

  fei::SharedPtr<feitester> fei_tester;

  FEI_OSTRINGSTREAM fullName;
  std::string fullInputFileName = fully_qualified_name(inputFileName);
  fullName<< fullInputFileName<<"."<< numProcs_<<"."<< localProc_;

  std::string fullName_str = fullName.str();
  const char* fullName_c_str = fullName_str.c_str();

  fei::SharedPtr<DataReader> data_reader(new DataReader);

  //let's read all the data out of the input file.
  CHK_ERR( data_reader->readData(fullName_c_str));

  //now we'll get the parameters out of another file.
  std::string fullParamFileName = fully_qualified_name(paramFile);
  const char* paramFile_c_str = fullParamFileName.c_str();

  std::vector<std::string> param_file_contents;
  fei_test_utils::read_input_file(paramFile_c_str, comm_,
                                  param_file_contents);

  data_reader->numParams_ = param_file_contents.size();
  data_reader->paramStrings_ = new char*[data_reader->numParams_];

  for(unsigned i=0; i<param_file_contents.size(); ++i) {
    std::string& str = param_file_contents[i];
    data_reader->paramStrings_[i] = new char[str.size()+1];
    for(unsigned j=0; j<str.size(); ++j) {
      data_reader->paramStrings_[i][j] = str[j];
    }
    data_reader->paramStrings_[i][str.size()] = '\0';
  }

  data_reader->solverLibraryName_ = solverName;

  std::string fullSolnFile = fully_qualified_name(solnFile);
  data_reader->solnFileName_ = fullSolnFile;

  std::string fullCheckFile = fully_qualified_name(checkFile);
  data_reader->checkFileName_ = fullCheckFile;

  //ok, all the data is in the 'data' object, so we're ready to start
  //handing it all over to an instantiation of the FEI.

  if (whichFEI == "OLDFEI") {
    fei_tester.reset(new FEI_tester(data_reader, comm_, localProc_, numProcs_));
  }
  if (whichFEI == "fei::FEI_Impl") {
    bool useNewFEI = true;
    fei_tester.reset(new FEI_tester(data_reader, comm_, localProc_, numProcs_,
      			      useNewFEI));
  }
  else if (whichFEI == "new_fei") {
    fei_tester.reset( new snl_fei_tester(data_reader, comm_, localProc_, numProcs_));
  }

  if (fei_tester.get() == NULL) {
    ERReturn(-1);
  }

  fei_tester->setPath(path_);

  int errCode = fei_tester->testInitialization();
  if (errCode < 0) {
    ERReturn(errCode);
  }
  if (errCode > 0) {
    FEI_COUT << "library " << solverName << " not available."<<FEI_ENDL;
    return(-1);
  }

  for(int solveCounter=1; solveCounter<=numSolves; solveCounter++) {

    CHK_ERR( fei_tester->testLoading() );

//    if (solverName == "Trilinos") {
//      fei_tester->dumpMatrixFiles();
//      fei_tester->setParameter("USE_FEI_MATRIX_LOCAL true");
//      fei_tester->testLoading();
//      fei_tester->dumpMatrixFiles();
//    }

    CHK_ERR( fei_tester->testSolve() );

    CHK_ERR( fei_tester->testCheckResult() );

  } //end of 'for solveCounter<=numSolves loop'

  return(0);
}

fei::SharedPtr<fei::ParameterSet>
test_FEI::get_test_parameters()
{
  std::vector<std::string> inputFileStrings;

  FEI_OSTRINGSTREAM filename;
  filename << path_;

  std::size_t length = path_.size();
  if (length > 0) {
    if (path_[length-1] != '/') filename<< "/";
  }

  filename << fileName_;

  std::string filename_str = filename.str();
  const char* filename_c_str = filename_str.c_str();
  fei_test_utils::read_input_file(filename_c_str, comm_, inputFileStrings);

  fei::SharedPtr<fei::ParameterSet> paramset(new fei::ParameterSet);
  fei::utils::parse_strings(inputFileStrings, " ", *paramset);

  return(paramset);
}

std::string test_FEI::fully_qualified_name(const std::string& fileName)
{
  FEI_OSTRINGSTREAM osstr;
  osstr << path_;

  std::size_t length = path_.size();
  if (length > 0) {
    if (path_[length-1] != '/') osstr<< "/";
  }

  osstr << fileName;

  return(osstr.str());
}
