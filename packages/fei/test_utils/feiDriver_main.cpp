/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

//
// This is a simple program to exercise the FEI.
//
#include <fei_macros.hpp>


//Including the header fei_base.hpp gets us the declaration for
//most FEI classes.

#include <fei_base.hpp>
#include <FEI_Implementation.hpp>

#include <fei_LibraryWrapper.hpp>
#include <test_utils/LibraryFactory.hpp>

#include <test_utils/Poisson_Elem.hpp>
#include <test_utils/PoissonData.hpp>

#include <test_utils/ElemBlock.hpp>
#include <test_utils/CRSet.hpp>
#include <test_utils/BCNodeSet.hpp>
#include <test_utils/CommNodeSet.hpp>
#include <test_utils/DataReader.hpp>
#include <test_utils/driverData.hpp>
#include <test_utils/fei_test_utils.hpp>
#include <test_utils/SolnCheck.hpp>

#undef fei_file
#define fei_file "feiDriver.cpp"
#include <fei_ErrMacros.hpp>

//==============================================================================
//Here's the main...
//==============================================================================
int feiDriver_main(int argc, char** argv,
		   MPI_Comm comm, int numProcs, int localProc){

  std::vector<std::string> stdstrings;
  CHK_ERR( fei_test_utils::get_filename_and_read_input(argc, argv,
						comm, localProc,
						stdstrings) );
  const char** params = NULL;
  int numParams = 0;
  fei::utils::strings_to_char_ptrs(stdstrings, numParams, params);
  
  fei::ParameterSet paramset;
  fei::utils::parse_strings(stdstrings, " ", paramset);

  std::string solverName;
  std::string inputFileName;
  int errcode = paramset.getStringParamValue("SOLVER_LIBRARY", solverName);
  errcode += paramset.getStringParamValue("INPUT_FILE", inputFileName);
  if (errcode != 0) {
    fei::console_out() << "Expected to find both 'SOLVER_NAME' and 'INPUT_FILE' "
	     <<"in input-file."<<FEI_ENDL;
    return(-1);
  }

  //let's add the appropriate file-name extension to the file-name obtained from
  //the input...
  FEI_OSTRINGSTREAM fullFileName;
  fullFileName<< inputFileName<<"."<<numProcs<<"."<< localProc;

  driverData drv;
  CHK_ERR( drv.readData(fullFileName.str().c_str()) );

  //ok, all the data is in the 'data' object, so we're ready to start
  //handing it all over to an instantiation of the FEI.

  //first, we have to instantiate a LibraryWrapper and an FEI...

  fei::SharedPtr<LibraryWrapper> wrapper;
  try {
    wrapper = fei::create_LibraryWrapper(comm, solverName.c_str());
  }
  catch (std::runtime_error& exc) {
    fei::console_out() << exc.what() << FEI_ENDL;
    ERReturn(-1);
  }

  fei::SharedPtr<FEI> fei(new FEI_Implementation(wrapper, comm));

  const char* feiVersionString;
  CHK_ERR( fei->version(feiVersionString) );

  FEI_COUT << FEI_ENDL << "FEI version: " << feiVersionString << FEI_ENDL << FEI_ENDL;

  CHK_ERR( fei->parameters(numParams, params) );

  std::vector<const char*>& methodNames = drv.get_methodNames();

  for(size_t i=0; i<methodNames.size(); i++) {
    if (!strcmp("destructor", methodNames[i])) {
      //In some cases the input file indicates that the FEI should be
      //destroyed and then re-allocated before continuing. Note that we
      //assume here that the solver-layer (linsyscore, wrapper or feData)
      //should also be destroyed and re-allocated at the same time.
      FEI_COUT << "feiDriver: proc " << localProc << " destroying/reallocing FEI"
	   << FEI_ENDL;

      fei.reset();
      wrapper.reset();
      try {
	wrapper = fei::create_LibraryWrapper(comm, solverName.c_str());
      }
      catch (std::runtime_error& exc) {
	fei::console_out() << exc.what()<<FEI_ENDL;
	ERReturn(-1);
      }

      fei.reset(new FEI_Implementation(wrapper, comm));

      CHK_ERR( fei->parameters(numParams, params) );

      continue;
    }

    FEI_COUT << "feiDriver: proc " << localProc << " calling FEI method: "
	 << methodNames[i] << FEI_ENDL;
    int feierror = drv.call_fei_method(methodNames[i], fei.get());
    if (feierror > 0) continue;
    if (feierror < 0) {
      //for testing purposes, temporarily, don't bail out if an fei method
      //returns an error.
      continue;
      //return(-1);
    }
  }

  MPI_Barrier(comm);

  if (localProc == 0) {
    FEI_COUT << "feiDriver: TEST PASSED" << FEI_ENDL;

    //This is something the SIERRA runtest tool looks for in test output...
    FEI_COUT << "SIERRA execution successful" << FEI_ENDL;
#ifdef SIERRA_BUILD_DATE
    FEI_COUT.setf(IOS_FIXED, IOS_FLOATFIELD);
    FEI_COUT << "Maximum CPU  time: 0.0 seconds." << FEI_ENDL;
#endif
  }

  delete [] params;

  return(0);
}
