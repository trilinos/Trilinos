/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <fei_mpi.h>

#include <test_utils/LibraryFactory.hpp>

#include <fei_LibraryWrapper.hpp>

#include <snl_fei_Factory.hpp>

#include <fei_Factory_Trilinos.hpp>
#ifdef HAVE_FEI_AZTECOO
#include <fei_Aztec_LinSysCore.hpp>
#endif

#ifdef HAVE_FEI_FETI
#include <FETI_DP_FiniteElementData.h>
#endif

//----------------------------------------------------------------------------
fei::SharedPtr<LibraryWrapper>
fei::create_LibraryWrapper(MPI_Comm comm,
			       const char* libraryName)
{
  std::string libname(libraryName);

  fei::SharedPtr<LinearSystemCore> lsc;
  fei::SharedPtr<FiniteElementData> fedata;
  fei::SharedPtr<LibraryWrapper> wrapper;

  if (libname == "Aztec") {
#ifdef HAVE_FEI_AZTECOO
    lsc.reset(new fei_trilinos::Aztec_LinSysCore(comm));
#else
    std::string msg("Aztec not available.");
    throw std::runtime_error(msg);
#endif
  }

  if (libname == "FETI") {
#ifdef HAVE_FEI_FETI
    fedata.reset(new FETI_DP_FiniteElementData(comm));
#endif
  }

  if (lsc.get() == NULL && fedata.get() == NULL) {
    //libraryName not found
    std::string msg("create_LibraryWrapper: ");
    msg += libraryName;
    msg += " not a valid name.";
    throw std::runtime_error(msg);
  }

  if (lsc.get() != NULL) {
    wrapper.reset(new LibraryWrapper(lsc));
    return(wrapper);
  }

  if (fedata.get() != NULL) {
    wrapper.reset(new LibraryWrapper(fedata));
    return(wrapper);
  }

  return(wrapper);
}

//----------------------------------------------------------------------------
fei::SharedPtr<fei::Factory>
fei::create_fei_Factory(MPI_Comm comm,
			    const char* libraryName)
{
  std::string libname(libraryName);

  if (libname.find("Trilinos") != std::string::npos) {
    fei::SharedPtr<fei::Factory> factory(new Factory_Trilinos(comm));

    if (libname.find("Amesos") != std::string::npos) {
      fei::ParameterSet paramset;
      paramset.add(fei::Param("Trilinos_Solver", "Amesos"));
      factory->parameters(paramset);
    }
    else if (libname.find("Aztec") != std::string::npos) {

      //if libname contains "AztecOO" we'll return the Trilinos factory
      //but if libname only contains "Aztec" then we want to skip on down
      //and return an snl_fei::Factory with a Aztec LibraryWrapper...

      if (libname.find("AztecOO") != std::string::npos) {
        return(factory);
      }
    }
    else {
      //This else handles the case where libname contains "Trilinos", but
      //doesn't contain "Aztec" or "Amesos"...
      return(factory);
    }
  }

  fei::SharedPtr<LibraryWrapper> wrapper;
  try {
    wrapper = fei::create_LibraryWrapper(comm, libraryName);
  }
  catch (std::runtime_error& exc) {
    std::string msg("create_fei_Factory: ");
    msg += exc.what();
    throw std::runtime_error(msg);
  }

  if (wrapper.get() != NULL) {
    fei::SharedPtr<fei::Factory> factory(new snl_fei::Factory(comm, wrapper));
    return(factory);
  }

  fei::SharedPtr<fei::Factory> empty;
  return(empty);
}

