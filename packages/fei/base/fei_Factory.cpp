/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <fei_Factory.hpp>
#include <fei_LogManager.hpp>
#include <fei_LogFile.hpp>
#include <fei_ParameterSet.hpp>

#include <FEI_Implementation.hpp>
#include <fei_FEI_Impl.hpp>

//----------------------------------------------------------------------------
fei::Factory::Factory(MPI_Comm comm)
{
  int numProcs = 1, localProc = 0;
#ifndef FEI_SER
  MPI_Comm_size(comm, &numProcs);
  MPI_Comm_rank(comm, &localProc);
#endif
  fei::LogManager::getLogManager().setNumProcs(numProcs, localProc);
}

//----------------------------------------------------------------------------
fei::Factory::~Factory()
{
  fei::LogFile::getLogFile().closeOutputStream();
  fei::LogManager::getLogManager().setOutputLevel(fei::NONE);
}

//----------------------------------------------------------------------------
void fei::Factory::parameters(const fei::ParameterSet& paramset)
{
  const fei::Param* param = paramset.get("FEI_OUTPUT_PATH");
  fei::Param::ParamType ptype = param != NULL ?
    param->getType() : fei::Param::BAD_TYPE;
  if (ptype == fei::Param::STRING) {
    fei::LogManager& log_manager = fei::LogManager::getLogManager();
    log_manager.setOutputPath(param->getStringValue().c_str());
  }

  param = paramset.get("debugOutput");
  ptype = param != NULL ? param->getType() : fei::Param::BAD_TYPE;
  if (ptype == fei::Param::STRING) {
    fei::LogManager& log_manager = fei::LogManager::getLogManager();
    log_manager.setOutputPath(param->getStringValue().c_str());
  }

  param = paramset.get("FEI_OUTPUT_LEVEL");
  ptype = param != NULL ? param->getType() : fei::Param::BAD_TYPE;
  if (ptype == fei::Param::STRING) {
    fei::LogManager& log_manager = fei::LogManager::getLogManager();
    log_manager.setOutputLevel(param->getStringValue().c_str());
  }
}

//----------------------------------------------------------------------------
fei::SharedPtr<FEI>
fei::Factory::createFEI(fei::SharedPtr<LibraryWrapper> wrapper,
                        MPI_Comm comm)
{
  //fei::SharedPtr<FEI> fei(new fei::FEI_Impl(wrapper, comm));
  fei::SharedPtr<FEI> fei(new FEI_Implementation(wrapper, comm));

  return(fei);
}

//----------------------------------------------------------------------------
fei::SharedPtr<FEI>
fei::Factory::createFEI(MPI_Comm comm)
{
  fei::SharedPtr<FEI> fei(new fei::FEI_Impl(this, comm));

  return(fei);
}

//----------------------------------------------------------------------------

