/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


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

