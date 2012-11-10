// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include "MueLu_VerboseObject.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <Teuchos_VerboseObject.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_VerbosityLevel.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

  VerboseObject::VerboseObject()
    : verbLevel_(NotSpecified) // = use global verbose level by default
  {
    // Note: using MPI_COMM_RANK is bad idea (because a subcommunicator may be used to run MueLu)
    // Belos have the same problem in the class BelosOutputManager.
    //
    // How to fix this: the FancyOStream provides setProcRankAndSize() to change the proc rank info used by setOutputToRootOnly().
    // Adding a method FancyOStream::getProcRank() would be enough. And it makes sense to use the info that come from the stream configuration.
    //
    // Documentation: after the patch, it migh be nice to add to the documentation (Teuchos and MueLu)
    //                that users of subcommunicators have to setup the output stream (of Teuchos::VerboseObject) separately.
    procRank_ = 0;
#ifdef HAVE_MPI
    int mpiStarted = 0; MPI_Initialized(&mpiStarted);
    if (mpiStarted)     MPI_Comm_rank(MPI_COMM_WORLD, &procRank_);
#endif
  }

  VerboseObject::~VerboseObject() { }

  VerbLevel VerboseObject::GetVerbLevel() const {
    if (verbLevel_ != NotSpecified)
      return verbLevel_;
    //     else if ()
    else
      return globalVerbLevel_;
  }

  void VerboseObject::SetVerbLevel(const VerbLevel verbLevel) { verbLevel_ = verbLevel; }

  int VerboseObject::GetProcRankVerbose() const { return procRank_; }

  bool VerboseObject::IsPrint(MsgType type, int thisProcRankOnly) const {
    return ((type & GetVerbLevel()) && (thisProcRankOnly < 0 || procRank_ == thisProcRankOnly));
  }

  Teuchos::FancyOStream & VerboseObject::GetOStream(MsgType type, int thisProcRankOnly) const {
    return (IsPrint(type, thisProcRankOnly)) ? *getOStream() : *blackHole_;
  }

  Teuchos::FancyOStream & VerboseObject::GetBlackHole() const { return *blackHole_; }

  RCP<Teuchos::FancyOStream> VerboseObject::blackHole_ = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));

  void VerboseObject::SetDefaultVerbLevel(const VerbLevel defaultVerbLevel) {
    TEUCHOS_TEST_FOR_EXCEPTION(defaultVerbLevel == NotSpecified, Exceptions::RuntimeError, "MueLu::VerboseObject::GetVerbLevel(): global verbose level cannot be 'NotSpecified'.");
    globalVerbLevel_ = defaultVerbLevel;
  }

  VerbLevel VerboseObject::GetDefaultVerbLevel() {
    return globalVerbLevel_;
  }

  VerbLevel VerboseObject::globalVerbLevel_ = High; // Default global verbose level.

} // namespace MueLu
