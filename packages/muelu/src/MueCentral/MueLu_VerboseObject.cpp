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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include "MueLu_VerboseObject.hpp"

#include <fstream>
#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <Teuchos_VerboseObject.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_VerbosityLevel.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

  VerboseObject::VerboseObject():
    verbLevel_(NotSpecified) // = use global verbose level by default
  { }

  VerboseObject::~VerboseObject() { }

  VerbLevel VerboseObject::GetVerbLevel() const {
    return (verbLevel_ != NotSpecified ?  verbLevel_ : globalVerbLevel_);
  }

  void VerboseObject::SetVerbLevel(const VerbLevel verbLevel) {
    verbLevel_ = verbLevel;
  }

  int VerboseObject::GetProcRankVerbose() const {
    return procRank_;
  }

  int VerboseObject::SetProcRankVerbose(int procRank) const {
    int oldRank = procRank_;
    procRank_ = procRank;
    return oldRank;
  }

  bool VerboseObject::IsPrint(MsgType type, int thisProcRankOnly) const {
    return ((type & GetVerbLevel()) && (thisProcRankOnly < 0 || procRank_ == thisProcRankOnly));
  }

  Teuchos::FancyOStream & VerboseObject::GetOStream(MsgType type, int thisProcRankOnly) const {
    if (!IsPrint(type, thisProcRankOnly))
      return *blackHole_;

    Teuchos::FancyOStream& os = *GetMueLuOStream();
    if (!(type & ((Extreme | Test) ^ Warnings)))
      os << "\n******* WARNING *******" << std::endl;

    return os;
  }

  Teuchos::FancyOStream& VerboseObject::GetBlackHole() const {
    return *blackHole_;
  }

  void VerboseObject::SetDefaultVerbLevel(const VerbLevel defaultVerbLevel) {
    TEUCHOS_TEST_FOR_EXCEPTION(defaultVerbLevel == NotSpecified, Exceptions::RuntimeError,
                               "MueLu::VerboseObject::GetVerbLevel(): global verbose level cannot be 'NotSpecified'.");
    globalVerbLevel_ = defaultVerbLevel;
  }

  VerbLevel VerboseObject::GetDefaultVerbLevel() {
    return globalVerbLevel_;
  }

  void VerboseObject::SetMueLuOStream(const Teuchos::RCP<Teuchos::FancyOStream>& mueluOStream) {
    mueluOStream->setOutputToRootOnly(-1);
    mueluOutputStream_ = mueluOStream;
  }

  void VerboseObject::SetMueLuOFileStream(const std::string& filename) {
    std::string fn;
#ifdef HAVE_MPI
    int mpiStarted = 0; MPI_Initialized(&mpiStarted);
    if (mpiStarted)     {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD , &rank); // CHECK: ALLOW MPI_COMM_WORLD
      fn = filename + "." + std::to_string(rank);
    } else
#endif
      fn = filename;
    RCP<std::ofstream> outFile(new std::ofstream(fn));
    Teuchos::RCP<Teuchos::FancyOStream> fancyOutFile = Teuchos::fancyOStream(Teuchos::rcp_implicit_cast<std::ostream>(outFile));
    SetMueLuOStream(fancyOutFile);
  }

  Teuchos::RCP<Teuchos::FancyOStream> VerboseObject::GetMueLuOStream() {
    if (mueluOutputStream_.get()==NULL) {
      mueluOutputStream_ = fancyOStream(rcpFromRef(std::cout));
      mueluOutputStream_->setOutputToRootOnly(-1);
    }
    return mueluOutputStream_;
  }


  VerbLevel VerboseObject::globalVerbLevel_ = High; // Default global verbose level.

  RCP<Teuchos::FancyOStream> VerboseObject::blackHole_ = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));

  RCP<Teuchos::FancyOStream> VerboseObject::mueluOutputStream_ = Teuchos::null;

} // namespace MueLu
