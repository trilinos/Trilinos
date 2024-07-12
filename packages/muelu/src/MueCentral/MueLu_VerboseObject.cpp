// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

VerboseObject::VerboseObject()
  : verbLevel_(NotSpecified)
  ,  // = use global verbose level by default
  numProcs_(0) {
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
  int mpiStarted = 0;
  MPI_Initialized(&mpiStarted);
  if (mpiStarted) {
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank_);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs_);
  }
#endif
}

VerboseObject::~VerboseObject() {}

VerbLevel VerboseObject::GetVerbLevel() const {
  return (verbLevel_ != NotSpecified ? verbLevel_ : globalVerbLevel_);
}

void VerboseObject::SetVerbLevel(const VerbLevel verbLevel) {
  verbLevel_ = verbLevel;
}

int VerboseObject::GetProcRankVerbose() const {
  return procRank_;
}

int VerboseObject::SetProcRankVerbose(int procRank) const {
  int oldRank = procRank_;
  procRank_   = procRank;
  return oldRank;
}

bool VerboseObject::IsPrint(MsgType type, int thisProcRankOnly) const {
  return ((type & GetVerbLevel()) && (thisProcRankOnly < 0 || procRank_ == thisProcRankOnly));
}

Teuchos::FancyOStream& VerboseObject::GetOStream(MsgType type, int thisProcRankOnly) const {
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
  int mpiStarted = 0;
  MPI_Initialized(&mpiStarted);
  if (mpiStarted) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    fn = filename + "." + std::to_string(procRank);
  } else
#endif
    fn = filename;
  RCP<std::ofstream> outFile(new std::ofstream(fn));
  Teuchos::RCP<Teuchos::FancyOStream> fancyOutFile = Teuchos::fancyOStream(Teuchos::rcp_implicit_cast<std::ostream>(outFile));
  SetMueLuOStream(fancyOutFile);
}

Teuchos::RCP<Teuchos::FancyOStream> VerboseObject::GetMueLuOStream() {
  if (mueluOutputStream_.get() == NULL) {
    mueluOutputStream_ = fancyOStream(rcpFromRef(std::cout));
    mueluOutputStream_->setOutputToRootOnly(-1);
  }
  return mueluOutputStream_;
}

VerbLevel VerboseObject::globalVerbLevel_ = High;  // Default global verbose level.

RCP<Teuchos::FancyOStream> VerboseObject::blackHole_ = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));

RCP<Teuchos::FancyOStream> VerboseObject::mueluOutputStream_ = Teuchos::null;

}  // namespace MueLu
