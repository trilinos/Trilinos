#include "MueLu_VerboseObject.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <Teuchos_VerboseObject.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_VerbosityLevel.hpp"

namespace MueLu {

  VerboseObject::VerboseObject()
    : verbLevel_(High) //TODO: Default
  { 
    // Note: using MPI_COMM_RANK is bad idea (because maybe a subcommunicator is used to run MueLu)
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

  VerbLevel VerboseObject::GetVerbLevel() const { return verbLevel_; }

  void VerboseObject::SetVerbLevel(const VerbLevel verbLevel) { verbLevel_ = verbLevel; }

  int VerboseObject::GetProcRankVerbose() const { return procRank_; }

  bool VerboseObject::IsPrint(MsgType type, int thisProcRankOnly) const { 
    return ((type & verbLevel_) && (thisProcRankOnly < 0 || procRank_ == thisProcRankOnly));
  }

  Teuchos::FancyOStream & VerboseObject::GetOStream(MsgType type, int thisProcRankOnly) const {
    return (IsPrint(type, thisProcRankOnly)) ? *getOStream() : *blackHole_;
  }

  Teuchos::FancyOStream & VerboseObject::GetBlackHole() const { return *blackHole_; }

  RCP<Teuchos::FancyOStream> VerboseObject::blackHole_ = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));

} // namespace MueLu
