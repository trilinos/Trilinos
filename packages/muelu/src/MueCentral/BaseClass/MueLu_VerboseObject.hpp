#ifndef MUELU_VERBOSECLASS_HPP
#define MUELU_VERBOSECLASS_HPP

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <Teuchos_VerboseObject.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_VerbosityLevel.hpp"

namespace MueLu {

  // This class is inspired by BelosOutputManager but:
  // - it uses Teuchos stream
  // - it allows to print message on proc0 or every processors
  // - it is a base class

  // TODO:
  // - Right now, you can't change the verbosity level globally.
  //   We have to take into account Teuchos::EVerbosityLevel from the base class Teuchos::VerboseObject and somehow convert it to a MsgType.
  // - Allow users to define the global verbosity level using both MueLu::MsgType and Teuchos::EVerbosityLevel?

  //! Verbose class for MueLu classes
  class VerboseObject
    : public Teuchos::VerboseObject<VerboseObject>
  {
    
  public:
    
    VerboseObject()
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

    //! Destructor.
    virtual ~VerboseObject() { }

    //! Get the verbosity level
    VerbLevel GetVerbLevel() const { return verbLevel_; }
  
    //! Set the verbosity level
    void SetVerbLevel(const VerbLevel verblevel) { verbLevel_ = verblevel; }

    //! Get proc rank used for printing (do not use this information for any other purpose)
    int GetProcRankVerbose() const { return procRank_; }

    //! Find out whether we need to print out information for a specific message type.
    /*! This method is used to determine whether computations are necessary for this message type. */
    bool IsPrint(MsgType type, int thisProcRankOnly = -1) const { 
      return ((type & verbLevel_) && (thisProcRankOnly < 0 || procRank_ == thisProcRankOnly));
    };

    //! Get an output stream for outputting the input message type.
    Teuchos::FancyOStream & GetOStream(MsgType type, int thisProcRankOnly = -1) const {
      return (IsPrint(type, thisProcRankOnly)) ? *getOStream() : *blackHole_;
    }

    Teuchos::FancyOStream & GetBlackHole() const { return *blackHole_; }

  private:
    VerbLevel verbLevel_;
    int procRank_;

    static RCP<Teuchos::FancyOStream> blackHole_;
  }; // class VerboseObject

} // namespace MueLu

#define MUELU_VERBOSECLASS_SHORT
#endif // ifndef MUELU_VERBOSECLASS_HPP
