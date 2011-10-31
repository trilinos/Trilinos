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
    ;

    //! Destructor.
    virtual ~VerboseObject() ;

    //! Get the verbosity level
    VerbLevel GetVerbLevel() const ;
  
    //! Set the verbosity level
    void SetVerbLevel(const VerbLevel verbLevel) ;

    //! Get proc rank used for printing (do not use this information for any other purpose)
    int GetProcRankVerbose() const ;

    //! Find out whether we need to print out information for a specific message type.
    /*! This method is used to determine whether computations are necessary for this message type. */
    bool IsPrint(MsgType type, int thisProcRankOnly = -1) const ;

    //! Get an output stream for outputting the input message type.
    Teuchos::FancyOStream & GetOStream(MsgType type, int thisProcRankOnly = -1) const ;

    Teuchos::FancyOStream & GetBlackHole() const ;

  private:
    VerbLevel verbLevel_;
    int procRank_;

    static RCP<Teuchos::FancyOStream> blackHole_;
  }; // class VerboseObject

} // namespace MueLu

#define MUELU_VERBOSECLASS_SHORT
#endif // ifndef MUELU_VERBOSECLASS_HPP
