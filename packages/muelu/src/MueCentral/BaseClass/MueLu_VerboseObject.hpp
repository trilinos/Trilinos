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
#ifndef MUELU_VERBOSEOBJECT_HPP
#define MUELU_VERBOSEOBJECT_HPP

#include "MueLu_ConfigDefs.hpp"
#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <Teuchos_VerboseObject.hpp>

#include "MueLu_VerbosityLevel.hpp"

namespace MueLu {

  /*!
    @class VerboseObject class.
    @brief Verbose class for MueLu classes

    This class is inspired by BelosOutputManager but:
      - it uses Teuchos stream
      - it allows to print message on proc0 or every processors
      - it is a base class

     @todo TODO
     - Right now, you can't change the verbosity level globally.
       We have to take into account Teuchos::EVerbosityLevel from the base class Teuchos::VerboseObject and somehow convert it to a MsgType.
     - Allow users to define the global verbosity level using both MueLu::MsgType and Teuchos::EVerbosityLevel?

     @todo TODO
     -# add a function to get the object 'local' verb level (GetVerbLevel does not do this)?
     -# interface without fancyOSstream but std::stream
     -# add 'const' and mutable

    @ingroup MueLuBaseClasses

  */
  class VerboseObject
    : public Teuchos::VerboseObject<VerboseObject>
  {

  public:

    //! @name Constructors/Destructors.
    //@{
    VerboseObject();

    //! Destructor.
    virtual ~VerboseObject();
    //@}

    /*! @brief Get the verbosity level

       If a verbosity level has not been specified for this object (using SetVerbLevel), this method returns the default/global verbose level.
    */
    VerbLevel GetVerbLevel() const;

    //! Set the verbosity level of this object
    void SetVerbLevel(const VerbLevel verbLevel);

    //! Get proc rank used for printing. <b>Do not use this information for any other purpose.</b>
    int GetProcRankVerbose() const;

    /*! @brief Find out whether we need to print out information for a specific message type.

        This method is used to determine whether computations are necessary for this message type.
    */
    bool IsPrint(MsgType type, int thisProcRankOnly = -1) const;

    //! Get an output stream for outputting the input message type.
    Teuchos::FancyOStream & GetOStream(MsgType type, int thisProcRankOnly = -1) const;

    Teuchos::FancyOStream & GetBlackHole() const;

    static void SetDefaultOStream(const RCP<Teuchos::FancyOStream> &defaultOStream);

    static RCP<Teuchos::FancyOStream> GetDefaultOStream();

    //! @name Public static member functions
    //@{

    //! Set the default (global) verbosity level.
    static void SetDefaultVerbLevel(const VerbLevel defaultVerbLevel);

    //! Get the default (global) verbosity level.
    static VerbLevel GetDefaultVerbLevel();

    //@}

  private:
    //! Verbose level specific to 'this'
    VerbLevel verbLevel_;
    int procRank_;

    static RCP<Teuchos::FancyOStream> blackHole_;

    //! Global verbose level. This verbose level is used when the verbose level of the object is not specified (verbLevel_ == NotSpecified)
    static VerbLevel globalVerbLevel_;
  }; // class VerboseObject

} // namespace MueLu

#define MUELU_VERBOSECLASS_SHORT
#endif // MUELU_VERBOSEOBJECT_HPP
