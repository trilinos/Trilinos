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
#ifndef MUELU_TWOLEVELFACTORY_HPP
#define MUELU_TWOLEVELFACTORY_HPP

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_Factory.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_TimeMonitor.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  /*!
    @class TwoLevelFactoryBase class.
    @brief Base class for factories that use two levels (fineLevel and coarseLevel).

    Examples of such factories are R, P, and A_coarse.

    @ingroup MueLuBaseClasses
  */

  class TwoLevelFactoryBase : public Factory {

  public:

    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    TwoLevelFactoryBase() { }

    //! Destructor.
    virtual ~TwoLevelFactoryBase()
    { }

    //@}

    //! Input
    //@{

    /*! @brief Specifies the data that this class needs, and the factories that generate that data.

        If the Build method of this class requires some data, but the generating factory is not specified in DeclareInput, then this class
        will fall back to the settings in FactoryManager.
    */
    virtual void DeclareInput(Level &fineLevel, Level &coarseLevel) const = 0;

    //!
    virtual void CallDeclareInput(Level & requestedLevel) const {
      if (requestedLevel.GetPreviousLevel() == Teuchos::null) {
        std::ostringstream errStr;
        errStr << "LevelID = " << requestedLevel.GetLevelID();
        throw Exceptions::DependencyError(errStr.str());
      }
      DeclareInput(*requestedLevel.GetPreviousLevel(), requestedLevel);
    }

    //@}

    //! @name Build methods.
    //@{


    //! Build an object with this factory.
    virtual void Build(Level & fineLevel, Level & coarseLevel) const = 0;

    //!
    virtual void CallBuild(Level& requestedLevel) const {
      int levelID = requestedLevel.GetLevelID();

#ifdef HAVE_MUELU_DEBUG
      // We cannot call Build method twice for the same level, but we can call it multiple times for different levels
      TEUCHOS_TEST_FOR_EXCEPTION((multipleCallCheck_ == ENABLED) && (multipleCallCheckGlobal_ == ENABLED) && (lastLevelID_ == levelID),
                                 Exceptions::RuntimeError,
                                 this->ShortClassName() << "::Build() called twice for the same level (levelID=" << levelID
                                 << "). This is likely due to a configuration error, or calling hierarchy setup multiple times "
                                 << "without resetting debug info through FactoryManager::ResetDebugData().");
      if (multipleCallCheck_ == FIRSTCALL)
        multipleCallCheck_ = ENABLED;

      lastLevelID_ = levelID;
#endif
      TEUCHOS_TEST_FOR_EXCEPTION(requestedLevel.GetPreviousLevel() == Teuchos::null, Exceptions::RuntimeError, "LevelID = " << levelID);

      RCP<const Teuchos::Comm<int> > comm = requestedLevel.GetComm();
      if (comm.is_null()) {
        // Some factories are called before we constructed Ac, and therefore,
        // before we set the level communicator. For such factories we can get
        // the comm from the previous level, as all processes go there
        RCP<Level>& prevLevel = requestedLevel.GetPreviousLevel();
        if (!prevLevel.is_null())
          comm = prevLevel->GetComm();
      }

      int oldRank = -1;
      if (!comm.is_null())
        oldRank = SetProcRankVerbose(comm->getRank());

      // Synchronization timer
      std::string syncTimer = this->ShortClassName() + ": Build sync (level=" + toString(requestedLevel.GetLevelID()) + ")";
      if (this->timerSync_ && !comm.is_null()) {
        TimeMonitor timer(*this, syncTimer);
        comm->barrier();
      }

      Build(*requestedLevel.GetPreviousLevel(), requestedLevel);

      if (this->timerSync_ && !comm.is_null()) {
        TimeMonitor timer(*this, syncTimer);
        comm->barrier();
      }

      GetOStream(Test) << *RemoveFactoriesFromList(GetParameterList()) << std::endl;

      if (oldRank != -1)
        SetProcRankVerbose(oldRank);
    }

    //@}

  }; //class TwoLevelFactoryBase


} //namespace MueLu

#define MUELU_TWOLEVELFACTORY_SHORT
#endif //ifndef MUELU_TWOLEVELFACTORY_HPP
