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
#ifndef MUELU_FACTORYMANAGERBASE_HPP
#define MUELU_FACTORYMANAGERBASE_HPP

#include <string>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"

#include "MueLu_FactoryBase_fwd.hpp"

namespace MueLu {

  /*!
     @class FactoryManagerBase
     @brief Class that provides default factories within Needs class.
     @ingroup MueLuBaseClasses
  */
  class FactoryManagerBase : public BaseClass {

  public:
    //@{ Constructors/Destructors.

    //! Destructor.
    virtual ~FactoryManagerBase() { }

    //@}

    //@{ Get/Set functions.

    //! Get
    // Return ref because user also give ref to the Hierarchy.
    const virtual RCP<const FactoryBase> & GetFactory(const std::string & varName) const = 0;
    //@}

    // Free temporarily hold data at the end of Hierarchy::Setup()
    // This method is const because the clean concerns only mutable data.
    virtual void Clean() const { } // TODO: should be used inside of MueLu::Hierarchy

    //! returns internal IgnoreUserData flag
    //! The FactoryManager has some control over the Level::GetFactory function
    //! If IgnoreUserData is set, the Level::GetFactory function always asks the factory manager for a valid factory
    //! Otherwise user-given data is preferred, if available (default behaviour)
    bool IgnoreUserData() const { return bIgnoreUserData_; }

    //! set IgnoreUserData flag
    void SetIgnoreUserData(bool bIgnoreUserData = false) { bIgnoreUserData_ = bIgnoreUserData; }

  private:
    //! boolean flag that controls behaviour of Level::GetFactory
    //! if bIgnoreUserData == true, the Level::GetFactory function always asks the Factory manager for a valid factory given a variable name
    //! if bIgnoreUserData == false, the Level::GetFactory prefers user-provided data for a variable name if available. Otherwise the factory manager is asked for a valid factory
    //! default: bIgnoreUserData = false;
    bool bIgnoreUserData_;

  }; // class FactoryManagerBase

} // namespace MueLu

#define MUELU_FACTORYMANAGERBASE_SHORT
#endif //ifndef MUELU_FACTORYMANAGERBASE_HPP
