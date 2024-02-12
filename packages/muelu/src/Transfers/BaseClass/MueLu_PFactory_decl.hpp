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
#ifndef MUELU_PFACTORY_DEF_HPP
#define MUELU_PFACTORY_DEF_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_PFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"

namespace MueLu {

/*!
  @class PFactory
  @brief Factory that provides an interface for a concrete implementation of a prolongation operator.

  For a concrete implementation the user has to overwrite the virtual Build method.
*/

class PFactory : public TwoLevelFactoryBase {
 protected:
  bool restrictionMode_;  //< true, if PFactory is used for generating the restriction operator

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  PFactory()
    : restrictionMode_(false) {}

  //! Destructor.
  virtual ~PFactory() {}
  //@}

  //! @name Build methods.
  //@{

  /*!
    @brief Abstract Build method.
  */
  virtual void BuildP(Level &fineLevel, Level &coarseLevel) const = 0;
  //@}

  //! @name Restriction mode
  //@{

  /// switch prolongator factory to restriction mode
  /// if set to true, the prolongation factory generates a restriciton operator instead of a prolongation operator
  void setRestrictionMode(bool bRestrictionMode = false) {
    restrictionMode_ = bRestrictionMode;
  }

  /// returns restrictionMode flag
  bool isRestrictionModeSet() { return restrictionMode_; }

  //@}

};  // class PFactory

}  // namespace MueLu

#define MUELU_PFACTORY_SHORT
#endif  // MUELU_PFACTORY_DEF_HPP
