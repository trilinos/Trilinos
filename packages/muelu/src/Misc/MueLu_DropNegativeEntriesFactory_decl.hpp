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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_DROPNEGATIVEENTRIESFACTORY_DECL_HPP
#define MUELU_DROPNEGATIVEENTRIESFACTORY_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_DropNegativeEntriesFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"

namespace MueLu {

/*!
  @class DropNegativeEntriesFactory class.
  @brief Application-specific filtering for A. Can be used in context of graph coarsening and aggregation.

  This factory drops all negative entries (or entries with a magnitude < 0). Only weak positive connections are kept.
  Do not use this kind of filtering for regular PDEs unless you have very good reasons.
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class DropNegativeEntriesFactory : public SingleLevelFactoryBase {
#undef MUELU_DROPNEGATIVEENTRIESFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  DropNegativeEntriesFactory() {}

  //! Destructor.
  virtual ~DropNegativeEntriesFactory() {}

  RCP<const ParameterList> GetValidParameterList() const;

  //@}

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const;

  //@}

  //! @name Build methods.
  //@{

  /*!
    @brief Build method.

    Builds filtered matrix and returns it in <tt>currentLevel</tt>.
    */
  void Build(Level& currentLevel) const;

  //@}

};  // class DropNegativeEntriesFactory

}  // namespace MueLu

#define MUELU_DROPNEGATIVEENTRIESFACTORY_SHORT
#endif  // MUELU_DROPNEGATIVEENTRIESFACTORY_DECL_HPP
