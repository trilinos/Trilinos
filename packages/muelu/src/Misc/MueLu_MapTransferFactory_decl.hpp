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
/*
 * MueLu_ContactMapTransferFactory_decl.hpp
 *
 *  Created on: Aug 2, 2012
 *      Author: wiesner
 */

#ifndef MUELU_MAPTRANSFERFACTORY_DECL_HPP_
#define MUELU_MAPTRANSFERFACTORY_DECL_HPP_

#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>


#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"

namespace MueLu {

  /*!
    @class MapTransferFactory class.
    @brief transfer factory for maps

    Factory that transfers a map (given by a variable name and a generating factory) for building
    a coarse version of the map. The coarse map is stored on the coarse level using the same variable name
    and generating factory than the original fine level map.

  */

  template <class Scalar = DefaultScalar,
            class LocalOrdinal = DefaultLocalOrdinal,
            class GlobalOrdinal = DefaultGlobalOrdinal,
            class Node = DefaultNode>
  class MapTransferFactory : public TwoLevelFactoryBase {
#undef MUELU_MAPTRANSFERFACTORY_SHORT
    #include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    MapTransferFactory();

    //! Destructor.
    virtual ~MapTransferFactory() {}

    RCP<const ParameterList> GetValidParameterList() const;
    //@}

    //! Input
    //@{

    void DeclareInput(Level& fineLevel, Level& coarseLevel) const;

    //@}

    //@{
    //! @name Build methods.

    //! Build an object with this factory.
    void Build(Level& fineLevel, Level& coarseLevel) const;

    //@}

  private:

    //std::string              mapName_;      ///< name of input and output variable
    mutable RCP<const FactoryBase>   mapFact_;      ///< generating factory of input variable

  }; // class MapTransferFactory

} // namespace MueLu

#define MUELU_MAPTRANSFERFACTORY_SHORT
#endif /* MUELU_MAPTRANSFERFACTORY_DECL_HPP_ */
