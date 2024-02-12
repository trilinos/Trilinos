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

#ifndef MUELU_ZEROSUBBLOCKAFACTORY_DECL_HPP_
#define MUELU_ZEROSUBBLOCKAFACTORY_DECL_HPP_

#include <Xpetra_Map_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_ZeroSubBlockAFactory_fwd.hpp"

namespace MueLu {

/*!
  @class ZeroSubBlockAFactory class.
  @brief Factory for extracting a zero block from a BlockedCrsMatrix.

  This is a very simple class to access a single matrix block in a blocked operator A,
  where the matrix block of interest is an actual zero block (as e.g. in saddle point systems).

*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class ZeroSubBlockAFactory : public SingleLevelFactoryBase {
#undef MUELU_ZEROSUBBLOCKAFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! Input
  //@{

  RCP<const ParameterList> GetValidParameterList() const override;

  void DeclareInput(Level &currentLevel) const override;

  //@}

  //@{
  //! @name Build methods

  /*!
    @brief Build a zero sub-block object with this factory

    Create a zero sub block matrix, that fits into a given blocked crs operator.

    Strided or block information is extracted in the following way:
      1) first check whether the corresponding sub maps are strided
         If yes, use the fixed block size and strided block id
      2) If no, get the full map of the map extractors.
         Check whether the full map is strided.
         If yes, use the strided information of the full map and build
         partial (strided) maps with it
         If no, throw an exception

    For blocked operators with block maps one should use the striding
    information from the sub maps. For strided operators, the striding
    information of the full map is the best choice.
  */
  void Build(Level &currentLevel) const override;

  //@}

};  // class ZeroSubBlockAFactory

}  // namespace MueLu

#define MUELU_ZEROSUBBLOCKAFACTORY_SHORT
#endif /* MUELU_ZEROSUBBLOCKAFACTORY_DECL_HPP_ */
