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
/*
 * MueLu_BlockedCoarseMapFactory_decl.hpp
 *
 *  Created on: Oct 16, 2012
 *      Author: tobias
 */

#ifndef MUELU_BLOCKEDCOARSEMAPFACTORY_DECL_HPP_
#define MUELU_BLOCKEDCOARSEMAPFACTORY_DECL_HPP_
#include "Xpetra_StridedMapFactory_fwd.hpp"

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_CoarseMapFactory.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

/*!
     @class BlockedCoarseMapFactory class.
     @brief Factory for generating coarse level map. Used by BlockedPFactory.

      overloads CoarseMapFactory. Uses a CoarseMapFactory as input parameter in the constructor and
      automatically calculates the domain offset using the max gid from the given coarse map.

     @ingroup MueLuTransferClasses
 */

template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
class BlockedCoarseMapFactory : public MueLu::CoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> { //SingleLevelFactoryBase {
#undef MUELU_BLOCKEDCOARSEMAPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

public:

  //! @name Constructors/Destructors.
  //@{

  //! Constructor
  BlockedCoarseMapFactory();

  //! Destructor
  virtual ~BlockedCoarseMapFactory();

  RCP<const ParameterList> GetValidParameterList(const ParameterList& paramList = ParameterList()) const;

  //@}

  //! @name Input
  //@{
  /*! @brief Specifies the data that this class needs, and the factories that generate that data.

        If the Build method of this class requires some data, but the generating factory is not specified in DeclareInput, then this class
        will fall back to the settings in FactoryManager.
   */

  void DeclareInput(Level &currentLevel) const;

  //@}

  //! @name Build methods.
  //@{

  //! Build an object with this factory.
  void Build(Level &currentLevel) const;

  //@}

  //! @name ovrloaded Get/Set functions

  /*! @brief setDomainMapOffset
     sets offset for domain map DOF Gids in tentative prolongation operator.
     offset must not be smaller than zero. Note: Direct solvers (Amesos/Amesos2) are not properly working with offset > 0.
   */
  void setDomainMapOffset(GlobalOrdinal offset);

  /*! @brief getDomainMapOffset
   * returns offset of the domain DOF map (=coarse map).
   */
  GlobalOrdinal getDomainMapOffset() const;

  //@}

private:

}; //class BlockedCoarseMapFactory

} //namespace MueLu

#define MUELU_BLOCKEDCOARSEMAPFACTORY_SHORT
#endif /* MUELU_BLOCKEDCOARSEMAPFACTORY_DECL_HPP_ */
