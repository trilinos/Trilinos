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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef PACKAGES_MUELU_SRC_GRAPH_MUELU_VARIABLEDOFLAPLACIANFACTORY_DECL_HPP_
#define PACKAGES_MUELU_SRC_GRAPH_MUELU_VARIABLEDOFLAPLACIANFACTORY_DECL_HPP_


#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_VariableDofLaplacianFactory_fwd.hpp"
#include "MueLu_Level_fwd.hpp"

#include "MueLu_Utilities_fwd.hpp"

namespace MueLu {

  template <class Scalar = double,
            class LocalOrdinal = int,
            class GlobalOrdinal = LocalOrdinal,
            class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class VariableDofLaplacianFactory : public SingleLevelFactoryBase {
#undef MUELU_VARIABLEDOFLAPLACIANFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors/Destructors.
    //@{

    //! Constructor
    VariableDofLaplacianFactory();

    //! Destructor
    virtual ~VariableDofLaplacianFactory() { }

    RCP<const ParameterList> GetValidParameterList() const;

    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const;

    //@}

    void Build(Level &currentLevel) const; // Build

  private:

    void buildPaddedMap(const Teuchos::ArrayRCP<const bool> & dofPresent, std::vector<LocalOrdinal> & map, size_t nDofs) const;
    void assignGhostLocalNodeIds(const Teuchos::RCP<const Map> & rowDofMap, const Teuchos::RCP<const Map> & colDofMap, std::vector<LocalOrdinal> & myLocalNodeIds, const std::vector<LocalOrdinal> & dofMap, size_t maxDofPerNode, size_t& nLocalNodes, size_t& nLocalPlusGhostNodes, Teuchos::RCP< const Teuchos::Comm< int > > comm) const;

    void MueLu_az_sort(size_t list[], size_t N, size_t list2[], Scalar list3[]) const;

  }; //class CoalesceDropFactory

} //namespace MueLu

#define MUELU_VARIABLEDOFLAPLACIANFACTORY_SHORT


#endif /* PACKAGES_MUELU_SRC_GRAPH_MUELU_VARIABLEDOFLAPLACIANFACTORY_DECL_HPP_ */
