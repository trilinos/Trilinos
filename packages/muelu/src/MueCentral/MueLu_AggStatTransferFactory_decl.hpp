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
 * MueLu_VariableTransferFactory_decl.hpp
 *
 *  Created on: Jul 30, 2012
 *      Author: wiesner
 */

#ifndef MUELU_AGGSTATTRANSFERFACTORY_DECL_HPP_
#define MUELU_AGGSTATTRANSFERFACTORY_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_AggStatTransferFactory_fwd.hpp"

namespace MueLu {

  /*!
    @class AggStatTransferFactory class.
    @brief Simple class for transferring aggregation status information to next coarser level

    This is to be used in conjunction with Muelu::RAPFactory::AddTransferFactory().
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class AggStatTransferFactory : public TwoLevelFactoryBase {
#undef MUELU_AGGSTATTRANSFERFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    /*! @brief Constructor.

       @param varName The name of the quantity to be restricted.
       @param genFact The generating factory
    */
    AggStatTransferFactory(std::string const & varName, Teuchos::RCP<const FactoryBase> const &genFact = Teuchos::null);

    //! Destructor.
    virtual ~AggStatTransferFactory();

    //@}

    //! @name Input
    //@{

    /*! @brief Specifies the data that this class needs, and the factories that generate that data.

        If the Build method of this class requires some data, but the generating factory is not specified in DeclareInput, then this class
        will fall back to the settings in FactoryManager.
    */
    void DeclareInput(Level &finelevel, Level &coarseLevel) const;

    //@}

    //! @name Build methods.
    //@{

    //! Build an object with this factory.
    void Build(Level & fineLevel, Level &coarseLevel) const;

    //@}

  private:
    //! name of variable to be transfered.
    std::string varName_;

    //! factory that generates variable
    Teuchos::RCP<const FactoryBase> genFact_;
  }; // class MultiVectorTransferFactory

} // namespace MueLu

#define MUELU_AGGSTATTRANSFERFACTORY_SHORT
#endif /* MUELU_AGGSTATTRANSFERFACTORY_DECL_HPP_ */
