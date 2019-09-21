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
#ifndef MUELU_REBALANCEACFACTORY_DECL_HPP
#define MUELU_REBALANCEACFACTORY_DECL_HPP

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_CrsMatrix_fwd.hpp>
#include <Xpetra_CrsMatrixWrap_fwd.hpp>
#include <Xpetra_MatrixFactory_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_RebalanceAcFactory_fwd.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_RAPFactory_fwd.hpp"
#include "MueLu_PerfUtils_fwd.hpp"

namespace MueLu {
  /*!
    @class RebalanceAcFactory
    @brief Factory for building coarse matrices.
  */
  template <class Scalar = DefaultScalar,
            class LocalOrdinal = DefaultLocalOrdinal,
            class GlobalOrdinal = DefaultGlobalOrdinal,
            class Node = DefaultNode>
  class RebalanceAcFactory : public TwoLevelFactoryBase {
#undef MUELU_REBALANCEACFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    RebalanceAcFactory() { }

    virtual ~RebalanceAcFactory() { }

    RCP<const ParameterList> GetValidParameterList() const;
    //@}

    //! @name Input
    //@{

    void DeclareInput(Level &fineLevel, Level &coarseLevel) const;

    //@}

    //! @name Build methods.
    //@{
    void Build(Level &fineLevel, Level &coarseLevel) const;
    //@}

    //@{
    /*! @brief Add rebalancing factory in the end of list of rebalancing factories in RebalanceAcFactory.

    Rebalancing factories are derived from SingleLevelFactoryBase and rebalance the underlaying object
    (e.g. map, vector,...) to fit to the rebalanced maps.
    */
    void AddRebalanceFactory(const RCP<const FactoryBase>& factory);

    //! Returns number of transfer factories.
    size_t NumRebalanceFactories() const { return rebalanceFacts_.size(); }

    //@}

  private:
    //! list of user-defined rebalancing Factories
    std::vector<RCP<const FactoryBase> > rebalanceFacts_;

  }; //class RebalanceAcFactory

} //namespace MueLu

#define MUELU_REBALANCEACFACTORY_SHORT
#endif // MUELU_REBALANCEACFACTORY_DECL_HPP
