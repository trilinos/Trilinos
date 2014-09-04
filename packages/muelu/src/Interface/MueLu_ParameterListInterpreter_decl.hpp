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
#ifndef MUELU_PARAMETERLISTINTERPRETER_DECL_HPP
#define MUELU_PARAMETERLISTINTERPRETER_DECL_HPP

#include <Teuchos_ParameterList.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_HierarchyManager.hpp"

#include "MueLu_AggregationExportFactory_fwd.hpp"
#include "MueLu_CoalesceDropFactory_fwd.hpp"
#include "MueLu_CoarseMapFactory_fwd.hpp"
#include "MueLu_ConstraintFactory_fwd.hpp"
#include "MueLu_CoordinatesTransferFactory_fwd.hpp"
#include "MueLu_CoupledAggregationFactory_fwd.hpp"
#include "MueLu_DirectSolver_fwd.hpp"
#include "MueLu_EminPFactory_fwd.hpp"
#include "MueLu_FactoryFactory_fwd.hpp"
#include "MueLu_FilteredAFactory_fwd.hpp"
#include "MueLu_GenericRFactory_fwd.hpp"
#include "MueLu_NullspaceFactory_fwd.hpp"
#include "MueLu_PatternFactory_fwd.hpp"
#include "MueLu_PgPFactory_fwd.hpp"
#include "MueLu_RAPFactory_fwd.hpp"
#include "MueLu_RebalanceAcFactory_fwd.hpp"
#include "MueLu_RebalanceTransferFactory_fwd.hpp"
#include "MueLu_RepartitionFactory_fwd.hpp"
#include "MueLu_SaPFactory_fwd.hpp"
#include "MueLu_SmootherFactory_fwd.hpp"
#include "MueLu_TentativePFactory_fwd.hpp"
#include "MueLu_TransPFactory_fwd.hpp"
#include "MueLu_UncoupledAggregationFactory_fwd.hpp"
#include "MueLu_ZoltanInterface_fwd.hpp"
#include "MueLu_Zoltan2Interface_fwd.hpp"

namespace MueLu {

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = KokkosClassic::DefaultNode::DefaultNodeType, class LocalMatOps = void>
  class ParameterListInterpreter :
    public HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_PARAMETERLISTINTERPRETER_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors
    //@{

    //! Constructor that accepts a user-provided ParameterList.
    ParameterListInterpreter(Teuchos::ParameterList& paramList);

    /*! @brief Constructor that reads parameters from an XML file.

        XML options are converted to ParameterList entries by Teuchos.
    */
    ParameterListInterpreter(const std::string& xmlFileName, const Teuchos::Comm<int>& comm);

    //! Destructor.
    virtual ~ParameterListInterpreter() { }

    //@}

    void SetParameterList(const Teuchos::ParameterList& paramList);

    void SetupHierarchy(Hierarchy& H) const;

  private:
    //! Setup Matrix object
    virtual void SetupMatrix(Matrix& A) const;

    int       blockSize_;
    CycleType Cycle_;

    //! Easy interpreter stuff
    //@{
    void SetEasyParameterList(const Teuchos::ParameterList& paramList);
    void Validate(const Teuchos::ParameterList& paramList) const;

    void UpdateFactoryManager(Teuchos::ParameterList& paramList, const Teuchos::ParameterList& defaultList, FactoryManager& manager) const;

    bool useCoordinates_;
    //@}

    //! Factory interpreter stuff
    // TODO:
    // - parameter list validator
    // - SetParameterList
    // - Set/Get directly Level manager
    // - build per level
    // - comments/docs
    // - use FactoryManager instead of FactoryMap
    //@{
    void SetFactoryParameterList(const Teuchos::ParameterList& paramList);

    typedef std::map<std::string, RCP<const FactoryBase>  > FactoryMap; //TODO: remove this line
    typedef std::map<std::string, RCP<FactoryManagerBase> > FactoryManagerMap;

    void BuildFactoryMap(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, FactoryMap& factoryMapOut, FactoryManagerMap& factoryManagers) const;
    //@}
  };

} // namespace MueLu

#define MUELU_PARAMETERLISTINTERPRETER_SHORT
#endif /* MUELU_PARAMETERLISTINTERPRETER_DECL_HPP */
