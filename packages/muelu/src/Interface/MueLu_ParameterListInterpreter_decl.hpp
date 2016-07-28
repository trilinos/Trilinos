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

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_Operator_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_HierarchyManager.hpp"

#include "MueLu_AggregationExportFactory_fwd.hpp"
#include "MueLu_BrickAggregationFactory_fwd.hpp"
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
#include "MueLu_LineDetectionFactory_fwd.hpp"
#include "MueLu_NullspaceFactory_fwd.hpp"
#include "MueLu_PatternFactory_fwd.hpp"
#include "MueLu_PgPFactory_fwd.hpp"
#include "MueLu_RAPFactory_fwd.hpp"
#include "MueLu_RebalanceAcFactory_fwd.hpp"
#include "MueLu_RebalanceTransferFactory_fwd.hpp"
#include "MueLu_RepartitionFactory_fwd.hpp"
#include "MueLu_SaPFactory_fwd.hpp"
#include "MueLu_SemiCoarsenPFactory_fwd.hpp"
#include "MueLu_SmootherFactory_fwd.hpp"
#include "MueLu_TentativePFactory_fwd.hpp"
#include "MueLu_TogglePFactory_fwd.hpp"
#include "MueLu_ToggleCoordinatesTransferFactory_fwd.hpp"
#include "MueLu_TransPFactory_fwd.hpp"
#include "MueLu_UncoupledAggregationFactory_fwd.hpp"
#include "MueLu_ZoltanInterface_fwd.hpp"
#include "MueLu_Zoltan2Interface_fwd.hpp"

#ifdef HAVE_MUELU_KOKKOS_REFACTOR
#include "MueLu_CoalesceDropFactory_kokkos_fwd.hpp"
#include "MueLu_CoarseMapFactory_kokkos_fwd.hpp"
#include "MueLu_CoordinatesTransferFactory_kokkos_fwd.hpp"
#include "MueLu_FilteredAFactory_kokkos_fwd.hpp"
#include "MueLu_NullspaceFactory_kokkos_fwd.hpp"
#include "MueLu_SaPFactory_kokkos_fwd.hpp"
#include "MueLu_TentativePFactory_kokkos_fwd.hpp"
#include "MueLu_UncoupledAggregationFactory_kokkos_fwd.hpp"
#endif

namespace MueLu {

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class ParameterListInterpreter :
    public HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_PARAMETERLISTINTERPRETER_SHORT
#include "MueLu_UseShortNames.hpp"
    typedef std::pair<std::string, const FactoryBase*> keep_pair;

  public:
    //! @name Constructors/Destructors
    //@{

  protected:
    /*! @brief Empty constructor
     *
     *  Constructor for derived classes
     */
    ParameterListInterpreter() {
      factFact_ = Teuchos::null;
    }

  public:
    /*! @brief Constructor that accepts a user-provided ParameterList.

        Constructor for parameter list interpreter which directly interprets Teuchos::ParameterLists

        @details The parameter list can be either in the easy parameter list format or in the factory driven parameter list format.

        @param[in] paramList (Teuchos::ParameterList): ParameterList containing the MueLu parameters
        @param[in] comm  (RCP<Teuchos::Comm<int> >): Optional RCP of a Teuchos communicator  (default: Teuchos::null)
        @param[in] factFact  (RCP<FactoryFactory>): Optional parameter allowing to define user-specific factory interpreters for user-specific extensions of the XML interface. (default: Teuchos::null)

     */
    ParameterListInterpreter(Teuchos::ParameterList& paramList, Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::null, Teuchos::RCP<FactoryFactory> factFact = Teuchos::null);

    /*! @brief Constructor that reads parameters from an XML file.

        XML options are converted to ParameterList entries by Teuchos.

        @param[in] xmlFileName (std::string): XML file to read
        @param[in] comm  (Teuchos::Comm<int>): Teuchos communicator
        @param[in] factFact  (RCP<FactoryFactory>): Optional parameter allowing to define user-specific factory interpreters for user-specific extensions of the XML interface. (default: Teuchos::null)

    */
    ParameterListInterpreter(const std::string& xmlFileName, const Teuchos::Comm<int>& comm, Teuchos::RCP<FactoryFactory> factFact = Teuchos::null);

    //! Destructor.
    virtual ~ParameterListInterpreter() { }

    //@}

    /*! @brief Set parameter list for Parameter list interpreter.

       The routine checks whether it is a parameter list in the easy parameter format or the more advanced factory-based parameter format and calls the corresponding interpreter routine.

       When finished, the parameter list is set that will used by the hierarchy build phase.

       This method includes validation and some pre-parsing of the list for:
           - verbosity level
           - data to export
           - cycle type
           - max coarse size
           - max levels
           - number of equations

       @param[in] paramList: ParameterList containing the MueLu parameters.
    */
    void SetParameterList(const Teuchos::ParameterList& paramList);

    //! Call the SetupHierarchy routine from the HiearchyManager object.
    void SetupHierarchy(Hierarchy& H) const;

  private:
    //! Setup Operator object
    virtual void SetupOperator(Operator& A) const;

    int       blockSize_;     ///< block size of matrix (fixed block size)
    CycleType Cycle_;         ///< multigrid cycle type (V-cycle or W-cycle)
    GlobalOrdinal dofOffset_; ///< global offset variable describing offset of DOFs in operator

    //! Easy interpreter stuff
    //@{
    // These two variables are only needed to print out proper [default]
    bool changedPRrebalance_;
    bool changedImplicitTranspose_;

    void SetEasyParameterList(const Teuchos::ParameterList& paramList);
    void Validate(const Teuchos::ParameterList& paramList) const;

    void UpdateFactoryManager(Teuchos::ParameterList& paramList, const Teuchos::ParameterList& defaultList, FactoryManager& manager,
                              int levelID, std::vector<keep_pair>& keeps) const;

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

    //! Internal factory for factories
    Teuchos::RCP<FactoryFactory> factFact_;
    //@}
  };

} // namespace MueLu

#define MUELU_PARAMETERLISTINTERPRETER_SHORT
#endif /* MUELU_PARAMETERLISTINTERPRETER_DECL_HPP */
