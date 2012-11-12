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
#ifndef MUELU_FACTORYMANAGER_DECL_HPP
#define MUELU_FACTORYMANAGER_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_FactoryManager_fwd.hpp"
#include "MueLu_FactoryManagerBase.hpp"

#include "MueLu_TentativePFactory_fwd.hpp"
#include "MueLu_SaPFactory_fwd.hpp"
#include "MueLu_RAPFactory_fwd.hpp"
#include "MueLu_NullspaceFactory_fwd.hpp"
#include "MueLu_TransPFactory_fwd.hpp"
#include "MueLu_SmootherFactory_fwd.hpp"
#include "MueLu_TrilinosSmoother_fwd.hpp"
#include "MueLu_DirectSolver_fwd.hpp"
#include "MueLu_UCAggregationFactory_fwd.hpp"
#include "MueLu_CoalesceDropFactory_fwd.hpp"
#include "MueLu_CoalesceDropFactory2_fwd.hpp"
#include "MueLu_RepartitionFactory_fwd.hpp"
#include "MueLu_ZoltanInterface_fwd.hpp"
#include "MueLu_AmalgamationFactory_fwd.hpp"
#include "MueLu_CoarseMapFactory_fwd.hpp"
#ifdef INCLUDE_MUELU_EXPERIMENTAL
#include "MueLu_PatternFactory_fwd.hpp"
#include "MueLu_ConstraintFactory_fwd.hpp"
#endif

namespace MueLu {

  /*!
    @class FactoryManager class.
    @brief This class specifies the default factory that should generate some data on a Level if the data does not exist and
    the generating factory has not been specified.

    Consider the following example.

    @code
      RCP<SingleLevelFactory> Afact;
      Level currentLevel;
      RCP<Matrix> thisLevelA;
      thisLevelA = currentLevel.Get<Matrix>("A",Afact.get());
    @endcode

    @todo If Afact is null (actually, Teuchos::null), then the FactoryManager associated with currentLevel will determine whether a default factory has
    been specified for creating A.  If "yes", then that factory will be called, A will be stored in currentLevel, and an RCP will be returned by
    the Get call.  If "no", then the FactoryManager will <b>throw an exception indicating that it does not know how to generate A</b>.
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class FactoryManager : public FactoryManagerBase {
#undef MUELU_FACTORYMANAGER_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructor/Destructors
    //@{

    /*! @brief Constructor.

            @param[in] PFact Factory to generate the prolongation operator.
            @param[in] RFact Factory to generate the restriction operator.
            @param[in] AcFact Factory to generate the coarse grid operator A.
    */
    FactoryManager(const RCP<const FactoryBase> PFact = Teuchos::null, const RCP<const FactoryBase> RFact = Teuchos::null, const RCP<const FactoryBase> AcFact = Teuchos::null);

    //! Constructor used by HierarchyFactory (temporary, will be removed)
    FactoryManager(const std::map<std::string, RCP<const FactoryBase> >& factoryTable)
    {
      factoryTable_ = factoryTable;
      SetIgnoreUserData(false); // set IgnorUserData flag to false (default behaviour) //TODO: use parent class constructor instead
    }

    //! Destructor.
    virtual ~FactoryManager();

    //@}

    //! @name Get/Set functions.
    //@{

    /*! @brief Set Factory

        Register the factory that should generate data if said factory is not specified in the request.

        @param[in] name of variable
        @param[in] factory that generates the data
    */
    void SetFactory(const std::string & varName, const RCP<const FactoryBase> & factory);

    /*! @brief Get factory associated with a particular data name.

       @param[in] varName name of variable.

    */
    const RCP<const FactoryBase> & GetFactory(const std::string & varName) const;

    //!
    const RCP<const FactoryBase> & GetDefaultFactory(const std::string & varName) const;

    //@}

    void Clean() const;

  private:

    //! @name Helper functions
    //@{

    /*! Add a factory to the default factory list and return it. This helper function is used by GetDefaultFactory()

     @todo TODO factory->setObjectLabel("Default " + varName + "Factory");
    */

    const RCP<const FactoryBase> & SetAndReturnDefaultFactory(const std::string & varName, const RCP<const FactoryBase> & factory) const;

    //! Test if factoryTable_[varName] exists
    static bool IsAvailable(const std::string & varName, const std::map<std::string, RCP<const FactoryBase> > & factoryTable);
    //@}

    /*! @brief User-defined factories.

      Note: we distinguish 'user defined factory' and 'default factory' to allow the deallocation of default factories separately.
    */
    std::map<std::string, RCP<const FactoryBase> > factoryTable_;

    /*! @brief Table that holds default factories.

      -# We distinguish 'user defined factory' and 'default factory' to allow the deallocation of default factories separately.
      -# <tt>defaultFactoryTable_</tt> is mutable because default factories are only added to the list when they are requested
      to avoid allocation of unused factories.
    */
    mutable
    std::map<std::string, RCP<const FactoryBase> > defaultFactoryTable_;

  }; // class

} // namespace MueLu

#define MUELU_FACTORYMANAGER_SHORT
#endif // MUELU_FACTORYMANAGER_DECL_HPP
