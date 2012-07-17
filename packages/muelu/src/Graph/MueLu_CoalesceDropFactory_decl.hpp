#ifndef MUELU_COALESCEDROPFACTORY_DECL_HPP
#define MUELU_COALESCEDROPFACTORY_DECL_HPP

#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_CrsGraph_fwd.hpp>
#include <Xpetra_CrsGraphFactory.hpp> //TODO
#include <Xpetra_Operator_fwd.hpp>
#include <Xpetra_StridedMap_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_CoalesceDropFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_Graph_fwd.hpp"
#include "MueLu_AmalgamationInfo_fwd.hpp"
#include "MueLu_SubBlockUnAmalgamationFactory_fwd.hpp"
#include "MueLu_PreDropFunctionBaseClass_fwd.hpp"

namespace MueLu {

  /*!
    @class CoalesceDropFactory
    @brief Factory for creating a graph base on a given matrix.

    Factory for creating graphs from matrices with entries selectively dropped.
  
    - TODO This factory is very incomplete.
    - TODO The Build method simply builds the matrix graph with no dropping.
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class CoalesceDropFactory : public SingleLevelFactoryBase {
#undef MUELU_COALESCEDROPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors/Destructors.
    //@{

    //! Constructor
    CoalesceDropFactory(RCP<const FactoryBase> AFact = Teuchos::null);

    //! Destructor
    virtual ~CoalesceDropFactory() { }

    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const;

    /// set predrop function
    void SetPreDropFunction(const RCP<MueLu::PreDropFunctionBaseClass<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &predrop);

    //@}

    void Build(Level &currentLevel) const; // Build



  private:

    //! A Factory
    RCP<const FactoryBase> AFact_;

    /// pre-drop function
    RCP<PreDropFunctionBaseClass> predrop_;

    /// SubBlockUnAmalgamationFactory (for generating the UnAmalgamation information)
    RCP<const FactoryBase> UnAmalgamationFact_;



  }; //class CoalesceDropFactory

} //namespace MueLu

#define MUELU_COALESCEDROPFACTORY_SHORT
#endif // MUELU_COALESCEDROPFACTORY_DECL_HPP
