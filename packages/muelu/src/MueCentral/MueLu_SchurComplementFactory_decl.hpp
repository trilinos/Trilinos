/*
 * MueLu_SchurComplementFactory_decl.hpp
 *
 *  Created on: Jun 18, 2012
 *      Author: wiesner
 */

#ifndef MUELU_SCHURCOMPLEMENTFACTORY_DECL_HPP_
#define MUELU_SCHURCOMPLEMENTFACTORY_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"

#include <Teuchos_ParameterList.hpp>


//Xpetra
#include <Xpetra_MapExtractor_fwd.hpp>
#include <Xpetra_CrsOperator_fwd.hpp>
#include <Xpetra_CrsMatrix_fwd.hpp>
//#include <Xpetra_Operator.hpp>
#include "Xpetra_Map_fwd.hpp"
#include "Xpetra_StridedMap_fwd.hpp"


#include "MueLu_SingleLevelFactoryBase.hpp"
//#include "MueLu_SchurComplementFactory_fwd.hpp"

#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_FactoryManagerBase_fwd.hpp"
#include "MueLu_SmootherBase_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"


namespace MueLu {

  /*!
    @class SchurComplementFactory class.
    @brief Factory for building the Schur Complement for a 2x2 block matrix.

    For a blocked matrix A = (F, G; D -Z) it computes the SchurComplement
    S = - 1/omega D \hat{F}^{-1} G + Z
    where omega is some scaling factor and \hat{F} an approximation of F (just the diagonal of F)
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class SchurComplementFactory : public SingleLevelFactoryBase {
#undef MUELU_SCHURCOMPLEMENTFACTORY_SHORT
    #include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    SchurComplementFactory(Teuchos::RCP<const FactoryBase> AFact,Scalar omega/*, size_t row, size_t col, LocalOrdinal blksize = 1*/);

    //! Destructor.
    virtual ~SchurComplementFactory();
    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const;

    //! Add a factory manager
    //void AddFactoryManager(RCP<const FactoryManagerBase> FactManager);

    //@}

    //@{
    //! @name Build methods.

    //! Build an object with this factory.
    void Build(Level & currentLevel) const;

    //@}


  private:
    Teuchos::RCP<const FactoryBase>     AFact_;                 ///< generating factory of input variable (blocked A operator)

    Scalar                              omega_;         ///< damping parameter

    bool switchRowOrder_;                                               ///< NOT USED YET

  }; // class SchurComplementFactory

} // namespace MueLu

#define MUELU_SCHURCOMPLEMENTFACTORY_SHORT
#endif /* MUELU_SCHURCOMPLEMENTFACTORY_DECL_HPP_ */
