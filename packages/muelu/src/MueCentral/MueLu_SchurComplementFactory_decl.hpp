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
#include <Xpetra_CrsOperator.hpp>
#include <Xpetra_Operator.hpp>
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

//TODO Update this!
    Example
    \code
    Teuchos::RCP<Xpetra::BlockedCrsOperator<Scalar,LO,GO,Node> > bOp = Teuchos::rcp(new Xpetra::BlockedCrsOperator<Scalar,LO,GO>(mapExtractor,mapExtractor,10));
    // ... let bOp be a 2x2 blocked operator ...
    bOp->fillComplete();

    // define factory for accessing block (0,0) in blocked operator A (assuming that the blocked operator is stored in Level class with NoFactory as generating factory)
    RCP<SubBlockAFactory> A11Fact = Teuchos::rcp(new SubBlockAFactory(MueLu::NoFactory::getRCP(), 0, 0));

    // define factory for accessing block (1,1) in blocked operator A
    RCP<SubBlockAFactory> A22Fact = Teuchos::rcp(new SubBlockAFactory(MueLu::NoFactory::getRCP(), 1, 1));

    RCP<Operator> A11 = level.Get<RCP<Operator> >("A", A11Fact); // extract (0,0) block from blocked operator A
    RCP<Operator> A22 = level.Get<RCP<Operator> >("A", A22Fact); // extract (1,1) block from blocked operator A
    \endcode
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class SchurComplementFactory : public SingleLevelFactoryBase {
#undef MUELU_SCHURCOMPLEMENTFACTORY_SHORT
    #include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    SchurComplementFactory(Teuchos::RCP<const FactoryBase> Afact,Scalar omega/*, size_t row, size_t col, LocalOrdinal blksize = 1*/);

    //! Destructor.
    virtual ~SchurComplementFactory();
    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const;

    //! Add a factory manager
    void AddFactoryManager(RCP<const FactoryManagerBase> FactManager);

    //@}

    //@{
    //! @name Build methods.

    //! Build an object with this factory.
    void Build(Level & currentLevel) const;

    //@}


  private:
    RCP<const FactoryManagerBase>       FactManager_;           //!< Factory manager for setting the schur complement
    Teuchos::RCP<const FactoryBase> AFact_;   ///< generating factory of input variable
    //const size_t                    row_;     ///< row id
    //const size_t                    col_;     ///< column id

    Scalar                                                      omega_;         ///< damping parameter

    bool switchRowOrder_;                                               ///< NOT USED YET

  }; // class SchurComplementFactory

} // namespace MueLu

#define MUELU_SCHURCOMPLEMENTFACTORY_SHORT
#endif /* MUELU_SCHURCOMPLEMENTFACTORY_DECL_HPP_ */
