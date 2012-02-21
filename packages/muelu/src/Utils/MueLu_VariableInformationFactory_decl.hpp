/*
 * MueLu_VariableInformationFactory_decl.hpp
 *
 *  Created on: 20.02.2012
 *      Author: tobias
 */

#ifndef MUELU_VARIABLEINFORMATIONFACTORY_DECL_HPP_
#define MUELU_VARIABLEINFORMATIONFACTORY_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_ThresholdAFilterFactory_fwd.hpp"

namespace MueLu {

  /*!
    @class VariableInformationFactory class.
    @brief Factory collects some information about given variable from Level

    VariableInformationFactory derives from TwoLevelFactoryBase such that it can be used by the RAPFactory class.
    It only uses the information from the fine level.
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class VariableInformationFactory : public TwoLevelFactoryBase {
#undef MUELU_VARIABLEINFORMATIONFACTORY_SHORT
    #include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    VariableInformationFactory(const std::string& ename, const FactoryBase* fac, bool bCoarseLevelInfo = false);

    //! Destructor.
    virtual ~VariableInformationFactory();
    //@}

    //! Input
    //@{

    void DeclareInput(Level &fineLevel, Level &coarseLevel) const;

    //@}

    //@{
    //! @name Build methods.

    //! Build an object with this factory.
    void Build(Level &fineLevel, Level &coarseLevel) const;

    //@}

  private:
    void Gershgorin(const RCP<Operator>& A, const RCP<Operator>& At = Teuchos::null) const;

    std::string        varName_;   ///< name of input and output variable
    const FactoryBase* factory_;   ///< generating factory of input variable
    bool               coarseLevel_; ///< if true, use coarseLevel information (otherwise fineLevel = default)

  }; // class VariableInformationFactory

} // namespace MueLu

#define MUELU_VARIABLEINFORMATIONFACTORY_SHORT
#endif /* MUELU_VARIABLEINFORMATIONFACTORY_DECL_HPP_ */
