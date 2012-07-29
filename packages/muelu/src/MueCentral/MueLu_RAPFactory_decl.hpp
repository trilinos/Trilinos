#ifndef MUELU_RAPFACTORY_DECL_HPP
#define MUELU_RAPFACTORY_DECL_HPP

#include <Xpetra_Operator_fwd.hpp>
#include <Xpetra_CrsOperator_fwd.hpp>
#include <Xpetra_OperatorFactory_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_RAPFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"

// MPI helper
#define sumAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out));

namespace MueLu {
  /*!
    @class RAPFactory
    @brief Factory for building coarse matrices.
  */
  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void, LocalOrdinal, Node>::SparseOps>
  class RAPFactory : public TwoLevelFactoryBase {
#undef MUELU_RAPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    RAPFactory(RCP<const FactoryBase> PFact = Teuchos::null, RCP<const FactoryBase> RFact = Teuchos::null, RCP<const FactoryBase> AFact = Teuchos::null);

    virtual ~RAPFactory();
    //@}

    //! @name Input
    //@{

    void DeclareInput(Level &fineLevel, Level &coarseLevel) const;

    //@}

    //! @name Build methods.
    //@{
    void Build(Level &fineLevel, Level &coarseLevel) const;
    //@}

    //! @name Handling of user-defined transfer factories
    //@{

    //! Indicate that the restriction operator action should be implicitly defined by the transpose of the prolongator.
    void SetImplicitTranspose(bool const &implicit);
    
    //! Indicate that zero entries on the diagonal of Ac shall be repaired (i.e. if A(i,i) == 0.0 set A(i,i) = 1.0)
    void SetRepairZeroDiagonal(bool const &repair) { 
      repairZeroDiagonals_ = repair; 
      if(repair) checkAc_ = true; // make sure that plausibility check is performed. Otherwise SetRepairZeroDiagonal(true) has no effect.
    }
    
    //! Indicate that a simple plausibility check shall be done for Ac after building RAP
    void SetPlausibilityCheck(bool const &check) {
      checkAc_ = check; 
    }

    /*! @brief Add transfer factory in the end of list of transfer factories in RAPFactory.

    Transfer factories are derived from TwoLevelFactoryBase and project some data from the fine level to
    the next coarser level.
    */
    void AddTransferFactory(const RCP<const FactoryBase>& factory);

    // TODO add a function to remove a specific transfer factory?

    //! Returns number of transfer factories.
    size_t NumTransferFactories() const;

    //! Returns factory that generates P.
    RCP<const FactoryBase> GetPFactory() const;

    //! Returns factory that generates R.
    RCP<const FactoryBase> GetRFactory() const;

    //@}
  private:

    //! @name internal Build methods.
    //@{
    RCP<Operator> BuildRAPExplicit(const RCP<Operator>& R, const RCP<Operator>& A, const RCP<Operator>& P, int levelID) const; //TODO: remove RCP for input args. TODO: static
    RCP<Operator> BuildRAPImplicit(const RCP<Operator>& A, const RCP<Operator>& P, int levelID) const;
    RCP<Operator> BuildRAPBlock   (const RCP<Operator>& R, const RCP<Operator>& A, const RCP<Operator>& P, int levelID) const;
    //@}

    //! @name internal print methods.
    void PrintMatrixInfo(const Operator & Ac, const std::string & msgTag) const; //TODO: static
    void PrintLoadBalancingInfo(const Operator & Ac, const std::string & msgTag) const;
    //@}
    
    //! @name internal plausibility check methods
    void CheckMainDiagonal(Operator & Ac) const;
    //@}

    //! P Factory
    RCP<const FactoryBase> PFact_;

    //! R Factory
    RCP<const FactoryBase> RFact_;
  
    //! A Factory
    RCP<const FactoryBase> AFact_;

    //! list of user-defined transfer Factories
    std::vector<RCP<const FactoryBase> > TransferFacts_;

    //! If true, the action of the restriction operator action is implicitly defined by the transpose of the prolongator.
    bool implicitTranspose_;
    
    //! If true, perform a basic plausibility check on Ac (default = false)
    //! note, that the repairZeroDiagonals_ flag only is valid for checkAc_ == true
    bool checkAc_;
    
    //! If true, the CheckMainDiagonal routine automatically repairs zero entries on main diagonal (default = false)
    //! i.e. if A(i,i) == 0.0 set A(i,i) = 1.0
    //! note, that the repairZeroDiagonals_ flag only is valid for checkAc_ == true
    bool repairZeroDiagonals_;
    
  }; //class RAPFactory

} //namespace MueLu

#define MUELU_RAPFACTORY_SHORT
#endif // MUELU_RAPFACTORY_DECL_HPP
