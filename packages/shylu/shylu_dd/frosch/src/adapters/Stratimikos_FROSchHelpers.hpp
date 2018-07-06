

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

#include "Thyra_FROSch_TwoLevelPreconditionerFactory_decl.hpp"


#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"

#include <string>
#include "Kokkos_DefaultNode.hpp"


namespace Stratimikos {
   
    template <typename LocalOrdinal = int, typename GlobalOrdinal = int, typename Node = KokkosClassic::DefaultNode::DefaultNodeType>
    void enableFROSchTwoLevel(DefaultLinearSolverBuilder& builder, const std::string& stratName = "FROSchTwoLevel")
    {
        const Teuchos::RCP<const Teuchos::ParameterList> precValidParams = Teuchos::sublist(builder.getValidParameters(), "Preconditioner Types");
        
        TEUCHOS_TEST_FOR_EXCEPTION(precValidParams->isParameter(stratName), std::logic_error,
                                   "Stratimikos::enableFROSch_TwoLevel cannot add \"" + stratName +"\" because it is already included in builder!");
        
        typedef Thyra::PreconditionerFactoryBase<double>                                     Base;
        typedef Thyra::FROSch_TwoLevelPreconditionerFactory<double, LocalOrdinal, GlobalOrdinal, Node> Impl;
        
        builder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), stratName);
    }
    
    
    
} // namespace Stratimikos


