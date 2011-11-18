#include <string>

#include <Teuchos_ParameterList.hpp>
#include <MueLu_VariableContainer.hpp>

#include "MueLu_ExplicitInstantiation.hpp"

#include "MueLu_TwoKeyMap_def.hpp"

#include "MueLu_FactoryBase_fwd.hpp"

template class MueLu::UTILS::TwoKeyMap<std::string, const MueLu::FactoryBase*, int>;
template class MueLu::UTILS::TwoKeyMap<std::string, const MueLu::FactoryBase*, Teuchos::ParameterEntry>;
template class MueLu::UTILS::TwoKeyMap<std::string, const MueLu::FactoryBase*, bool>;

// experimental
template class MueLu::UTILS::TwoKeyMap<const MueLu::FactoryBase*, std::string, std::map<const MueLu::FactoryBase*, int> >;
template class MueLu::UTILS::TwoKeyMap<const MueLu::FactoryBase*, std::string, Teuchos::RCP<MueLu::VariableContainer> >;
