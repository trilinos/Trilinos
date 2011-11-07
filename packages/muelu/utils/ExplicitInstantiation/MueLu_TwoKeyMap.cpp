#include "MueLu_ExplicitInstantiation.hpp"

#include "MueLu_TwoKeyMap_def.hpp"

#include "MueLu_FactoryBase.hpp"
#include <string>

template class MueLu::UTILS::TwoKeyMap<std::string, const MueLu::FactoryBase*, int>;
template class MueLu::UTILS::TwoKeyMap<std::string, const MueLu::FactoryBase*, Teuchos::ParameterEntry>;
template class MueLu::UTILS::TwoKeyMap<std::string, const MueLu::FactoryBase*, bool>;
