#include "Panzer_config.hpp"

#ifdef PANZER_EXPLICIT_TEMPLATE_INSTANTIATION

#include "Panzer_Traits.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_EpetraLinearObjFactoryT.hpp"

template class panzer::EpetraLinearObjFactory<panzer::Traits,int>;
template class panzer::EpetraLinearObjFactory<panzer::Traits,short>;

#endif
