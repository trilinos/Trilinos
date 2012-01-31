#include "Panzer_config.hpp"

#ifdef HAVE_PANZER_EXPLICIT_INSTANTIATION

#include "Panzer_Traits.hpp"
#include "Panzer_EpetraLinearObjFactory_decl.hpp"
#include "Panzer_EpetraLinearObjFactory_impl.hpp"

template class panzer::EpetraLinearObjFactory<panzer::Traits,int>;
template class panzer::EpetraLinearObjFactory<panzer::Traits,short>;

#endif
