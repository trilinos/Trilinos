#include "Panzer_config.hpp"

#ifdef HAVE_STOKHOS

#ifdef HAVE_PANZER_EXPLICIT_INSTANTIATION

#include "Panzer_Traits.hpp"
#include "Panzer_SGEpetraLinearObjFactory_decl.hpp"
#include "Panzer_SGEpetraLinearObjFactory_impl.hpp"

template class panzer::SGEpetraLinearObjFactory<panzer::Traits,int>;
template class panzer::SGEpetraLinearObjFactory<panzer::Traits,short>;

#endif

#endif
