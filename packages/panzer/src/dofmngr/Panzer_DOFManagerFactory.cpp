#include "Panzer_config.hpp"

#ifdef HAVE_PANZER_EXPLICIT_INSTANTIATION

#include "Panzer_DOFManagerFactory_decl.hpp"
#include "Panzer_DOFManagerFactory_impl.hpp"

template class panzer::DOFManagerFactory<int,long int>;
template class panzer::DOFManagerFactory<int,int>;
template class panzer::DOFManagerFactory<short,int>;
template class panzer::DOFManagerFactory<char,long int>;

#endif
