#include "Panzer_config.hpp"

#ifdef HAVE_PANZER_EXPLICIT_INSTANTIATION

#include "Panzer_BlockedDOFManager_decl.hpp"
#include "Panzer_BlockedDOFManager_impl.hpp"

template class panzer::BlockedDOFManager<int,long int>;
template class panzer::BlockedDOFManager<int,int>;
template class panzer::BlockedDOFManager<short,int>;
template class panzer::BlockedDOFManager<char,long int>;

#endif
