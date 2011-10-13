#include "Panzer_config.hpp"

// build it just for fun
#ifdef PANZER_EXPLICIT_TEMPLATE_INSTANTIATION

#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_BlockedDOFManagerT.hpp"

template class panzer::BlockedDOFManager<int,long int>;
template class panzer::BlockedDOFManager<int,int>;
template class panzer::BlockedDOFManager<short,int>;
template class panzer::BlockedDOFManager<char,long int>;

#endif
