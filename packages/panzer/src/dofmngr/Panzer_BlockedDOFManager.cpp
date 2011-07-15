#include "Panzer_BlockedDOFManager.hpp"

// FEI includes
#include "fei_Factory_Trilinos.hpp"

// build it just for fun
#ifndef NO_EXPLICIT_TEMPLATE_INSTANTIATION

template class panzer::BlockedDOFManager<int,long int>;
template class panzer::BlockedDOFManager<int,int>;
template class panzer::BlockedDOFManager<char,long int>;

#endif
