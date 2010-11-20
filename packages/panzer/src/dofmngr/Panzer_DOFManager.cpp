#ifndef NO_EXPLICIT_TEMPLATE_INSTANTIATION

#include "Panzer_DOFManager.hpp"

template class panzer::DOFManager<int,long int>;
template class panzer::DOFManager<int,int>;
template class panzer::DOFManager<char,long int>;

#endif
