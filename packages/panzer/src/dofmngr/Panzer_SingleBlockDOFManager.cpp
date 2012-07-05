#include "Panzer_config.hpp"

#include "Panzer_SingleBlockDOFManager_impl.hpp"
#include "Panzer_SingleBlockDOFManager_decl.hpp"

#ifdef HAVE_PANZER_EXPLICIT_INSTANTIATION

template class panzer::SingleBlockDOFManager<int,long int>;
template class panzer::SingleBlockDOFManager<int,int>;
//template class panzer::SingleBlockDOFManager<short,int>;
//template class panzer::SingleBlockDOFManager<char,long int>;

#endif
//Explicit Instantiation stuff. Later.
namespace panzer {

} /* panzer */
