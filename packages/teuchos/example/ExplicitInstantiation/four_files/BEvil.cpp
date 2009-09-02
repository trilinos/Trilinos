#include "BEvil_decl.hpp"

#ifdef DO_EXPLICIT_INSTANTIATION

#include "BEvil_def.hpp"

namespace EvilPack {

template class BEvil<double>;
template class BEvil<int>;

} // namespace EvilPack

#endif // DO_EXPLICIT_INSTANTIATION
