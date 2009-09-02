#include "EvilBase_decl.hpp"

#ifdef DO_EXPLICIT_INSTANTIATION

#include "EvilBase_def.hpp"

namespace EvilPack {

template class EvilBase<double>;
template class EvilBase<int>;

} // namespace EvilPack

#endif // DO_EXPLICIT_INSTANTIATION
