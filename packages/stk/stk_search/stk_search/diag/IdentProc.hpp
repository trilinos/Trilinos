#ifndef stk_search_diag_IdentProc_hpp
#define stk_search_diag_IdentProc_hpp

#include <stk_util/diag/Writer.hpp>
#include <stk_search/IdentProc.hpp>


namespace stk {
namespace search {
namespace ident {

template <class K, class P>
stk::diag::Writer &operator<<(stk::diag::Writer &dout, const IdentProc<K, P> &ident_proc) {
  if (dout.shouldPrint()) {
    dout << "id " << std::hex << ident_proc.ident << std::dec << ", proc " << ident_proc.proc;
  }

  return dout;
}

} // namespace ident
} // namespace search
} // namespace stk

#endif // stk_search_diag_IdentProc_hpp
