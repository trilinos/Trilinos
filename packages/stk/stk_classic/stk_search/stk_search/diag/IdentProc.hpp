/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_search_diag_IdentProc_hpp
#define stk_search_diag_IdentProc_hpp

#include <stk_util/diag/Writer.hpp>
#include <stk_search/IdentProc.hpp>


namespace stk_classic {
namespace search {
namespace ident {

template <class K, class P>
stk_classic::diag::Writer &operator<<(stk_classic::diag::Writer &dout, const IdentProc<K, P> &ident_proc) {
  if (dout.shouldPrint()) {
    dout << "id " << std::hex << ident_proc.ident << std::dec << ", proc " << ident_proc.proc;
  }

  return dout;
}

} // namespace ident
} // namespace search
} // namespace stk_classic

#endif // stk_search_diag_IdentProc_hpp
