#ifndef STK_SEARCH_UTIL_INTREPID2_HasParamDist_hpp
#define STK_SEARCH_UTIL_INTREPID2_HasParamDist_hpp

#include <stk_util/stk_config.h>

#ifdef STK_HAVE_INTREPID2

#include "Intrepid2_CellTools.hpp"
#include <type_traits>

//-- hacky?? way to detect whether Intrepid2 has the ParametricDistance struct
//
//First, create a forward-declaration of ParametricDistance:
namespace Intrepid2 {
template<unsigned CellTopoKey>
struct ParametricDistance;
}//namespace Intrepid2


namespace stk::search {

// Base case: Assume ParametricDistance doesn't exist:
template <typename T, typename = void>
struct is_complete_type : std::false_type {};

// SFINAE case: If sizeof(Intrepid2::ParametricDistance<shards::Line<>::key>) is valid, the type is complete
template <typename T>
struct is_complete_type<T, std::void_t<decltype(sizeof(T))>> : std::true_type {};

inline constexpr bool intrepid2_has_param_dist = is_complete_type<Intrepid2::ParametricDistance<shards::Line<>::key>>::value;

}

#endif
#endif
