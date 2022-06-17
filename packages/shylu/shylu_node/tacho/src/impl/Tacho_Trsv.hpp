#ifndef __TACHO_TRSV_HPP__
#define __TACHO_TRSV_HPP__

/// \file Tacho_Trsv.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

///
/// Trsv:
///

/// various implementation for different uplo and algo parameters
template <typename ArgUplo, typename ArgTrans, typename ArgAlgo> struct Trsv;

struct TrsvAlgorithm {
  using type = ActiveAlgorithm::type;
};

} // namespace Tacho

#endif
