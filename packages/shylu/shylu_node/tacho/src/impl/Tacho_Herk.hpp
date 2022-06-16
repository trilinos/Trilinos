#ifndef __TACHO_HERK_HPP__
#define __TACHO_HERK_HPP__

/// \file Tacho_Herk.hpp
/// \brief Front interface for Herk operators
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

///
/// Herk:
///

/// various implementation for different uplo and algo parameters
template <typename ArgUplo, typename ArgTrans, typename ArgAlgo> struct Herk;

struct HerkAlgorithm {
  using type = ActiveAlgorithm::type;
};

} // namespace Tacho

#endif
