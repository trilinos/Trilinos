#ifndef __TACHO_MODIFY_DIAGONALS_HPP__
#define __TACHO_MODIFY_DIAGONALS_HPP__

/// \file Tacho_ModifyDiagonals.hpp
/// \brief Front interface for ModifyDiagonals
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

    ///
    /// ModifyDiagonals
    ///

    /// various implementation for different uplo and algo parameters
    template<typename ArgAlgo>
    struct ModifyDiagonals;

}

#endif
