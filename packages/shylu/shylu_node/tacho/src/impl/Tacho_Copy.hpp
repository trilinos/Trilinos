#ifndef __TACHO_COPY_HPP__
#define __TACHO_COPY_HPP__

/// \file Tacho_Copy.hpp
/// \brief Front interface for Copy
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

    ///
    /// Copy
    ///

    /// various implementation for different uplo and algo parameters
    template<typename ArgAlgo>
    struct Copy;
}

#endif
