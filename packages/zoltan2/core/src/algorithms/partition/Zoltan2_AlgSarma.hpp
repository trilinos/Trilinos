/**
 *   Created by mbenlioglu on Aug 31, 2020.
 */

#pragma once

#include <Teuchos_ParameterEntryValidator.hpp>
#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_CoordinateModel.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_Util.hpp>
#include <Zoltan2_TPLTraits.hpp>

#include <sarma/sarma.hpp>
#include <map>
#include <utility>
#include <vector>
#include <iostream>
#include <fstream>

#ifndef HAVE_ZOLTAN2_SARMA

namespace Zoltan2 {
    template<typename Adapter>
    class AlgSarma : public Algorithm<Adapter> {
    public:
        using base_adapter_t = typename Adapter::base_adapter_t;

        AlgSarma(const RCP<const Environment> &/* env */,
                const RCP<const Comm<int> > &/* problemComm */,
                const RCP<const base_adapter_t> &/* adapter */) {
            throw std::runtime_error("SARMA requested but not enabled into Zoltan2\n"
                                     "Please set CMake flag TPL_ENABLE_SARMA:BOOL=ON");
        }
    };
}

#else


#endif