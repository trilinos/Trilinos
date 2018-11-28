//@HEADER
// ************************************************************************
//
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Alexander Heinlein (alexander.heinlein@uni-koeln.de)
//
// ************************************************************************
//@HEADER

#ifndef THYRA_FROSCH_XPETRA_FACTORY_HPP
#define THYRA_FROSCH_XPETRA_FACTORY_HPP

#include <ShyLU_DDFROSch_config.h>

#ifdef HAVE_SHYLU_DDFROSCH_STRATIMIKOS
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

#include "Thyra_FROSchFactory_def.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"

#include <string>
#include "Kokkos_DefaultNode.hpp"


namespace Stratimikos {
    
    using namespace Teuchos;
    using namespace Thyra;
    
    template <typename LO = int, typename GO = int, typename NO = KokkosClassic::DefaultNode::DefaultNodeType>
    void enableFROSch (DefaultLinearSolverBuilder& builder,
                       const std::string& stratName = "FROSch")
    {
        const RCP<const ParameterList> precValidParams = sublist(builder.getValidParameters(), "Preconditioner Types");
        
        TEUCHOS_TEST_FOR_EXCEPTION(precValidParams->isParameter(stratName), std::logic_error,
                                   "Stratimikos::enableFROSch cannot add \"" + stratName +"\" because it is already included in builder!");
        
        typedef PreconditionerFactoryBase<double>  Base;
        if (!stratName.compare("FROSch")) {
            typedef FROSchFactory<double, LO, GO, NO> Impl;
            builder.setPreconditioningStrategyFactory(abstractFactoryStd<Base, Impl>(), stratName);
        }
    }
    
} // namespace Stratimikos

#endif

#endif
