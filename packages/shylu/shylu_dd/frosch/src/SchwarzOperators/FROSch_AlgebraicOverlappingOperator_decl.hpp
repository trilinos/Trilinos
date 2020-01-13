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

#ifndef _FROSCH_ALGEBRAICOVERLAPPINGOPERATOR_DECL_HPP
#define _FROSCH_ALGEBRAICOVERLAPPINGOPERATOR_DECL_HPP

#include <FROSch_OverlappingOperator_def.hpp>


namespace FROSch {
    
    using namespace Teuchos;
    using namespace Xpetra;

    enum AddingLayersStrategy {LayersFromMatrix=0,LayersFromGraph=1,LayersOld=2};

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = KokkosClassic::DefaultNode::DefaultNodeType>
    class AlgebraicOverlappingOperator : public OverlappingOperator<SC,LO,GO,NO> {

    protected:

        using CommPtr               = typename SchwarzOperator<SC,LO,GO,NO>::CommPtr;

        using XMapPtr               = typename SchwarzOperator<SC,LO,GO,NO>::XMapPtr;
        using ConstXMapPtr          = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtr;

        using XMatrixPtr            = typename SchwarzOperator<SC,LO,GO,NO>::XMatrixPtr;
        using ConstXMatrixPtr       = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMatrixPtr;

        using ConstXCrsGraphPtr     = typename SchwarzOperator<SC,LO,GO,NO>::ConstXCrsGraphPtr;

        using ParameterListPtr      = typename SchwarzOperator<SC,LO,GO,NO>::ParameterListPtr;

    public:

        AlgebraicOverlappingOperator(ConstXMatrixPtr k,
                                     ParameterListPtr parameterList);

        virtual int initialize()
        {
            FROSCH_ASSERT(false,"AlgebraicOverlappingOperator cannot be built without input parameters.");
            return 0;
        };

        int initialize(int overlap,
                       ConstXMapPtr repeatedMap = null);

        int compute();

        void describe(FancyOStream &out,
                      const EVerbosityLevel verbLevel=Describable::verbLevel_default) const;

        std::string description() const;

    protected:

        int buildOverlappingMatrices(int overlap,
                                     ConstXMapPtr repeatedMap);
        
        virtual int updateLocalOverlappingMatrices();
        

        AddingLayersStrategy AddingLayersStrategy_ = LayersFromGraph;
    };

}

#endif
