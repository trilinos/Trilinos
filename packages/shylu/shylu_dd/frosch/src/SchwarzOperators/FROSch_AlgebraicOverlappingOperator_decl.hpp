// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_ALGEBRAICOVERLAPPINGOPERATOR_DECL_HPP
#define _FROSCH_ALGEBRAICOVERLAPPINGOPERATOR_DECL_HPP

#include <FROSch_OverlappingOperator_def.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    enum AddingLayersStrategy {LayersFromMatrix=0,LayersFromGraph=1,LayersOld=2};

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
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
        };

        int initialize(int overlap,
                       ConstXMapPtr repeatedMap = null);

        int compute();

        void describe(FancyOStream &out,
                      const EVerbosityLevel verbLevel=Describable::verbLevel_default) const;

        string description() const;

    protected:

        int buildOverlappingMatrices(int overlap,
                                     ConstXMapPtr repeatedMap);

        virtual int updateLocalOverlappingMatrices();
        virtual int updateLocalOverlappingMatrices_Symbolic();

        virtual void extractLocalSubdomainMatrix_Symbolic();

        AddingLayersStrategy AddingLayersStrategy_ = LayersFromGraph;
    };

}

#endif
