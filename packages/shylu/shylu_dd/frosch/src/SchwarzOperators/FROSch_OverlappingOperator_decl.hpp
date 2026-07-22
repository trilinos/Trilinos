// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_OVERLAPPINGOPERATOR_DECL_HPP
#define _FROSCH_OVERLAPPINGOPERATOR_DECL_HPP

#include <FROSch_SchwarzOperator_def.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class OverlappingOperator : public SchwarzOperator<SC,LO,GO,NO> {

    protected:

        using CommPtr               = typename SchwarzOperator<SC,LO,GO,NO>::CommPtr;

        using XMapPtr               = typename SchwarzOperator<SC,LO,GO,NO>::XMapPtr;
        using ConstXMapPtr          = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtr;

        using XMatrixPtr            = typename SchwarzOperator<SC,LO,GO,NO>::XMatrixPtr;
        using ConstXMatrixPtr       = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMatrixPtr;

        using XMultiVector          = typename SchwarzOperator<SC,LO,GO,NO>::XMultiVector;
        using XMultiVectorPtr       = typename SchwarzOperator<SC,LO,GO,NO>::XMultiVectorPtr;

        using XImportPtr            = typename SchwarzOperator<SC,LO,GO,NO>::XImportPtr;
        using XExportPtr            = typename SchwarzOperator<SC,LO,GO,NO>::XExportPtr;

        using ParameterListPtr      = typename SchwarzOperator<SC,LO,GO,NO>::ParameterListPtr;

        using SolverPtr             = typename SchwarzOperator<SC,LO,GO,NO>::SolverPtr;
        using SolverFactoryPtr      = typename SchwarzOperator<SC,LO,GO,NO>::SolverFactoryPtr;

        using SCVecPtr              = typename SchwarzOperator<SC,LO,GO,NO>::SCVecPtr;
        using ConstSCVecPtr         = typename SchwarzOperator<SC,LO,GO,NO>::ConstSCVecPtr;

        using UN                    = typename SchwarzOperator<SC,LO,GO,NO>::UN;

    public:

        OverlappingOperator(ConstXMatrixPtr k,
                            ParameterListPtr parameterList);

        ~OverlappingOperator();

        virtual int initialize() = 0;

        virtual int compute() = 0;

        virtual void apply(const XMultiVector &x,
                           XMultiVector &y,
                           bool usePreconditionerOnly,
                           ETransp mode=NO_TRANS,
                           SC alpha=ScalarTraits<SC>::one(),
                           SC beta=ScalarTraits<SC>::zero()) const;

    protected:

        enum CombinationType {Averaging,Full,Restricted};

        virtual int initializeOverlappingOperator();

        virtual int initializeSubdomainSolver(ConstXMatrixPtr localMat);

        virtual int computeOverlappingOperator();

        virtual int updateLocalOverlappingMatrices() = 0;


        ConstXMatrixPtr OverlappingMatrix_;

        ConstXMapPtr OverlappingMap_;

        // Temp Vectors for apply()
        mutable XMultiVectorPtr XTmp_;
        mutable XMultiVectorPtr XOverlap_;
        mutable XMultiVectorPtr XOverlapTmp_;
        mutable XMultiVectorPtr YOverlap_;

        XImportPtr Scatter_;

        SolverPtr SubdomainSolver_;

        XMultiVectorPtr Multiplicity_;

        CombinationType Combine_ = Averaging;
    };

}

#endif
