// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_COARSEOPERATOR_DECL_HPP
#define _FROSCH_COARSEOPERATOR_DECL_HPP

#include <FROSch_SchwarzOperator_def.hpp>

// #define FROSCH_COARSEOPERATOR_DETAIL_TIMERS
// #define FROSCH_COARSEOPERATOR_EXPORT_AND_IMPORT

// TODO: Member sortieren!?

#ifdef HAVE_SHYLU_DDFROSCH_EPETRAEXT
#include "EpetraExt_OperatorOut.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_RowMatrixOut.h"
#endif

#ifdef HAVE_SHYLU_DDFROSCH_ZOLTAN2
#include <Zoltan2_MatrixAdapter.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#endif


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class CoarseOperator : public SchwarzOperator<SC,LO,GO,NO> {

    protected:

        using CommPtr                       = typename SchwarzOperator<SC,LO,GO,NO>::CommPtr;

        using XMap                          = typename SchwarzOperator<SC,LO,GO,NO>::XMap;
        using XMapPtr                       = typename SchwarzOperator<SC,LO,GO,NO>::XMapPtr;
        using ConstXMapPtr                  = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtr;
        using XMapPtrVecPtr                 = typename SchwarzOperator<SC,LO,GO,NO>::XMapPtrVecPtr;
        using ConstXMapPtrVecPtr            = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtrVecPtr;
        using ConstXMapPtrVecPtr2D          = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtrVecPtr2D;

        using XMatrixPtr                    = typename SchwarzOperator<SC,LO,GO,NO>::XMatrixPtr;
        using ConstXMatrixPtr               = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMatrixPtr;

        using XCrsGraph                     = typename SchwarzOperator<SC,LO,GO,NO>::XCrsGraph;
        using GraphPtr                      = typename SchwarzOperator<SC,LO,GO,NO>::GraphPtr;
        using ConstXCrsGraphPtr             = typename SchwarzOperator<SC,LO,GO,NO>::ConstXCrsGraphPtr;

        using XMultiVector                  = typename SchwarzOperator<SC,LO,GO,NO>::XMultiVector;
        using XMultiVectorPtr               = typename SchwarzOperator<SC,LO,GO,NO>::XMultiVectorPtr;
        using XMultiVectorPtrVecPtr         = typename SchwarzOperator<SC,LO,GO,NO>::XMultiVectorPtrVecPtr;
        using ConstXMultiVectorPtr          = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMultiVectorPtr;
        using ConstXMultiVectorPtrVecPtr    = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMultiVectorPtrVecPtr;

        using XImport                       = typename SchwarzOperator<SC,LO,GO,NO>::XImport;
        using XImportPtr                    = typename SchwarzOperator<SC,LO,GO,NO>::XImportPtr;
        using XImportPtrVecPtr              = typename SchwarzOperator<SC,LO,GO,NO>::XImportPtrVecPtr;

        using XExport                       = typename SchwarzOperator<SC,LO,GO,NO>::XExport;
        using XExportPtr                    = typename SchwarzOperator<SC,LO,GO,NO>::XExportPtr;
        using XExportPtrVecPtr              = typename SchwarzOperator<SC,LO,GO,NO>::XExportPtrVecPtr;

        using ParameterListPtr              = typename SchwarzOperator<SC,LO,GO,NO>::ParameterListPtr;

        using CoarseSpacePtr                = typename SchwarzOperator<SC,LO,GO,NO>::CoarseSpacePtr;

        using SolverPtr                     = typename SchwarzOperator<SC,LO,GO,NO>::SolverPtr;
        using SolverFactoryPtr              = typename SchwarzOperator<SC,LO,GO,NO>::SolverFactoryPtr;

        using UN                            = typename SchwarzOperator<SC,LO,GO,NO>::UN;

        using IntVec                        = typename SchwarzOperator<SC,LO,GO,NO>::IntVec;
        using IntVec2D                      = typename SchwarzOperator<SC,LO,GO,NO>::IntVec2D;

        using GOVec                         = typename SchwarzOperator<SC,LO,GO,NO>::GOVec;
        using GOVecPtr                      = typename SchwarzOperator<SC,LO,GO,NO>::GOVecPtr;

        using LOVec                         = typename SchwarzOperator<SC,LO,GO,NO>::LOVec;
        using LOVecPtr2D                    = typename SchwarzOperator<SC,LO,GO,NO>::LOVecPtr2D;

        using SCVec                         = typename SchwarzOperator<SC,LO,GO,NO>::SCVec;

        using ConstLOVecView                = typename SchwarzOperator<SC,LO,GO,NO>::ConstLOVecView;

        using ConstGOVecView                = typename SchwarzOperator<SC,LO,GO,NO>::ConstGOVecView;

        using ConstSCVecView                = typename SchwarzOperator<SC,LO,GO,NO>::ConstSCVecView;

    public:

        CoarseOperator(ConstXMatrixPtr k,
                       ParameterListPtr parameterList);

        ~CoarseOperator();

        virtual int initialize() = 0;

        virtual int compute();

        virtual ConstXMapPtr computeCoarseSpace(CoarseSpacePtr coarseSpace) = 0;

        virtual int clearCoarseSpace();

        virtual void apply(const XMultiVector &x,
                           XMultiVector &y,
                           bool usePreconditionerOnly,
                           ETransp mode=NO_TRANS,
                           SC alpha=ScalarTraits<SC>::one(),
                           SC beta=ScalarTraits<SC>::zero()) const;

        virtual void applyPhiT(const XMultiVector& x,
                               XMultiVector& y) const;

        virtual void applyCoarseSolve(XMultiVector& x,
                                      XMultiVector& y,
                                      ETransp mode=NO_TRANS) const;

        virtual void applyPhi(const XMultiVector& x,
                              XMultiVector& y) const;

        virtual CoarseSpacePtr getCoarseSpace() const;

        //Repeated Coarse map
        virtual int buildElementNodeList() = 0;

        virtual int buildGlobalGraph(Teuchos::RCP<DDInterface<SC,LO,GO,NO> > theDDInterface_) = 0;

        virtual int buildCoarseGraph() = 0;

        // AH: Can this be moved to protected?
        virtual XMapPtr BuildRepeatedMapCoarseLevel(ConstXMapPtr &nodesMap,
                                                    UN dofsPerNode,
                                                    ConstXMapPtrVecPtr dofsMaps,
                                                    UN partition) = 0;

    protected:

        virtual int setUpCoarseOperator();

        XMatrixPtr buildCoarseMatrix();

        int buildCoarseSolveMap(ConstXMapPtr coarseMapUnique);

        virtual void extractLocalSubdomainMatrix_Symbolic() = 0;

        CommPtr CoarseSolveComm_;

        bool OnCoarseSolveComm_ = false;

        int NumProcsCoarseSolve_ = 0;

        CoarseSpacePtr CoarseSpace_;

        XMatrixPtr Phi_;
        XMatrixPtr CoarseMatrix_;

        // result of ExtractLocalSubdomainMatrix_Symbolic
        bool coarseExtractLocalSubdomainMatrix_Symbolic_Done_ = false;
        XMatrixPtr coarseSubdomainMatrix_;
        XMatrixPtr coarseLocalSubdomainMatrix_;
        XImportPtr coarseScatter_;

        // Temp Vectors for apply()
        mutable XMultiVectorPtr XTmp_;
        mutable XMultiVectorPtr XCoarse_;
        mutable XMultiVectorPtr XCoarseSolve_;
        mutable XMultiVectorPtr XCoarseSolveTmp_;
        mutable XMultiVectorPtr YTmp_;
        mutable XMultiVectorPtr YCoarse_;
        mutable XMultiVectorPtr YCoarseSolve_;
        mutable XMultiVectorPtr YCoarseSolveTmp_;

        ConstXMapPtrVecPtr GatheringMaps_ = ConstXMapPtrVecPtr(0);
        XMapPtrVecPtr MLGatheringMaps_ = XMapPtrVecPtr(0); // AH: Could this be const as well?

        // AH: Improve the naming
        XMapPtr CoarseMap_;
        XMapPtr CoarseSolveMap_;
        XMapPtr CoarseSolveRepeatedMap_;
        XMapPtr RepMapCoarse_;
        XMapPtr MLCoarseMap_;

        GraphPtr SubdomainConnectGraph_;
        GraphPtr ElementNodeList_;

        // What are the default values?
        UN CoarseDofsPerNode_;

        UN PartitionType_;

        SolverPtr CoarseSolver_;

        ParameterListPtr DistributionList_;

        XExportPtrVecPtr CoarseSolveExporters_ = XExportPtrVecPtr(0);
        XExportPtrVecPtr MLCoarseSolveExporters_;

#ifdef FROSCH_COARSEOPERATOR_EXPORT_AND_IMPORT
        XImportPtrVecPtr CoarseSolveImporters_ = XImportPtrVecPtr(0);
#endif
    };

}

#endif
